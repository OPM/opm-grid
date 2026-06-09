/*
  Copyright 2026 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "config.h"

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/DefaultGeometryPolicy.hpp>
#include <opm/grid/cpgrid/Entity.hpp>
#include <opm/grid/cpgrid/LgrFaultHelpers.hpp>
#include <opm/grid/cpgrid/LgrHelpers.hpp>
#include <opm/grid/cpgrid/LevelCartesianIndexMapper.hpp>
#include <opm/grid/cpgrid/ParentToChildCellToPointGlobalIdHandle.hpp>
#include <opm/grid/cpgrid/OrientedEntityTable.hpp>
#include <opm/grid/utility/OpmLog.hpp>

#include <algorithm>    // for std::max
#include <array>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>  // for std::integral_constant
#include <unordered_map>
#include <utility> // for std::pair
#include <vector>

namespace Opm
{
namespace Lgr
{

std::array<std::vector<int>,6> groupFaceIndicesByType(const Dune::cpgrid::CpGridData& gridData,
                                                       const Dune::cpgrid::Entity<0>& element)
{
    std::array<std::vector<int>,6> faceIdxsByType{};
    const auto& cellToFace = gridData.cellToFace(element.index());
    
    for (const auto& face : cellToFace) {
        
        const auto faceIdx = face.index();
        const auto tag = gridData.faceTag(faceIdx);
        
        if ((tag == I_FACE) || (tag == J_FACE) || (tag == K_FACE)) {
            faceIdxsByType[faceGroupIndex(tag, face.orientation())].push_back(faceIdx);  
        }
        else {
            OPM_THROW(std::logic_error, "Face " + std::to_string(faceIdx) +
                      " has an invalid classification: expected one of +/-I, +/-J, or +/-K");
        }
    }
    return faceIdxsByType;
}

bool hasAtMostOneFacePerGroup(const std::array<std::vector<int>,6>& classifiedFaces)
{
    for (const auto& faces : classifiedFaces) {
        if (faces.size()> 1)
            return false;
    }
    return true;
}

std::tuple<Dune::FieldVector<double,3>, double, Dune::FieldVector<double,3>>
computeFaceCenterAreaNormal(const std::vector<Dune::FieldVector<double,3>>& faceToCoord)
{
    assert(faceToCoord.size() == 4); // for now, only quadrilateral face 

    Dune::FieldVector<double,3> faceCenter = {0., 0.,0.};
    for (const auto& coord : faceToCoord) {
        faceCenter += coord;
    }
    faceCenter /= 4.;
    
    // Calculate face area by adding the 4 areas of the triangles partitioning the face.
    double faceArea = 0.0;
    for (std::size_t i = 0; i < 4; ++i) {
        const auto& p0 = faceToCoord[i];
        const auto& p1 = faceToCoord[(i + 1) % 4];
        
        Dune::cpgrid::Geometry<0,3>::GlobalCoordinate trianCorners[3] = {p0, p1, faceCenter};
        faceArea += std::abs(Dune::area(trianCorners));
    }
    
    Dune::FieldVector<double,3> faceNormal = {0., 0.,0.};
    // Construct two vectors that lie on the face, e.g. difference of two conners with the face center. 
    // then obtain an orthogonal vector to both of them. Finally, normalize.
    Dune::cpgrid::Geometry<0,3>::GlobalCoordinate v0 = faceToCoord[0] - faceCenter;
    Dune::cpgrid::Geometry<0,3>::GlobalCoordinate v1 = faceToCoord[1] - faceCenter;
    faceNormal = Dune::FMatrixHelp::Impl::crossProduct(v0, v1);

    assert(faceNormal.two_norm()>0);
    faceNormal /= faceNormal.two_norm();
    
    return {faceCenter, faceArea, faceNormal};
}

std::optional<std::pair<Dune::FieldVector<double,3>, Dune::FieldVector<double,3>>>
computeSegmentIntersection(const Dune::FieldVector<double,3>& startA, const Dune::FieldVector<double,3>& endA,
                           const Dune::FieldVector<double,3>& startB, const Dune::FieldVector<double,3>& endB,
                           bool& isInteriorInA,
                           bool& isInteriorInB)
{
    const auto directionA = endA - startA;
    const auto directionB = endB - startB;
    constexpr double eps = 1e-8; // tolerance

    const auto c = Dune::FMatrixHelp::Impl::crossProduct(directionA, directionB);
    // dot(c, directionA) = 0, dot(c, directionB) = 0

    // Segment A: startA + t*directionA, t in [0,1]
    // Segment B: startB + s*directionB, s in [0,1]

    // Discard skew segments
    const auto coplanar = Dune::dot(startB - startA, c);
    if (std::abs(coplanar) > eps){ // segments are skew (no intersection)
        return std::nullopt;
    }

    // Parallel segments
    if (c.two_norm() < eps) {
        // Check if colinear
        const auto crossOffset = Dune::FMatrixHelp::Impl::crossProduct(startB - startA, directionA);

        if (crossOffset.two_norm() >  eps) {
            return std::nullopt; // parallel but not colinear
        }

        // Parallel and colinear: project onto the largest segment
        const auto lenA = Dune::dot(directionA, directionA);
        const auto lenB = Dune::dot(directionB, directionB);

        if (lenA <  eps  || lenB < eps) {
            return std::nullopt; // degenerate segment
        }
        if (lenA > lenB) {
            auto project= [&](const Dune::FieldVector<double,3>& p) {
                return  Dune::dot(p - startA, directionA) / lenA;
            };

            double t0 = project(startB);
            double t1 = project(endB);

            if (t0 > t1) std::swap(t0, t1);

            double t_min = std::max(0.0, t0);
            double t_max = std::min(1.0, t1);

            if (t_min > t_max /*+ 1e-8*/) {
                return std::nullopt; // disjoint
            }

            // Single-point intersection
            if (std::abs(t_min - t_max) < 1e-8) {
                auto p = startA + t_min * directionA;

                isInteriorInA = (t_min > 0  && t_min < 1.0);
                isInteriorInB = true; // approximate

                return std::make_optional<std::pair<Dune::FieldVector<double,3>,Dune::FieldVector<double,3>>>(p,p);
            }

            // Overlapping segment
            Dune::FieldVector<double,3> p0 = startA + t_min * directionA;
            Dune::FieldVector<double,3> p1 = startA + t_max * directionA;

            isInteriorInA =  (t_min > 0  && t_min < 1.0);
            isInteriorInB = true;// (t_max > 0  && t_max < 1.0); // recall segmentB overlap in segmentA

            return std::make_optional<std::pair<Dune::FieldVector<double,3>,Dune::FieldVector<double,3>>>(p0,p1);
        }
        else {

            auto project= [&](const Dune::FieldVector<double,3>& p) {
                return  Dune::dot(p - startB, directionB) / lenB;
            };

            double t0 = project(startA);
            double t1 = project(endA);

            if (t0 > t1) std::swap(t0, t1);

            double t_min = std::max(0.0, t0);
            double t_max = std::min(1.0, t1);

            if (t_min > t_max /*+ 1e-8*/) {
                return std::nullopt; // disjoint
            }

            // Single-point intersection
            if (std::abs(t_min - t_max) < 1e-8) {
                auto p = startB + t_min * directionB;

                isInteriorInA = true; // // recall segmentB overlap in segmentA
                isInteriorInB = (t_min > 0  && t_min < 1.0);

                return std::make_optional<std::pair<Dune::FieldVector<double,3>,Dune::FieldVector<double,3>>>(p,p);
            }

            // Overlapping segment
            Dune::FieldVector<double,3> p0 = startB + t_min * directionB;
            Dune::FieldVector<double,3> p1 = startB + t_max * directionB;

            isInteriorInA = true;
            isInteriorInB = true;

            return std::make_optional<std::pair<Dune::FieldVector<double,3>,Dune::FieldVector<double,3>>>(p0,p1);

        }
    }
    else {
        // Non-parallel, compute intersection point
        //
        // X denotes the crossProduct and <,> the dot:
        // startA + t*dirA = startB + s*dirB
        //          t*dirA = startB - startA + s*dirB
        //  (t*dirA)X dirB = (startB - startA  + s*dirB)X dirB
        // t*(dirA X dirB) = (startSegmentB - startSegmentA)X dirB + s*(dirB X dirB)
        //             t*c = (startSegmentB - startSegmentA)X dirB + s*0
        //        <t*c, c> = <(startSegmentB - startSegmentA)X dirB, c>
        //         t <c,c> = <(startSegmentB - startSegmentA)X dirB, c>
        double t = Dune::dot(Dune::FMatrixHelp::Impl::crossProduct(startB - startA, directionB), c) / Dune::dot(c,c);
        // Analogously, for s:
        double s = Dune::dot(Dune::FMatrixHelp::Impl::crossProduct(startB - startA, directionA), c) / Dune::dot(c,c);

        if ((t >= 0) && (t <= 1) && (s >= 0) && (s <= 1)) { // segments intersect
            isInteriorInA = (t > 0) && (t < 1);
            isInteriorInB = (s > 0) && (s < 1);
            return std::make_optional<std::pair<Dune::FieldVector<double,3>,Dune::FieldVector<double,3>>>(startA + (t*directionA), startA + (t*directionA));
        } else { // lines intersect, but not the segments
            return std::nullopt;
        }
    }
}

std::vector<std::array<int,2>> createEdges(const auto& faceToPoint)
{
    std::vector<std::array<int,2>> edges{};
    edges.reserve(faceToPoint.size());
    for (std::size_t i = 0; i < faceToPoint.size(); ++i) {
        edges.push_back(std::array<int,2>{ faceToPoint[i], faceToPoint[(i+1)%faceToPoint.size()]});
    }
    return edges;
}

bool inSemiplane(const Dune::FieldVector<double,3>& vertex,
                 const Dune::FieldVector<double,3>& faceVertex,
                 const Dune::FieldVector<double,3>& faceNormal,
                 const Dune::FieldVector<double,3>& directionEdge)
{

    // 1. Check if vertex lies on the plane
    const auto v = vertex - faceVertex;
    if (std::fabs(Dune::dot(faceNormal, v)) > 0)
        return false;

    // 2. Compute in-plane perpendicular vector
    const auto u = Dune::FMatrixHelp::Impl::crossProduct(faceNormal, directionEdge);

    // 3. Check side
    return Dune::dot(u, v) >= 0;
}

bool isVertexInsideFace(const Dune::FieldVector<double,3>& vertex,
                        const Dune::cpgrid::CpGridData& face_gridData,
                        const std::vector<std::array<int,2>>& face_edges,
                        const Dune::FieldVector<double,3>& face_normal)
{
    bool vertexIsInsideFace = true;
    for (const auto& edge : face_edges) {
        const auto edge0 = Dune::cpgrid::Entity<3>(face_gridData, edge[0], true).geometry().center();
        const auto edge1 = Dune::cpgrid::Entity<3>(face_gridData, edge[1], true).geometry().center();
            
        vertexIsInsideFace = vertexIsInsideFace && inSemiplane(vertex, edge0, face_normal, edge1-edge0);
        if (!vertexIsInsideFace)
            return false;
    } 
    return  vertexIsInsideFace;
}


bool isVertexInSegmentInterior(const Dune::FieldVector<double,3>& vertex,
                               const Dune::FieldVector<double,3>& startSegment,
                               const Dune::FieldVector<double,3>& endSegment)
{
    constexpr double tol = 1e-12;
    // segment(t) = (endSegment - startSegment)*t + startSegment, t in [0,1] 
    const double segLen2 = (endSegment - startSegment).two_norm2();

    // Degenerate segment
    if (segLen2 < tol)
        return false;

    // Check collinearity
    const auto cp = Dune::FMatrixHelp::Impl::crossProduct(vertex - startSegment,
                                                          endSegment - startSegment);

    // Scale tolerance with segment length (if collinear, d*vertex ~ cp -> d^2 vertex.two_norm ~)
    if (cp.two_norm2() > tol * tol * segLen2)
        return false;

    // Position along the segment
    const double t = Dune::dot(vertex - startSegment, endSegment - startSegment) / segLen2;

    return (t > tol) && (t < 1.0 - tol);
}

std::vector<int> isVertexInElementFaceEdge(const Dune::FieldVector<double,3>& vertex,
                                                           int elemIdx,
                                                           const Dune::cpgrid::CpGridData& gridData)
{
    std::vector<int> isInFace{};

    bool isInAtLeastOneEdge = false;
    for (const auto& face : gridData.cellToFace(elemIdx)) {
        const auto edges = createEdges(gridData.faceToPoint(face.index()));
        for (const auto& edge : edges) {
            const auto edge0 = Dune::cpgrid::Entity<3>(gridData, edge[0], true).geometry().center();
            const auto edge1 = Dune::cpgrid::Entity<3>(gridData, edge[1], true).geometry().center();
            isInAtLeastOneEdge = isInAtLeastOneEdge || isVertexInSegmentInterior(vertex, edge0, edge1);
            if (isVertexInSegmentInterior(vertex, edge0, edge1))
                isInFace.push_back(face.index());
        }
    }
    return isInFace;
}

std::vector<int> isVertexInElementEdge(const Dune::FieldVector<double,3>& vertex,
                                       int elemIdx,
                                       const Dune::cpgrid::CpGridData& gridData)
{
    std::vector<int> isInEdge{};
    
    const auto& cellToPoint = gridData.cellToPoint(elemIdx);
    
    bool isInAtLeastOneEdge = false;
    for (const auto& face : gridData.cellToFace(elemIdx)) {
        const auto edges = createEdges(gridData.faceToPoint(face.index()));
        for (const auto& edge : edges) {
            
            bool edgeHasNoCellVertex = std::find(cellToPoint.begin(), cellToPoint.end(), edge[0]) == cellToPoint.end() &&
                std::find(cellToPoint.begin(), cellToPoint.end(), edge[1]) == cellToPoint.end();
            
            if (edgeHasNoCellVertex) 
                continue;
            
            const auto edge0 = Dune::cpgrid::Entity<3>(gridData, edge[0], true).geometry().center();
            const auto edge1 = Dune::cpgrid::Entity<3>(gridData, edge[1], true).geometry().center();
            isInAtLeastOneEdge = isInAtLeastOneEdge || isVertexInSegmentInterior(vertex, edge0, edge1);
            if (isVertexInSegmentInterior(vertex, edge0, edge1))
                isInEdge.push_back(face.index());
        }
    }
    return isInEdge;
}

std::set<Dune::FieldVector<double,3>,FieldVectorLess>
collectNewVertices(const Dune::cpgrid::CpGridData&                                                    parentGridData,
                   const Dune::cpgrid::Entity<0>&                                                     parentElem,
                   std::map<Dune::FieldVector<double,3>, int, FieldVectorLess>&                       parentFaceExistingVertex_to_parentGridCornerIdx,
                   const std::array<std::vector<int>,6>&                                              classifiedParentCellFaces,
                   const Dune::cpgrid::CpGridData&                                                    singleCellRefinementData, 
                   const Dune::cpgrid::Entity<0>&                                                     refinedElem,
                   std::map<Dune::FieldVector<double,3>, int, FieldVectorLess>&                       refinedFaceExistingVertex_to_refinedGridCornerIdx,
                   std::map<int,std::vector<std::vector<Dune::FieldVector<double,3>>>>&               overlapFaces,
                   std::vector<std::vector<std::pair<int,std::vector<Dune::FieldVector<double,3>>>>>& vanishedRefinedFaceToNewRefinedFaces,
                   std::vector<int>&                                                                  fullyContainedInParentFace)
{
    std::set<Dune::FieldVector<double,3>,FieldVectorLess> foundNewVertices{};

    const auto& refinedCellToFace = singleCellRefinementData.cellToFace(refinedElem.index());
    const auto& parentCellToFace = parentGridData.cellToFace(parentElem.index());
    
    for (const auto& parentFace : parentCellToFace) {
        // skip if face type is not repeated (i.e. there is only one face of type {face_tag, face_orientation})
        const auto parentFaceIdx = parentFace.index();
        if (classifiedParentCellFaces[faceGroupIndex(parentGridData.faceTag(parentFaceIdx), parentFace.orientation())].size() == 1)
            continue;
        
        for (const auto& refinedFace : refinedCellToFace) {
            const auto refinedFaceIdx = refinedFace.index();
            // Skip face if it is not on the boundary of the single cell refinement grid
            if (singleCellRefinementData.faceToCellSize(refinedFaceIdx) != 1)
                continue;

            bool refinedFaceFullyContainedInParentFace = false;
            const auto newFace = computeFaceOverlapVertices(parentFace,                                       // face1
                                                            parentGridData,                                   // face1_gridData
                                                            parentFaceExistingVertex_to_parentGridCornerIdx,  // face1ExistingVertex_to_face1GridCornerIdx
                                                            refinedFace,                                      // face2
                                                            singleCellRefinementData,                         // face2_gridData       
                                                            refinedFaceExistingVertex_to_refinedGridCornerIdx,// face2ExistingVertex_to_face2GridCornerIdx 
                                                            foundNewVertices,              
                                                            refinedFaceFullyContainedInParentFace);           // face2FullyContainedInFace1
            
            if (refinedFaceFullyContainedInParentFace){
                fullyContainedInParentFace[refinedFaceIdx] = parentFaceIdx;
            }
            
            if (newFace.has_value() && !newFace.value().empty()) {
                
                // Save face vertices/coordinates as one of the overlapping refined faces with the parent cell
                auto it  = overlapFaces.find(parentFaceIdx);
                if (it != overlapFaces.end()) {
                    auto& collectedFaces = *it; 
                    collectedFaces.second.push_back(newFace.value());
                }
                else {
                    auto& collectedFaces = overlapFaces[parentFaceIdx];
                    collectedFaces.push_back(newFace.value());
                }
                 
                if (!refinedFaceFullyContainedInParentFace) {
                    vanishedRefinedFaceToNewRefinedFaces[refinedFaceIdx].push_back(std::make_pair(parentFaceIdx, newFace.value()));
                }
            }
        }
    }
    return foundNewVertices;
}

std::vector<Dune::FieldVector<double,3>>
orderFaceVertices(int faceTag,
                  const std::set<Dune::FieldVector<double,3>,FieldVectorLess>& faceVertices)
{
    assert( (faceTag >= 0) && (faceTag < 3));
    std::vector<Dune::FieldVector<double,3>> orderedVertices{};
    
    if (!faceVertices.empty()) {
        constexpr double eps = 1e-8; // tolerance
        assert(faceVertices.size() == 4); // for general face, do sth else...
        orderedVertices.resize(faceVertices.size());
        
        // Box bounds (I -> yz box, J-> xz box, K-> xy box)
        double minA = std::numeric_limits<double>::max();
        double maxA = std::numeric_limits<double>::lowest();

        double minB = std::numeric_limits<double>::max();
        double maxB = std::numeric_limits<double>::lowest();
        
        for (const auto& vertex : faceVertices) { 
            minA = std::min(minA, vertex[(faceTag+1)%3]);  
            maxA = std::max(maxA, vertex[(faceTag+1)%3]);  
            minB = std::min(minB, vertex[(faceTag+2)%3]);  
            maxB = std::max(maxB, vertex[(faceTag+2)%3]);  
        }
        if (faceTag < 2) { //  Vertex order
            // 0-> I_FACE: jk, (j+1)k, (j+1)(k+1), j(k+1)
            // 1-> J_FACE: (i+1)k, ik, i(k+1), (i+1)(k+1)
            for (const auto& vertex : faceVertices) {
                const auto a = vertex[(faceTag+1)%3];
                const auto b = vertex[(faceTag+2)%3];
                if (b-minB<eps) { 
                    if (a-minA<eps) 
                        orderedVertices[faceTag] = vertex;
                    else 
                        orderedVertices[faceTag+1] = vertex;
                }
                else if (maxB-b<eps) {
                    if (a-minA<eps)
                        orderedVertices[(faceTag+3)%4] = vertex;
                    else 
                        orderedVertices[(faceTag+2)%4] = vertex;
                }
            }
        }
        else { //  Vertex order (no faults in z-direction, maybe remove this case?)      
            // K_FACE: ij, (i+1)j, (i+1)(j+1), (i+1)(j+1)
            for (const auto& vertex : faceVertices) {
                const auto a = vertex[0]; // x
                const auto b = vertex[1]; // y
                if (b-minB<eps) {
                    if (a-minA<eps)
                        orderedVertices[0] = vertex;
                    else
                        orderedVertices[1] = vertex;
                }
                else if (maxB-b<eps) {
                    if (a-minA<eps)
                        orderedVertices[3] = vertex;
                    else
                        orderedVertices[2] = vertex;
                }
            }
        }
    }
    return orderedVertices;
}


} // namespace Lgr
} // namespace Opm
