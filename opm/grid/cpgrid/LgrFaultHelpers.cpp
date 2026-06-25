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
    constexpr double eps = 1e-10;

    isInteriorInA = false;
    isInteriorInB = false;

    const auto dA = endA - startA; // directionA
    const auto dB = endB - startB; // directionB

    const double lenA2 = Dune::dot(dA,dA);
    const double lenB2 = Dune::dot(dB,dB);
    // Degenerate A: startA + t*directionA, t in [0,1]
    if (lenA2 < eps) {
        if (isVertexInSegmentInterior(startA, startB, endB)) { 
            isInteriorInB = Dune::dot(startA-startB,dB)/lenB2 > eps &&
                Dune::dot(startA-startB,dB)/lenB2 < 1-eps;
            return {{startA,startA}};
        }
        return std::nullopt;
    }

    // Degenerate B: startB + s*directionB, s in [0,1]
    if (lenB2 < eps) {
        if (isVertexInSegmentInterior(startB, startA, endA)) { 
            isInteriorInA = Dune::dot(startB-startA,dA)/lenA2 > eps &&
                Dune::dot(startB-startA,dA)/lenA2 < 1-eps;
            return {{startB,startB}};
        }
        return std::nullopt;
    }

    const auto c = Dune::FMatrixHelp::Impl::crossProduct(dA,dB);
    const double c2 = Dune::dot(c,c);
    
    // Parallel
    if (c2 < eps) {// not collinear
        const auto offset = Dune::FMatrixHelp::Impl::crossProduct(startB-startA,dA);

        if (offset.two_norm() > eps)
            return std::nullopt;
        
        // project B onto A
        auto paramA = [&](const auto& p) {
            return Dune::dot(p-startA,dA)/lenA2;
        };

        double t0 = paramA(startB);
        double t1 = paramA(endB);

        if (t0 > t1)
            std::swap(t0,t1);
        
        const double lo = std::max(0.0,t0);
        const double hi = std::min(1.0,t1);

        if (lo > hi + eps)
            return std::nullopt;
        
        auto p0 = startA + lo*dA;
        auto p1 = startA + hi*dA;

        isInteriorInA = (lo > eps && lo < 1-eps) || (hi > eps && hi < 1-eps);

        double s0 = Dune::dot(p0-startB,dB)/lenB2;
        double s1 = Dune::dot(p1-startB,dB)/lenB2;
        
        isInteriorInB = (s0 > eps && s0 < 1-eps) || (s1 > eps && s1 < 1-eps);

        return {{p0,p1}};
    }
    
    // Non-parallel
    const auto r = startB-startA;
    // Compute intersection point
    //
    // X denotes the crossProduct and <,> the dot:
    // startA + t*dirA = startB + s*dirB
    //          t*dirA = startB - startA + s*dirB
    //  (t*dirA)X dirB = (startB - startA  + s*dirB)X dirB
    // t*(dirA X dirB) = (startSegmentB - startSegmentA)X dirB + s*(dirB X dirB)
    //             t*c = (startSegmentB - startSegmentA)X dirB + s*0
    //        <t*c, c> = <(startSegmentB - startSegmentA)X dirB, c>
    //         t <c,c> = <(startSegmentB - startSegmentA)X dirB, c>
    const double t = Dune::dot(Dune::FMatrixHelp::Impl::crossProduct(r,dB), c) / c2;
    const double s = Dune::dot(Dune::FMatrixHelp::Impl::crossProduct(r,dA), c) / c2;

    if ((t < -eps) || (t > 1+eps) || (s < -eps) || (s > 1+eps))
        return std::nullopt;

    isInteriorInA = t > eps && t < 1-eps;
    isInteriorInB = s > eps && s < 1-eps;
    
    const auto p = startA + t*dA;

    return {{p,p}};
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
    const double eps =  1e-8;
    // 1. Check if vertex lies on the plane
    const auto v = vertex - faceVertex;
    if (std::fabs(Dune::dot(faceNormal, v)) > eps)
        return false;

    // 2. Compute in-plane perpendicular vector
    const auto u = Dune::FMatrixHelp::Impl::crossProduct(faceNormal, directionEdge);

    // 3. Check side
    return Dune::dot(u, v) >= -eps*directionEdge.two_norm();
}

bool isVertexInsideFace(const Dune::FieldVector<double,3>& vertex,
                        const Dune::cpgrid::CpGridData& face_gridData,
                        const std::vector<std::array<int,2>>& face_edges,
                        const Dune::FieldVector<double,3>& face_normal)
{
    assert(face_edges.size());
    for (const auto& edge : face_edges) {
        const auto edge0 = Dune::cpgrid::Entity<3>(face_gridData, edge[0], true).geometry().center();
        const auto edge1 = Dune::cpgrid::Entity<3>(face_gridData, edge[1], true).geometry().center();
        
        if (!inSemiplane(vertex, edge0, face_normal, edge1-edge0))
            return false;
    }
    return  true; 
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

std::vector<int> isVertexOnBoundaryExcludingCellEdges(const Dune::FieldVector<double,3>& vertex,
                                                        int elemIdx,
                                                        const Dune::cpgrid::CpGridData& gridData)
{
    std::vector<int> faces{};
    
    for (const auto& face : gridData.cellToFace(elemIdx)) {
        
        const auto edges = createEdges(gridData.faceToPoint(face.index()));
        const auto& faceNormal = gridData.faceNormals(face.index());
        if (!isVertexInsideFace(vertex, gridData, edges, faceNormal))
            continue;

        if (isVertexInElementEdge(vertex, elemIdx, gridData).empty())
            faces.push_back(face.index());
    }
    return faces;
}


std::vector<int> isVertexInElementFaceEdge(const Dune::FieldVector<double,3>& vertex,
                                           int elemIdx,
                                           const Dune::cpgrid::CpGridData& gridData)
{
    std::vector<int> isInFace{};

    const auto& cellToPoint = gridData.cellToPoint(elemIdx);
    
    bool isInAtLeastOneEdge = false;
    for (const auto& face : gridData.cellToFace(elemIdx)) {
        const auto edges = createEdges(gridData.faceToPoint(face.index()));
        for (const auto& edge : edges) {
            
            bool isNotCellEdge = std::find(cellToPoint.begin(), cellToPoint.end(), edge[0]) == cellToPoint.end() ||
                std::find(cellToPoint.begin(), cellToPoint.end(), edge[1]) == cellToPoint.end();
            
            if (!isNotCellEdge) 
                continue;
            
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
                   std::map<Dune::FieldVector<double,3>, int, FieldVectorLess>&                       parentGridData_vertex_to_vertexIdx,
                   const std::array<std::vector<int>,6>&                                              classifiedParentCellFaces,
                   const Dune::cpgrid::CpGridData&                                                    cellRefinementData, 
                   const Dune::cpgrid::Entity<0>&                                                     refinedElem,
                   std::map<Dune::FieldVector<double,3>, int, FieldVectorLess>&                       cellRefinementData_vertex_to_vertexIdx,
                   std::map<int,std::vector<std::vector<Dune::FieldVector<double,3>>>>&               overlapFaces,
                   std::vector<std::vector<std::pair<int,std::vector<Dune::FieldVector<double,3>>>>>& vanishedRefinedFaceToNewRefinedFaces,
                   std::vector<int>&                                                                  fullyContainedInParentFace)
{
    std::set<Dune::FieldVector<double,3>,FieldVectorLess> foundNewVertices{};

    const auto& refinedCellToFace = cellRefinementData.cellToFace(refinedElem.index());
    const auto& parentCellToFace = parentGridData.cellToFace(parentElem.index());
    
    for (const auto& parentFace : parentCellToFace) {
        // skip if face type is not repeated (i.e. there is only one face of type {face_tag, face_orientation})
        const auto parentFaceIdx = parentFace.index();
        if (classifiedParentCellFaces[faceGroupIndex(parentGridData.faceTag(parentFaceIdx), parentFace.orientation())].size() == 1)
            continue;
        
        for (const auto& refinedFace : refinedCellToFace) {
            const auto refinedFaceIdx = refinedFace.index();
            // Skip face if it is not on the boundary of the single cell refinement grid
            if (cellRefinementData.faceToCellSize(refinedFaceIdx) != 1)
                continue;

            bool refinedFaceFullyContainedInParentFace = false;
            const auto newFace = computeFaceOverlapVertices(parentFace,                             // face1
                                                            parentGridData,                         // gridData1
                                                            parentGridData_vertex_to_vertexIdx,     // gridData1_vertex_to_vertexIdx
                                                            refinedFace,                            // face2
                                                            cellRefinementData,                     // gridData2       
                                                            cellRefinementData_vertex_to_vertexIdx, // gridData2_vertex_to_vertexIdx 
                                                            foundNewVertices,                                
                                                            refinedFaceFullyContainedInParentFace); // face2FullyContainedInFace1

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

void assignCorrectChildFace(bool                                             hasOnlyOneFacePerType,
                            int                                              faceIdx,
                            int                                              faceType_count,
                            int                                              child_face,
                            std::vector<int>&                                children_faces,
                            std::vector<int>&                                correctRefinedFace_to_parentFace,
                            const std::vector<std::vector<int>>&             oldToNewFaceIdx,
                            const std::unordered_map<int, std::vector<int>>& vanishedFaceToNewFaces,
                            const std::vector<std::vector<int>>&             parentFace_to_correctRefinedFaces,
                            const std::vector<int>&                          fullyContainedInParentFace)
{
    auto addChildFace = [&](int idx) {
        children_faces.push_back(idx);
        correctRefinedFace_to_parentFace[idx] = faceIdx;
    };

    if (hasOnlyOneFacePerType) {
        addChildFace(child_face);
        return;
    }

    const auto& mappedFaces = oldToNewFaceIdx[child_face];
    
    if (faceType_count == 1) {
        for (int idx : mappedFaces)
            addChildFace(idx);
        return;
    }
    
    auto it = vanishedFaceToNewFaces.find(child_face);

    if (it != vanishedFaceToNewFaces.end()) {
        const auto& newFaces = it->second;
        std::unordered_set<int> newSet(newFaces.begin(), newFaces.end());

        for (int idx : parentFace_to_correctRefinedFaces[faceIdx]) {
            if (newSet.contains(idx))
                addChildFace(idx);
        }
    }
    else if (fullyContainedInParentFace[child_face] == faceIdx) {
        for (int idx : mappedFaces) {
            addChildFace(idx);
        }
    }
}

void provideRefinementParentFaceIdxRelations(bool                                                        hasOnlyOneFacePerType,
                                             const std::array<std::vector<int>, 6>&                      classifiedParentCellFaces,
                                             const Dune::cpgrid::CpGridData&                             parentGridData,
                                             const Dune::cpgrid::Entity<0>&                              parentElem,
                                             const std::array<int,3>&                                    cells_per_dim,
                                             std::vector<int>&                                           correctRefinedFace_to_parentFace,
                                             std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                             const std::vector<std::vector<int>>&                        oldToNewFaceIdx,
                                             const std::unordered_map<int, std::vector<int>>&            vanishedFaceToNewFaces,
                                             const std::vector<std::vector<int>>&                        parentFace_to_correctRefinedFaces,
                                             const std::vector<int>&                                     fullyContainedInParentFace)
{
    // Auxiliary integers to simplify new-born-face-index notation.
    const int k_faces = cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1);
    const int i_faces = (cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2];

    for (const auto& face : parentGridData.cellToFace(parentElem.index())) {
        
        const auto& parent_face_tag = parentGridData.faceTag(face.index());
        const auto faceType_count = classifiedParentCellFaces[faceGroupIndex(parent_face_tag, face.orientation())].size();

        std::vector<int> children_faces; // Cannot reserve/resize "now", it depends of the type of face.
        
        if (parent_face_tag == face_tag::K_FACE) {
#ifndef NDEBUG   
            assert(faceType_count == 1);
#endif    
            children_faces.reserve(cells_per_dim[0]*cells_per_dim[1]);
            
            for (int j = 0; j < cells_per_dim[1]; ++j) {
                for (int i = 0; i < cells_per_dim[0]; ++i) {
                    int child_face;
                    if (!face.orientation()) // false -> BOTTOM FACE -> k=0
                        child_face = (j*cells_per_dim[0]) + i;
                    else // true -> TOP FACE -> k=cells_per_dim[2]
                        child_face = (cells_per_dim[2]*cells_per_dim[0]*cells_per_dim[1]) +(j*cells_per_dim[0]) + i;

                    if (hasOnlyOneFacePerType) {
                        children_faces.push_back(child_face);
                        correctRefinedFace_to_parentFace[child_face] = face.index();
                    }
                    else {
                        for (const auto& newIdx : oldToNewFaceIdx[child_face]) {
                            children_faces.push_back(newIdx);
                            correctRefinedFace_to_parentFace[newIdx] = face.index();
                        }
#ifndef NDEBUG                     
                        assert(oldToNewFaceIdx[child_face].size()==1);
                        auto it = vanishedFaceToNewFaces.find(child_face);
                        assert(it == vanishedFaceToNewFaces.end());
#endif
                    }
                } // i-for-lopp
            } //j-for-loop
        } // if-K_FACE
        else if (parent_face_tag == face_tag::I_FACE) {
            for (int k = 0; k < cells_per_dim[2]; ++k) {
                for (int j = 0; j < cells_per_dim[1]; ++j) {
                    int child_face;
                    if (!face.orientation()) // false -> LEFT FACE -> i=0
                        child_face = k_faces + (k*cells_per_dim[1]) + j;
                    else // true -> RIGHT FACE -> i=cells_per_dim[0]
                        child_face = k_faces + (cells_per_dim[0]*cells_per_dim[1]*cells_per_dim[2]) + (k*cells_per_dim[1]) + j;

                    assignCorrectChildFace(hasOnlyOneFacePerType,
                                           face.index(),
                                           faceType_count,
                                           child_face,
                                           children_faces,
                                           correctRefinedFace_to_parentFace,
                                           oldToNewFaceIdx,
                                           vanishedFaceToNewFaces,
                                           parentFace_to_correctRefinedFaces,
                                           fullyContainedInParentFace);
                } // j-for-loop
            } // k-for-loop
        } // if-I_FACE
        else if (parent_face_tag == face_tag::J_FACE) {
            for (int i = 0; i < cells_per_dim[0]; ++i) {
                for (int k = 0; k < cells_per_dim[2]; ++k) {
                    int child_face;
                    if (!face.orientation()) // false -> FRONT FACE -> j=0
                        child_face = k_faces + i_faces + (i*cells_per_dim[2]) + k;
                    else  // true -> BACK FACE -> j=cells_per_dim[1]
                        child_face = k_faces + i_faces  + (cells_per_dim[1]*cells_per_dim[0]*cells_per_dim[2])
                            + (i*cells_per_dim[2]) + k;

                    assignCorrectChildFace(hasOnlyOneFacePerType,
                                           face.index(),
                                           faceType_count,
                                           child_face,
                                           children_faces,
                                           correctRefinedFace_to_parentFace,
                                           oldToNewFaceIdx,
                                           vanishedFaceToNewFaces,
                                           parentFace_to_correctRefinedFaces,
                                           fullyContainedInParentFace);
                } // k-for-loop
            } // i-for-loop
        } // if-J_FACE
        children_faces.shrink_to_fit();
        faceInMarkedElemAndRefinedFaces[face.index()].push_back(std::make_pair(parentElem.index(), children_faces));  
    }
}

void addVertices(int&                                                               grid1_numVertices, 
                 const Dune::cpgrid::CpGridData&                                    grid1,
                 const std::map<Dune::FieldVector<double,3>, int, FieldVectorLess>& grid2_vertex_to_vertexIdx, 
                 const std::set<Dune::FieldVector<double,3>, FieldVectorLess>&      foundNewVertices, // new for grid1 or grid 2?? 
                 GeomData&                                                          correctedGrid1GeomData,
                 std::map<int,int>&                                                 correctedGrid1VertexIdx_to_grid2VertexIdx,
                 std::map<Dune::FieldVector<double,3>, int, FieldVectorLess>&       newVertex_to_correctedGrid1VertexIdx)
{
    const auto& grid_vertices = *(grid1.getGeometry().geomVector(std::integral_constant<int,3>()));
    Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<0,3>>& correctedGrid1_vertices =
        *(correctedGrid1GeomData.geometries.geomVector(std::integral_constant<int,3>()));
    
    correctedGrid1_vertices.resize(grid1_numVertices + foundNewVertices.size()); 
    
    for (int i = 0; i < grid1.size(3); ++i) { /** Check here if we're adding too many?? */
        correctedGrid1_vertices[i] = grid_vertices.get(i);
    }
  
    for (const auto& vertex : foundNewVertices) { /**Do we also need to check grid1_vertex_to_vertexIdx?? */
        auto it =  grid2_vertex_to_vertexIdx.find(vertex);
        
        assert(it !=  grid2_vertex_to_vertexIdx.end()); 
        auto grid2VertexIdx = it->second;
        
        correctedGrid1VertexIdx_to_grid2VertexIdx[grid1_numVertices] = grid2VertexIdx;
        
        newVertex_to_correctedGrid1VertexIdx[vertex] = grid1_numVertices; 
        correctedGrid1_vertices[grid1_numVertices] = Dune::cpgrid::Geometry<0, 3>(vertex);
        
        ++grid1_numVertices;
    }
}

void addVertices(int&                                                               grid1_numVertices,
                 const Dune::cpgrid::CpGridData&                                    grid1,
                 const std::map<Dune::FieldVector<double,3>, int, FieldVectorLess>& parentGrid_vertex_to_vertexIdx, // grid2_vertex_to_vertexIdx
                 const std::set<Dune::FieldVector<double,3>, FieldVectorLess>&      foundNewVertices, // new for grid1 or grid 2?? 
                 GeomData&                                                          correctedGrid1GeomData,
                 std::map<int,int>&                                                 correctedGrid1VertexIdx_to_parentGridVertexIdx, // correctedGrid1VertexIdx_to_grid2VertexIdx,
                 std::map<Dune::FieldVector<double,3>, int, FieldVectorLess>&       newVertex_to_correctedGrid1VertexIdx,
                 const std::array<int,8>&                                           parentCellToPoint,
                 std::vector<std::array<int,2>>&                                    parentBoundaryVertexIdx_to_correctedGrid1BoundaryVertexIdx)
{
    const auto& grid1_vertices = *(grid1.getGeometry().geomVector(std::integral_constant<int,3>()));
    Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<0,3>>& correctedGrid1_vertices =
        *(correctedGrid1GeomData.geometries.geomVector(std::integral_constant<int,3>()));
    correctedGrid1_vertices.resize(grid1_numVertices + foundNewVertices.size());
    
    for (int i = 0; i < grid1.size(3); ++i) {
        correctedGrid1_vertices[i] = grid1_vertices.get(i);
    }
        
    for (const auto& vertex : foundNewVertices) {
        auto it = parentGrid_vertex_to_vertexIdx.find(vertex);
        
        if (it != parentGrid_vertex_to_vertexIdx.end()) {
            auto parentGridVertexIdx = it->second;
            if (std::find(parentCellToPoint.begin(), parentCellToPoint.end(), parentGridVertexIdx) == parentCellToPoint.end()) {
                
                parentBoundaryVertexIdx_to_correctedGrid1BoundaryVertexIdx.push_back(std::array<int,2>{parentGridVertexIdx, grid1_numVertices});
                correctedGrid1VertexIdx_to_parentGridVertexIdx[grid1_numVertices] = parentGridVertexIdx;

                newVertex_to_correctedGrid1VertexIdx[vertex] = grid1_numVertices;
                correctedGrid1_vertices[grid1_numVertices] = Dune::cpgrid::Geometry<0, 3>(vertex);
            }
        }

        newVertex_to_correctedGrid1VertexIdx[vertex] = grid1_numVertices;
        correctedGrid1_vertices[grid1_numVertices] = Dune::cpgrid::Geometry<0, 3>(vertex);

        ++grid1_numVertices;
    }
}

bool isAtGridBoundary(const Dune::cpgrid::CpGridData& gridData,
                      const Dune::cpgrid::Entity<0>& element)
{
    const auto& cellToFace = gridData.cellToFace(element.index());
    for (const auto& face : cellToFace) {
        if (gridData.faceToCellSize(face.index()) == 1)
            return true;
    }
    return false;
}


void addFaces(int                                                                                      grid1_numFaces,
              int                                                                                      grid2_numFaces,
              const std::map<int,std::vector<std::vector<Dune::FieldVector<double,3>>>>&               grid1_overlapFaces,
              const Dune::cpgrid::CpGridData&                                                          grid1,
              const std::vector<std::vector<std::pair<int,std::vector<Dune::FieldVector<double,3>>>>>& grid1_vanishedFaceIdx_to_newFaces,
              std::vector<std::vector<int>>&                                                           grid1FaceIdx_to_correctedGrid1FaceIdx, 
              std::unordered_map<int, int>&                                                            correctedGrid1FaceIdx_to_grid1FaceIdx,
              GeomData&                                                                                correctedGrid1GeomData,
              std::vector<std::vector<int>>&                                                           grid2VanishedFaceIdx_to_correctedGrid1FaceIndices,
              const  std::map<Dune::FieldVector<double,3>, int, FieldVectorLess>&                      newVertex_to_correctedGrid1VertexIdx,
              const std::map<Dune::FieldVector<double,3>, int, FieldVectorLess>&                       grid1_vertex_to_vertexIdx,
              std::unordered_map<int, std::vector<int>>&                                               grid1VanishedFaceIdx_to_correctedGrid1FaceIndices,
              std::unordered_map<int,int>&                                                             correctedGrid1BoundaryVertexIdx_to_grid2FaceIdx,
              const std::vector<int>&                                                                  grid1Face_fullyContainedIn_grid2Face)
{
    grid1FaceIdx_to_correctedGrid1FaceIdx.resize(grid1_numFaces);
    grid2VanishedFaceIdx_to_correctedGrid1FaceIndices.resize(grid2_numFaces);
        
    int invalidIdx = -1;
    //  grid1VertexIdx_to_grid2FaceIdx.resize(grid1_numVertices, invalidIdx);
    // correctedGrid1BoundaryVertexIdx_to_grid2FaceIdx,

    Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<2,3>>& correctedGrid1_faces =
        *(correctedGrid1GeomData.geometries.geomVector(std::integral_constant<int,1>()));
    Dune::cpgrid::EntityVariableBase<enum face_tag>& mutable_correctedGrid1_face_tags = correctedGrid1GeomData.face_tags;
    Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>& mutable_correctedGrid1_face_normals = correctedGrid1GeomData.face_normals;

    // Estimate new size of total number of faces in correctedGridData1
    int upperBoundFaceSize = grid1_numFaces; // total faces before correcting the data
    for (const auto& [grid2_faceIdx, overlapFacesInfo] : grid1_overlapFaces) {
        upperBoundFaceSize += overlapFacesInfo.size();
    }
    correctedGrid1_faces.resize(upperBoundFaceSize); // maybe larger than needed
    mutable_correctedGrid1_face_tags.resize(upperBoundFaceSize);
    mutable_correctedGrid1_face_normals.resize(upperBoundFaceSize);

    std::vector<std::vector<int>> aux_correctedGrid1_face_to_point{};
    aux_correctedGrid1_face_to_point.resize(upperBoundFaceSize);
    
    int correctedGrid1_num_points = 0;
    int correctedGrid1_face_count = 0;
    
    for (int i = 0; i < grid1.numFaces(); ++i) {
        const auto& newFacesInfo = grid1_vanishedFaceIdx_to_newFaces[i];
        const int grid2_containerFaceIdx = grid1Face_fullyContainedIn_grid2Face[i];
        
        if ( newFacesInfo.empty() || (grid2_containerFaceIdx != invalidIdx)) { // face is not "affected by gridData2", then store it "like it is" geometrically

            const auto face =  Dune::cpgrid::EntityRep<1>(i, true);

            grid1FaceIdx_to_correctedGrid1FaceIdx[i] = std::vector<int>{correctedGrid1_face_count};
            correctedGrid1FaceIdx_to_grid1FaceIdx[correctedGrid1_face_count] = i;

            correctedGrid1_faces[correctedGrid1_face_count] = (*grid1.getGeometry().geomVector(std::integral_constant<int,1>()))[face];
            mutable_correctedGrid1_face_tags[correctedGrid1_face_count] = grid1.faceTag(i);
            mutable_correctedGrid1_face_normals[correctedGrid1_face_count] = grid1.faceNormals(i);

            std::vector<int> faceToPoint{};
            faceToPoint.reserve(grid1.faceToPoint(i).size());
            for (const auto& vertexIdx : grid1.faceToPoint(i)) {
                faceToPoint.push_back(vertexIdx);
            }
            correctedGrid1_num_points += faceToPoint.size();
            aux_correctedGrid1_face_to_point[correctedGrid1_face_count] = faceToPoint;

            const auto& faceToCell = grid1.faceToCell(i);
            correctedGrid1GeomData.face_to_cell.appendRow(faceToCell.begin(), faceToCell.end());

            ++correctedGrid1_face_count;
        }
        else {
            std::vector<int> newFaceIndices{};

            for (const auto& [grid2_faceIdx, newFaceToCoord] : newFacesInfo) {

                const auto overlapIt = grid1_overlapFaces.find(grid2_faceIdx);
                assert(overlapIt != grid1_overlapFaces.end());
                const auto& overlapFacesInfo = overlapIt->second;

                newFaceIndices.push_back(correctedGrid1_face_count);
                grid1FaceIdx_to_correctedGrid1FaceIdx[i].push_back(correctedGrid1_face_count);
                correctedGrid1FaceIdx_to_grid1FaceIdx[correctedGrid1_face_count] = i;
                
                grid2VanishedFaceIdx_to_correctedGrid1FaceIndices[grid2_faceIdx].push_back(correctedGrid1_face_count);

                const auto [faceCenter, faceArea, faceNormal] = computeFaceCenterAreaNormal(newFaceToCoord);

                correctedGrid1_faces[correctedGrid1_face_count] = Dune::cpgrid::Geometry<2,3>(faceCenter, faceArea);
                mutable_correctedGrid1_face_tags[correctedGrid1_face_count] = grid1.faceTag(i); // shared tag
                mutable_correctedGrid1_face_normals[correctedGrid1_face_count] = faceNormal;

                std::vector<int> faceToPoint{};
                faceToPoint.reserve(newFaceToCoord.size());

                for (const auto& vertex : newFaceToCoord) {
                    int vertexIdx = -1; // invalid to be rewritten
                    // check if vertex (is new) and has already been stored 
                    auto it = newVertex_to_correctedGrid1VertexIdx.find(vertex);
                    if (it != newVertex_to_correctedGrid1VertexIdx.end()) {
                        vertexIdx = it->second;
                    }
                    else { // otherwise, vertex existed already in gridData1 (before "correction")
                        auto iit = grid1_vertex_to_vertexIdx.find(vertex);
                        assert( iit != grid1_vertex_to_vertexIdx.end());
                        vertexIdx = iit->second;

                        // correctedGrid1BoundaryVertexIdx_to_grid2BoundaryVertexIdx[vertexIdx] = grid2_faceIdx;
                    }
                    assert(vertexIdx>=0);
                    faceToPoint.push_back(vertexIdx);

                    correctedGrid1BoundaryVertexIdx_to_grid2FaceIdx[vertexIdx] = grid2_faceIdx;

                    // grid1VertexIdx_to_grid2FaceIdx[vertexIdx] = grid2_faceIdx;
                }
                correctedGrid1_num_points += faceToPoint.size();
                aux_correctedGrid1_face_to_point[correctedGrid1_face_count] = faceToPoint;

                const auto& face_to_cell = grid1.faceToCell(i);
                correctedGrid1GeomData.face_to_cell.appendRow(face_to_cell.begin(), face_to_cell.end());

                /*for (const auto& setCoord : overlapFacesInfo) {
                    if (setCoord == newFaceToCoord)
                        grid2VanishedFaceIdx_to_correctedGrid1FaceIndices[grid2_faceIdx].push_back(correctedGrid1_face_count);
                        }*/
                ++correctedGrid1_face_count;
            }
            grid1VanishedFaceIdx_to_correctedGrid1FaceIndices[i] = newFaceIndices;
        }
    }
    correctedGrid1GeomData.face_to_point.reserve(correctedGrid1_face_count, correctedGrid1_num_points);
    for (int face = 0; face < correctedGrid1_face_count; ++face) {
        correctedGrid1GeomData.face_to_point.appendRow(aux_correctedGrid1_face_to_point[face].begin(),
                                                       aux_correctedGrid1_face_to_point[face].end());
    }

}

void addCells(GeomData& geomData,
              const Dune::cpgrid::CpGridData& gridData)
{
    Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<3,3>>& cells = *(geomData.geometries.geomVector(std::integral_constant<int,0>()));
    cells.resize(gridData.size(0));
    geomData.cell_to_point.resize(gridData.size(0));

    for (int i = 0; i < gridData.size(0); ++i) {

        const auto elem = Dune::cpgrid::Entity<0>(gridData, i, true);
        geomData.cell_to_point[i] = gridData.cellToPoint(i);

        int* indices_storage_ptr = geomData.cell_to_point[i].data();
        cells[i] = Dune::cpgrid::Geometry<3,3>(elem.geometry().center(),
                                               elem.geometry().volume(),
                                               geomData.geometries.geomVector(std::integral_constant<int,3>()),
                                               indices_storage_ptr);
    }
    geomData.face_to_cell.makeInverseRelation(geomData.cell_to_face);
}

void makeCellRefinementParentFaceAware(bool                                                        hasOnlyOneFacePerType,
                                       const std::array<std::vector<int>, 6>&                      classifiedParentCellFaces,
                                       const Dune::cpgrid::CpGridData&                             parentGrid,
                                       const Dune::cpgrid::Entity<0>&                              parentCellElem,
                                       const Dune::cpgrid::CpGridData&                             cellRefGrid,
                                       Dune::cpgrid::CpGridData&                                   correctedCellRefGrid,
                                       GeomData&                                                   correctedCellRefGeomData,
                                       std::unordered_map<int,int>&                                correctedCellRefBoundaryVertexIdx_to_parentGridFaceIdx,
                                       std::vector<int>&                                           correctedCellRefBoundaryFaceIdx_to_parentGridFaceIdx,
                                       std::vector<std::array<int,2>>&                             parentBoundaryVertexIdx_to_correctedCellRefBoundaryVertexIdx,
                                       std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                       const std::array<int,3>&                                    cells_per_dim)
{
    int numFaces = cellRefGrid.numFaces(); // to be increased

    std::set<Dune::FieldVector<double,3>, FieldVectorLess>      foundNewVertices{};
    std::map<Dune::FieldVector<double,3>, int, FieldVectorLess> newVertex_to_correctedCellRefinementDataVertexIdx{};
    std::map<Dune::FieldVector<double,3>, int, FieldVectorLess> cellRefGridBoundaryVertex_to_cellRefGridVertexIdx{};
    std::map<Dune::FieldVector<double,3>, int, FieldVectorLess> parentGrid_vertex_to_vertexIdx{};
    

    std::map<int, std::vector<std::vector<Dune::FieldVector<double,3>>>> allOverlapFaces{};
    // allOverlapFaces[ parent cell face idx ] = { {overlapFace0 set of Coord}, ..., {overlapFaceN set of Coord}}
    std::vector<std::vector<std::pair<int,std::vector<Dune::FieldVector<double,3>>>>> vanishedRefFace_to_parentFaceAndOverlapCoords{};
    vanishedRefFace_to_parentFaceAndOverlapCoords.resize(numFaces);
    // vanishedRefFace_to_parentFacaAndOverlapCoords[ cellRefGrid face index ] = { { parent face index, {set of coordinates of overlap area} }_0, ...}

    int invalidIdx = -1;
    std::vector<int> fullyContainedInParentFace{}; // fullyContainedInParentFace [ refined face idx ] = parent face index
    fullyContainedInParentFace.resize(numFaces, -1);

    for (int i = 0; i < cellRefGrid.size(0); ++i) {

        const auto refinedElem = Dune::cpgrid::Entity<0>(cellRefGrid, i, true);
        if (!isAtGridBoundary(cellRefGrid, refinedElem))
            continue;

        const auto collectedVertices = collectNewVertices(parentGrid,
                                                          parentCellElem,
                                                          parentGrid_vertex_to_vertexIdx,
                                                          classifiedParentCellFaces,
                                                          cellRefGrid,
                                                          refinedElem,
                                                          cellRefGridBoundaryVertex_to_cellRefGridVertexIdx,
                                                          allOverlapFaces,
                                                          vanishedRefFace_to_parentFaceAndOverlapCoords,
                                                          fullyContainedInParentFace);

        foundNewVertices.insert(collectedVertices.begin(), collectedVertices.end());
    }

    int numVertices = cellRefGrid.size(3);
    const auto& parentCellToPoint = parentGrid.cellToPoint(parentCellElem.index());
    std::map<int,int> correctedCellRefDataVertexIdx_to_parentGridVertexIdx{};

    addVertices(numVertices,
                cellRefGrid,
                parentGrid_vertex_to_vertexIdx,
                foundNewVertices,
                correctedCellRefGeomData,
                correctedCellRefDataVertexIdx_to_parentGridVertexIdx,
                newVertex_to_correctedCellRefinementDataVertexIdx,
                parentCellToPoint, 
                parentBoundaryVertexIdx_to_correctedCellRefBoundaryVertexIdx);
    
    std::vector<std::vector<int>> cellRefGridFace_to_correctedCellRefGridFaces{};
    cellRefGridFace_to_correctedCellRefGridFaces.resize(numFaces);
    std::unordered_map<int, int> correctedCellRefGridFace_to_cellRefGridFace{};

    std::vector<std::vector<int>> parentFace_to_correctRefinedFaces{};
    parentFace_to_correctRefinedFaces.resize(parentGrid.numFaces());
    
    std::unordered_map<int, std::vector<int>> cellRefGridVanishedFaceIdx_to_correctedCellRefGridFaceIndices{};
    
    addFaces(numFaces,                                                       // grid1_numFaces,
             parentGrid.numFaces(),                                          // grid2_numFaces,
             allOverlapFaces,                                                // grid1_overlapFaces,
             cellRefGrid,                                                    // grid1,
             vanishedRefFace_to_parentFaceAndOverlapCoords,                  // grid1_vanishedFaceIdx_to_newFaces,
             cellRefGridFace_to_correctedCellRefGridFaces,                   // grid1FaceIdx_to_correctedGrid1FaceIdx, 
             correctedCellRefGridFace_to_cellRefGridFace,                    // correctedGrid1FaceIdx_to_grid1FaceIdx,
             correctedCellRefGeomData,                                       // correctedGrid1GeomData,
             parentFace_to_correctRefinedFaces,                              // grid2VanishedFaceIdx_to_correctedGrid1FaceIndices,
             newVertex_to_correctedCellRefinementDataVertexIdx,              // newVertex_to_correctedGrid1VertexIdx,
             cellRefGridBoundaryVertex_to_cellRefGridVertexIdx,              // grid1_vertex_to_vertexIdx,
             cellRefGridVanishedFaceIdx_to_correctedCellRefGridFaceIndices,  // grid1VanishedFaceIdx_to_correctedGrid1FaceIndices,
             correctedCellRefBoundaryVertexIdx_to_parentGridFaceIdx,         // correctedGrid1VertexIdx_to_grid2FaceIdx,
             fullyContainedInParentFace);                                    // grid1Face_fullyContainedIn_grid2Face);
    
    addCells(correctedCellRefGeomData,
             cellRefGrid);
    
    correctedCellRefBoundaryFaceIdx_to_parentGridFaceIdx.resize(correctedCellRefGrid.numFaces(), invalidIdx);

    provideRefinementParentFaceIdxRelations(hasOnlyOneFacePerType,
                                            classifiedParentCellFaces,
                                            parentGrid,
                                            parentCellElem,
                                            cells_per_dim,
                                            correctedCellRefBoundaryFaceIdx_to_parentGridFaceIdx,
                                            faceInMarkedElemAndRefinedFaces,
                                            cellRefGridFace_to_correctedCellRefGridFaces, 
                                            cellRefGridVanishedFaceIdx_to_correctedCellRefGridFaceIndices,
                                            parentFace_to_correctRefinedFaces,
                                            fullyContainedInParentFace);
    
}


void computeNewRefinedGeometriesOnSharedCoarseFace(int numFaces1,
                                                   const std::vector<int>& refinedFaces1,
                                                   const Dune::cpgrid::CpGridData& parentFaceAwareCellRefinement1,
                                                   std::set<Dune::FieldVector<double,3>, FieldVectorLess>& foundNewVertices1,
                                                   std::map<Dune::FieldVector<double,3>, int, FieldVectorLess>& faceExistingVertex1_to_cellRef1CornerIdx,
                                                   std::map<int,std::vector<std::vector<Dune::FieldVector<double,3>>>>& overlapFaces1,
                                                   std::vector<std::vector<std::pair<int,std::vector<Dune::FieldVector<double,3>>>>>& vanishedCellRef1Face_to_newRefinedFaces,
                                                   std::vector<int>& face2FullyContainedInNeighborFace1,
                                                   int numFaces2,
                                                   const std::vector<int>& refinedFaces2,
                                                   const Dune::cpgrid::CpGridData& parentFaceAwareCellRefinement2,
                                                   std::set<Dune::FieldVector<double,3>, FieldVectorLess>& foundNewVertices2,
                                                   std::map<Dune::FieldVector<double,3>, int, FieldVectorLess>& faceExistingVertex2_to_cellRef2CornerIdx,
                                                   std::map<int,std::vector<std::vector<Dune::FieldVector<double,3>>>>& overlapFaces2,
                                                   std::vector<std::vector<std::pair<int,std::vector<Dune::FieldVector<double,3>>>>>& vanishedCellRef2Face_to_newRefinedFaces,
                                                   std::vector<int>&  face1FullyContainedInNeighborFace2)
{
    vanishedCellRef1Face_to_newRefinedFaces.resize(numFaces1);
    vanishedCellRef2Face_to_newRefinedFaces.resize(numFaces2);

    int invalidIdx = -1;
    face2FullyContainedInNeighborFace1.resize(numFaces2, invalidIdx);   
    face1FullyContainedInNeighborFace2.resize(numFaces1, invalidIdx);

    for (const auto& refinedFace1 : refinedFaces1) {
        const auto face1 = Dune::cpgrid::EntityRep<1>(refinedFace1, true);
        const auto& cellToFace1 = parentFaceAwareCellRefinement1.faceToCell(refinedFace1);
        assert(cellToFace1.size() == 1); // it's at bounday cell-refinement
        const auto cell1 = cellToFace1[0];

        for (const auto& refinedFace2 : refinedFaces2) {
            const auto face2 = Dune::cpgrid::EntityRep<1>(refinedFace2, true);
            const auto& cellToFace2 = parentFaceAwareCellRefinement2.faceToCell(refinedFace2);
            assert(cellToFace2.size() == 1); // it's at bounday cell-refinement
            const auto cell2 = cellToFace2[0];

            bool face2FullyContainedInFace1 = false;
            const auto newFaceA = computeFaceOverlapVertices(face1,
                                                             parentFaceAwareCellRefinement1,
                                                             faceExistingVertex1_to_cellRef1CornerIdx,
                                                             face2,
                                                             parentFaceAwareCellRefinement2,
                                                             faceExistingVertex2_to_cellRef2CornerIdx,
                                                             foundNewVertices2,
                                                             face2FullyContainedInFace1);

            bool face1FullyContainedInFace2 = false;
            const auto newFaceB = computeFaceOverlapVertices(face2,
                                                             parentFaceAwareCellRefinement2,
                                                             faceExistingVertex2_to_cellRef2CornerIdx,
                                                             face1,
                                                             parentFaceAwareCellRefinement1,
                                                             faceExistingVertex1_to_cellRef1CornerIdx,
                                                             foundNewVertices1,
                                                             face1FullyContainedInFace2);

            assert(newFaceA.has_value() == newFaceB.has_value());

            if (face2FullyContainedInFace1){
                face2FullyContainedInNeighborFace1[refinedFace2] = refinedFace1;
            }

            if (face1FullyContainedInFace2){
                face1FullyContainedInNeighborFace2[refinedFace1] = refinedFace2;
            }

            if (newFaceA.has_value() && !newFaceA.value().empty()) {
                assert( newFaceA.value() == newFaceB.value() );
                // Save face vertices/coordinates as one of the overlapping refined faces with the parent cell
                auto it1  = overlapFaces1.find(refinedFace2);
                if (it1 != overlapFaces1.end()) {
                    auto& collectedFaces = *it1;
                    collectedFaces.second.push_back(newFaceA.value());
                }
                else {
                    auto& collectedFaces = overlapFaces1[refinedFace2];
                    collectedFaces.push_back(newFaceA.value());
                }

                auto it2  = overlapFaces2.find(refinedFace1);
                if (it2 != overlapFaces2.end()) {
                    auto& collectedFaces = *it2;
                    collectedFaces.second.push_back(newFaceA.value());
                }
                else {
                    auto& collectedFaces = overlapFaces2[refinedFace1];
                    collectedFaces.push_back(newFaceA.value());
                }

                if (!face2FullyContainedInFace1) {
                    vanishedCellRef2Face_to_newRefinedFaces[refinedFace2].push_back(std::make_pair(refinedFace1, newFaceA.value()));
                }
                if (!face1FullyContainedInFace2) {
                    vanishedCellRef1Face_to_newRefinedFaces[refinedFace1].push_back(std::make_pair(refinedFace2, newFaceA.value()));
                }
            }
        }
    }
}

void computeNewGeometries(const Dune::CpGrid& grid,
                          const Dune::cpgrid::CpGridData& parentFaceAwareCellRefinement1,
                          const Dune::cpgrid::CpGridData& parentFaceAwareCellRefinement2,
                          const Dune::cpgrid::Entity<0>& parentCell1,
                          const Dune::cpgrid::Entity<0>& parentCell2,
                          std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces)
{
    const auto& leafView = grid.leafGridView();
    bool parentsShareFace = false;

    int shared_face_count = 0;
    int shared_face_idx = -1;
    for (const auto& intersection : Dune::intersections(leafView, parentCell1)) {
        if (!intersection.neighbor())
            continue;

        const auto& neighborIdx = intersection.outside().index();
        if (neighborIdx == parentCell2.index()) {
            parentsShareFace = true;
            shared_face_idx = intersection.id();
            ++shared_face_count;
            break;
        }
    }

    if (!parentsShareFace) {
        std::cout<< " do nothing " << std::endl;
    }

    assert(shared_face_count == 1); // parent cells share exactly one face
    assert(shared_face_idx >= 0); // valid index

    assert(faceInMarkedElemAndRefinedFaces[shared_face_idx].size() == 2);

    const auto& [p1, refinedFaces1] = faceInMarkedElemAndRefinedFaces[shared_face_idx][0];
    const auto& [p2, refinedFaces2] = faceInMarkedElemAndRefinedFaces[shared_face_idx][1];

    std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>> neighborAwareCellRef1Data;
    std::shared_ptr<Dune::cpgrid::CpGridData> neighborAwareCellRef1_ptr = std::make_shared<Dune::cpgrid::CpGridData>(neighborAwareCellRef1Data); // ccobj_
    auto& neighborAwareCellRef1 = *neighborAwareCellRef1_ptr;
    GeomData neighborAwareGeomData1(neighborAwareCellRef1);

    std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>> neighborAwareCellRef2Data;
    std::shared_ptr<Dune::cpgrid::CpGridData> neighborAwareCellRef2_ptr = std::make_shared<Dune::cpgrid::CpGridData>(neighborAwareCellRef2Data); // ccobj_
    auto& neighborAwareCellRef2 = *neighborAwareCellRef2_ptr;
    GeomData neighborAwareGeomData2(neighborAwareCellRef2);

    int numFaces1 = parentFaceAwareCellRefinement1.numFaces();
    int numFaces2 = parentFaceAwareCellRefinement2.numFaces();

    std::set<Dune::FieldVector<double,3>, FieldVectorLess> foundNewVertices1{};
    std::map<Dune::FieldVector<double,3>, int, FieldVectorLess> faceExistingVertex1_to_cellRef1CornerIdx{};
    std::map<int,std::vector<std::vector<Dune::FieldVector<double,3>>>> overlapFaces1{};
    std::vector<std::vector<std::pair<int,std::vector<Dune::FieldVector<double,3>>>>> vanishedCellRef1Face_to_newRefinedFaces{};
    std::vector<int> face2FullyContainedInNeighborFace1{}; // face1 index in celRef1 that contains face2

    std::set<Dune::FieldVector<double,3>, FieldVectorLess> foundNewVertices2{};
    std::map<Dune::FieldVector<double,3>, int, FieldVectorLess> faceExistingVertex2_to_cellRef2CornerIdx{};
    std::map<int,std::vector<std::vector<Dune::FieldVector<double,3>>>> overlapFaces2{};
    std::vector<std::vector<std::pair<int,std::vector<Dune::FieldVector<double,3>>>>> vanishedCellRef2Face_to_newRefinedFaces{};
    std::vector<int> face1FullyContainedInNeighborFace2{}; // face2 index in celRef2 that contains face1

    computeNewRefinedGeometriesOnSharedCoarseFace(numFaces1,
                                                  refinedFaces1,
                                                  parentFaceAwareCellRefinement1,
                                                  foundNewVertices1,
                                                  faceExistingVertex1_to_cellRef1CornerIdx,
                                                  overlapFaces1,
                                                  vanishedCellRef1Face_to_newRefinedFaces,
                                                  face2FullyContainedInNeighborFace1,
                                                  numFaces2,
                                                  refinedFaces2,
                                                  parentFaceAwareCellRefinement2,
                                                  foundNewVertices2,
                                                  faceExistingVertex2_to_cellRef2CornerIdx,
                                                  overlapFaces2,
                                                  vanishedCellRef2Face_to_newRefinedFaces,
                                                  face1FullyContainedInNeighborFace2);

    int numVertices1 = parentFaceAwareCellRefinement1.size(3);
    std::map<Dune::FieldVector<double,3>, int, FieldVectorLess> newVertex_to_neighborAwareCellRef1CornerIdx{};
    std::map<int,int> newVertexIdx1_to_vertexIdx2{}; // neighborAwareCellRef1CornerIdx-> cellRef2CornerIdx

    addVertices(numVertices1,
                parentFaceAwareCellRefinement1, // assumed to be aware of parent-cell-faces!!!
                faceExistingVertex2_to_cellRef2CornerIdx,
                foundNewVertices1,
                neighborAwareGeomData1,// neighborAwareCellRef1_corners,
                newVertexIdx1_to_vertexIdx2,
                newVertex_to_neighborAwareCellRef1CornerIdx);

    int numVertices2 = parentFaceAwareCellRefinement2.size(3);
    std::map<Dune::FieldVector<double,3>, int, FieldVectorLess> newVertex_to_neighborAwareCellRef2CornerIdx{};
    std::map<int,int> newVertexIdx2_to_vertexIdx1{}; // neighborAwareCellRef2CornerIdx-> cellRef1CornerIdx

    addVertices(numVertices2,
                parentFaceAwareCellRefinement2, // assumed to be aware of parent-cell-faces!!!
                faceExistingVertex1_to_cellRef1CornerIdx,
                foundNewVertices2,
                neighborAwareGeomData2, //neighborAwareCellRef2_corners,
                newVertexIdx2_to_vertexIdx1,
                newVertex_to_neighborAwareCellRef2CornerIdx);

    std::unordered_map<int, std::vector<int>> vanishedFace1ToNewFaces{};
    std::unordered_map<int, int> newToOldFaceIdx1{};// the vanished faces are not included??
    std::vector<std::vector<int>> oldToNewFaceIdx1{};

    std::vector<std::vector<int>> vanishedFace2_to_correctFaces1{};
    //std::vector<int> corn1Idx_to_face2Idx{};

    std::unordered_map<int,int> corn1Idx_to_face2Idx; //correctedGrid1BoundaryVertexIdx_to_grid2FaceIdx{};

    addFaces(//numVertices1,
             numFaces1,
             numFaces2,
             overlapFaces1,
             parentFaceAwareCellRefinement1,
             vanishedCellRef1Face_to_newRefinedFaces,
             oldToNewFaceIdx1,
             newToOldFaceIdx1,
             neighborAwareGeomData1,
             vanishedFace2_to_correctFaces1,
             newVertex_to_neighborAwareCellRef1CornerIdx,
             faceExistingVertex1_to_cellRef1CornerIdx,
             vanishedFace1ToNewFaces,
             corn1Idx_to_face2Idx,
             face1FullyContainedInNeighborFace2);

    std::unordered_map<int, std::vector<int>> vanishedFace2ToNewFaces{};
    std::unordered_map<int, int> newToOldFaceIdx2{};// the vanished faces are not included??
    std::vector<std::vector<int>> oldToNewFaceIdx2{};

    std::vector<std::vector<int>> vanishedFace1_to_correctFaces2{};
    // std::vector<int> corn2Idx_to_face1Idx{};

    std::unordered_map<int,int> corn2Idx_to_face1Idx; //correctedGrid2BoundaryVertexIdx_to_grid1FaceIdx{};
    addFaces(//numVertices2,
             numFaces2,
             numFaces1,
             overlapFaces2,
             parentFaceAwareCellRefinement2,
             vanishedCellRef2Face_to_newRefinedFaces,
             oldToNewFaceIdx2,
             newToOldFaceIdx2,
             neighborAwareGeomData2,
             vanishedFace1_to_correctFaces2,
             newVertex_to_neighborAwareCellRef2CornerIdx,
             faceExistingVertex2_to_cellRef2CornerIdx,
             vanishedFace2ToNewFaces,
             corn2Idx_to_face1Idx,
             face2FullyContainedInNeighborFace1);

    addCells(neighborAwareGeomData1,
             parentFaceAwareCellRefinement1);

    addCells(neighborAwareGeomData2,
             parentFaceAwareCellRefinement2);
}


/*void makeCellRefinementNeighboringCellRefinementAware(const Dune::CpGrid&                                                grid,
                                                      const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>&      parentFaceAwareCellRefinements,
                                                      const Dune::cpgrid::Entity<0>&                                     parentCell,
                                                      std::vector<std::array<int,2>>&                                    extended_parent_to_refined_corners,
                                                      std::vector<std::vector<std::pair<int, std::vector<int>>>>&        faceInMarkedElemAndRefinedFaces,
                                                      const Dune::cpgrid::DefaultGeometryPolicy&                         beforeCorrection_geometries,
                                                      std::unordered_map<int,int>&                                       refinedCornIdx_to_parentFaceIdx,
                                                      std::vector<int>&                                                  correctRefinedFace_to_parentFace,
                                                      Dune::cpgrid::CpGridData&                                          corrected_cellRef_data,
                                                      Dune::cpgrid::DefaultGeometryPolicy&                               corrected_cellRef_geometries,
                                                      std::vector<std::array<int,8>>&                                    corrected_cellRef_cell_to_point,
                                                      Dune::cpgrid::OrientedEntityTable<0,1>&                            corrected_cellRef_cell_to_face,
                                                      Opm::SparseTable<int>&                                             corrected_cellRef_face_to_point,
                                                      Dune::cpgrid::OrientedEntityTable<1,0>&                            corrected_cellRef_face_to_cell,
                                                      Dune::cpgrid::EntityVariable<enum face_tag,1>&                     corrected_cellRef_face_tags,
                                                      Dune::cpgrid::SignedEntityVariable<Dune::FieldVector<double,3>,1>& corrected_cellRef_face_normals)
{
}*/



} // namespace Lgr
} // namespace Opm
