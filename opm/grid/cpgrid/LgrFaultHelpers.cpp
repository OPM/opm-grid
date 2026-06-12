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

void assignCorrectChildFace(bool                                             noCorrectionNeeded,
                            int                                              faceIdx,
                            int                                              faceType_count,
                            int                                              child_face,
                            std::vector<int>&                                children_faces,
                            std::vector<std::vector<int>>&                   refinedFace_to_parentFaces,
                            const std::vector<std::vector<int>>&             oldToNewFaceIdx,
                            const std::unordered_map<int, std::vector<int>>& vanishedFaceToNewFaces,
                            const std::vector<std::vector<int>>&             parentFace_to_correctRefinedFaces,
                            const std::vector<int>&                          fullyContainedInParentFace)
{
    auto addChildFace = [&](int idx) {
        children_faces.push_back(idx);
        refinedFace_to_parentFaces[idx].push_back(faceIdx);
    };

    if (noCorrectionNeeded) {
        addChildFace(child_face);
        return;
    }

    const auto& mappedFaces = oldToNewFaceIdx[child_face];

    if (faceType_count <= 1) {
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
        for (int idx : mappedFaces)
            addChildFace(idx);
    }
}


void provideRefinementParentFaceIdxRelations(bool                                                        noCorrectionNeeded,
                                             const std::array<std::vector<int>, 6>&                      classifiedFaces,
                                             const Dune::cpgrid::CpGridData&                             parentGridData,
                                             const Dune::cpgrid::Entity<0>&                              parentElem,
                                             const std::array<int,3>&                                    cells_per_dim,
                                             // std::vector<std::vector<int>>&                              parentFace_to_refinedFaces,
                                             std::vector<std::vector<int>>&                              refinedFace_to_parentFaces,
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
        const auto faceType_count = classifiedFaces[faceGroupIndex(parent_face_tag, face.orientation())].size();

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

                    if (noCorrectionNeeded) {
                        children_faces.push_back(child_face);
                        refinedFace_to_parentFaces[child_face] = std::vector{face.index()};
                    }
                    else {
                        for (const auto& newIdx : oldToNewFaceIdx[child_face]) {
                            children_faces.push_back(newIdx);
                            refinedFace_to_parentFaces[newIdx].push_back(face.index());
                            // parentFace_to_refinedFaces[face.index()].push_back(newIdx);
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

                    assignCorrectChildFace(noCorrectionNeeded,
                                           face.index(),
                                           faceType_count,
                                           child_face,
                                           children_faces,
                                           refinedFace_to_parentFaces,
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

                    assignCorrectChildFace(noCorrectionNeeded,
                                           face.index(),
                                           faceType_count,
                                           child_face,
                                           children_faces,
                                           refinedFace_to_parentFaces,
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

/*std::map<Dune::FieldVector<double,3>, int, FieldVectorLess> getCoordToPoint(const Dune::cpgrid::CpGridData& gridData,
                                                                            const Dune::cpgrid::Entity<0>& element)
{
    std::map<Dune::FieldVector<double,3>, int, FieldVectorLess> coordToVertexIdx{};
    
    for (const auto& pointIdx : gridData.cellToPoint(element.index())) {
        const auto point = Dune::cpgrid::Entity<3>(gridData, pointIdx, true).geometry().center();
        coordToVertexIdx[point] = pointIdx;
        }
        return coordToVertexIdx;
        }
*/

void makeCellRefinementParentFaceAware(bool                                                               noCorrectionNeeded,
                                       const std::array<std::vector<int>, 6>&                             classifiedParentFaces,
                                       const Dune::cpgrid::CpGridData&                                    singleCellRefinementData,
                                       const Dune::cpgrid::CpGridData&                                    parentGridData,
                                       const Dune::cpgrid::Entity<0>&                                     parentElem,
                                       std::vector<std::array<int,2>>&                                    extended_parent_to_refined_corners,
                                       Dune::cpgrid::DefaultGeometryPolicy&                               refined_geometries,
                                       std::vector<std::vector<std::pair<int, std::vector<int>>>>&        faceInMarkedElemAndRefinedFaces,
                                       std::unordered_map<int,int>&                                       refinedCornIdx_to_parentFaceIdx,
                                       std::vector<std::vector<int>>&                                     refinedFace_to_parentFaces,
                                       Dune::cpgrid::CpGridData&                                          correctedRefinementData,
                                       Dune::cpgrid::DefaultGeometryPolicy&                               corrected_refined_geometries,
                                       std::vector<std::array<int,8>>&                                    corrected_refined_cell_to_point,
                                       Dune::cpgrid::OrientedEntityTable<0,1>&                            corrected_refined_cell_to_face,
                                       Opm::SparseTable<int>&                                             corrected_refined_face_to_point,
                                       Dune::cpgrid::OrientedEntityTable<1,0>&                            corrected_refined_face_to_cell,
                                       Dune::cpgrid::EntityVariable<enum face_tag,1>&                     corrected_refined_face_tags,
                                       Dune::cpgrid::SignedEntityVariable<Dune::FieldVector<double,3>,1>& corrected_refined_face_normals,
                                       const std::array<int,3>& cells_per_dim)
{
    int numVertices = singleCellRefinementData.size(3);
    int numFaces = singleCellRefinementData.numFaces(); // to be increased

    std::set<Dune::FieldVector<double,3>, FieldVectorLess> missingVertices{};
    std::map<Dune::FieldVector<double,3>, int, FieldVectorLess> newVertexToIdx{};
    std::map<Dune::FieldVector<double,3>, int, FieldVectorLess> refinedFaceExistingVertex_to_refinedGridCornerIdx{};
    std::map<Dune::FieldVector<double,3>, int, FieldVectorLess> parentFaceExistingVertex_to_parentGridCornerIdx{};
    
    Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<0,3>>& refined_corners =
                    *(refined_geometries.geomVector(std::integral_constant<int,3>()));

    Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<0,3>>& corrected_corners =
        *(corrected_refined_geometries.geomVector(std::integral_constant<int,3>()));
    Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<2,3>>& corrected_faces =
        *(corrected_refined_geometries.geomVector(std::integral_constant<int,1>()));
    Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<3,3>>& corrected_cells =
        *(corrected_refined_geometries.geomVector(std::integral_constant<int,0>()));
    Dune::cpgrid::EntityVariableBase<enum face_tag>& mutable_corrected_face_tags = corrected_refined_face_tags;
    Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>& mutable_corrected_face_normals = corrected_refined_face_normals;
   

    // overlapFaces[parent cell face idx] = { refined face index, {Coord0, ..., Coord3} }
    // std::map<int,std::vector<std::pair<int,std::vector<Coordinate>>>> overlapFaces{};

    std::map<int, std::vector<std::vector<Dune::FieldVector<double,3>>>> allOverlapFaces{};
    // allOverlapFaces[ parent cell face ] = { {overlapFaceIdx0, {its set of Coord}}, ..., {overlapFaceN, {its set of Coord}}}
    std::vector<std::vector<std::pair<int,std::vector<Dune::FieldVector<double,3>>>>> vanishedRefinedFaceToNewRefinedFaces{};
    vanishedRefinedFaceToNewRefinedFaces.resize(numFaces);

    // std::vector<std::vector<int>> vanishedFace_to_parentFaces{};
    //vanishedFace_to_parentFaces.resize(singleCellRefinementData.numFaces());
    
    std::vector<int> fullyContainedInParentFace{};
    fullyContainedInParentFace.resize(singleCellRefinementData.numFaces());

    for (int i = 0; i < singleCellRefinementData.size(0); ++i) {

        const auto refinedElem = Dune::cpgrid::Entity<0>(singleCellRefinementData, i, true);
        const auto collectedVertices = collectNewVertices(parentGridData,
                                                          parentElem,
                                                          parentFaceExistingVertex_to_parentGridCornerIdx, 
                                                          classifiedParentFaces,
                                                          singleCellRefinementData,
                                                          refinedElem,
                                                          refinedFaceExistingVertex_to_refinedGridCornerIdx,
                                                          //  vanishedFace_to_parentFaces,
                                                          allOverlapFaces,
                                                          vanishedRefinedFaceToNewRefinedFaces,
                                                          fullyContainedInParentFace);
        
        missingVertices.insert(collectedVertices.begin(), collectedVertices.end());
    }
    corrected_corners.resize(singleCellRefinementData.size(3) + missingVertices.size());
    
    // Copy all vertices from single cell refinement (the missing vertices arising from
    // edge-intersections with parent cell faces when the parent cell has more than one
    // face of the same face tag and face oriantation are not including here, and will be
    // added later).
    for (int i = 0; i < singleCellRefinementData.size(3); ++i) {
        corrected_corners[i] = refined_corners.get(i);
    }

    const auto& parentCellToPoint = parentGridData.cellToPoint(parentElem.index());
    
    // Assign vertex index in the single-cell-refinement to new vertices
    // that come from intersection of edges with parent cell faces.
    for (const auto& vertex : missingVertices) {
        auto it = parentFaceExistingVertex_to_parentGridCornerIdx.find(vertex); 
        if (it != parentFaceExistingVertex_to_parentGridCornerIdx.end()) {
            auto coarseIdx = it->second;
            if (std::find(parentCellToPoint.begin(),
                          parentCellToPoint.end(), coarseIdx) == parentCellToPoint.end()) {
                // extended_parent_to_refined_corners.push_back(entry);
                extended_parent_to_refined_corners.push_back(std::array<int,2>{coarseIdx, numVertices});
                newVertexToIdx[vertex] = numVertices;
                refinedFaceExistingVertex_to_refinedGridCornerIdx[vertex] = numVertices;
                
                corrected_corners[numVertices] = Dune::cpgrid::Geometry<0, 3>(vertex);
                ++numVertices;
            }
        }
        else {
        
            newVertexToIdx[vertex] = numVertices;
            refinedFaceExistingVertex_to_refinedGridCornerIdx[vertex] = numVertices;
            
            corrected_corners[numVertices] = Dune::cpgrid::Geometry<0, 3>(vertex);
            ++numVertices;
        }
    }

    // Estimate new size of total number of refined face
    int updateFaceSize = singleCellRefinementData.numFaces(); // total faces before correcting the data
    for (const auto& [parentFaceIdx, overlapFacesInfo] : allOverlapFaces) {
        updateFaceSize += overlapFacesInfo.size();
    }
    corrected_faces.resize(updateFaceSize); // maybe larger than needed
    mutable_corrected_face_tags.resize(updateFaceSize);
    mutable_corrected_face_normals.resize(updateFaceSize);


    std::vector<std::vector<int>> aux_face_to_point{};
    aux_face_to_point.resize(updateFaceSize);

    std::vector<std::vector<int>> parentFace_to_refinedFaces{};
    parentFace_to_refinedFaces.resize(parentGridData.numFaces());

    std::unordered_map<int, std::vector<int>> vanishedFaceToNewFaces{};
    
    std::unordered_map<int, int> newToOldFaceIdx{};// the vanished faces are not included??
    std::vector<std::vector<int>> oldToNewFaceIdx{};
    oldToNewFaceIdx.resize(singleCellRefinementData.numFaces());


    std::vector<std::vector<int>> parentFace_to_correctRefinedFaces{};
    parentFace_to_correctRefinedFaces.resize(parentGridData.numFaces());
    
    int corrected_num_points = 0;
    int corrected_face_count = 0;
    for (int i = 0; i < singleCellRefinementData.numFaces(); ++i) {
       
        if ( vanishedRefinedFaceToNewRefinedFaces[i].empty()) {
        
            const auto face =  Dune::cpgrid::EntityRep<1>(i, true);

            oldToNewFaceIdx[i] = std::vector<int>{corrected_face_count};
            newToOldFaceIdx[corrected_face_count] = i;
            
            corrected_faces[corrected_face_count] = (*singleCellRefinementData.getGeometry().geomVector(std::integral_constant<int,1>()))[face];
            mutable_corrected_face_tags[corrected_face_count] = singleCellRefinementData.faceTag(i);
            mutable_corrected_face_normals[corrected_face_count] = singleCellRefinementData.faceNormals(i);

            std::vector<int> faceToPoint{};
            faceToPoint.reserve(singleCellRefinementData.faceToPoint(i).size());
            for (const auto& vIdx : singleCellRefinementData.faceToPoint(i)) {
                faceToPoint.push_back(vIdx);
            }
            corrected_num_points += faceToPoint.size();
            aux_face_to_point[corrected_face_count] = faceToPoint;

            const auto& face_to_cell = singleCellRefinementData.faceToCell(i);
            corrected_refined_face_to_cell.appendRow(face_to_cell.begin(), face_to_cell.end());
            
            ++corrected_face_count;
        }
        else {
            const auto& newFacesInfo = vanishedRefinedFaceToNewRefinedFaces[i];
            assert(newFacesInfo.size()>1);
            
            std::vector<int> newFaces{};
            
            for (const auto& [parentFaceIdx, newRefinedFaceToCoord] : newFacesInfo) {

                newFaces.push_back(corrected_face_count);
                oldToNewFaceIdx[i].push_back(corrected_face_count);
                newToOldFaceIdx[corrected_face_count] = i;
                parentFace_to_correctRefinedFaces[parentFaceIdx].push_back(corrected_face_count);

                const auto [faceCenter, faceArea, faceNormal] = computeFaceCenterAreaNormal(newRefinedFaceToCoord);
                
                corrected_faces[corrected_face_count] = Dune::cpgrid::Geometry<2,3>(faceCenter, faceArea);
                mutable_corrected_face_tags[corrected_face_count] = singleCellRefinementData.faceTag(i); // shared tag
                mutable_corrected_face_normals[corrected_face_count] = faceNormal; 
                
                std::vector<int> faceToPoint{};
                faceToPoint.reserve(4);
                
                for (const auto& v : newRefinedFaceToCoord) {
                    faceToPoint.push_back(refinedFaceExistingVertex_to_refinedGridCornerIdx[v]);
                    refinedCornIdx_to_parentFaceIdx[refinedFaceExistingVertex_to_refinedGridCornerIdx[v]] = parentFaceIdx;
                }

                // Add the amount of points to the count num_points.
                corrected_num_points += faceToPoint.size();
                aux_face_to_point[corrected_face_count] = faceToPoint;
                
                const auto& face_to_cell = singleCellRefinementData.faceToCell(i);
                corrected_refined_face_to_cell.appendRow(face_to_cell.begin(), face_to_cell.end());

                const auto& overlapFacesInfo = allOverlapFaces.at(parentFaceIdx);

                for (const auto& setCoord : overlapFacesInfo) {
                    if (setCoord == newRefinedFaceToCoord)
                        parentFace_to_refinedFaces[parentFaceIdx].push_back(corrected_face_count);
                }
                    ++numFaces;
                    ++corrected_face_count;
            }
            vanishedFaceToNewFaces[i] = newFaces;
        }
    }
    corrected_refined_face_to_point.reserve(corrected_face_count, corrected_num_points);
    for (int face = 0; face < corrected_face_count; ++face) {
        corrected_refined_face_to_point.appendRow(aux_face_to_point[face].begin(), aux_face_to_point[face].end());
    }
    
    corrected_cells.resize(singleCellRefinementData.size(0));
    corrected_refined_cell_to_point.resize(singleCellRefinementData.size(0));

    for (int i = 0; i < singleCellRefinementData.size(0); ++i) {
        
        const auto refinedElem = Dune::cpgrid::Entity<0>(singleCellRefinementData, i, true);
        corrected_refined_cell_to_point[i] = singleCellRefinementData.cellToPoint(i);

        int* indices_storage_ptr = corrected_refined_cell_to_point[i].data();
        corrected_cells[i] = Dune::cpgrid::Geometry<3,3>(refinedElem.geometry().center(),
                                                         refinedElem.geometry().volume(),
                                                         corrected_refined_geometries.geomVector(std::integral_constant<int,3>()),
                                                         indices_storage_ptr);
    }
    corrected_refined_face_to_cell.makeInverseRelation(corrected_refined_cell_to_face);

    // refinedFace_to_parentFaces[ refinedFaceIdx ] = empty vector for refined faces
    // interior to the single-cell-refinement; { parentFaceIdx0, parentFaceIdx1, ...}
    // otherwise, where the parent cell faces overlap the refined face. 
    refinedFace_to_parentFaces.resize(correctedRefinementData.numFaces());

    provideRefinementParentFaceIdxRelations(noCorrectionNeeded,
                                            classifiedParentFaces,
                                            parentGridData,
                                            parentElem,
                                            cells_per_dim,
                                            // parentFace_to_refinedFaces,
                                            refinedFace_to_parentFaces,
                                            faceInMarkedElemAndRefinedFaces,
                                            oldToNewFaceIdx,
                                            vanishedFaceToNewFaces,
                                            parentFace_to_correctRefinedFaces,
                                            fullyContainedInParentFace);
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

    for (const auto& intersection : Dune::intersections(leafView, parentCell1)) {
        if (!intersection.neighbor())
            continue;
            
        const auto& neighborIdx = intersection.outside().index();
        if (neighborIdx == parentCell2.index()) {
            parentsShareFace = true;
            break;
        }
    }

    if (!parentsShareFace) {
        std::cout<< " do nothing " << std::endl;
    }


    std::vector<int> cellRef1CornerIdx_to_cellRef2CornerIdx{};
    cellRef1CornerIdx_to_cellRef2CornerIdx.resize(parentFaceAwareCellRefinement1.size(3), -1/* invalid index */);
    
    std::vector<int>  cellRef2CornerIdx_to_cellRef1CornerIdx{};
    cellRef2CornerIdx_to_cellRef1CornerIdx.resize(parentFaceAwareCellRefinement2.size(3), -1 /*invalid index */);

    std::vector<int> cellRef1FaceIdx_to_cellRef2FaceIdx{};
    cellRef1FaceIdx_to_cellRef2FaceIdx.resize(parentFaceAwareCellRefinement1.numFaces(), -1/* invalid index */);
    
    std::vector<int> cellRef2FaceIdx_to_cellRef1FaceIdx{};
    cellRef2FaceIdx_to_cellRef1FaceIdx.resize(parentFaceAwareCellRefinement2.numFaces(), -1 /*invalid index */);


    std::set<Dune::FieldVector<double,3>, FieldVectorLess> foundNewVertices{};
    std::map<Dune::FieldVector<double,3>, int, FieldVectorLess> face1ExistingVertex_to_cellRef1CornerIdx{}; 
    std::map<Dune::FieldVector<double,3>, int, FieldVectorLess> face2ExistingVertex_to_cellRef2CornerIdx{};

    for (std::size_t i = 0; i < faceInMarkedElemAndRefinedFaces.size(); ++i) {
        if (faceInMarkedElemAndRefinedFaces[i].size() == 1)
            continue;
        
        assert(faceInMarkedElemAndRefinedFaces[i].size() == 2);

        const auto& [p1, refinedFaces1] = faceInMarkedElemAndRefinedFaces[i][0];
        const auto& [p2, refinedFaces2] = faceInMarkedElemAndRefinedFaces[i][1];
        
      
        for (const auto& refinedFace1 : refinedFaces1) {  
            const auto face1 = Dune::cpgrid::EntityRep<1>(refinedFace1, true);
            
            for (const auto& refinedFace2 : refinedFaces2) {
                const auto face2 = Dune::cpgrid::EntityRep<1>(refinedFace2, true);
                 
                bool face2FullyContainedInFace1 = false;
                const auto newFace = Opm::Lgr::computeFaceOverlapVertices(face1,
                                                                                 parentFaceAwareCellRefinement1,
                                                                                 face1ExistingVertex_to_cellRef1CornerIdx,
                                                                                 face2,
                                                                                 parentFaceAwareCellRefinement2,
                                                                                 face2ExistingVertex_to_cellRef2CornerIdx,
                                                                                 foundNewVertices,
                                                                                 face2FullyContainedInFace1);

                if (face2FullyContainedInFace1){
                    std::cout<< "face1 " << refinedFace1 << " contains in face2: " << refinedFace2 << std::endl;
                    std::cout<< " save it" <<std::endl;
                }    
            
                if (newFace.has_value() && !newFace.value().empty()) {
                    std::cout<< "Begin new face in face1 and face2! "<< refinedFace1 << std::endl;
                    for (const auto& coord : newFace.value()) {
                        std::cout<< coord[0] << " " << coord[1] << " " << coord[2] <<std::endl;
                    }
                    std::cout << std::endl;
                }
                // Create a method to track/relate refined corners at the boundary of neighboring
                // single cell refinements
            }
        }
    }   
}


void makeCellRefinementNeighboringCellRefinementAware(const Dune::CpGrid& grid,
                                                      const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>&      parentFaceAwareCellRefinements,
                                                      const Dune::cpgrid::Entity<0>&                                     parentCell,
                                                      std::vector<std::array<int,2>>&                                    extended_parent_to_refined_corners,
                                                      std::vector<std::vector<std::pair<int, std::vector<int>>>>&        faceInMarkedElemAndRefinedFaces,
                                                      const Dune::cpgrid::DefaultGeometryPolicy&                         beforeCorrection_geometries,
                                                      std::unordered_map<int,int>&                                       refinedCornIdx_to_parentFaceIdx,
                                                      std::vector<std::vector<int>>&                                     refinedFace_to_parentFaces,
                                                      Dune::cpgrid::CpGridData&                                          corrected_cellRef_data,
                                                      Dune::cpgrid::DefaultGeometryPolicy&                               corrected_cellRef_geometries,
                                                      std::vector<std::array<int,8>>&                                    corrected_cellRef_cell_to_point,
                                                      Dune::cpgrid::OrientedEntityTable<0,1>&                            corrected_cellRef_cell_to_face,
                                                      Opm::SparseTable<int>&                                             corrected_cellRef_face_to_point,
                                                      Dune::cpgrid::OrientedEntityTable<1,0>&                            corrected_cellRef_face_to_cell,
                                                      Dune::cpgrid::EntityVariable<enum face_tag,1>&                     corrected_cellRef_face_tags,
                                                      Dune::cpgrid::SignedEntityVariable<Dune::FieldVector<double,3>,1>& corrected_cellRef_face_normals)
{
    /*// To store the corrected single-cell-refinement data
      std::vector<std::shared_ptr<CpGridData>>                     corrected_data;
      std::shared_ptr<CpGridData>                                  corrected_cellRef_ptr = std::make_shared<CpGridData>(corrected_refined_data); // ccobj_
        auto&                                                        corrected_cellRef_grid = *corrected_refined_grid_ptr;
        DefaultGeometryPolicy&                                       corrected_cellRef_geometries = corrected_refined_grid.geometry_;
        std::vector<std::array<int,8>>&                              corrected_cellRef_cell_to_point = corrected_refined_grid.cell_to_point_;
        cpgrid::OrientedEntityTable<0,1>&                            corrected_cellRef_cell_to_face = corrected_refined_grid.cell_to_face_;
        Opm::SparseTable<int>&                                       corrected_cellRef_face_to_point = corrected_refined_grid.face_to_point_;
        cpgrid::OrientedEntityTable<1,0>&                            corrected_cellRef_face_to_cell = corrected_refined_grid.face_to_cell_;
        cpgrid::EntityVariable<enum face_tag,1>&                     corrected_cellRef_face_tags = corrected_refined_grid.face_tag_;
        cpgrid::SignedEntityVariable<Dune::FieldVector<double,3>,1>& corrected_cellRef_face_normals = corrected_refined_grid.face_normals_;

    std::vector<bool> coarseFaceHasBeenCorrected{};
    coarseFaceHasBeenCorrected.resize(grid.numFaces(), false); // define this where this method is called

    assert(grid.numFace() == faceInMarkedElemAndRefinedFaces.size());
    
    for (int i = 0; i < grid.numFace(); ++i) {
        if ((grid.currentLeafData().faceToCellSize(face) == 1) || faceInMarkedElemAndRefinedFaces[i].empty())
            coarseFaceHasBeenCorrected[i] = true; // no correction needed at grid boundary, or at face not involved in refinement     
    }

    const auto& beforeCorrectionCellRefinement = parentFaceAwareCellsRefinements[parentCell.index()];
    assert(beforeCorrectionCellRefinement != nullptr);


    int numVertices = beforeCorrectionCellRefinement->size(3);
    int numFaces = beforeCorrectionCellRefinement->numFaces(); 

    std::set<Dune::FieldVector<double,3>, FieldVectorLess> foundNewVertices{};
    std::map<Dune::FieldVector<double,3>, int, FieldVectorLess> newVertex_to_correctedCellRefCornerIdx{};
    std::map<Dune::FieldVector<double,3>, int, FieldVectorLess> refinedFaceExistingVertex_to_refinedGridCornerIdx{};
    std::map<Dune::FieldVector<double,3>, int, FieldVectorLess> neihborFaceExistingVertex_to_neighborGridCornerIdx{};
    
    Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<0,3>>& beforeCorrectionCorners =
                    *(beforeCorrection_geometries.geomVector(std::integral_constant<int,3>()));

    Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<0,3>>& corrected_corners =
        *(corrected_cellRef_geometries.geomVector(std::integral_constant<int,3>()));
    Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<2,3>>& corrected_faces =
        *(corrected_cellRef_geometries.geomVector(std::integral_constant<int,1>()));
    Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<3,3>>& corrected_cells =
        *(corrected_cellRef_geometries.geomVector(std::integral_constant<int,0>()));
    Dune::cpgrid::EntityVariableBase<enum face_tag>& mutable_corrected_face_tags = corrected_cellRef_face_tags;
    Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>& mutable_corrected_face_normals = corrected_cellRef_face_normals;
    
    std::map<int, std::vector<std::vector<Dune::FieldVector<double,3>>>> allOverlapFaces{};
    // allOverlapFaces[ cell face ] = { {overlapFaceIdx0, {its set of Coord}}, ..., {overlapFaceN, {its set of Coord}}}
    std::vector<std::vector<std::pair<int,std::vector<Dune::FieldVector<double,3>>>>> vanishedRefinedFaceToNewRefinedFaces{};
    vanishedRefinedFaceToNewRefinedFaces.resize(numFaces);
    
    std::vector<int> fullyContainedInNeighborFace{};
    fullyContainedInParentFace.resize(beforeCorrectionCellRefinement->numFaces());


     const auto& leafView = grid.leafGridView();
   
        for (const auto& intersection : Dune::intersections(leafView, parentCell)) {
            if (!intersection.neighbor())
                continue;
            
            const auto& neighbor = intersection.outside();
            const auto& neighborRefinement = parentFaceAwareCellRefinements[neighbor.index()];
            if ( neighborRefinement == nullptr)
                continue; // neighbor not involved in refinement

            // Correct

            for (int i = 0; i < beforeCorrectionCellRefinement->size(0); ++i) {

                const auto refinedElem = Dune::cpgrid::Entity<0>(*beforeCorrectionCellRefinement, i, true);
                const auto collectedVertices = collectNewVertices(parentGridData,
                                                                  parentCell,
                                                                  parentFaceExistingVertex_to_parentGridCornerIdx, 
                                                                  classifiedParentFaces,
                                                                  singleCellRefinementData,
                                                                  refinedElem,
                                                                  refinedFaceExistingVertex_to_refinedGridCornerIdx,
                                                                  allOverlapFaces,
                                                                  vanishedRefinedFaceToNewRefinedFaces,
                                                                  fullyContainedInParentFace);
        
                missingVertices.insert(collectedVertices.begin(), collectedVertices.end());
            }
            corrected_corners.resize(singleCellRefinementData.size(3) + missingVertices.size());



            
            coarseFaceHasBeenCorrected[intersection.id()] = true; // check!!!
        }
        

   
    
    // Copy all vertices from single cell refinement (the missing vertices arising from
    // edge-intersections with parent cell faces when the parent cell has more than one
    // face of the same face tag and face oriantation are not including here, and will be
    // added later).
    for (int i = 0; i < singleCellRefinementData.size(3); ++i) {
        corrected_corners[i] = refined_corners.get(i);
    }

    const auto& parentCellToPoint = parentGridData.cellToPoint(parentElem.index());
    
    // Assign vertex index in the single-cell-refinement to new vertices
    // that come from intersection of edges with parent cell faces.
    for (const auto& vertex : missingVertices) {
        auto it = parentFaceExistingVertex_to_parentGridCornerIdx.find(vertex); 
        if (it != parentFaceExistingVertex_to_parentGridCornerIdx.end()) {
            auto coarseIdx = it->second;
            if (std::find(parentCellToPoint.begin(),
                          parentCellToPoint.end(), coarseIdx) == parentCellToPoint.end()) {
                // extended_parent_to_refined_corners.push_back(entry);
                extended_parent_to_refined_corners.push_back(std::array<int,2>{coarseIdx, numVertices});
                newVertexToIdx[vertex] = numVertices;
                refinedFaceExistingVertex_to_refinedGridCornerIdx[vertex] = numVertices;
                
                corrected_corners[numVertices] = Dune::cpgrid::Geometry<0, 3>(vertex);
                ++numVertices;
            }
        }
        else {
        
            newVertexToIdx[vertex] = numVertices;
            refinedFaceExistingVertex_to_refinedGridCornerIdx[vertex] = numVertices;
            
            corrected_corners[numVertices] = Dune::cpgrid::Geometry<0, 3>(vertex);
            ++numVertices;
        }
    } */
    
}



} // namespace Lgr
} // namespace Opm
