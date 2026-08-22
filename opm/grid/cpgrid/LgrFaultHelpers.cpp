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

void addFaceVerticesCoordinatesToIndex(const Dune::cpgrid::CpGridData& grid,
                                       int faceIndex,
                                       std::map<Dune::FieldVector<double,3>, int, FieldVectorLess>& pointMap)
{
    for (const auto& pointIndex : grid.faceToPoint(faceIndex)) {
        const auto pointPosition = Dune::cpgrid::Entity<3>(grid, pointIndex, true).geometry().center();
        pointMap[pointPosition] = pointIndex;
    }
}

std::map<Dune::FieldVector<double,3>, int, FieldVectorLess>
buildExtendedCellPointVertexMap(const Dune::cpgrid::CpGridData& grid,
                                int cellIndex)
{
    std::map<Dune::FieldVector<double,3>, int, FieldVectorLess> pointMap;

    for (const auto& face : grid.cellToFace(cellIndex)) {
        addFaceVerticesCoordinatesToIndex(grid, face.index(), pointMap);
    }

    return pointMap;
}

void collectNewVertices(const Dune::cpgrid::CpGridData&                              parentGridData,
                        const Dune::cpgrid::Entity<0>&                               parentElem,
                        const std::array<std::vector<int>,6>&                        classifiedParentCellFaces,
                        const Dune::cpgrid::CpGridData&                              cellRefinementData, 
                        const Dune::cpgrid::Entity<0>&                               refinedElem,
                        BoundaryFaceInfo&                                            boundaryFaceInfo)
{
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

            addFaceVerticesCoordinatesToIndex(cellRefinementData, refinedFaceIdx,
                                              boundaryFaceInfo.boundaryVertexCoordinates_to_vertexIdx);

            bool refinedFaceFullyContainedInParentFace = false;
            const auto newFace = computeFaceOverlapVertices(parentFace,                             // face1
                                                            parentGridData,                         // gridData1
                                                            refinedFace,                            // face2
                                                            cellRefinementData,                     // gridData2
                                                            boundaryFaceInfo.foundNewVertices,
                                                            refinedFaceFullyContainedInParentFace); // face2FullyContainedInFace1

            if (refinedFaceFullyContainedInParentFace){
                boundaryFaceInfo.face_fullyContainedIn_otherGridFaceIdx[refinedFaceIdx] = parentFaceIdx;
            }
            
            if (newFace.has_value() && !newFace.value().empty()) {
                // Save face vertices/coordinates as one of the overlapping refined faces with the parent cell
                auto it  = boundaryFaceInfo.otherGridFaceIdx_to_newFacesToCoords.find(parentFaceIdx);
                if (it != boundaryFaceInfo.otherGridFaceIdx_to_newFacesToCoords.end()) {
                    auto& collectedFaces = *it; 
                    collectedFaces.second.push_back(newFace.value());
                }
                else {
                    auto& collectedFaces = boundaryFaceInfo.otherGridFaceIdx_to_newFacesToCoords[parentFaceIdx];
                    collectedFaces.push_back(newFace.value());
                }
                 
                if (!refinedFaceFullyContainedInParentFace) {
                    boundaryFaceInfo.vanishedFaceIdx_to_otherGridFaceIdxAndNewFaceToCoords[refinedFaceIdx].push_back(std::make_pair(parentFaceIdx, newFace.value()));
                }
            }
        }
    }
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

void provideRefinementParentFaceIdxRelations(const std::array<std::vector<int>, 6>&                      classifiedParentCellFaces,
                                             const Dune::cpgrid::CpGridData&                             parentGridData,
                                             const Dune::cpgrid::Entity<0>&                              parentElem,
                                             const std::array<int,3>&                                    cells_per_dim,
                                             CellRefinementBoundaryInfo&                                 cellRefBoundaryInfo,
                                             std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                             const std::vector<std::vector<int>>&                        oldGridFaceIdx_to_newGridFaceIdxList,
                                             const std::vector<std::vector<int>>&                        parentFaceIdx_to_newGridFaceIdxList,
                                             const std::vector<int>&                                     oldGridFaceIdx_fullyContainedIn_parentFaceIdx)
{
    bool hasOnlyOneFacePerType = cellRefBoundaryInfo.parentCellhasSingleFacePerType;

    auto addChildFace = [&](int newChildFaceIdx, int parentFaceIdx, std::vector<int>& children_faces)
    {
        children_faces.push_back(newChildFaceIdx);
        cellRefBoundaryInfo.boundaryRefinedFaceIdx_to_parentFaceIdx[newChildFaceIdx] = parentFaceIdx;
    };

    auto assignCorrectChildFace = [&](int oldChildFaceIdx, int parentFaceIdx, int faceType_count, std::vector<int>& children_faces)
    {
        if (hasOnlyOneFacePerType) {
            addChildFace(oldChildFaceIdx, parentFaceIdx, children_faces);
            return;
        }

        const auto& newFaceIdxList = oldGridFaceIdx_to_newGridFaceIdxList[oldChildFaceIdx];

        if ((faceType_count == 1) || (oldGridFaceIdx_fullyContainedIn_parentFaceIdx[oldChildFaceIdx] == parentFaceIdx)) {
            for (int newChildFaceIdx : newFaceIdxList)
                addChildFace(newChildFaceIdx, parentFaceIdx, children_faces);
            return;
        }

        if (newFaceIdxList.size()>1) {
            for (int newChildFaceIdx : parentFaceIdx_to_newGridFaceIdxList[parentFaceIdx]) {
                if (std::find(newFaceIdxList.begin(), newFaceIdxList.end(), newChildFaceIdx) != newFaceIdxList.end())
                    addChildFace(newChildFaceIdx, parentFaceIdx, children_faces);
            }
        }
    };

    const int k_faces = cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1);
    const int i_faces = (cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2];
    const int offSet = cells_per_dim[0]*cells_per_dim[1]*cells_per_dim[2];
    
    for (const auto& face : parentGridData.cellToFace(parentElem.index())) {

        const int parentFaceIdx = face.index();
        const auto parentFaceTag = parentGridData.faceTag(parentFaceIdx);
        const auto faceType_count = classifiedParentCellFaces[faceGroupIndex(parentFaceTag, face.orientation())].size();

        std::vector<int> newFaceIdxList{};

        if (parentFaceTag == face_tag::K_FACE) {
            for (int j = 0; j < cells_per_dim[1]; ++j) {
                for (int i = 0; i < cells_per_dim[0]; ++i) {
                    
                    int oldChildFaceIdx =  (j*cells_per_dim[0]) + i + (face.orientation()*offSet);
                    assignCorrectChildFace(oldChildFaceIdx, parentFaceIdx, faceType_count, newFaceIdxList);
                } 
            }
        }
        else if (parentFaceTag == face_tag::I_FACE) {  
            for (int k = 0; k < cells_per_dim[2]; ++k) {
                for (int j = 0; j < cells_per_dim[1]; ++j) {
                    
                    int oldChildFaceIdx = k_faces + (k*cells_per_dim[1]) + j + (face.orientation()*offSet);
                    assignCorrectChildFace(oldChildFaceIdx, parentFaceIdx, faceType_count, newFaceIdxList);
                } 
            }
        } 
        else if (parentFaceTag == face_tag::J_FACE) {
            for (int i = 0; i < cells_per_dim[0]; ++i) {
                for (int k = 0; k < cells_per_dim[2]; ++k) {
                    
                    int oldChildFaceIdx = k_faces + i_faces + (i*cells_per_dim[2]) + k + (face.orientation()*offSet);
                    assignCorrectChildFace(oldChildFaceIdx, parentFaceIdx, faceType_count, newFaceIdxList);
                } 
            } 
        } 
        newFaceIdxList.shrink_to_fit();
        faceInMarkedElemAndRefinedFaces[face.index()].push_back(std::make_pair(parentElem.index(), newFaceIdxList));
    }
}

void populateVertexData(const Dune::cpgrid::CpGridData&                                    oldGrid,
                        const std::map<int,BoundaryFaceInfo>&                              parentCellFaceIdx_to_boundaryFaceInfo,
                        GeomData&                                                          newGridGeomData,
                        std::map<Dune::FieldVector<double,3>, int, FieldVectorLess>&       newGridVertexCoordinates_to_newGridVertexIdx,
                        std::unordered_map<int, int>&                                      newGridBoundaryVertexIdx_to_parentFaceIdx,
                        std::vector<bool>&                                                 newGridBoundaryVertexCoincidesWithParentVertex)
{
    int numVertices = oldGrid.size(3);
    const auto& oldGridVertices = *(oldGrid.getGeometry().geomVector(std::integral_constant<int,3>()));
    Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<0,3>>& newGridVertices = *(newGridGeomData.geometries.geomVector(std::integral_constant<int,3>()));

    int newVerticesCount = 0;
    for (const auto& [_, boundInfo] : parentCellFaceIdx_to_boundaryFaceInfo) {
        newVerticesCount += boundInfo.foundNewVertices.size();
    }
    newGridVertices.resize(numVertices + newVerticesCount);
      
    // add the "old" vertices (the parent-cell-face-aware are includeed here)
    for (int i = 0; i < oldGrid.size(3); ++i) {
        newGridVertices[i] = oldGridVertices.get(i);
    }

   for (const auto& [parentFaceIdx, boundFaceInfo] : parentCellFaceIdx_to_boundaryFaceInfo) {
        for (const auto& vertex : boundFaceInfo.foundNewVertices) {
            
            newGridVertexCoordinates_to_newGridVertexIdx[vertex] = numVertices; 
            newGridVertices[numVertices] = Dune::cpgrid::Geometry<0, 3>(vertex);
            newGridBoundaryVertexIdx_to_parentFaceIdx[numVertices] = parentFaceIdx;
            
            ++numVertices;
        }
   }
   newGridBoundaryVertexCoincidesWithParentVertex.resize(numVertices, false); // false is used to initialized the new entries, old entries stay the same
}

void populateVertexData(const Dune::cpgrid::CpGridData&                                    oldGrid,
                        const std::map<Dune::FieldVector<double,3>, int, FieldVectorLess>& extendedParentCellToPointVertexMap, 
                        const std::set<Dune::FieldVector<double,3>, FieldVectorLess>&      foundNewVertices,                             // boundaryFaceInfo
                        GeomData&                                                          newGridGeomData,
                        std::map<Dune::FieldVector<double,3>, int, FieldVectorLess>&       newGridVertexCoordinates_to_newGridVertexIdx, // boundaryFaceInfo
                        std::vector<std::array<int,2>>&                                    parentBoundaryVertexIdx_to_newGridBoundaryVertexIdx)
{
    int numVertices = oldGrid.size(3);
    const auto& oldGridVertices = *(oldGrid.getGeometry().geomVector(std::integral_constant<int,3>()));
    Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<0,3>>& newGridVertices = *(newGridGeomData.geometries.geomVector(std::integral_constant<int,3>()));
    newGridVertices.resize(numVertices + foundNewVertices.size());

    auto addParentBoundaryVertexMapping = [&](Dune::FieldVector<double,3> vertex, int i)
    {
        auto it = extendedParentCellToPointVertexMap.find(vertex);

        if (it == extendedParentCellToPointVertexMap.end()) {
            return;
        }

        const auto parentGridVertexIdx = it->second;

        auto iit = std::find_if(parentBoundaryVertexIdx_to_newGridBoundaryVertexIdx.begin(),
                                parentBoundaryVertexIdx_to_newGridBoundaryVertexIdx.end(),
                                [parentGridVertexIdx](const std::array<int, 2>& arr) {
                                    return arr[0] == parentGridVertexIdx;
                                });

        if (iit == parentBoundaryVertexIdx_to_newGridBoundaryVertexIdx.end()) {
            parentBoundaryVertexIdx_to_newGridBoundaryVertexIdx.push_back( {parentGridVertexIdx, i} );
        }
    };
    
    for (int i = 0; i < oldGrid.size(3); ++i) {
        addParentBoundaryVertexMapping(oldGridVertices.get(i).center(), i);
        
        newGridVertices[i] = oldGridVertices.get(i);
    }
        
    for (const auto& vertex : foundNewVertices) {
        addParentBoundaryVertexMapping(vertex, numVertices);
        
        newGridVertexCoordinates_to_newGridVertexIdx[vertex] = numVertices;
        newGridVertices[numVertices] = Dune::cpgrid::Geometry<0, 3>(vertex);

        ++numVertices;
    }
}

bool isAtGridBoundary(const Dune::cpgrid::CpGridData& grid,
                      const Dune::cpgrid::Entity<0>&  element)
{
    for (const auto& face : grid.cellToFace(element.index())) {
        if (grid.faceToCellSize(face.index()) == 1)
            return true;
    }
    return false;
}

void populateFaceData(int                                        oldGridNumFaces,
                      const Dune::cpgrid::CpGridData&            oldGrid,
                      GridModificationMapping&                   modificationMapps,
                      GeomData&                                  newGridGeomData,
                      const BoundaryFaceInfo&                    boundaryFaceInfo,
                      std::vector<std::vector<int>>&             parentGridFaceIdx_to_newGridFaceIdxList,
                      std::unordered_map<int,int>&               newGridBoundaryVertexIdx_to_parentGridFaceIdx)
{
    modificationMapps.oldGridFaceIdx_to_newGridFaceIdxList.resize(oldGridNumFaces);

    Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<2,3>>& newGridFaces =
        *(newGridGeomData.geometries.geomVector(std::integral_constant<int,1>()));
    Dune::cpgrid::EntityVariableBase<enum face_tag>& mutable_newGrid_face_tags = newGridGeomData.face_tags;
    Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>& mutable_newGrid_face_normals = newGridGeomData.face_normals;
    
    int upperBoundFaceSize = oldGridNumFaces; // total faces before correcting the data
    for (const auto& [_, overlapFacesInfo] : boundaryFaceInfo.otherGridFaceIdx_to_newFacesToCoords) {
        upperBoundFaceSize += overlapFacesInfo.size();
    }
    newGridFaces.resize(upperBoundFaceSize); // maybe larger than needed
    mutable_newGrid_face_tags.resize(upperBoundFaceSize);
    mutable_newGrid_face_normals.resize(upperBoundFaceSize);

    std::vector<std::vector<int>> aux_newGrid_face_to_point{};
    aux_newGrid_face_to_point.resize(upperBoundFaceSize);

    int newGrid_num_points = 0;
    int newGrid_face_count = 0;

    auto registerFace = [&](int oldFaceIdx, int newFaceIdx, int parentGridFaceIdx,
                            auto&& faceToPoint, auto&& geometry, auto&& normal, auto tag)
    {
        modificationMapps.oldGridFaceIdx_to_newGridFaceIdxList[oldFaceIdx].push_back(newFaceIdx);

        if (parentGridFaceIdx != -1) {
            parentGridFaceIdx_to_newGridFaceIdxList[parentGridFaceIdx].push_back(newGrid_face_count);
        }

        newGridFaces[newFaceIdx] = geometry;
        mutable_newGrid_face_tags[newFaceIdx] = tag;
        mutable_newGrid_face_normals[newFaceIdx] = normal;

        aux_newGrid_face_to_point[newGrid_face_count] = std::vector<int>(faceToPoint.begin(), faceToPoint.end());
        newGrid_num_points += faceToPoint.size();

        const auto faceToCell = oldGrid.faceToCell(oldFaceIdx);
        newGridGeomData.face_to_cell.appendRow(faceToCell.begin(), faceToCell.end());

        ++newGrid_face_count;
    };

    for (int i = 0; i < oldGrid.numFaces(); ++i) {
        const auto& newFacesInfo =  boundaryFaceInfo.vanishedFaceIdx_to_otherGridFaceIdxAndNewFaceToCoords[i];
        const int parentGrid_containerFaceIdx = boundaryFaceInfo.face_fullyContainedIn_otherGridFaceIdx[i];

        if ( newFacesInfo.empty() || (parentGrid_containerFaceIdx != -1 /* invalidIdx */)) {

            const auto face =  Dune::cpgrid::EntityRep<1>(i, true);
            const auto faceToPoint = oldGrid.faceToPoint(i);

            registerFace(i, newGrid_face_count, parentGrid_containerFaceIdx, faceToPoint,
                         (*oldGrid.getGeometry().geomVector(std::integral_constant<int,1>()))[face],
                         oldGrid.faceNormals(i), oldGrid.faceTag(i));
        }
        else {
            std::vector<int> newFaceIndices{};

            for (const auto& [parentGridFaceIdx, newFaceToCoord] : newFacesInfo) {

                newFaceIndices.push_back(newGrid_face_count);
                
                const auto [faceCenter, faceArea, faceNormal] = computeFaceCenterAreaNormal(newFaceToCoord);
                const auto faceToPoint =  modifyFaceToPoint(newFaceToCoord,
                                                            modificationMapps,
                                                            boundaryFaceInfo,
                                                            parentGridFaceIdx,
                                                            &newGridBoundaryVertexIdx_to_parentGridFaceIdx);

                registerFace(i, newGrid_face_count, parentGridFaceIdx, faceToPoint,
                             Dune::cpgrid::Geometry<2,3>(faceCenter, faceArea),
                             faceNormal, oldGrid.faceTag(i));
            }
        }
    }
    populateFaceToPoint(newGrid_face_count, newGrid_num_points, aux_newGrid_face_to_point,  newGridGeomData);
}

void populateFaceToPoint(int face_count,
                         int num_points,
                         const std::vector<std::vector<int>>& face_to_point,
                         GeomData& geomData)
{
    geomData.face_to_point.reserve(face_count, num_points);
    for (int face = 0; face < face_count; ++face) {
        geomData.face_to_point.appendRow(face_to_point[face].begin(), face_to_point[face].end());
    }
}

std::vector<int> modifyFaceToPoint(const std::vector<Dune::FieldVector<double,3>>& newFaceToCoord,
                                   const GridModificationMapping&                  modificationMap,
                                   const BoundaryFaceInfo&                         boundaryFaceInfo,
                                   int                                             parentGridFaceIdx,
                                   std::unordered_map<int,int>*                    newGridBoundaryVertexIdx_to_parentGridFaceIdx)
{
    std::vector<int> faceToPoint;

    for (const auto& vertex : newFaceToCoord) {
        int vertexIdx = -1;
        // Check if the vertex is new (e.g., added when modifying the grid) 
        auto newIt = modificationMap.newGridVertexCoordinates_to_newGridVertexIdx.find(vertex);
        if (newIt != modificationMap.newGridVertexCoordinates_to_newGridVertexIdx.end()) {
            vertexIdx = newIt->second;
        }
        else { // Vertex already existed in the original grid ()
            auto oldIt = boundaryFaceInfo.boundaryVertexCoordinates_to_vertexIdx.find(vertex);
            assert(oldIt != boundaryFaceInfo.boundaryVertexCoordinates_to_vertexIdx.end());

            vertexIdx = oldIt->second;
        }
        if ((parentGridFaceIdx != -1 /* invalid-index */) &&  (newGridBoundaryVertexIdx_to_parentGridFaceIdx != nullptr)) {
            (*newGridBoundaryVertexIdx_to_parentGridFaceIdx)[vertexIdx] = parentGridFaceIdx;
        }
        assert(vertexIdx >= 0);
        faceToPoint.push_back(vertexIdx);
    }

    return faceToPoint;
}

void populateFaceData(int                                   oldGridNumFaces,
                      const Dune::cpgrid::CpGridData&       oldGrid,
                      GridModificationMapping&              modificationMap,
                      GeomData&                             newGridGeomData,
                      const std::map<int,BoundaryFaceInfo>& boundaryFacesInfo,
                      std::vector<int>&                     newGridFaceIdx_to_parentFaceIdx)
{
    int invalidIdx = -1;
    modificationMap.oldGridFaceIdx_to_newGridFaceIdxList.resize(oldGridNumFaces);

    Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<2,3>>& newGridFaces = *(newGridGeomData.geometries.geomVector(std::integral_constant<int,1>()));
    Dune::cpgrid::EntityVariableBase<enum face_tag>& mutable_newGridFaceTags = newGridGeomData.face_tags;
    Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>& mutable_newGridFaceNormals = newGridGeomData.face_normals;

    // Estimate new size of total number of faces in correctedGridData1
    int upperBoundFaceSize = oldGridNumFaces; // total faces before correcting the data
    for (const auto& [unused_parentFaceIdx, boundaryFaceInfo] : boundaryFacesInfo) {
        for (const auto& [unsused_neighborGrid_faceIdx, overlapFacesInfo] : boundaryFaceInfo.otherGridFaceIdx_to_newFacesToCoords) {
            upperBoundFaceSize += overlapFacesInfo.size();
        }
    }
    newGridFaces.resize(upperBoundFaceSize); // maybe larger than needed
    mutable_newGridFaceTags.resize(upperBoundFaceSize);
    mutable_newGridFaceNormals.resize(upperBoundFaceSize);
    
    modificationMap.newGridFaceIdx_to_oldGridFaceIdx.resize(upperBoundFaceSize, invalidIdx);

    // intentional copy
    const auto hasParent = newGridFaceIdx_to_parentFaceIdx; 
    newGridFaceIdx_to_parentFaceIdx.resize(upperBoundFaceSize, invalidIdx);

    std::vector<std::vector<int>> aux_newGrid_face_to_point{};
    aux_newGrid_face_to_point.resize(upperBoundFaceSize);
    
    int newGrid_num_points = 0;
    int newGrid_face_count = 0;
    
    auto registerFace = [&](int oldFaceIdx, int newFaceIdx, int parentFaceIdx,
                            auto&& faceToPoint, auto&& geometry, auto&& normal, auto tag)
    {
        modificationMap.oldGridFaceIdx_to_newGridFaceIdxList[oldFaceIdx].push_back(newFaceIdx);
        modificationMap.newGridFaceIdx_to_oldGridFaceIdx[newFaceIdx] = oldFaceIdx;

        if (parentFaceIdx != -1) {
            newGridFaceIdx_to_parentFaceIdx[newFaceIdx] = parentFaceIdx;
        }

        newGridFaces[newFaceIdx] = geometry;
        mutable_newGridFaceTags[newFaceIdx] = tag;
        mutable_newGridFaceNormals[newFaceIdx] = normal;

        aux_newGrid_face_to_point[newGrid_face_count] = std::vector<int>(faceToPoint.begin(), faceToPoint.end());
        newGrid_num_points += faceToPoint.size();

        const auto faceToCell = oldGrid.faceToCell(oldFaceIdx);
        newGridGeomData.face_to_cell.appendRow(faceToCell.begin(), faceToCell.end());

        ++newGrid_face_count;
    };
    
    for (int i = 0; i < oldGrid.numFaces(); ++i) {
        
        const auto face =  Dune::cpgrid::EntityRep<1>(i, true); 
        int parentFaceIdx = hasParent[i]; 
        auto itP = boundaryFacesInfo.find(parentFaceIdx);
            
        if ((parentFaceIdx == -1) || (itP == boundaryFacesInfo.end()) ||
            (itP != boundaryFacesInfo.end() && (itP->second.faceFullyContainedInNeighbor[i]))) {

            const auto faceToPoint = oldGrid.faceToPoint(i);

            registerFace(i, newGrid_face_count, parentFaceIdx, faceToPoint,
                         (*oldGrid.getGeometry().geomVector(std::integral_constant<int,1>()))[face],
                         oldGrid.faceNormals(i), oldGrid.faceTag(i));
        }
        else {
            const auto& boundaryFaceInfo = itP->second;
            const auto& newFacesInfo = boundaryFaceInfo.vanishedFaceIdx_to_otherGridFaceIdxAndNewFaceToCoords[i];
            
            if (!newFacesInfo.empty()) {
                for (const auto& [grid2_faceIdx, newFaceToCoord] : newFacesInfo) {
                    
                    const auto faceToPoint =  modifyFaceToPoint(newFaceToCoord,
                                                                modificationMap,
                                                                boundaryFaceInfo);

                    const auto [faceCenter, faceArea, faceNormal] = computeFaceCenterAreaNormal(newFaceToCoord);
                    
                    registerFace(i, newGrid_face_count, parentFaceIdx, faceToPoint,
                                 Dune::cpgrid::Geometry<2,3>(faceCenter, faceArea),
                                 faceNormal, oldGrid.faceTag(i));
                }
            }
        }
    }
    populateFaceToPoint(newGrid_face_count, newGrid_num_points, aux_newGrid_face_to_point, newGridGeomData);
}

void populateCellData(const Dune::cpgrid::CpGridData& sourceGrid,
                      GeomData&                       geomData)
{
    Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<3,3>>& cells = *(geomData.geometries.geomVector(std::integral_constant<int,0>()));
    cells.resize(sourceGrid.size(0));
    geomData.cell_to_point.resize(sourceGrid.size(0));

    for (int i = 0; i < sourceGrid.size(0); ++i) {
        
        geomData.cell_to_point[i] = sourceGrid.cellToPoint(i);
        int* indices_storage_ptr = geomData.cell_to_point[i].data();

        const auto element = Dune::cpgrid::Entity<0>(sourceGrid, i, true);
        
        cells[i] = Dune::cpgrid::Geometry<3,3>(element.geometry().center(), element.geometry().volume(),
                                               geomData.geometries.geomVector(std::integral_constant<int,3>()),
                                               indices_storage_ptr);
    }
    geomData.face_to_cell.makeInverseRelation(geomData.cell_to_face);
}

void populateBoundaryRefVertexCoincidesWithParentVertex(int numVertices,
                                                        CellRefinementBoundaryInfo& cellRefBoundaryInfo)
{
    cellRefBoundaryInfo.boundaryRefinedVertexCoincidesWithParentVertex.resize(numVertices, false);
    for (const auto& [_, equivRefVertexIdx] : cellRefBoundaryInfo.parentVertexIdx_to_boundaryRefinedVertexIdx) {
        cellRefBoundaryInfo.boundaryRefinedVertexCoincidesWithParentVertex[equivRefVertexIdx] = true;
    }
}

void makeCellRefinementParentFaceAware(const std::array<std::vector<int>,6>&                       classifiedParentCellFaces,
                                       const Dune::cpgrid::CpGridData&                             parentGrid,
                                       const Dune::cpgrid::Entity<0>&                              parentCellElem,
                                       const Dune::cpgrid::CpGridData&                             cellRefGrid,
                                       Dune::cpgrid::CpGridData&                                   correctedCellRefGrid,
                                       GeomData&                                                   correctedCellRefGeomData,
                                       CellRefinementBoundaryInfo&                                 cellRefBoundaryInfo,
                                       std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                       const std::array<int,3>&                                    cells_per_dim)
{
    int numFaces = cellRefGrid.numFaces();
    BoundaryFaceInfo boundaryFaceInfo{};
    boundaryFaceInfo.vanishedFaceIdx_to_otherGridFaceIdxAndNewFaceToCoords.resize(numFaces);
    boundaryFaceInfo.face_fullyContainedIn_otherGridFaceIdx.resize(numFaces, /*invalidIdx */-1);

    for (int i = 0; i < cellRefGrid.size(0); ++i) {

        const auto refinedElem = Dune::cpgrid::Entity<0>(cellRefGrid, i, true);
        if (!isAtGridBoundary(cellRefGrid, refinedElem))
            continue;
        
        collectNewVertices(parentGrid,
                           parentCellElem,
                           classifiedParentCellFaces,
                           cellRefGrid,
                           refinedElem,
                           boundaryFaceInfo);
    }
    
    GridModificationMapping modificationMapps{};
     
    const auto& extendedParentCellToPointVertexMap = buildExtendedCellPointVertexMap(parentGrid, parentCellElem.index());;
    populateVertexData(cellRefGrid,
                       extendedParentCellToPointVertexMap, 
                       boundaryFaceInfo.foundNewVertices,
                       correctedCellRefGeomData,
                       modificationMapps.newGridVertexCoordinates_to_newGridVertexIdx,
                       cellRefBoundaryInfo.parentVertexIdx_to_boundaryRefinedVertexIdx);
    
    std::vector<std::vector<int>> parentFaceIdx_to_correctedFaceIdxList{}; // empty entry for face not affected, is it better std::unordered_map? only few entries are populated...
    parentFaceIdx_to_correctedFaceIdxList.resize(parentGrid.numFaces());
    
    populateFaceData(numFaces,
                     cellRefGrid,
                     modificationMapps,      
                     correctedCellRefGeomData,
                     boundaryFaceInfo,
                     parentFaceIdx_to_correctedFaceIdxList,                            
                     cellRefBoundaryInfo.boundaryRefinedVertexIdx_to_parentFaceIdx);
    
    populateCellData(cellRefGrid, correctedCellRefGeomData);
    
    cellRefBoundaryInfo.boundaryRefinedFaceIdx_to_parentFaceIdx.resize(correctedCellRefGrid.numFaces(), /* invalidIdx */ -1);
    provideRefinementParentFaceIdxRelations(classifiedParentCellFaces,
                                            parentGrid,
                                            parentCellElem,
                                            cells_per_dim,
                                            cellRefBoundaryInfo, 
                                            faceInMarkedElemAndRefinedFaces,
                                            modificationMapps.oldGridFaceIdx_to_newGridFaceIdxList,
                                            parentFaceIdx_to_correctedFaceIdxList,
                                            boundaryFaceInfo.face_fullyContainedIn_otherGridFaceIdx);

    populateBoundaryRefVertexCoincidesWithParentVertex(correctedCellRefGrid.size(3), cellRefBoundaryInfo);
}


void computeNewRefinedGeometriesOnSharedCoarseFace(const std::vector<int>& refinedFaces1,
                                                   const Dune::cpgrid::CpGridData& parentFaceAwareCellRefinement1,
                                                   BoundaryFaceInfo& boundaryFaceInfo,
                                                   const std::vector<int>& refinedFaces2,
                                                   const Dune::cpgrid::CpGridData& parentFaceAwareCellRefinement2)
{
    for (const auto& refinedFace1 : refinedFaces1) {
        const auto face1 = Dune::cpgrid::EntityRep<1>(refinedFace1, true);
        
        addFaceVerticesCoordinatesToIndex(parentFaceAwareCellRefinement1, refinedFace1,
                                          boundaryFaceInfo.boundaryVertexCoordinates_to_vertexIdx);
            
        for (const auto& refinedFace2 : refinedFaces2) {
            const auto face2 = Dune::cpgrid::EntityRep<1>(refinedFace2, true);

            bool face1FullyContainedInFace2 = false;
            const auto newFace = computeFaceOverlapVertices(face2,
                                                            parentFaceAwareCellRefinement2,
                                                            face1,
                                                            parentFaceAwareCellRefinement1,
                                                            boundaryFaceInfo.foundNewVertices,
                                                            face1FullyContainedInFace2);

            boundaryFaceInfo.faceFullyContainedInNeighbor[refinedFace1] = boundaryFaceInfo.faceFullyContainedInNeighbor[refinedFace1] || face1FullyContainedInFace2;

            if (face1FullyContainedInFace2){
                boundaryFaceInfo.face_fullyContainedIn_otherGridFaceIdx[refinedFace1] = refinedFace2;
            }

            if (newFace.has_value() && !newFace.value().empty()) {
                // Save face vertices/coordinates as one of the overlapping refined faces with the parent cell
                auto it1  = boundaryFaceInfo.otherGridFaceIdx_to_newFacesToCoords.find(refinedFace2);
                if (it1 != boundaryFaceInfo.otherGridFaceIdx_to_newFacesToCoords.end()) {
                    auto& collectedFaces = *it1;
                    collectedFaces.second.push_back(newFace.value());
                }
                else {
                    auto& collectedFaces = boundaryFaceInfo.otherGridFaceIdx_to_newFacesToCoords[refinedFace2];
                    collectedFaces.push_back(newFace.value());
                }
                if (!face1FullyContainedInFace2) {
                    boundaryFaceInfo.vanishedFaceIdx_to_otherGridFaceIdxAndNewFaceToCoords[refinedFace1].push_back(std::make_pair(refinedFace2, newFace.value()));
                }
            }
        }
    }
}

std::pair<std::shared_ptr<Dune::cpgrid::CpGridData>, GridModificationMapping>
makeCellRefinementNeighborsAware(const Dune::cpgrid::Entity<0>& parentCell,
                                 const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& cellRefs, // to get neighbor information
                                 const Dune::cpgrid::CpGridData& parentGrid,
                                 const Dune::cpgrid::CpGridData& cellRefGrid,
                                 CellRefinementBoundaryInfo& cellRefGridBoundaryInfo,
                                 std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces) // needs to be updated
{
    std::map<int,BoundaryFaceInfo> boundaryFacesInfo{};
    int numFaces = cellRefGrid.numFaces();

    for (const auto& parentFace : parentGrid.cellToFace(parentCell.index())) {
      
        const auto faceToCell = parentGrid.faceToCell(parentFace.index());

        if ((faceToCell.size() == 1) || (faceInMarkedElemAndRefinedFaces[parentFace.index()].size()<=1)){
            continue;  // face at parent grid boundary does need to be corrected (cell refinement is already parentCellFaces aware).
        }
        
        assert(faceToCell.size() == 2);
        assert(faceInMarkedElemAndRefinedFaces[parentFace.index()].size() == 2);

        const auto& [p1, refinedFaces1] = faceInMarkedElemAndRefinedFaces[parentFace.index()][0];
        const auto& [p2, refinedFaces2] = faceInMarkedElemAndRefinedFaces[parentFace.index()][1];
        
        int neighborCellRefIdx = (p1 == parentCell.index())? p2 : p1;

        const auto& cellRefFaces = (p1 == parentCell.index())? refinedFaces1 : refinedFaces2;
        const auto& neighborCellRefFaces = (p1 == parentCell.index())? refinedFaces2 : refinedFaces1;

        auto& boundaryFaceInfo = boundaryFacesInfo[parentFace.index()];
        boundaryFaceInfo.faceFullyContainedInNeighbor.resize(numFaces, false);
        boundaryFaceInfo.face_fullyContainedIn_otherGridFaceIdx.resize(numFaces, -1);
        boundaryFaceInfo.vanishedFaceIdx_to_otherGridFaceIdxAndNewFaceToCoords.resize(numFaces);

        computeNewRefinedGeometriesOnSharedCoarseFace(cellRefFaces,
                                                      cellRefGrid,
                                                      boundaryFaceInfo,
                                                      neighborCellRefFaces,
                                                      *cellRefs[neighborCellRefIdx]);
    }
    
    std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>> neighborAwareCellRefData;
    std::shared_ptr<Dune::cpgrid::CpGridData> neighborAwareCellRef_ptr = std::make_shared<Dune::cpgrid::CpGridData>(neighborAwareCellRefData); // ccobj_
    auto& neighborAwareCellRef = *neighborAwareCellRef_ptr;
    GeomData neighborAwareGeomData(neighborAwareCellRef);

    GridModificationMapping modificationMap{};
    
    populateVertexData(cellRefGrid,
                       boundaryFacesInfo,
                       neighborAwareGeomData,
                       modificationMap.newGridVertexCoordinates_to_newGridVertexIdx, 
                       cellRefGridBoundaryInfo.boundaryRefinedVertexIdx_to_parentFaceIdx,
                       cellRefGridBoundaryInfo.boundaryRefinedVertexCoincidesWithParentVertex);
    
    populateFaceData(numFaces,
                     cellRefGrid,
                     modificationMap,
                     neighborAwareGeomData,
                     boundaryFacesInfo,
                     cellRefGridBoundaryInfo.boundaryRefinedFaceIdx_to_parentFaceIdx);

    populateCellData(cellRefGrid, neighborAwareGeomData);

    return {neighborAwareCellRef_ptr, modificationMap};
}

void makeCellRefinementsNeighborsAware(std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& cellRefs, // to get neighbor information
                                       std::vector<CellRefinementBoundaryInfo>& cellRefsBoundaryInfo,
                                       const Dune::cpgrid::CpGridData& parentGrid,
                                       std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                       std::vector<GridModificationMapping>& modificationMapps) // needs to be updated
{
    const auto cellRefsBeforeCorrection = cellRefs; // copy on propose
    
    for (int elemIdx = 0; elemIdx < parentGrid.size(0); ++elemIdx) {
        auto& cellRef = cellRefs[elemIdx];
        if (cellRef == nullptr)
            continue;
        
        const auto& parentCell = Dune::cpgrid::Entity<0>(parentGrid, elemIdx, true);
        
        auto [correctedCellRef, correctedModificationMap] =  makeCellRefinementNeighborsAware(parentCell,
                                                                                            cellRefsBeforeCorrection, // to get neighbor information             
                                                                                            parentGrid,
                                                                                            *cellRef,
                                                                                            cellRefsBoundaryInfo[elemIdx],       
                                                                                            faceInMarkedElemAndRefinedFaces);
        cellRefs[elemIdx] = correctedCellRef;
        modificationMapps[elemIdx] = correctedModificationMap;
    }

    for (std::size_t faceIdx = 0; faceIdx < faceInMarkedElemAndRefinedFaces.size(); ++faceIdx) {
        
        for (std::size_t cellRefCount = 0; cellRefCount < faceInMarkedElemAndRefinedFaces[faceIdx].size(); ++cellRefCount) {

            auto& [parentCell, refinedFaces] = faceInMarkedElemAndRefinedFaces[faceIdx][cellRefCount];  
            const auto& oldToNewList = modificationMapps[parentCell].oldGridFaceIdx_to_newGridFaceIdxList;

            if (oldToNewList.empty()) 
                continue;
              
            std::vector<int> updatedFaceIndices{};

            for (const auto& oldIdx : refinedFaces) {
                for (const auto& newIdx : oldToNewList[oldIdx]) {
                    assert(newIdx>=0);
                    updatedFaceIndices.push_back(newIdx);
                }
            }
            faceInMarkedElemAndRefinedFaces[faceIdx][cellRefCount] = {parentCell, updatedFaceIndices};
        }
    }   
}

bool areClose(const Dune::FieldVector<double,3>& v,
              const Dune::FieldVector<double,3>& w)
{
    return (std::abs(v[0] - w[0]) < 1e-12) && (std::abs(v[1] - w[1]) < 1e-12) && (std::abs(v[2] - w[2])< 1e-12);
}

std::optional<int> findVertexIdxInFaces(const Dune::cpgrid::CpGridData&    grid,
                                        const std::vector<int>&            candidateFaces,
                                        const Dune::FieldVector<double,3>& targetVertex)
{
    for (const auto& faceIdx : candidateFaces) {
        const auto matchingVertexIdx = findMatchingVertexIdx(grid, faceIdx, targetVertex);
        
        if (matchingVertexIdx)
            return matchingVertexIdx;
    }
    return std::nullopt;
}

std::optional<int> findMatchingVertexIdx(const Dune::cpgrid::CpGridData&    grid,
                                         int                                faceIndex,
                                         const Dune::FieldVector<double,3>& targetVertex)
{
    const auto& facePoints = grid.faceToPoint(faceIndex);
    
    for (const auto& pointIndex : facePoints) {
        const auto& pointPosition = Dune::cpgrid::Entity<3>(grid, pointIndex, true).geometry().center();

        if (areClose(pointPosition, targetVertex))
            return pointIndex;
    }
    return std::nullopt;
}

std::vector<Dune::FieldVector<double,3>> getFaceVertexCoordinates(const Dune::cpgrid::CpGridData& grid,
                                                                  int                             faceIdx)
{
    std::vector<Dune::FieldVector<double,3>> faceCoords{};
    
    for (const auto& vertexIdx : grid.faceToPoint(faceIdx)){
        faceCoords.push_back(Dune::cpgrid::Entity<3>(grid, vertexIdx, true).geometry().center());
    }
    return faceCoords;
}

std::optional<int> findMatchingFaceIdx(const Dune::cpgrid::CpGridData& targetGrid,
                                       const std::vector<int>&         candidateFaces,
                                       const Dune::cpgrid::CpGridData& sourceGrid,
                                       int                             sourceFaceIdx)
{
    const auto& sourceFaceCoords = getFaceVertexCoordinates(sourceGrid, sourceFaceIdx);
    
    for (const auto& candidateFaceIdx : candidateFaces) {
        
        const auto& candidateFaceCoords = getFaceVertexCoordinates(targetGrid, candidateFaceIdx);

        if (sourceFaceCoords.size() != candidateFaceCoords.size())
            continue;

        bool matchingFaceFound = std::all_of(sourceFaceCoords.begin(), sourceFaceCoords.end(),
                                             [&](const auto& v) {
                                                 return std::any_of(candidateFaceCoords.begin(), candidateFaceCoords.end(),
                                                                    [&](const auto& w) {
                                                                        return areClose(v, w);
                                                                    });
                                             });
        if (matchingFaceFound) {
            return candidateFaceIdx;
        }
    }
    return std::nullopt;
}

} // namespace Lgr
} // namespace Opm
