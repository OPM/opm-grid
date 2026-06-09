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

#ifndef OPM_GRID_CPGRID_LGRFAULTHELPERS_HEADER_INCLUDED
#define OPM_GRID_CPGRID_LGRFAULTHELPERS_HEADER_INCLUDED


#include <dune/common/dotproduct.hh>
#include <dune/common/fmatrixev.hh>
#include <dune/common/fvector.hh>
#include <dune/common/version.hh>
#include <dune/grid/common/mcmgmapper.hh>


#include <opm/grid/CpGrid.hpp>
#include <opm/grid/common/Volumes.hpp>
#include <opm/grid/cpgrid/LgrHelpers.hpp>


#include <array>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility> // for std::pair
#include <vector>

namespace Dune
{
namespace cpgrid
{
class CpGridData;
class DefaultGeometryPolicy;
}
}

namespace Opm
{

template<class Grid> class LevelCartesianIndexMapper;

namespace Lgr
{

struct ParentAwareCellRefinement {

    std::vector<std::unordered_map<int,int>> extraRefinedCornIdx_to_parentFaceIdx{};
    // extraRefinedCornIdx_to_parentFaceIdx.resize(currentLeafData().size(0));

    std::vector<std::vector<std::vector<int>>> refinedFace_to_parentFaces{};
    // refinedFace_to_parentFaces.resize(currentLeafData().size(0));

    std::vector<std::vector<bool>> singleCellRef_coincideWithCoarseCorner{};
    // singleCellRef_coincideWithCoarseCorner.resize(currentLeafData().size(0));

    std::vector<bool> singleCellRef_hasOnlyOneFacePerType{};
   // singleCellRef_hasOnlyOneFacePerType.resize(currentLeafData().size(0), true);
};


struct FieldVectorLess {
    bool operator()(const Dune::FieldVector<double,3>& v,
                    const Dune::FieldVector<double,3>& w) const
    {
        for (int i = 0; i < 3; ++i) {
            if (v[i] < w[i]) return true;
            if (v[i] > w[i]) return false;
        }
        return false;
    }
};


/// @brief Computes the index of a face group from its tag and orientation.
///
/// Maps a face tag and orientation pair to one of six face-group indices:
/// - 0: I-face, orientation == false
/// - 1: I-face, orientation == true
/// - 2: J-face, orientation == false
/// - 3: J-face, orientation == true
/// - 4: K-face, orientation == false
/// - 5: K-face, orientation == true
///
/// @param [in] tag         Face tag (I_FACE, J_FACE, or K_FACE).
/// @param [in] orientation Face orientation.
///
/// @return Index in the range [0, 5] corresponding to the face group.
///
/// @note Assumes that the face tags are encoded as
///       I_FACE == 0, J_FACE == 1, and K_FACE == 2.
constexpr int faceGroupIndex(int tag, bool orientation)
{
    return 2*tag + static_cast<int>(orientation);
}

/// @brief Groups the face indices of an element by face type and orientation.
///
/// Collects the indices of all faces incident to the given element and
/// groups them into six categories corresponding to the face tag
/// (I, J, or K) and orientation (false or true).
///
/// The returned array is organized as follows:
/// - [0]: I-face, orientation == false
/// - [1]: I-face, orientation == true
/// - [2]: J-face, orientation == false
/// - [3]: J-face, orientation == true
/// - [4]: K-face, orientation == false
/// - [5]: K-face, orientation == true
///
/// @param [in] gridData Grid data containing face connectivity and face tags.
/// @param [in] element  Element whose incident faces are to be grouped.
///
/// @return An array containing the face indices grouped by face tag and orientation.
///
/// @throws std::logic_error If a face has an invalid tag and cannot be
///         classified as an I-, J-, or K-face.
std::array<std::vector<int>,6> groupFaceIndicesByType(const Dune::cpgrid::CpGridData& gridData,
                                                      const Dune::cpgrid::Entity<0>& element);

/// @brief Checks whether each face group contains at most one face.
///
/// The face groups correspond to the six combinations of face tag
/// (I, J, K) and orientation (false, true). The function returns
/// true if no group contains more than one face index.
///
/// @param classifiedFaces Face indices grouped by face tag and orientation.
///
/// @return true if every face group contains at most one face; otherwise false.
bool hasAtMostOneFacePerGroup(const std::array<std::vector<int>,6>& classifiedFaces);

/// @brief Computes the center, area, and unit normal of a quadrilateral face.
///
/// The face center is computed as the arithmetic mean of the four face
/// vertices. The face area is obtained by partitioning the face into four
/// triangles formed by the face center and each face edge and summing their
/// areas. The face normal is computed as the normalized cross product of two
/// vectors spanning the face.
///
/// @param [in] faceToCoord Coordinates of the four face vertices.
///
/// @return A tuple containing:
/// - std::get<0>(result): face center.
/// - std::get<1>(result): face area.
/// - std::get<2>(result): unit face normal.
///
/// @note The function currently assumes a quadrilateral face (faceToCoord.size() == 4).
/// @note The orientation of the returned normal depends on the ordering of the face vertices.
std::tuple<Dune::FieldVector<double,3>, double, Dune::FieldVector<double,3>>
computeFaceCenterAreaNormal(const std::vector<Dune::FieldVector<double,3>>& faceToCoord);

/// @brief Computes the intersection of two 3D line segments.
///
/// Determines whether two line segments intersect and returns the
/// intersection geometry if it exists.
///
/// The function handles the following cases:
/// - Non-parallel segments intersecting at a single point.
/// - Parallel but disjoint segments.
/// - Colinear segments with a single-point contact.
/// - Colinear segments with an overlapping segment.
/// - Skew (non-coplanar) segments.
///
/// The returned pair represents the intersection geometry:
/// - (p, p) for a single-point intersection.
/// - (p0, p1) for an overlapping segment.
/// - std::nullopt if the segments do not intersect.
///
/// @param [in] startA   Start point of segment A.
/// @param [in] endA     End point of segment A.
/// @param [in] startB   Start point of segment B.
/// @param [in] endB     End point of segment B.
/// @param[out] isInteriorInA
///     Set to true if the intersection lies in the interior of
///     segment A (i.e. not only at an endpoint). For overlapping
///     segments, indicates whether part of the overlap is interior
///     to segment A.
/// @param[out] isInteriorInB
///     Set to true if the intersection lies in the interior of
///     segment B (i.e. not only at an endpoint). For overlapping
///     segments, indicates whether part of the overlap is interior
///     to segment B.
///
/// @return
/// - std::nullopt if the segments do not intersect.
/// - std::pair{p, p} if the intersection consists of a single point.
/// - std::pair{p0, p1} if the intersection is an overlapping segment
///   with endpoints p0 and p1.
///
/// @note Segments are considered parallel or coplanar using a tolerance
///       of approximately 1e-8.
/// @note Degenerate segments (near-zero length) are treated as
///       non-intersecting.
std::optional<std::pair<Dune::FieldVector<double,3>, Dune::FieldVector<double,3>>>
computeSegmentIntersection(const Dune::FieldVector<double,3>& startA, const Dune::FieldVector<double,3>& endA,
                           const Dune::FieldVector<double,3>& startB, const Dune::FieldVector<double,3>& endB,
                           bool& isInteriorInA,
                           bool& isInteriorInB);

/// @brief Creates the edge list of a polygon from its ordered vertex indices.
///
/// Given a sequence of point/vertex indices describing a face boundary,
/// this function generates one edge for each consecutive pair of vertices.
/// The last vertex is automatically connected back to the first vertex,
/// producing a closed polygonal loop.
///
/// @param faceToPoint Container of vertex indices defining the face in
///                    boundary order. Must support size() and operator[].
///
/// @return A vector of edges, where each edge is represented as
///         std::array<int, 2>{startVertexIdx, endVertexIdx}.
std::vector<std::array<int,2>> createEdges(const auto& faceToPoint);

/// @brief Determines whether a point lies in the half-plane defined by a face edge.
///
/// The half-plane is bounded by the line passing through faceVertex in the
/// direction directionEdge and oriented according to faceNormal. The point
/// must lie on the face plane and on the non-negative side of the boundary
/// line to be considered inside the half-plane.
///
/// @param [in] vertex        Point to test.
/// @param [in] faceVertex    A point on the boundary line.
/// @param [in] faceNormal    Normal vector of the face plane.
/// @param [in] directionEdge Direction vector of the boundary line.
///
/// @return true if vertex lies in the half-plane; false otherwise.
bool inSemiplane(const Dune::FieldVector<double,3>& vertex,
                 const Dune::FieldVector<double,3>& faceVertex,
                 const Dune::FieldVector<double,3>& faceNormal,
                 const Dune::FieldVector<double,3>& directionEdge);

/// @brief Determines whether a point lies inside or on the boundary of a face.
///
/// The test is performed by constructing a half-plane for each face edge
/// using the face normal and verifying that the point lies within all
/// corresponding half-planes.
/// Assumes that the face is planar and convex. For more general faces,
/// an alternative approach may be required, for example by mapping to a
/// reference element, performing the computation in reference
/// coordinates, and transforming the result back to global coordinates.
///
/// @param [in] vertex
///        Point to test. Typically a Dune::FieldVector<double, 3>.
///
/// @param [in] face_gridData
///        CpGridData instance to which the face belongs.
///
/// @param [in] face_edges
///        Vector of corner-index pairs defining the edges of the face.
///        The indices refer to corners stored in face_gridData.
///
/// @param [in] face_normal
///        Normal vector of the face. Used to construct the edge half-planes.
///
/// @return true if the point lies inside the face or on its boundary; false otherwise.
bool isVertexInsideFace(const Dune::FieldVector<double,3>& vertex,
                        const Dune::cpgrid::CpGridData& face_gridData,
                        const std::vector<std::array<int,2>>& face_edges,
                        const Dune::FieldVector<double,3>& face_normal);

bool isVertexInSegmentInterior(const Dune::FieldVector<double,3>& vertex,
                               const Dune::FieldVector<double,3>& startSegment,
                               const Dune::FieldVector<double,3>& endSegment);

// face indices where the vertex appear interior of an edge
std::vector<int> isVertexInElementFaceEdge(const Dune::FieldVector<double,3>& vertex,
                                           int elemIdx,
                                           const Dune::cpgrid::CpGridData& gridData);

// edges form only by the 8 vertices of the element
std::vector<int> isVertexInElementEdge(const Dune::FieldVector<double,3>& vertex,
                                       int elemIdx,
                                       const Dune::cpgrid::CpGridData& gridData);


/// @brief Orders the vertices of a (quadrilateral) face according to its face tag.
///
/// Returns the four face vertices in a "canonical" ordering that depends on
/// the face tag (I_FACE, J_FACE, or K_FACE). The ordering is chosen to be
/// consistent with the logical grid directions associated with the face.
///
/// @param [in] faceTag      Face tag (I_FACE, J_FACE, or K_FACE).
/// @param [in] faceVertices Set containing the four face vertices.
///
/// @return The face vertices ordered according to the convention of the
///         specified face tag. 
///
/// @note The function currently assumes a quadrilateral face and requires
///       exactly four vertices.
/// @note The returned ordering differs for I-, J-, and K-faces.
std::vector<Dune::FieldVector<double,3>>
orderFaceVertices(int faceTag,
                  const std::set<Dune::FieldVector<double,3>,FieldVectorLess>& faceVertices);

/// @brief Computes the vertices defining the overlap polygon between two
/// geometrically coincident faces.
///
/// The function identifies existing and newly created vertices that
/// contribute to the overlap region between face1 and face2, including
/// vertices generated by edge-edge intersections. It also detects the
/// special case where face2 is fully contained in face1.
///
/// Only faces with matching face tags and orientations are considered.
/// If no valid overlap polygon can be constructed, std::nullopt is
/// returned.
///
/// @param [in] face1
/// @param [in] face2
/// @param [in] face1_gridData
///        CpGridData instance to which face1 belongs.
///
/// @param [in] face2_gridData
///        CpGridData instance to which face2 belongs.
///
/// @param [out] foundNewVertices
///        Set used to collect newly created vertices that do not already exist
///        in either face1_gridData or face2_gridData. Such vertices may be
///        generated, for example, when two edges intersect at an interior point.
///
/// @param [out] face1ExistingVertex_to_face1GridCornerIdx
///        Maps the coordinates of an existing vertex in face1 to its corner
///        index in face1_gridData.
///
/// @param [out] face2ExistingVertex_to_face2GridCornerIdx
///        Maps the coordinates of an existing vertex in face2 to its corner
///        index in face2_gridData.
///
/// @param [out] face2FullyContainedInFace1
///        Set to true if all vertices of face2 lie on face1.
///
/// @return
///        An optional containing the coordinates of the vertices that delimit
///        the overlap region between the two faces. Returns std::nullopt if
///        no overlap exists.
template<typename Face>
std::optional<std::vector<Dune::FieldVector<double,3>>>
computeFaceOverlapVertices(const Face&                                                  face1,
                           const Dune::cpgrid::CpGridData&                              face1_gridData,
                           std::map<Dune::FieldVector<double,3>, int, FieldVectorLess>& face1ExistingVertex_to_face1GridCornerIdx,
                           const Face&                                                  face2,
                           const Dune::cpgrid::CpGridData&                              face2_gridData,
                           std::map<Dune::FieldVector<double,3>, int, FieldVectorLess>& face2ExistingVertex_to_face2GridCornerIdx,
                           std::set<Dune::FieldVector<double,3>,FieldVectorLess>&       foundNewVertices,
                           bool&                                                        face2FullyContainedInFace1)
{
    const auto face1_tag = face1_gridData.faceTag(face1.index());
    const auto face2_tag = face2_gridData.faceTag(face2.index());
    
    // The orientation may differ when the grid data objects correspond
    // to neighboring single-cell refinements rather than a parent grid
    // and its associated single-cell refinement. Therefore, we only verify
    // that the face tags are identical, regardless of their orientation.
    if (face1_tag != face2_tag)
        return std::nullopt;
    
    const auto& face1_to_point = face1_gridData.faceToPoint(face1.index());
    const auto& face2_to_point = face2_gridData.faceToPoint(face2.index());

    // Used to create semi-planes to determine whether a vertex of face2
    // lies inside face1.
    const auto face1_edges = createEdges(face1_to_point);
    const auto& face1_normal = face1_gridData.faceNormals(face1.index());
    
    std::set<Dune::FieldVector<double,3>,FieldVectorLess> overlapFaceVertices{};
    
    // Case 1: face2 is fully contained in face1
    for (std::size_t i = 0; i < face2_to_point.size(); ++i) {
        const auto currentPoint = Dune::cpgrid::Entity<3>(face2_gridData, face2_to_point[i], true).geometry().center();
        if (isVertexInsideFace(currentPoint, face1_gridData, face1_edges, face1_normal)) {
            overlapFaceVertices.insert(currentPoint);
            face2ExistingVertex_to_face2GridCornerIdx[currentPoint] = face2_to_point[i];
        }
    }
    if (overlapFaceVertices.size() == face2_to_point.size()) {
        face2FullyContainedInFace1 = true;
    }

    // Case 2: Intersection/overlap area between face1 and face2 is not trivial (new vertices!) 
    for (const auto& face1_edge : face1_edges) {
        
        const auto edge0 = Dune::cpgrid::Entity<3>(face1_gridData, face1_edge[0], true).geometry().center();
        const auto edge1 = Dune::cpgrid::Entity<3>(face1_gridData, face1_edge[1], true).geometry().center();
        
        for (std::size_t i = 0; i < face2_to_point.size(); ++i) {

            const auto currentPoint = Dune::cpgrid::Entity<3>(face2_gridData, face2_to_point[i], true).geometry().center();
            const auto previousPoint = Dune::cpgrid::Entity<3>(face2_gridData, face2_to_point[(i-1)%face2_to_point.size()], true).geometry().center();

            bool isInteriorInFace1Edge{};
            bool isInteriorInFace2Edge{};
            
            const auto segmentInter = computeSegmentIntersection(edge0, edge1,
                                                                 previousPoint, currentPoint,
                                                                 isInteriorInFace1Edge, isInteriorInFace2Edge);
            
            bool currentPointInSemiplaneEdge = inSemiplane(currentPoint, edge0, face1_normal, edge1-edge0);
            bool previousPointInSemiplaneEdge = inSemiplane(previousPoint, edge0, face1_normal, edge1-edge0);
            
            if (currentPointInSemiplaneEdge) {
                if (!previousPointInSemiplaneEdge) {
                    if (segmentInter.has_value()) {
                        const auto& [p0,p1] = segmentInter.value();
                        if (isInteriorInFace2Edge) {
                            foundNewVertices.insert(p0); 
                            foundNewVertices.insert(p1); // if p0 and p1 coincide, no problem->it's a set
                            if (!isInteriorInFace1Edge) {
                                if (Dune::dot(p0-edge0,p0-edge0) < 1e-8) {
                                    face1ExistingVertex_to_face1GridCornerIdx[edge0] = face1_edge[0];
                                }
                                else if (Dune::dot(p0-edge1,p0-edge1) < 1e-8) {
                                    face1ExistingVertex_to_face1GridCornerIdx[edge1] = face1_edge[1];
                                }  
                            }
                        }
                        overlapFaceVertices.insert(p0);
                        overlapFaceVertices.insert(p1);
                        face2FullyContainedInFace1 = false;
                    }
                    overlapFaceVertices.insert(currentPoint);
                    face2ExistingVertex_to_face2GridCornerIdx[currentPoint] = face2_to_point[i];
                    face2FullyContainedInFace1 = false;
                }
            }
            else if (previousPointInSemiplaneEdge){
                if (segmentInter.has_value()) {
                    const auto& [p0,p1] = segmentInter.value();
                    if (isInteriorInFace2Edge) {
                        foundNewVertices.insert(p0);
                        foundNewVertices.insert(p1);
                    }
                    overlapFaceVertices.insert(p0);
                    overlapFaceVertices.insert(p1);
                }
                overlapFaceVertices.insert(previousPoint);
                face2ExistingVertex_to_face2GridCornerIdx[previousPoint] = face2_to_point[(i-1)%face2_to_point.size()];
                face2FullyContainedInFace1 = false;
            }
        } // end-for-refinedFaceToPoint-loop
    } // end-for-edges-loop
    if (overlapFaceVertices.size() != 4) { // faces might share only an edge, which does not create a new face
        return std::nullopt;
    }
    // Make sure all points forming the new face are actually ("inside") lying on face1
    bool allPointsInsideFace1 = true;
    for (const auto& vertex : overlapFaceVertices) {
        allPointsInsideFace1 = allPointsInsideFace1 && isVertexInsideFace(vertex, face1_gridData, face1_edges, face1_normal);
        if (!allPointsInsideFace1) {
            return std::nullopt;
        }
    }
    const auto ordered = orderFaceVertices(face1_tag, overlapFaceVertices);
    return std::make_optional<std::vector<Dune::FieldVector<double,3>>>(ordered);
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
                   std::vector<int>&                                                                  fullyContainedInParentFace);


} // namespace Lgr
} // namespace Opm


#endif // OPM_GRID_CPGRID_LGRFAUTLHELPERS_HEADER_INCLUDED
