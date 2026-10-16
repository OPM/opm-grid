//===========================================================================
//
// File: Geometry.hpp
//
// Created: Fri May 29 23:29:24 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            B�rd Skaflestad     <bard.skaflestad@sintef.no>
//            Antonella Ritorto   <antonella.ritorto@opm-op.com>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2009, 2010, 2011 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010, 2011, 2022 Equinor ASA.

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

#ifndef OPM_GEOMETRY_REFINEMENT_HEADER
#define OPM_GEOMETRY_REFINEMENT_HEADER

// Definitions of Geometry<3, cdim>'s refinement members; included by Geometry.hpp.

namespace Dune
{
    namespace cpgrid
    {

template <int cdim>
void Geometry<3, cdim>::refineCellifiedPatch(const std::array<int,3>& cells_per_dim,
DefaultGeometryPolicy& all_geom,
std::vector<std::array<int,8>>&  refined_cell_to_point,
cpgrid::OrientedEntityTable<0,1>& refined_cell_to_face,
Opm::SparseTable<int>& refined_face_to_point,
cpgrid::OrientedEntityTable<1,0>& refined_face_to_cell,
cpgrid::EntityVariable<enum face_tag, 1>& refined_face_tags,
cpgrid::SignedEntityVariable<typename Geometry<3, cdim>::PointType, 1>& refined_face_normals,
const std::array<int,3>& patch_dim,
const std::vector<double>& widthsX,
const std::vector<double>& lengthsY,
const std::vector<double>& heightsZ) const
            {
                EntityVariableBase<cpgrid::Geometry<0,3>>& refined_corners =
                    all_geom.geomVector(std::integral_constant<int,3>());
                EntityVariableBase<cpgrid::Geometry<2,3>>& refined_faces =
                    all_geom.geomVector(std::integral_constant<int,1>());
                EntityVariableBase<cpgrid::Geometry<3,3>>& refined_cells =
                    all_geom.geomVector(std::integral_constant<int,0>());
                EntityVariableBase<enum face_tag>& mutable_face_tags = refined_face_tags;
                EntityVariableBase<PointType>& mutable_face_normals = refined_face_normals;

                /// --- REFINED CORNERS ---
                // The strategy is to compute the local refined corners
                // of the unit/reference cube, and then apply the map global().
                // Determine the size of the vector containing all the corners
                // of all the global refined cells (children cells).
                // For easier notation:
                const std::array<int,3>& refined_dim = { cells_per_dim[0]*patch_dim[0],
                                                         cells_per_dim[1]*patch_dim[1],
                                                         cells_per_dim[2]*patch_dim[2]};
                refined_corners.resize((refined_dim[0] + 1)*(refined_dim[1] + 1)*(refined_dim[2] + 1));
                // The nummbering starts at the botton, so k=0 (z-axis), and j=0 (y-axis), i=0 (x-axis).
                // Then, increasing k ('going up'), followed by increasing i ('going right->'),
                // and finally, increasing j ('going back'). This order criteria for corners
                // 'Up [increasing k]- Right [incresing i]- Back [increasing j]'
                // is consistant with cpgrid numbering.
                //
                assert(static_cast<int>(widthsX.size()) == patch_dim[0]);
                assert(static_cast<int>(lengthsY.size()) == patch_dim[1]);
                assert(static_cast<int>(heightsZ.size()) == patch_dim[2]);
                const auto localCoordNumerator = []( const std::vector<double>& vec, int sumLimit, double multiplier) {
                    double lcn = 0;
                    assert(!vec.empty());
                    assert(sumLimit < static_cast<int>(vec.size()));
                    lcn += multiplier*vec[sumLimit];
                    for (int m = 0; m < sumLimit; ++m) {
                        lcn += vec[m];
                    }
                    return lcn;
                };
                // E.g. localCoordNumerator( dx, 3, 0.25) =  x0 + x1 + x2 + 0.25.x3
                //
                const double sumWidths = std::accumulate(widthsX.begin(), widthsX.end(), double(0));
                // x0 + x1 + ... + xL, if dx = {x0, x1, ..., xL}
                const double sumLengths = std::accumulate(lengthsY.begin(), lengthsY.end(), double(0));
                // y0 + y1 + ... + yM, if dy = {y0, y1, ..., yM}
                const double sumHeights = std::accumulate(heightsZ.begin(), heightsZ.end(), double(0));
                // z0 + z1 + ... + zN, if dz = {z0, z1, ..., zN}

                for (int j = 0; j < refined_dim[1] +1; ++j) {
                    double local_y = 0;
                    for (int i = 0; i < refined_dim[0] +1; ++i) {
                        double local_x = 0.;
                        for (int k = 0; k < refined_dim[2] +1; ++k) {
                            double local_z = 0.;

                            // Compute the index of each global refined corner associated with 'jik'.
                            int refined_corner_idx =
                                (j*(refined_dim[0]+1)*(refined_dim[2]+1)) + (i*(refined_dim[2]+1)) + k;

                            // Compute the local refined corner of the unit/reference cube associated with 'jik'.
                            if ( i == refined_dim[0]) { // last corner in the x-direction
                                local_x = sumWidths;
                            } else {
                                local_x = localCoordNumerator(widthsX, i/cells_per_dim[0], double((i % cells_per_dim[0])) / cells_per_dim[0]);
                            }
                            if ( j == refined_dim[1]) { // last corner in the y-direction
                                local_y = sumLengths;
                            } else {
                                local_y = localCoordNumerator(lengthsY, j/cells_per_dim[1], double((j % cells_per_dim[1])) / cells_per_dim[1]);
                            }
                            if ( k == refined_dim[2]) { // last corner in the z-direction
                                local_z = sumHeights;
                            } else {
                                local_z = localCoordNumerator(heightsZ, k/cells_per_dim[2], double((k % cells_per_dim[2])) /cells_per_dim[2]);
                            }

                            const LocalCoordinate& local_refined_corner = { local_x/sumWidths, local_y/sumLengths, local_z/sumHeights };
                            assert(local_x/sumWidths <= 1.);
                            assert(local_y/sumLengths <= 1.);
                            assert(local_z/sumHeights <= 1.);

                            // Compute the global refined corner 'jik' and add it in its corresponfing entry in "refined_corners".
                            refined_corners[refined_corner_idx] = Geometry<0, 3>(this->global(local_refined_corner));
                        } // end k-for-loop
                    } // end i-for-loop
                } // end j-for-loop
                /// --- END REFINED CORNERS ---
                //
                /// --- REFINED FACES ---
                // We want to populate "refined_faces". The size of "refined_faces" is
                const int refined_faces_size =
                    (refined_dim[0]*refined_dim[1]*(refined_dim[2]+1)) // 'bottom/top faces'
                    + ((refined_dim[0]+1)*refined_dim[1]*refined_dim[2]) // 'left/right faces'
                    + (refined_dim[0]*(refined_dim[1]+1)*refined_dim[2]); // 'front/back faces'
                refined_faces.resize(refined_faces_size);
                refined_face_tags.resize(refined_faces_size);
                refined_face_normals.resize(refined_faces_size);
                //
                // To create a face as a Geometry<2,3> type object we need its CENTROID and its VOLUME(area).
                // We store the centroids/areas  in the following order:
                // - Bottom-top faces -> 3rd coordinate constant in each face.
                // - Left-right faces -> 1st coordinate constant in each face.
                // - Front-back faces -> 2nd coordinate constant in each face.
                //
                // REFINED FACE AREAS
                // To compute the area of each face, we divide it in 4 triangles,
                // compute the area of those with "area()", where the arguments
                // are the 3 corners of each triangle. Then, sum them up to get the area
                // of the global refined face.
                // For each face, we construct 4 triangles with
                // (1) centroid of the face,
                // (2) one of the edges of the face.
                //
                // A triangle has 3 edges. Once we choose a face to base a triangle on,
                // we choose an edge of that face as one of the edges of the triangle.
                // The other 2 edges are fixed, since the centroid of the face the triangle
                // is based on is one of its corners. That's why to identify
                // a triangle we only need two things:
                // (1) the face it's based on and
                // (2) one of the four edges of that face.
                //
                // For each face, we need
                // 1. index of the refined face
                //    [needed to access indices of the 4 edges of the face in "refined_face_to_edges"]
                // 2. centroid of the face (common corner of the 4 triangles based on that face).
                //    [available via "['face'].center()"
                // 3. container of 4 entries (the 4 edges of the face).
                //    Each entry consists in the 2 indices defining each edge of the face.
                //    [available in "refined_face_to_edges"].
                //
                // Populate "mutable_face_tags/normals", "refined_face_to_point/cell",
                // "refined_faces".
                //
                for (int constant_direction = 0; constant_direction < 3; ++constant_direction){
                    // adding %3 and constant_direction, we go through the 3 type of faces.
                    // 0 -> 3rd coordinate constant: l('k') < cells_per_dim[2]+1, m('j') < cells_per_dim[1], n('i') < cells_per_dim[0]
                    // 1 -> 1rt coordinate constant: l('i') < cells_per_dim[0]+1, m('k') < cells_per_dim[2], n('j') < cells_per_dim[1]
                    // 2 -> 2nd coordinate constant: l('j') < cells_per_dim[1]+1, m('i') < cells_per_dim[0], n('k') < cells_per_dim[2]
                    std::array<int,3> refined_dim_mixed = {
                        refined_dim[(2+constant_direction)%3],
                        refined_dim[(1+constant_direction)%3],
                        refined_dim[constant_direction % 3] };
                    for (int l = 0; l < refined_dim_mixed[0] + 1; ++l) {
                        for (int m = 0; m < refined_dim_mixed[1]; ++m) {
                            for (int n = 0; n < refined_dim_mixed[2]; ++n) {
                                // Compute the face data.
                                auto [face_tag, idx, face_to_point, face_to_cell, wrong_local_centroid] =
                                    getIndicesFace(l, m, n, constant_direction, refined_dim);
                                // Add the tag to "refined_face_tags".
                                mutable_face_tags[idx]= face_tag;
                                // Add the 4 corners of the face to "refined_face_to_point".
                                refined_face_to_point.appendRow(face_to_point.begin(), face_to_point.end());
                                // Add the neighboring cells of the face to "refined_face_to_cell".
                                refined_face_to_cell.appendRow(face_to_cell.begin(), face_to_cell.end());
                                // Compute the centroid as the average of the 4 corners of the face
                                GlobalCoordinate face_center = { 0., 0., 0.};
                                for (int corn = 0; corn < 4; ++corn){
                                    face_center += refined_corners[face_to_point[corn]].center();
                                }
                                face_center /= 4.;
                                // Construct global face normal(s) (only one 'needed') and add it to "mutable_face_normals"
                                // Construct two vectors in the face, e.g. difference of two conners with the centroid,
                                // then obtain an orthogonal vector to both of them. Finally, normalize.
                                // Auxuliary vectors on the face:
                                GlobalCoordinate face_vector0 = refined_corners[face_to_point[0]].center() - face_center;
                                GlobalCoordinate face_vector1 = refined_corners[face_to_point[1]].center() - face_center;
                                mutable_face_normals[idx] = {
                                    (face_vector0[1]*face_vector1[2]) -  (face_vector0[2]*face_vector1[1]),
                                    (face_vector0[2]*face_vector1[0]) -  (face_vector0[0]*face_vector1[2]),
                                    (face_vector0[0]*face_vector1[1]) -  (face_vector0[1]*face_vector1[0])};
                                mutable_face_normals[idx] /= mutable_face_normals[idx].two_norm();
                                if (face_tag == J_FACE) {
                                    mutable_face_normals[idx] *= -1;
                                }
                                // Construct "refined_face_to_edges"
                                // with the {edge_indix[0], edge_index[1]} for each edge of the refined face.
                                std::vector<std::array<int,2>> refined_face_to_edges = {
                                    { face_to_point[0], face_to_point[1] },
                                    { face_to_point[1], face_to_point[2] },
                                    { face_to_point[2], face_to_point[3] },
                                    { face_to_point[3], face_to_point[0] }
                                };
                                // Calculate the AREA of each face of a global refined face,
                                // by adding the 4 areas of the triangles partitioning each face.
                                double refined_face_area = 0.0;
                                for (int edge = 0; edge < 4; ++edge) {
                                    // Construction of each triangle on the current face with one
                                    // of its edges equal to "edge".
                                    Geometry<0,3>::GlobalCoordinate trian_corners[3] = {
                                        refined_corners[refined_face_to_edges[edge][0]].center(),
                                        refined_corners[refined_face_to_edges[edge][1]].center(),
                                        face_center };
                                    refined_face_area += std::fabs(area(trian_corners));
                                } // end edge-for-loop
                                //
                                //
                                // Construct the Geometry<2,3> of the global refined face.
                                refined_faces[idx] = Geometry<2,cdim>(face_center, refined_face_area);
                            } // end n-for-loop
                        } // end m-for-loop
                    } // end l-for-loop
                } // end r-for-loop
                /// --- END REFINED FACES ---
                //
                /// --- REFINED CELLS ---
                // We need to populate "refined_cells"
                // "refined_cells"'s size is cells_per_dim[0] * cells_per_dim[1] * cells_per_dim[2].
                // To build each global refined cell, we need
                // 1. its global refined CENTER
                // 2. its VOLUME
                // 3. all global refined corners [available in "refined_corners"]
                // 4. indices of its 8 corners.
                //
                refined_cells.resize(refined_dim[0] * refined_dim[1] * refined_dim[2]);
                // Vector to store, in each entry, the 8 indices of the 8 corners
                // of each global refined cell. Determine its size.
                refined_cell_to_point.resize(refined_dim[0] * refined_dim[1] * refined_dim[2]);
                // The numembering starts with index 0 for the refined cell with corners
                // {0,0,0}, ...,{1/cells_per_dim[0], 1/cells_per_dim[1], 1/cells_per_dim[2]},
                // then the indices grow first picking the cells in the x-axis (Right, i), then y-axis (Back, j), and
                // finally, z-axis (Up, k).
                //
                // CENTERS
                // REFINED CELL CENTERS
                // The strategy is to compute the centers of the refined local
                // unit/reference cube, and then apply the map global().
                //
                // VOLUMES OF THE REFINED CELLS
                // REMARK: Each global refined 'cell' is a hexahedron since it may not be cube-shaped
                // since its a 'deformation' of unit/reference cube. We may use 'hexahedron' to refer
                // to the global refined cell in the computation of its volume.
                //
                // The strategy is to construct 24 tetrahedra in each hexahedron.
                // Each tetrahedron is built with
                // (1) the center of the hexahedron,
                // (2) the middle point of the face the tetrahedron is based on, and
                // (3) one of the edges of the face mentioned in 2.
                // Each face 'supports' 4 tetrahedra, and we have 6 faces per hexahedron, which
                // gives us the 24 tetrahedra per 'cell' (hexahedron).
                //
                // To compute the volume of each tetrahedron, we use "simplex_volume()" with
                // the 6 corners of the tetrahedron as arguments. Summing up the 24 volumes,
                // we get the volumne of the hexahedorn (global refined 'cell').
                //
                // Sum of all the volumes of all the (children) global refined cells.
                double sum_all_refined_cell_volumes = 0.0;
                //
                // For each (global refined 'cell') hexahedron, to create 24 tetrahedra and their volumes,
                // we introduce
                // Vol1. "hexa_to_face" (needed to access face centroids).
                // Vol2. "hexa_face_centroids" (one of the 6 corners of all 4 tetrahedra based on that face).
                // Vol3.  the center of the global refined 'cell' (hexahedron)
                //       (common corner of the 24 tetrahedra).
                // Vol4. "tetra_edge_indices" indices of the 4x6 tetrahedra per 'cell',
                //        grouped by the face they are based on.
                // Then we construct and compute the volume of the 24 tetrahedra with mainly
                // "hexa_face_centroids" (Vol2.), global refined cell center (Vol3.), and "tetra_edge_indices" (Vol4.).
                //
                for (int k = 0; k < refined_dim[2]; ++k) {
                    for (int j = 0; j < refined_dim[1]; ++j) {
                        for (int i = 0; i < refined_dim[0]; ++i) {
                            // INDEX of the global refined cell associated with 'kji'.
                            int refined_cell_idx = (k*refined_dim[0]*refined_dim[1]) + (j*refined_dim[0]) +i;
                            // Obtain the global refined center with 'this->global(local_refined_cell_center)'.
                            // 2. VOLUME of the global refined 'kji' cell
                            double refined_cell_volume = 0.0; // (computed below!)
                            // 3. All Global refined corners ("refined_corners")
                            // 4. Indices of the 8 corners of the global refined cell associated with 'kji'.
                            std::array<int,8> cell_to_point = { //
                                (j*(refined_dim[0]+1)*(refined_dim[2]+1))     + (i*(refined_dim[2]+1))      +k, // fake '0' {0,0,0}
                                (j*(refined_dim[0]+1)*(refined_dim[2]+1))     + ((i+1)*(refined_dim[2]+1))  +k, // fake '1' {1,0,0}
                                ((j+1)*(refined_dim[0]+1)*(refined_dim[2]+1)) + (i*(refined_dim[2]+1))      +k, // fake '2' {0,1,0}
                                ((j+1)*(refined_dim[0]+1)*(refined_dim[2]+1)) + ((i+1)*(refined_dim[2]+1))  +k, // fake '3' {1,1,0}
                                (j*(refined_dim[0]+1)*(refined_dim[2]+1))     + (i*(refined_dim[2]+1))      +k+1, // fake '4' {0,0,1}
                                (j*(refined_dim[0]+1)*(refined_dim[2]+1))     + ((i+1)*(refined_dim[2]+1))  +k+1, // fake '5' {1,0,1}
                                ((j+1)*(refined_dim[0]+1)*(refined_dim[2]+1)) + (i*(refined_dim[2]+1))      +k+1, // fake '6' {0,1,1}
                                ((j+1)*(refined_dim[0]+1)*(refined_dim[2]+1)) + ((i+1)*(refined_dim[2]+1))  +k+1 // fake '7' {1,1,1}
                            };
                            // Add this 8 corners to the corresponding entry of "refined_cell_to_point".
                            refined_cell_to_point[refined_cell_idx] = cell_to_point;
                            // 1. CENTER of the global refined cell associated with 'kji' (Vol3.)
                            // Compute the center of the local refined unit/reference cube associated with 'kji'.
                            GlobalCoordinate refined_cell_center = {0., 0., 0.};
                            for (int corn = 0; corn < 8; ++corn) {
                                refined_cell_center += refined_corners[cell_to_point[corn]].center();
                            }
                            refined_cell_center /= 8.;
                            //
                            // VOLUME HEXAHEDRON (GLOBAL REFINED 'CELL')
                            // Vol1. INDICES ('from 0 to 5') of the faces of the hexahedron (needed to access face centroids).
                            std::vector<int> hexa_to_face = { //hexa_face_0to5_indices = {
                                // index face '0' bottom
                                (k*refined_dim[0]*refined_dim[1]) + (j*refined_dim[0]) + i,
                                // index face '1' front
                                (refined_dim[0]*refined_dim[1]*(refined_dim[2]+1))
                                +  ((refined_dim[0]+1)*refined_dim[1]*refined_dim[2])
                                + (j*refined_dim[0]*refined_dim[2]) + (i*refined_dim[2]) + k,
                                // index face '2' left
                                (refined_dim[0]*refined_dim[1]*(refined_dim[2]+1))
                                + (i*refined_dim[1]*refined_dim[2]) + (k*refined_dim[1]) + j,
                                // index face '3' right
                                (refined_dim[0]*refined_dim[1]*(refined_dim[2]+1))
                                + ((i+1)*refined_dim[1]*refined_dim[2]) + (k*refined_dim[1]) + j,
                                // index face '4' back
                                (refined_dim[0]*refined_dim[1]*(refined_dim[2]+1)) +
                                ((refined_dim[0]+1)*refined_dim[1]*refined_dim[2])
                                + ((j+1)*refined_dim[0]*refined_dim[2]) + (i*refined_dim[2]) + k,
                                // index face '5' top
                                ((k+1)*refined_dim[0]*refined_dim[1]) + (j*refined_dim[0]) + i};
                            //
                            //  We add the 6 faces of the cell into "refined_cell_to_face".
                            using cpgrid::EntityRep;
                            // First value is index, Second is orientation.
                            // Still have to find out what the orientation should be.
                            // right face ('3') outer normal points 'from left to right' -> orientation true
                            // back face ('4') outer normal points 'from front to back' -> orientation true
                            // top face ('5') outer normal points 'from bottom to top' -> orientation true
                            // (the other cases are false)
                            std::vector<cpgrid::EntityRep<1>> cell_to_face = {
                                { hexa_to_face[0], false}, {hexa_to_face[1], false},
                                { hexa_to_face[2], false}, {hexa_to_face[3], true},
                                { hexa_to_face[4], true},  {hexa_to_face[5], true} };
                            refined_cell_to_face.appendRow(cell_to_face.begin(), cell_to_face.end());
                            //
                            // Vol2. CENTROIDS of the faces of the hexahedron.
                            // (one of the 6 corners of all 4 tetrahedra based on that face).
                            std::vector<Geometry<0,3>::GlobalCoordinate> hexa_face_centroids;
                            for (auto& idx : hexa_to_face) {
                                hexa_face_centroids.push_back(refined_faces[idx].center());
                            }
                            // Indices of the 4 edges of each face of the hexahedron.
                            // A tetrahedron has six edges. Once we choose a face to base a
                            // tetrahedron on, we choose an edge of that face as one of the
                            // edges of the tetrahedron. The other five edges are fixed, since
                            // the center of the hexahedron and the center of the face are
                            // the reminder 2 corners of the tetrahedron. That's why to identify
                            // a tetrahedron we only need two things:
                            // (1) the face it's based on and
                            // (2) one of the four edges of that face.
                            //
                            // Container with 6 entries, one per face. Each entry has the
                            // 4 indices of the 4 corners of each face.
                            std::vector<std::array<int,4>> cell_face4corners;
                            cell_face4corners.reserve(6);
                            for (int face = 0; face < 6;  ++face) {
                                cell_face4corners.push_back({
                                        refined_face_to_point[hexa_to_face[face]][0],
                                        refined_face_to_point[hexa_to_face[face]][1],
                                        refined_face_to_point[hexa_to_face[face]][2],
                                        refined_face_to_point[hexa_to_face[face]][3] });
                            }
                            // Vol4. Container with indices of the edges of the 4 tetrahedra per face
                            // [according to description above]
                            std::vector<std::vector<std::array<int,2>>> tetra_edge_indices;
                            tetra_edge_indices.reserve(6);
                            for (auto& face_indices : cell_face4corners)
                            {
                                std::vector<std::array<int,2>> face4edges_indices = {
                                    { face_indices[0], face_indices[1]}, // fake '{0,1}'/'{4,5}'
                                    { face_indices[0], face_indices[2]}, // fake '{0,2}'/'{4,6}'
                                    { face_indices[1], face_indices[3]}, // fake '{1,3}'/'{5,7}'
                                    { face_indices[2], face_indices[3]} }; // fake '{2,3}'/'{6,7}'
                                tetra_edge_indices.push_back(face4edges_indices);
                            }
                            // Sum of the 24 volumes to get the volume of the hexahedron,
                            // stored in "refined_cell_volume".
                            // Calculate the volume of each hexahedron, by adding
                            // the 4 tetrahedra at each face (4x6 = 24 tetrahedra).
                            for (int face = 0; face < 6; ++face) {
                                for (int edge = 0; edge < 4; ++edge) {
                                    // Construction of each tetrahedron based on "face" with one
                                    // of its edges equal to "edge".
                                    const Geometry<0, 3>::GlobalCoordinate tetra_corners[4] = {
                                        refined_corners[tetra_edge_indices[face][edge][0]].center(),  // (see Vol4.)
                                        refined_corners[tetra_edge_indices[face][edge][1]].center(),  // (see Vol4.)
                                        hexa_face_centroids[face],  // (see Vol2.)
                                        refined_cell_center };  // (see Vol3.)
                                    refined_cell_volume += std::fabs(simplex_volume(tetra_corners));
                                } // end edge-for-loop
                            } // end face-for-loop
                            // Add the volume of the hexahedron (global refined 'cell')
                            // to the container with of all volumes of all the refined cells.
                            sum_all_refined_cell_volumes += refined_cell_volume;
                            // Create a pointer to the first element of "refined_cell_to_point"
                            // (required as the fourth argement to construct a Geometry<3,3> type object).
                            int* indices_storage_ptr = refined_cell_to_point[refined_cell_idx].data();
                            // Construct the Geometry of the refined cell associated with 'kji'.
                            refined_cells[refined_cell_idx] =
                                Geometry<3,cdim>(refined_cell_center,
                                                 refined_cell_volume,
                                                 all_geom.pointGeometryPtr(),
                                                 indices_storage_ptr);
                        } // end i-for-loop
                    }  // end j-for-loop
                } // end k-for-loop
                // Rescale all volumes if the sum of volume of all the global refined 'cells' does not match the
                // volume of the 'parent cell'.
                // Compare the sum of all the volumes of all refined cells with 'parent cell' volume.
                if (std::fabs(sum_all_refined_cell_volumes - this->volume())
                    > std::numeric_limits<Geometry<3, cdim>::ctype>::epsilon()) {
                    Geometry<3, cdim>::ctype correction = this->volume() / sum_all_refined_cell_volumes;
                    for(auto& cell: refined_cells){
                        cell.vol_ *= correction;
                    }
                } // end if-statement
                /// --- END REFINED CELLS ---
            } /// --- END of refine(dx, dy, dz)

template <int cdim>
const std::tuple< enum face_tag, int,
std::array<int, 4>, std::vector<cpgrid::EntityRep<0>>,
typename Geometry<3, cdim>::LocalCoordinate>
Geometry<3, cdim>::getIndicesFace(int l, int m, int n, int constant_direction, const std::array<int, 3>& cells_per_dim) const
            {
                using cpgrid::EntityRep;
                std::vector<cpgrid::EntityRep<0>> neighboring_cells_of_one_face; // {index, orientation}
                switch(constant_direction) {
                case 0:  // {l,m,n} = {k,j,i}, constant in z-direction
                    // Orientation true when outer normal points 'from bottom to top'
                    // Orientation false when outer normal points 'from top to bottom'
                    if (l != 0) {
                        neighboring_cells_of_one_face.push_back({((l-1)*cells_per_dim[0]*cells_per_dim[1])
                                + (m*cells_per_dim[0]) + n, true});
                    }
                    if (l != cells_per_dim[2]) {
                        neighboring_cells_of_one_face.push_back({ (l*cells_per_dim[0]*cells_per_dim[1])
                                + (m*cells_per_dim[0]) + n, false});
                    }
                    return { face_tag::K_FACE,
                             (l*cells_per_dim[0]*cells_per_dim[1]) + (m*cells_per_dim[0]) + n,
                             {(m*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (n*(cells_per_dim[2]+1)) +l,
                              (m*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + ((n+1)*(cells_per_dim[2]+1)) +l,
                              ((m+1)*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + ((n+1)*(cells_per_dim[2]+1)) +l,
                               ((m+1)*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (n*(cells_per_dim[2]+1)) +l},
                             neighboring_cells_of_one_face,
                             {(.5 + n)/cells_per_dim[0], (.5 + m)/cells_per_dim[1], double(l)/cells_per_dim[2]}};
                case 1:  // {l,m,n} = {i,k,j}, constant in the x-direction
                    // Orientation true when outer normal points 'from left to right'
                    // Orientation false when outer normal points 'from right to left'
                    if (l != 0) {
                        neighboring_cells_of_one_face.push_back({(m*cells_per_dim[0]*cells_per_dim[1])
                                + (n*cells_per_dim[0]) +l-1, true});
                    }
                    if (l != cells_per_dim[0]) {
                        neighboring_cells_of_one_face.push_back({ (m*cells_per_dim[0]*cells_per_dim[1])
                                + (n*cells_per_dim[0]) + l, false});
                    }
                    return { face_tag::I_FACE,
                             (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1))
                        + (l*cells_per_dim[1]*cells_per_dim[2]) + (m*cells_per_dim[1]) + n,
                             {(n*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (l*(cells_per_dim[2]+1)) +m,
                              ((n+1)*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (l*(cells_per_dim[2]+1)) +m,
                              ((n+1)*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (l*(cells_per_dim[2]+1)) +m+1,
                              (n*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (l*(cells_per_dim[2]+1)) +m+1},
                             neighboring_cells_of_one_face,
                             { double(l)/cells_per_dim[0], (.5 + n)/cells_per_dim[1], (.5 + m)/cells_per_dim[2]}};
                case 2: // {l,m,n} = {j,i,k}, constant in the y-direction
                    // Orientation true when outer normal points 'from front to back'
                    // Orientation false when outer normal points 'from back to front'
                    if (l != 0) {
                        neighboring_cells_of_one_face.push_back({(n*cells_per_dim[0]*cells_per_dim[1])
                                + ((l-1)*cells_per_dim[0]) +m, true});
                    }
                    if (l != cells_per_dim[1]) {
                        neighboring_cells_of_one_face.push_back({(n*cells_per_dim[0]*cells_per_dim[1])
                                + (l*cells_per_dim[0]) + m, false});
                    }
                    return { face_tag::J_FACE,
                             (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2] +1))
                        + ((cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2])
                        + (l*cells_per_dim[0]*cells_per_dim[2]) + (m*cells_per_dim[2]) + n,
                             {(l*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (m*(cells_per_dim[2]+1)) +n,
                              (l*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + ((m+1)*(cells_per_dim[2]+1)) +n,
                              (l*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + ((m+1)*(cells_per_dim[2]+1)) +n+1,
                              (l*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (m*(cells_per_dim[2]+1)) +n+1},
                             neighboring_cells_of_one_face,
                             {(.5 + m)/cells_per_dim[0], double(l)/cells_per_dim[1], (.5 + n)/cells_per_dim[2]}};
                default:
                    // Should never be reached, but prevents compiler warning
                    OPM_THROW(std::logic_error, "Unhandled dimension. This should never happen!");
                }
            }

    } // namespace cpgrid
} // namespace Dune

#endif // OPM_GEOMETRY_REFINEMENT_HEADER
