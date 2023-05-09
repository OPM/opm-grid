#include <config.h>

#if HAVE_ECL_INPUT

#include <algorithm>
#include <cstddef>
#include <type_traits>

#include <opm/common/utility/numeric/cmp.hpp>

#include <opm/grid/UnstructuredGrid.h>

/*
   The grid_equal() function is separated out into a separate file to
   ensure that it is compiled with a C++ compiler, so that the
   namespace features used in the opm/common/util/numeric/cmp.cpp
   implementation compiles.
*/

template <typename T>
bool array_equal_helper(const T*          x,
                        const T*          y,
                        const std::size_t n,
                        std::true_type)
{
    return Opm::cmp::array_equal(x, y, n);
}

template <typename T>
bool array_equal_helper(const T*          x,
                        const T*          y,
                        const std::size_t n,
                        std::false_type)
{
    return std::equal(x, x + n, y);
}

template <typename T>
bool array_equal(const T* x, const T* y, const std::size_t n)
{
    using Float_P = typename std::is_floating_point<T>::type;

    return array_equal_helper(x, y, n, Float_P());
}

template <typename T>
bool equal_alloc_status(const T* x, const T* y)
{
    return (x == nullptr) == (y == nullptr);
}

template <typename T>
bool equal_vector(const T* x, const T* y, const std::size_t n)
{
    return equal_alloc_status(x, y) &&
           ((x == nullptr) || array_equal(x, y, n));
}

// ----------------------------------------------------------------------
bool grid_equal(const UnstructuredGrid* g1, const UnstructuredGrid* g2)
// ----------------------------------------------------------------------
{
    if ((g1 == nullptr) || (g2 == nullptr)) {
        // Grids must be allocated to qualify as equal.
        //
        // Note: We also return false if *neither* are allocated even
        // though nullptr==nullptr.

        return false;
    }

    bool        eq;
    std::size_t n;

    // Basic sanity check.  Array dimensions.
    {
        eq = ((g1->dimensions      == g2->dimensions)      &&
              (g1->number_of_nodes == g2->number_of_nodes) &&
              (g1->number_of_faces == g2->number_of_faces) &&
              (g1->number_of_cells == g2->number_of_cells));
    }

    // ============================================================

    // Topology checks.  Exact integer equality.
    {
        // 1) Face->node topology
        n  =       g1->number_of_faces + 1;
        eq = eq && equal_vector(g1->face_nodepos, g2->face_nodepos, n);

        n  =       g1->face_nodepos[ g1->number_of_faces ];
        eq = eq && equal_vector(g1->face_nodes, g2->face_nodes, n);

        // 2) Cell->face topology
        n  =       g1->number_of_cells + 1;
        eq = eq && equal_vector(g1->cell_facepos, g2->cell_facepos, n);

        n  =       g1->cell_facepos[ g1->number_of_cells ];
        eq = eq && equal_vector(g1->cell_faces, g2->cell_faces, n);

        // 3) Face->cell topology
        n  =       2 * g1->number_of_faces;
        eq = eq && equal_vector(g1->face_cells, g2->face_cells, n);

        // 4) Local-to-global cell map
        n  =       1 * g1->number_of_cells;
        eq = eq && equal_vector(g1->global_cell, g2->global_cell, n);

        // 5) Connection's cardinal directions per cell
        n  =       g1->cell_facepos[ g1->number_of_cells ];
        eq = eq && equal_vector(g1->cell_facetag, g2->cell_facetag, n);
    }

    // ============================================================

    // Geometry checks.  Approximate, floating point comparisons.
    {
        // 1) Node coordinates.
        n  =       g1->dimensions * g1->number_of_nodes;
        eq = eq && equal_vector(g1->node_coordinates,
                                g2->node_coordinates, n);

        // 2) Face centroids
        n  =       g1->dimensions * g1->number_of_faces;
        eq = eq && equal_vector(g1->face_centroids,
                                g2->face_centroids, n);

        // 3) Face normals
        n  =       g1->dimensions * g1->number_of_faces;
        eq = eq && equal_vector(g1->face_normals,
                                g2->face_normals, n);

        // 4) Face areas
        n  =       1 * g1->number_of_faces;
        eq = eq && equal_vector(g1->face_areas,
                                g2->face_areas, n);

        // 5) Cell centroids
        n  =       g1->dimensions * g1->number_of_cells;
        eq = eq && equal_vector(g1->cell_centroids,
                                g2->cell_centroids, n);

        // 6) Cell volumes
        n  =       1 * g1->number_of_cells;
        eq = eq && equal_vector(g1->cell_volumes,
                                g2->cell_volumes, n);
    }

    return eq;
}

#endif // #if HAVE_ECL_INPUT
