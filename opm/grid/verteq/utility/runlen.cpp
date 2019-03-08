#include "runlen.hpp"
#include <opm/grid/UnstructuredGrid.h>
using namespace Opm;

rlw_int Opm::grid_cell_facetag (const UnstructuredGrid& g) {
    // tag for faces in cells
    return rlw_int (g.number_of_cells,
                    g.cell_facepos,
                    g.cell_facetag);
}

rlw_int Opm::grid_cell_faces (const UnstructuredGrid& g) {
    // id for faces in cells
    return rlw_int (g.number_of_cells,
                    g.cell_facepos,
                    g.cell_faces);
}
