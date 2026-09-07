#include <config.h>

// Warning suppression for Dune includes.
#include <opm/grid/utility/platform_dependent/disable_warnings.h>

#include <dune/common/unused.hh>
#include <opm/grid/CpGrid.hpp>
#include <opm/grid/UnstructuredGrid.h>
#include <opm/grid/cornerpoint_grid.h>
#include <opm/grid/cpgrid/GridHelpers.hpp>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <opm/grid/cpgrid/dgfparser.hh>

#if HAVE_OPM_COMMON
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#endif

#define DISABLE_DEPRECATED_METHOD_CHECK 1
using Dune::referenceElement; //grid check assume usage of Dune::Geometry
#include <dune/grid/test/gridcheck.hh>


// Re-enable warnings.
#include <opm/grid/utility/platform_dependent/reenable_warnings.h>

#include <cstdio>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>

namespace {

void serializeUnstructuredGridToFile(const UnstructuredGrid& grid, const std::string& filename)
{
    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("Failed to open output file for serialized UnstructuredGrid: " + filename);
    }

    const int ndims = grid.dimensions;
    const int ncells = grid.number_of_cells;
    const int nfaces = grid.number_of_faces;
    const int nnodes = grid.number_of_nodes;
    const int nfacenodes = grid.face_nodepos[nfaces];
    const int ncellfaces = grid.cell_facepos[ncells];
    const int has_tag = (grid.cell_facetag != nullptr) ? 1 : 0;
    const int has_indexmap = (grid.global_cell != nullptr) ? 1 : 0;

    out << ndims << ' ' << ncells << ' ' << nfaces << ' ' << nnodes << ' '
        << nfacenodes << ' ' << ncellfaces << ' ' << has_tag << ' ' << has_indexmap << '\n';

    for (int d = 0; d < ndims; ++d) {
        out << grid.cartdims[d] << (d + 1 == ndims ? '\n' : ' ');
    }

    for (int i = 0; i < nnodes * ndims; ++i) {
        out << grid.node_coordinates[i] << (i + 1 == nnodes * ndims ? '\n' : ' ');
    }

    for (int i = 0; i < nfaces + 1; ++i) {
        out << grid.face_nodepos[i] << (i + 1 == nfaces + 1 ? '\n' : ' ');
    }
    for (int i = 0; i < nfacenodes; ++i) {
        out << grid.face_nodes[i] << (i + 1 == nfacenodes ? '\n' : ' ');
    }
    for (int i = 0; i < 2 * nfaces; ++i) {
        out << grid.face_cells[i] << (i + 1 == 2 * nfaces ? '\n' : ' ');
    }
    for (int i = 0; i < nfaces; ++i) {
        out << grid.face_areas[i] << (i + 1 == nfaces ? '\n' : ' ');
    }
    for (int i = 0; i < nfaces * ndims; ++i) {
        out << grid.face_centroids[i] << (i + 1 == nfaces * ndims ? '\n' : ' ');
    }
    for (int i = 0; i < nfaces * ndims; ++i) {
        out << grid.face_normals[i] << (i + 1 == nfaces * ndims ? '\n' : ' ');
    }

    for (int i = 0; i < ncells + 1; ++i) {
        out << grid.cell_facepos[i] << (i + 1 == ncells + 1 ? '\n' : ' ');
    }

    for (int i = 0; i < ncellfaces; ++i) {
        if (has_tag) {
            out << grid.cell_faces[i] << ' ' << grid.cell_facetag[i];
        } else {
            out << grid.cell_faces[i];
        }
        out << (i + 1 == ncellfaces ? '\n' : ' ');
    }

    if (has_indexmap) {
        for (int i = 0; i < ncells; ++i) {
            out << grid.global_cell[i] << (i + 1 == ncells ? '\n' : ' ');
        }
    }

    for (int i = 0; i < ncells; ++i) {
        out << grid.cell_volumes[i] << (i + 1 == ncells ? '\n' : ' ');
    }
    for (int i = 0; i < ncells * ndims; ++i) {
        out << grid.cell_centroids[i] << (i + 1 == ncells * ndims ? '\n' : ' ');
    }
}

void compareCpGridCellsByGlobalId(const Dune::CpGrid& lhs,
                                  const Dune::CpGrid& rhs,
                                  const double tol = 1e-12)
{
    // Compare top-level observable grid metadata.
    assert(lhs.logicalCartesianSize() == rhs.logicalCartesianSize());
    assert(lhs.numCells() == rhs.numCells());
    assert(lhs.numFaces() == rhs.numFaces());
    assert(lhs.numVertices() == rhs.numVertices());
    assert(lhs.numCellFaces() == rhs.numCellFaces());
    assert(lhs.globalCell() == rhs.globalCell());
    assert(lhs.uniqueBoundaryIds() == rhs.uniqueBoundaryIds());
    assert(lhs.numBoundarySegments() == rhs.numBoundarySegments());

    // NOTE: zcornData is not part of the serialized UnstructuredGrid file format,
    // so we cannot require exact equivalence for this member here.

    // Vertex geometry.
    for (int v = 0; v < lhs.numVertices(); ++v) {
        const auto& lp = lhs.vertexPosition(v);
        const auto& rp = rhs.vertexPosition(v);
        for (int d = 0; d < 3; ++d) {
            assert(std::abs(lp[d] - rp[d]) < tol);
        }
    }

    // Face topology and geometry.
    for (int f = 0; f < lhs.numFaces(); ++f) {
        assert(lhs.numFaceVertices(f) == rhs.numFaceVertices(f));
        for (int lv = 0; lv < lhs.numFaceVertices(f); ++lv) {
            assert(lhs.faceVertex(f, lv) == rhs.faceVertex(f, lv));
        }

        assert(lhs.faceCell(f, 0) == rhs.faceCell(f, 0));
        assert(lhs.faceCell(f, 1) == rhs.faceCell(f, 1));
        assert(std::abs(lhs.faceArea(f) - rhs.faceArea(f)) < tol);

        const auto& lfc = lhs.faceCentroid(f);
        const auto& rfc = rhs.faceCentroid(f);
        const auto& lfn = lhs.faceNormal(f);
        const auto& rfn = rhs.faceNormal(f);
        for (int d = 0; d < 3; ++d) {
            assert(std::abs(lfc[d] - rfc[d]) < tol);
            assert(std::abs(lfn[d] - rfn[d]) < tol);
        }
        assert(lhs.boundaryId(f) == rhs.boundaryId(f));
    }

    // Cell topology and geometry.
    for (int c = 0; c < lhs.numCells(); ++c) {
        assert(lhs.numCellFaces(c) == rhs.numCellFaces(c));
        for (int lf = 0; lf < lhs.numCellFaces(c); ++lf) {
            assert(lhs.cellFace(c, lf) == rhs.cellFace(c, lf));
        }

        assert(std::abs(lhs.cellVolume(c) - rhs.cellVolume(c)) < tol);
        assert(std::abs(lhs.cellCenterDepth(c) - rhs.cellCenterDepth(c)) < tol);

        const auto& lcc = lhs.cellCentroid(c);
        const auto& rcc = rhs.cellCentroid(c);
        for (int d = 0; d < 3; ++d) {
            assert(std::abs(lcc[d] - rcc[d]) < tol);
        }

        std::array<int, 3> lijk{};
        std::array<int, 3> rijk{};
        lhs.getIJK(c, lijk);
        rhs.getIJK(c, rijk);
        assert(lijk == rijk);

        const auto lecl = lhs.getEclCentroid(c);
        const auto recl = rhs.getEclCentroid(c);
        for (int d = 0; d < 3; ++d) {
            assert(std::abs(lecl[d] - recl[d]) < tol);
        }
    }

    // Compare Dune ids for codim 0 and 3 (supported codims).
    const auto& lhsGlobalIdSet = lhs.globalIdSet();
    const auto& rhsGlobalIdSet = rhs.globalIdSet();
    const auto lhsLeaf = lhs.leafGridView();
    const auto rhsLeaf = rhs.leafGridView();

    std::unordered_map<int, std::size_t> rhsCellIdByIndex;
    std::unordered_map<int, std::size_t> rhsVertexIdByIndex;
    rhsCellIdByIndex.reserve(static_cast<std::size_t>(rhs.numCells()));
    rhsVertexIdByIndex.reserve(static_cast<std::size_t>(rhs.numVertices()));

    for (const auto& elem : Dune::elements(rhsLeaf)) {
        rhsCellIdByIndex.emplace(elem.index(), rhsGlobalIdSet.id(elem));
    }
    for (const auto& vertex : Dune::vertices(rhsLeaf)) {
        rhsVertexIdByIndex.emplace(vertex.index(), rhsGlobalIdSet.id(vertex));
    }

    for (const auto& elem : Dune::elements(lhsLeaf)) {
        const auto it = rhsCellIdByIndex.find(elem.index());
        assert(it != rhsCellIdByIndex.end());
        assert(static_cast<std::int64_t>(lhsGlobalIdSet.id(elem)) == static_cast<std::int64_t>(it->second));
    }
    for (const auto& vertex : Dune::vertices(lhsLeaf)) {
        const auto it = rhsVertexIdByIndex.find(vertex.index());
        assert(it != rhsVertexIdByIndex.end());
        assert(static_cast<std::int64_t>(lhsGlobalIdSet.id(vertex)) == static_cast<std::int64_t>(it->second));
    }

    const auto& lhsGlobal = lhs.globalCell();
    const auto& rhsGlobal = rhs.globalCell();
    assert(lhsGlobal.size() == rhsGlobal.size());

    std::unordered_map<int, int> lhsG2L;
    std::unordered_map<int, int> rhsG2L;
    lhsG2L.reserve(lhsGlobal.size());
    rhsG2L.reserve(rhsGlobal.size());

    for (int i = 0; i < static_cast<int>(lhsGlobal.size()); ++i) {
        lhsG2L[lhsGlobal[i]] = i;
    }
    for (int i = 0; i < static_cast<int>(rhsGlobal.size()); ++i) {
        rhsG2L[rhsGlobal[i]] = i;
    }

    for (const auto& gidLocal : lhsG2L) {
        const int gid = gidLocal.first;
        const int li = gidLocal.second;
        const auto riIt = rhsG2L.find(gid);
        assert(riIt != rhsG2L.end());
        const int ri = riIt->second;

        assert(std::abs(lhs.cellVolume(li) - rhs.cellVolume(ri)) < tol);
        const auto& lc = lhs.cellCentroid(li);
        const auto& rc = rhs.cellCentroid(ri);
        for (int d = 0; d < 3; ++d) {
            assert(std::abs(lc[d] - rc[d]) < tol);
        }

        const auto lecl = lhs.getEclCentroid(li);
        const auto recl = rhs.getEclCentroid(ri);
        for (int d = 0; d < 3; ++d) {
            assert(std::abs(lecl[d] - recl[d]) < tol);
        }
    }
}

} // namespace

template <class GridView>
void testGridIteration( const GridView& gridView, const int nElem )
{
    typedef typename GridView::template Codim<0>::Iterator ElemIterator;
    typedef typename GridView::IntersectionIterator IsIt;
    typedef typename GridView::template Codim<0>::Geometry Geometry;

    int numElem = 0;
    ElemIterator elemIt = gridView.template begin<0>();
    ElemIterator elemEndIt = gridView.template end<0>();
    for (; elemIt != elemEndIt; ++elemIt) {
        const Geometry& elemGeom = elemIt->geometry();
        if (std::abs(elemGeom.volume() - 1.0) > 1e-8)
            std::cout << "element's " << numElem << " volume is wrong:"<<elemGeom.volume()<<"\n";

        typename Geometry::LocalCoordinate local( 0.5 );
        typename Geometry::GlobalCoordinate global = elemGeom.global( local );
        typename Geometry::GlobalCoordinate center = elemGeom.center();
        if( (center - global).two_norm() > 1e-6 )
        {
          std::cout << "center = " << center << " global( localCenter ) = " << global << std::endl;
        }


        int numIs = 0;
        IsIt isIt = gridView.ibegin(*elemIt);
        IsIt isEndIt = gridView.iend(*elemIt);
        for (; isIt != isEndIt; ++isIt, ++ numIs)
        {
            const auto& intersection = *isIt;
            const auto& isGeom = intersection.geometry();
            //std::cout << "Checking intersection id = " << localIdSet.id( intersection ) << std::endl;
            if (std::abs(isGeom.volume() - 1.0) > 1e-8)
                std::cout << "volume of intersection " << numIs << " of element " << numElem << " volume is wrong: " << isGeom.volume() << "\n";

            if (intersection.neighbor())
            {
              if( numIs != intersection.indexInInside() )
                  std::cout << "num iit = " << numIs << " indexInInside " << intersection.indexInInside() << std::endl;

              if (std::abs(intersection.outside().geometry().volume() - 1.0) > 1e-8)
                  std::cout << "outside element volume of intersection " << numIs << " of element " << numElem
                            << " volume is wrong: " << intersection.outside().geometry().volume() << std::endl;

              if (std::abs(intersection.inside().geometry().volume() - 1.0) > 1e-8)
                  std::cout << "inside element volume of intersection " << numIs << " of element " << numElem
                            << " volume is wrong: " << intersection.inside().geometry().volume() << std::endl;
            }
        }

        if (numIs != 2 * GridView::dimension )
            std::cout << "number of intersections is wrong for element " << numElem << "\n";

        ++ numElem;
    }

    if (numElem != nElem )
        std::cout << "number of elements is wrong: " << numElem << ", expected " << nElem << std::endl;
}


template <class Grid>
void testGrid(Grid& grid, const std::string& name, const size_t nElem, const size_t nVertices)
{
    typedef typename Grid::LeafGridView GridView;
    /*

    try {
      gridcheck( grid );
    }
    catch ( const Dune::Exception& e)
    {
      std::cerr << "Warning: " << e.what() << std::endl;
    }
*/
    std::cout << name << std::endl;

    testGridIteration( grid.leafGridView(), nElem );

    std::cout << "create vertex mapper\n";
  Dune::MultipleCodimMultipleGeomTypeMapper<GridView> mapper(grid.leafGridView(), Dune::mcmgVertexLayout());

    std::cout << "VertexMapper.size(): " << mapper.size() << "\n";
    if (static_cast<size_t>(mapper.size()) != nVertices ) {
        std::cout << "Wrong size of vertex mapper. Expected " << nVertices << "!" << std::endl;
        //std::abort();
    }

    // VTKWriter does not work with geometry type none at the moment
    if( true || grid.geomTypes( 0 )[ 0 ].isCube() )
    {
      std::cout << "create vtkWriter\n";
      typedef Dune::VTKWriter<GridView> VtkWriter;
      VtkWriter vtkWriter(grid.leafGridView());

      std::cout << "create cellData\n";
      int numElems = grid.size(0);
      std::vector<double> tmpData(numElems, 0.0);

      std::cout << "add cellData\n";
      vtkWriter.addCellData(tmpData, name);

      std::cout << "write data\n";
      vtkWriter.write(name, Dune::VTK::ascii);
    }

}

int main(int argc, char** argv )
{
    // initialize MPI
    Dune::MPIHelper::instance( argc, argv );

    // test CpGrid
    typedef Dune::CpGrid Grid;

#if HAVE_OPM_COMMON
    const char *deckString =
        "RUNSPEC\n"
        "METRIC\n"
        "DIMENS\n"
        "2 2 2 /\n"
        "GRID\n"
        "DXV\n"
        "2*1 /\n"
        "DYV\n"
        "2*1 /\n"
        "DZ\n"
        "8*1 /\n"
        "TOPS\n"
        "8*100.0 /\n";

    Opm::Parser parser;
    const auto deck = parser.parseString(deckString);

    Grid grid;
    Opm::EclipseGrid ecl_grid(deck);

    grid.processEclipseFormat(&ecl_grid, nullptr, false, false, false);
    testGrid( grid, "CpGrid_ecl", 8, 27 );

    const auto& grid_leafView = grid.leafGridView();
    Dune::CartesianIndexMapper<Grid> grid_cartMapper =  Dune::CartesianIndexMapper<Grid>(grid);
    for (const auto& element: Dune::elements(grid_leafView)){
        const auto& elemEclCentroid = ecl_grid.getCellCenter(grid_cartMapper.cartesianIndex(element.index()));
        const auto& elemCpGridEclCentroid_Entity = grid.getEclCentroid(element);
        const auto& elemCpGridEclCentroid_Index = grid.getEclCentroid(element.index());
        for (int coord = 0; coord < 3; ++coord)
        {
            assert(elemEclCentroid[coord] == elemCpGridEclCentroid_Entity[coord]);
            assert(elemEclCentroid[coord] == elemCpGridEclCentroid_Index[coord]);
            std::cout << "From Eclipse: " << elemEclCentroid[coord]
                      << " From CpGrid (Entity): " << elemCpGridEclCentroid_Entity[coord]
                      << " From CpGrid (Index): " << elemCpGridEclCentroid_Index[coord]<< '\n';
        }
        std::cout << " " << '\n';
    }

    {
        // Roundtrip: cornerpoint -> UnstructuredGrid file -> CpGrid,
        // and compare with direct CpGrid cornerpoint construction.
        std::vector<double> coord = ecl_grid.getCOORD();
        std::vector<double> zcorn = ecl_grid.getZCORN();
        std::vector<int> actnum = ecl_grid.getACTNUM();

        grdecl g{};
        g.dims[0] = static_cast<int>(ecl_grid.getNX());
        g.dims[1] = static_cast<int>(ecl_grid.getNY());
        g.dims[2] = static_cast<int>(ecl_grid.getNZ());
        g.coord = coord.data();
        g.zcorn = zcorn.data();
        g.actnum = actnum.empty() ? nullptr : actnum.data();

        std::unique_ptr<UnstructuredGrid, decltype(&destroy_grid)> ug(
            create_grid_cornerpoint(&g, 0.0, 0),
            &destroy_grid);
        assert(ug);

        const std::string roundtripFile = "cpgrid_cornerpoint_roundtrip.txt";
        serializeUnstructuredGridToFile(*ug, roundtripFile);

        Grid roundtripGrid(roundtripFile);
        compareCpGridCellsByGlobalId(grid, roundtripGrid);

        std::remove(roundtripFile.c_str());
    }

#endif

    std::stringstream dgfFile;
    // create unit cube with 8 cells in each direction
    dgfFile << "DGF" << std::endl;
    dgfFile << "Interval" << std::endl;
    dgfFile << "0 0 0" << std::endl;
    dgfFile << "4 4 4" << std::endl;
    dgfFile << "4 4 4" << std::endl;
    dgfFile << "#" << std::endl;

    Dune::GridPtr< Grid > gridPtr( dgfFile );
    testGrid( *gridPtr, "CpGrid_dgf", 64, 125 );

    {
        // Test unstructured grid read from file, constructed from filename.
        const std::string gridFileName = "cpgrid_test_from_filename.txt";
        std::stringstream gridFile;
        std::ofstream out(gridFileName);

        // UnstructuredGrid file format:
        // ndims ncells nfaces nnodes nfacenodes ncellfaces has_tag has_indexmap
        gridFile << "3 1 6 8 24 6 1 1" << std::endl;
        gridFile << "1 1 1" << std::endl;

        // node coordinates (8 * 3)
        gridFile << "0 0 0  1 0 0  0 1 0  1 1 0  0 0 1  1 0 1  0 1 1  1 1 1" << std::endl;

        // face_nodepos (nfaces + 1)
        gridFile << "0 4 8 12 16 20 24" << std::endl;
        // face_nodes (6 faces * 4 nodes)
        gridFile << "0 2 6 4  1 5 7 3  0 4 5 1  2 3 7 6  0 1 3 2  4 6 7 5" << std::endl;

        // face_cells (2 * nfaces)
        gridFile << "0 -1  0 -1  0 -1  0 -1  0 -1  0 -1" << std::endl;

        // face_areas (nfaces)
        gridFile << "1 1 1 1 1 1" << std::endl;
        // face_centroids (nfaces * 3)
        gridFile << "0 0.5 0.5  1 0.5 0.5  0.5 0 0.5  0.5 1 0.5  0.5 0.5 0  0.5 0.5 1" << std::endl;
        // face_normals (nfaces * 3)
        gridFile << "-1 0 0  1 0 0  0 -1 0  0 1 0  0 0 -1  0 0 1" << std::endl;

        // cell_facepos (ncells + 1)
        gridFile << "0 6" << std::endl;
        // cell_faces + facetags (ncellfaces pairs)
        // facetags: 0/1 I-/I+, 2/3 J-/J+, 4/5 K-/K+
        gridFile << "0 0  1 1  2 2  3 3  4 4  5 5" << std::endl;

        // global_cell map (ncells)
        gridFile << "0" << std::endl;
        // cell_volumes (ncells)
        gridFile << "1" << std::endl;
        // cell_centroids (ncells * 3)
        gridFile << "0.5 0.5 0.5" << std::endl;

        out << gridFile.str();
        out.close();

        Grid fileGrid(gridFileName);
        assert(fileGrid.size(0) == 1);
        assert(fileGrid.size(3) == 8);

        std::remove(gridFileName.c_str());
    }

    return 0;
}
