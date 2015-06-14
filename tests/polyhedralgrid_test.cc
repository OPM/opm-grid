#include "config.h"

#include <dune/grid/polyhedralgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>

#include <iostream>

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
    "4*100.0 /\n";

template <class GridView>
void testGrid( const GridView& gridView )
{
    typedef typename GridView::template Codim<0>::Iterator ElemIterator;
    typedef typename GridView::IntersectionIterator IsIt;

    int numElem = 0;
    ElemIterator elemIt = gridView.template begin<0>();
    ElemIterator elemEndIt = gridView.template end<0>();
    for (; elemIt != elemEndIt; ++elemIt) {
        const auto& elem = *elemIt;
        const auto& elemGeom = elem.geometry();
        if (std::abs(elemGeom.volume() - 1.0) > 1e-8)
            std::cout << "element's " << numElem << " volume is wrong:"<<elemGeom.volume()<<"\n";

        int numIs = 0;
        IsIt isIt = gridView.ibegin(elem);
        IsIt isEndIt = gridView.iend(elem);
        for (; isIt != isEndIt; ++isIt, ++ numIs) {
            const auto& isGeom = (*isIt).geometry();
            if (std::abs(isGeom.volume() - 1.0) > 1e-8)
                std::cout << "volume of intersection " << numIs << " of element " << numElem << " volume is wrong: " << isGeom.volume() << "\n";

            if (isIt->neighbor())
                if (std::abs(isIt->outside().geometry().volume() - 1.0) > 1e-8)
                    std::cout << "outside element volume of intersection " << numIs << " of element " << numElem << " volume is wrong: " << isIt->outside().geometry().volume() << "\n";

            if (std::abs(isIt->inside().geometry().volume() - 1.0) > 1e-8)
                std::cout << "inside element volume of intersection " << numIs << " of element " << numElem << " volume is wrong: " << isIt->inside().geometry().volume() << "\n";
        }

        if (numIs != 6)
            std::cout << "number of intersections is wrong for element " << numElem << "\n";

        ++ numElem;
    }

    if (numElem != 2*2*2)
        std::cout << "number of elements is wrong: " << numElem << "\n";
}

int main()
{
    Opm::Parser parser;
    const auto deck = parser.parseString(deckString);

    std::vector<double> porv;
    typedef Dune::PolyhedralGrid< 3, 3 > Grid;
    typedef Grid::LeafGridView GridView;
    Grid grid(deck, porv);

    testGrid( grid.leafGridView() );

    std::cout << "create vertex mapper\n";
    Dune::MultipleCodimMultipleGeomTypeMapper<GridView,
                                              Dune::MCMGVertexLayout> mapper(grid.leafGridView());

    std::cout << "VertexMapper.size(): " << mapper.size() << "\n";
    if (mapper.size() != 27) {
        std::cout << "Wrong size of vertex mapper. Expected 27!\n";
        std::abort();
    }

    std::cout << "create vtkWriter\n";
    typedef Dune::VTKWriter<GridView> VtkWriter;
    VtkWriter vtkWriter(grid.leafGridView());

    std::cout << "create cellData\n";
    int numElems = grid.size(0);
    std::vector<double> tmpData(numElems, 0.0);

    std::cout << "add cellData\n";
    vtkWriter.addCellData(tmpData, "testdata");

    std::cout << "write\n";
    vtkWriter.write("polyhedralgrid_test", Dune::VTK::ascii);

    return 0;
}
