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
        for (; isIt != isEndIt; ++isIt) {
            const auto& isGeom = (*isIt).geometry();
            if (std::abs(isGeom.volume() - 1.0) > 1e-8)
                std::cout << "volume of intersection " << numIs << " of element " << numElem << " volume is wrong: " << isGeom.volume() << "\n";
            ++ numIs;
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

    typedef Dune::VTKWriter<GridView> VtkWriter;
    VtkWriter vtkWriter(grid.leafGridView());
    vtkWriter.write("polyhedralgrid_test");
    return 0;
}
