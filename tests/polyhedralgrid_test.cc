#include "config.h"

#include <dune/grid/polyhedralgrid.hh>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>

#include <iostream>

const char *deckString =
    "RUNSPEC\n"
    "DIMENS\n"
    "10 10 10 /\n"
    "GRID\n"
    "DXV\n"
    "10*1 /\n"
    "DYV\n"
    "10*1 /\n"
    "DZV\n"
    "10*1 /\n";

template <class GridView>
void testGrid( const GridView& gridView )
{
  typedef typename GridView :: template Codim< 0 > :: Iterator Iterator;
  typedef typename GridView :: template Codim< 0 > :: Entity   Entity;
  typedef typename GridView :: IndexSet IndexSet;
  const IndexSet& indexSet = gridView.indexSet();
  for( Iterator it = gridView.template begin< 0 > (),
       end = gridView.template begin< 0 > (); it != end; ++it )
  {
    const Entity& entity = *it ;
    std::cout << "Entity[ " << indexSet.index( entity ) << " ] = "
              << entity.geometry().center() << std::endl;
  }
}

int main()
{
    Opm::Parser parser;
    const auto deck = parser.parseString(deckString);

    std::vector<double> porv;
    Grid grid(deck, porv);

    typedef Grid::template Codim<0>::EntityIterator ElemIterator;
    ElemIterator elemIt = grid.template begin<0>();
    ElemIterator elemEndIt = grid.template end<0>();
    for (; elemIt != elemEndIt; ++elemIt) {
    }

    testGrid( grid.leafGridView() );

    return 0;
}
