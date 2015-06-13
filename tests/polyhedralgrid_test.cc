#include "config.h"

#include <dune/grid/polyhedralgrid.hh>

#include <opm/parser/eclipse/Parser/Parser.hh>
#include <opm/parser/eclipse/Deck/Deck.hh>

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

int main()
{
    Opm::Parser parser;
    const auto deck = parser.parseString(deckString);

    PolyhedralGrid<dim, dimworld> grid(deck);

    return 0;
}
