#include <config.h>

// Warning suppression for Dune includes.
#include <opm/grid/utility/platform_dependent/disable_warnings.h>

#include <dune/common/unused.hh>
#include <opm/grid/polyhedralgrid.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <opm/grid/cpgrid/dgfparser.hh>
#include <opm/grid/polyhedralgrid/dgfparser.hh>


#define DISABLE_DEPRECATED_METHOD_CHECK 1
#include <dune/grid/test/gridcheck.hh>

// Re-enable warnings.
#include <opm/grid/utility/platform_dependent/reenable_warnings.h>

#include <opm/grid/utility/OpmParserIncludes.hpp>

#include <iostream>

// two hexahedrons using polygon/polyhedron format
static const char* hexaPoly = "\
DGF\n \
vertex\n \
0.0   0.0   0.0\n \
0.5   0.0   0.0\n \
1.0   0.0   0.0\n \
0.0   1.0   0.0\n \
0.5   1.0   0.0\n \
1.0   1.0   0.0\n \
0.0   0.0   1.0\n \
0.5   0.0   1.0\n \
1.0   0.0   1.0\n \
0.0   1.0   1.0\n \
0.5   1.0   1.0\n \
1.0   1.0   1.0\n \
#\n \
Polygon\n \
0  3   9   6\n \
1  4  10   7\n \
0  1   7   6\n \
3  4  10   9\n \
0  1   4   3\n \
6  7  10   9\n \
2  5  11   8\n \
1  2   8   7\n \
4  5  11  10\n \
1  2   5   4\n \
7  8  11  10\n \
#\n \
Polyhedron\n \
0  1  2  3  4   5\n \
1  6  7  8  9  10\n \
#";

static const char* tetraPoly = "\
DGF\n \
VERTEX\n \
firstindex 1\n \
0.0 0.0 0.0\n \
1.0 0.0 0.0\n \
1.0 1.0 0.0\n \
0.0 1.0 0.0\n \
0.0 0.0 1.0\n \
1.0 0.0 1.0\n \
1.0 1.0 1.0\n \
0.0 1.0 1.0\n \
0.5 0.0 0.0\n \
0.0 0.5 0.0\n \
0.0 0.5 1.0\n \
5.0 0.0 1.0\n \
1.0 0.0 5.0\n \
1.0 5.0 0.0\n \
1.0 1.0 5.0\n \
0.5 1.0 0.0\n \
0.5 0.5 0.0\n \
1.0 0.5 0.5\n \
1.0 0.5 1.0\n \
0.5 1.0 1.0\n \
0.0 0.0 0.5\n \
0.0 1.0 0.5\n \
0.5 0.0 0.5\n \
0.5 1.0 0.5\n \
0.5 0.5 1.0\n \
0.0 0.5 0.5\n \
#\n \
Polyhedron\n \
39 40 41 83\n \
8 10 11 51\n \
81 82 83 96\n \
35 36 38 70\n \
12 13 14 78\n \
71 72 74 98\n \
86 87 88 110\n \
5 6 7 92\n \
9 10 14 53\n \
84 85 88 95\n \
106 107 108 112\n \
68 69 71 75\n \
19 20 21 84\n \
22 23 24 57\n \
57 58 62 87\n \
34 35 37 68\n \
104 105 107 110\n \
96 97 100 103\n \
1 2 6 54\n \
56 58 60 85\n \
53 54 55 89\n \
46 47 49 64\n \
76 77 80 82\n \
70 72 73 97\n \
78 79 80 90\n \
16 17 21 79\n \
15 17 18 77\n \
43 44 45 106\n \
59 60 61 93\n \
3 4 5 59\n \
51 52 55 75\n \
28 29 32 74\n \
65 66 67 111\n \
30 31 33 109\n \
41 42 45 103\n \
25 26 32 65\n \
92 93 94 109\n \
25 27 33 66\n \
48 49 50 105\n \
63 64 67 108\n \
90 91 95 101\n \
100 101 102 112\n \
89 91 94 99\n \
98 99 102 111\n \
#\n \
Polygon\n \
1 9 17\n \
1 9 23\n \
1 10 17\n \
1 10 21\n \
1 17 21\n \
1 17 23\n \
1 21 23\n \
2 9 13\n \
2 9 17\n \
2 9 18\n \
2 13 18\n \
2 14 17\n \
2 14 18\n \
2 17 18\n \
3 14 15\n \
3 14 17\n \
3 14 24\n \
3 15 24\n \
3 16 17\n \
3 16 24\n \
3 17 24\n \
4 10 16\n \
4 10 22\n \
4 16 22\n \
5 11 23\n \
5 11 25\n \
5 11 26\n \
5 12 23\n \
5 12 25\n \
5 21 23\n \
5 21 26\n \
5 23 25\n \
5 23 26\n \
6 12 13\n \
6 12 18\n \
6 12 19\n \
6 13 18\n \
6 18 19\n \
7 15 19\n \
7 15 24\n \
7 19 24\n \
7 19 25\n \
7 20 24\n \
7 20 25\n \
7 24 25\n \
8 11 20\n \
8 11 26\n \
8 20 22\n \
8 20 26\n \
8 22 26\n \
9 13 18\n \
9 13 23\n \
9 17 18\n \
9 17 23\n \
9 18 23\n \
10 16 17\n \
10 16 22\n \
10 16 26\n \
10 17 21\n \
10 17 26\n \
10 21 26\n \
10 22 26\n \
11 20 25\n \
11 20 26\n \
11 23 25\n \
11 23 26\n \
11 25 26\n \
12 13 18\n \
12 13 23\n \
12 18 19\n \
12 18 23\n \
12 18 25\n \
12 19 25\n \
12 23 25\n \
13 18 23\n \
14 15 18\n \
14 15 24\n \
14 17 18\n \
14 17 24\n \
14 18 24\n \
15 18 19\n \
15 18 24\n \
15 19 24\n \
16 17 24\n \
16 17 26\n \
16 22 24\n \
16 22 26\n \
16 24 26\n \
17 18 23\n \
17 18 24\n \
17 18 26\n \
17 21 23\n \
17 21 26\n \
17 23 26\n \
17 24 26\n \
18 19 24\n \
18 19 25\n \
18 23 25\n \
18 23 26\n \
18 24 25\n \
18 24 26\n \
18 25 26\n \
19 24 25\n \
20 22 24\n \
20 22 26\n \
20 24 25\n \
20 24 26\n \
20 25 26\n \
21 23 26\n \
22 24 26\n \
23 25 26\n \
24 25 26\n \
#\n \
BOUNDARYDOMAIN\n \
default 1\n \
#\n \
GRIDPARAMETER\n \
closure none\n \
#";

int main(int argc, char** argv )
{
    // initialize MPI
    Dune::MPIHelper::instance( argc, argv );

    // test PolyhedralGrid
    {
        typedef Dune::PolyhedralGrid< 3, 3 > Grid;

#if HAVE_ECL_INPUT
        const char *deckString =
            "RUNSPEC\n"
            "METRIC\n"
            "DIMENS\n"
            "4 4 4 /\n"
            "GRID\n"
            "DXV\n"
            "4*1 /\n"
            "DYV\n"
            "4*1 /\n"
            "DZ\n"
            "16*1 /\n"
            "TOPS\n"
            "16*100.0 /\n";

        Opm::Parser parser;
        const auto deck = parser.parseString(deckString);
        Opm::EclipseGrid eclgrid( deck);
        std::vector<double> porv;

        std::cout <<"Check 3d grid created from deck" << std::endl << std::endl;
        Grid grid(eclgrid, porv);
        gridcheck( grid );
        std::cout << std::endl;
#endif
        // test DGF grid creation capabilities
        std::stringstream dgfFile;
        // create unit cube with 8 cells in each direction
        dgfFile << "DGF" << std::endl;
        dgfFile << "Interval" << std::endl;
        dgfFile << "0 0 0" << std::endl;
        dgfFile << "1 1 1" << std::endl;
        dgfFile << "2 1 1" << std::endl;
        dgfFile << "#" << std::endl;

        std::cout <<"Check 3d Cartesian grid created from DGF file" << std::endl << std::endl;
        Dune::GridPtr< Grid > gridPtr( dgfFile );
        gridcheck( *gridPtr );
        std::cout << std::endl;

        {
            std::cout <<"Check 3d hexahedral grid created from DGF file" << std::endl << std::endl;
            std::stringstream poly;
            poly << hexaPoly;
            Dune::GridPtr< Grid > gridPoly( poly );
            gridcheck( *gridPoly );
            std::cout << std::endl;
        }

        {
            std::cout <<"Check 3d tetrahedral grid created from DGF file" << std::endl << std::endl;
            std::stringstream poly;
            poly << tetraPoly;
            Dune::GridPtr< Grid > gridPoly( poly );
            gridcheck( *gridPoly );
            std::cout << std::endl;
        }

    }

    {
        std::cout <<"Check 2d Cartesian grid created from DGF file" << std::endl << std::endl;
        std::stringstream dgfFile;
        // create unit cube with 8 cells in each direction
        dgfFile << "DGF" << std::endl;
        dgfFile << "Interval" << std::endl;
        dgfFile << "0 0" << std::endl;
        dgfFile << "1 1" << std::endl;
        dgfFile << "8 8" << std::endl;
        dgfFile << "#" << std::endl;
        typedef Dune::PolyhedralGrid< 2, 2 > Grid;
        Dune::GridPtr< Grid > gridPtr( dgfFile );
        gridcheck( *gridPtr );
    }
    return 0;
}
