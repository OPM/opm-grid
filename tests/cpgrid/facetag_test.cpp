/*
  Copyright 2016 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2016 Statoil ASA.

  This file is part of The Open Porous Media project  (OPM).

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
#include <config.h>

#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif
#define NVERBOSE // to suppress our messages when throwing


#define BOOST_TEST_MODULE FaceTagTests
#include <boost/test/unit_test.hpp>
#include <dune/grid/CpGrid.hpp>
#include <dune/grid/cpgrid/GridHelpers.hpp>

#include <array>

BOOST_AUTO_TEST_CASE(facetag)
{
    int m_argc = boost::unit_test::framework::master_test_suite().argc;
    char** m_argv = boost::unit_test::framework::master_test_suite().argv;
    Dune::MPIHelper::instance(m_argc, m_argv);
    Dune::CpGrid grid;
    std::array<int, 3>    dims     = { 3, 3, 3 };
    std::array<double, 3> cellsize = { 1., 1., 1. };
    grid.createCartesian(dims, cellsize);
    Dune::cpgrid::Cell2FacesContainer c2f(&grid);
    
    for( int cell=0; cell < grid.numCells(); ++cell)
    {
        std::cout<<"cell="<<cell;
        auto cell_faces = c2f[cell];
        for(auto face = cell_faces.begin(), endFace = cell_faces.end();
            face != endFace; ++face)
        {
            auto tag = grid.faceTag(face);
            auto c0 = grid.faceCell(*face, 0), c1 = grid.faceCell(*face, 1);
            std::cout<<"   face="<<*face<<" c0="<<c0<<" c1="<<c1<<std::endl;
            if ( c1 < 0 || c0 < 0)
            {
                // boundary face
                std::array<int, 3> ijk = {{-1, -1, -1}};
                grid.getIJK(c0<0?c1:c0, ijk);
                bool valid_tag = false;
                std::cout<<"      ijk="<<ijk[0]<<" "<<ijk[1]<<" "<<ijk[2]<<" tag="<<tag<<std::endl;
                for ( int dim = 0; dim < 3; ++dim)
                {
                    if ( ijk[dim] == 0 )
                    {
                        valid_tag = valid_tag || ( tag == 2*dim );
                    }
                    if  ( ijk[dim] == 2 )
                    {
                        valid_tag = valid_tag || ( tag == 2*dim+1 );
                    }
                }
                BOOST_CHECK( valid_tag );
            }
            else
            {
                std::array<int, 3> ijk0 = {{-1, -1, -1}}, ijk1 = {{-1, -1, -1}};
                grid.getIJK(c0, ijk0);
                grid.getIJK(c1, ijk1);
                std::cout<<"      ijk0="<<ijk0[0]<<" "<<ijk0[1]<<" "<<ijk0[2]
                         <<" ijk1="<<ijk1[0]<<" "<<ijk1[1]<<" "<<ijk1[2]<<" tag="<<tag<<std::endl;
                BOOST_ASSERT( ijk0[0] <= ijk1[0] &&  ijk0[1] <= ijk1[1] &&  ijk0[2] <= ijk1[2]);
                int firstInside = (cell==c0) ? 1 : 0;
                
                for ( int dim = 0; dim < 3; ++dim)
                {
                    if ( ijk0[dim] < ijk1[dim] )
                    {
                        BOOST_CHECK( tag == 2 * dim + firstInside);
                    }
                }
            }
        }
    }
}
