/*
  Contributed by Martin Nolte, University of Freiburg.

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

#ifndef DUNE_CPGRID_DGFPARSER_HH
#define DUNE_CPGRID_DGFPARSER_HH

#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <opm/grid/CpGrid.hpp>

namespace Dune
{

  // DGFGridFactory< CpGrid >
  // ------------------------

  template<>
  struct DGFGridFactory< CpGrid >
  {
    typedef CpGrid Grid;

    typedef MPIHelper::MPICommunicator MPICommunicatorType;

    static const int dimension = Grid::dimension;
    typedef Grid::Codim< 0 >::Entity Element;
    typedef Grid::Codim< dimension >::Entity Vertex;

    explicit DGFGridFactory ( std::istream &input,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() );

    explicit DGFGridFactory ( const std::string &filename,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() );

    Grid *grid () const
    {
      return grid_;
    }

    template< class Intersection >
    bool wasInserted ( const Intersection &intersection ) const
    {
      return false;
    }

    template< class Intersection >
    int boundaryId ( const Intersection &intersection ) const
    {
      return intersection.boundaryId();
    }

    template< int codim >
    int numParameters () const
    {
      return 0;
    }

    template< class Entity >
    std::vector< double > &parameter ( const Entity& )
    {
      DUNE_THROW( InvalidStateException,
                  "Calling DGFGridFactory::parameter is only allowed if there are parameters." );
    }

    bool haveBoundaryParameters () const { return false; }

    template< class Intersection >
    const typename DGFBoundaryParameter::type &
    boundaryParameter ( const Intersection& ) const
    {
      return DGFBoundaryParameter::defaultValue();
    }

  private:
    void generate ( std::istream &input );

    Grid *grid_;
  };



  // Implementation of DGFGridFactory< CpGrid >
  // ------------------------------------------

  inline DGFGridFactory< CpGrid >
    ::DGFGridFactory ( std::istream &input, MPICommunicatorType )
  : grid_( 0 )
  {
    generate( input );
  }


  inline DGFGridFactory< CpGrid >
    ::DGFGridFactory ( const std::string &filename, MPICommunicatorType )
  : grid_( 0 )
  {
    std::ifstream input( filename.c_str() );
    if( !input )
      DUNE_THROW( DGFException, "Unable to open file: " << filename << "." );
    generate( input );
    input.close();
  }


  inline void DGFGridFactory< CpGrid >::generate ( std::istream &input )
  {
    dgf::IntervalBlock intervalBlock( input );

    if( !intervalBlock.isactive() )
      DUNE_THROW( DGFException, "DGF stream must contain an interval block to be used with CpGrid." );
    if( intervalBlock.numIntervals() != 1 )
      DUNE_THROW( DGFException, "Currently, CpGrid can only handle 1 interval block." );

    if( intervalBlock.dimw() != dimension )
      DUNE_THROW( DGFException, "CpGrid cannot handle an interval of dimension " << intervalBlock.dimw() << "." );
    const dgf::IntervalBlock::Interval &interval = intervalBlock.get( 0 );

    // compute some values

    const double dx = (interval.p[ 1 ][ 0 ] - interval.p[ 0 ][ 0 ]) / interval.n[ 0 ];
    const double dy = (interval.p[ 1 ][ 1 ] - interval.p[ 0 ][ 1 ]) / interval.n[ 1 ];
    const double dz = (interval.p[ 1 ][ 2 ] - interval.p[ 0 ][ 2 ]) / interval.n[ 2 ];

    const double bottom = interval.p[ 0 ][ 2 ];
    const double top = interval.p[ 1 ][ 2 ];

    // compute pillars

    std::vector< double > coord;
    coord.reserve( 6*(interval.n[ 0 ] + 1)*(interval.n[ 1 ] + 1) );
    for( int j = 0; j <= interval.n[ 1 ]; ++j )
    {
      const double y = j*dy;
      for( int i = 0; i <= interval.n[ 0 ]; ++i )
      {
        const double x = i*dx;
        const double pillar[ 6 ] = { x, y, bottom, x, y, top };
        coord.insert( coord.end(), pillar, pillar + 6 );
      }
    }

    // create cubes
    const int num_per_layer = 4*interval.n[ 0 ]*interval.n[ 1 ];
    std::vector< double > zcorn( 2*num_per_layer*interval.n[ 2 ] );
    double *offset = &zcorn[ 0 ];
    for( int k = 0; k < interval.n[ 2 ]; ++k )
    {
      const double zlow = k*dz;
      std::fill( offset, offset + num_per_layer, zlow );
      offset += num_per_layer;
      const double zhigh = (k+1)*dz;
      std::fill( offset, offset + num_per_layer, zhigh );
      offset += num_per_layer;
    }

    // ???
    std::vector< int > actnum( interval.n[ 0 ]*interval.n[ 1 ]*interval.n[ 2 ], 1 );

    // create eclipse format description
    grdecl g;
    g.dims[ 0 ] = interval.n[ 0 ];
    g.dims[ 1 ] = interval.n[ 1 ];
    g.dims[ 2 ] = interval.n[ 2 ];
    g.coord = &coord[ 0 ];
    g.zcorn = &zcorn[ 0 ];
    g.actnum = &actnum[ 0 ];

    // create grid
    grid_ = new Grid;
    grid_->processEclipseFormat( g, 0.0, false, false );
  }



  // DGFGridInfo< CpGrid >
  // ---------------------

  template<>
  struct DGFGridInfo< CpGrid >
  {
    typedef CpGrid Grid;

    static int refineStepsForHalf ()
    {
      return 1;
    }

    static double refineWeight ()
    {
      return 1.0 / 8.0;
    }
  };

}

#endif // #ifndef DUNE_CPGRID_DGFPARSER_HH
