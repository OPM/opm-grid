/*
  Copyright 2015 IRIS AS

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
#ifndef DUNE_POLYHEDRALGRID_GRIDHELPERS_HEADER_INCLUDED
#define DUNE_POLYHEDRALGRID_GRIDHELPERS_HEADER_INCLUDED

#include <opm/grid/GridHelpers.hpp>
#include <dune/grid/polyhedralgrid.hh>

namespace Opm
{
namespace UgGridHelpers
{

template<int dim, int dimworld>
struct CellCentroidTraits< Dune::PolyhedralGrid< dim, dimworld > >
  : public CellCentroidTraits<UnstructuredGrid>
{
};

template<int dim, int dimworld>
struct CellVolumeIteratorTraits< Dune::PolyhedralGrid< dim, dimworld > >
  : public CellVolumeIteratorTraits<UnstructuredGrid>
{
};

template<int dim, int dimworld>
struct FaceCentroidTraits< Dune::PolyhedralGrid< dim, dimworld > >
 : public FaceCentroidTraits< UnstructuredGrid >
{
};

template<int dim, int dimworld>
struct Cell2FacesTraits< Dune::PolyhedralGrid< dim, dimworld > >
 : public Cell2FacesTraits<UnstructuredGrid>
{
};

template<int dim, int dimworld>
struct Face2VerticesTraits< Dune::PolyhedralGrid< dim, dimworld > >
 : public Face2VerticesTraits<UnstructuredGrid>
{
};

template<int dim, int dimworld>
struct FaceCellTraits< Dune::PolyhedralGrid< dim, dimworld > >
 : public FaceCellTraits<UnstructuredGrid>
{
};

} // end namespace UGGridHelpers
} // end namespace OPM
#endif
