/*
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services.
  Copyright 2015 NTNU
  Copyright 2015 Statoil AS

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
#ifndef DUNE_CPGRID_ZOLTANPARTITION_HEADER
#define DUNE_CPGRID_ZOLTANPARTITION_HEADER

#include <dune/grid/CpGrid.hpp>
#include <dune/grid/common/ZoltanGraphFunctions.hpp>

#if HAVE_ZOLTAN
namespace Dune
{
namespace cpgrid
{
std::vector<int> zoltanGraphPartitionGrid(const CpGrid& grid,
                                          const CollectiveCommunication<MPI_Comm>& cc,
                                          bool globalGridOnAllProcs);
}
}
#endif // HAVE_ZOLTAN
#endif // header guard
