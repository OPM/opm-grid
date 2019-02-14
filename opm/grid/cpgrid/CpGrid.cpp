//===========================================================================
//
// File: CpGrid.cpp
//
// Created: Thu Jun  4 12:55:28 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            B�rd Skaflestad     <bard.skaflestad@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.

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
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#if HAVE_MPI
#include <opm/grid/utility/platform_dependent/disable_warnings.h>
#include "mpi.h"
#include <opm/grid/utility/platform_dependent/reenable_warnings.h>
#endif

#include "../CpGrid.hpp"
#include "CpGridData.hpp"
#include <opm/grid/common/ZoltanPartition.hpp>
#include <opm/grid/common/GridPartitioning.hpp>
#include <opm/grid/common/WellConnections.hpp>

#include <fstream>
#include <iostream>

namespace Dune
{

    CpGrid::CpGrid()
        : data_( new cpgrid::CpGridData(*this)),
          current_view_data_(data_.get()),
          distributed_data_(),
          cell_scatter_gather_interfaces_(new InterfaceMap)
    {}





std::pair<bool, std::unordered_set<std::string> >
CpGrid::scatterGrid(const std::vector<const cpgrid::OpmWellType *> * wells,
                    const double* transmissibilities, int overlapLayers)
{
    // Silence any unused argument warnings that could occur with various configurations.
    static_cast<void>(wells);
    static_cast<void>(transmissibilities);
    static_cast<void>(overlapLayers);
#if HAVE_MPI
    if(distributed_data_)
    {
        std::cerr<<"There is already a distributed version of the grid."
                 << " Maybe scatterGrid was called before?"<<std::endl;
        return std::make_pair(false, std::unordered_set<std::string>());
    }

    CollectiveCommunication cc(MPI_COMM_WORLD);

    int my_num=cc.rank();
#ifdef HAVE_ZOLTAN
    auto part_and_wells =
        cpgrid::zoltanGraphPartitionGridOnRoot(*this, wells, transmissibilities, cc, 0);
    int num_parts = cc.size();
    using std::get;
    auto cell_part = std::get<0>(part_and_wells);
    auto defunct_wells = std::get<1>(part_and_wells);
#else
    std::vector<int> cell_part(current_view_data_->global_cell_.size());
    int  num_parts=-1;
    std::array<int, 3> initial_split;
    initial_split[1]=initial_split[2]=std::pow(cc.size(), 1.0/3.0);
    initial_split[0]=cc.size()/(initial_split[1]*initial_split[2]);
    partition(*this, initial_split, num_parts, cell_part, false, false);
    const auto& cpgdim =  logicalCartesianSize();
    std::vector<int> cartesian_to_compressed(cpgdim[0]*cpgdim[1]*cpgdim[2], -1);
    for( int i=0; i < numCells(); ++i )
    {
        cartesian_to_compressed[globalCell()[i]] = i;
    }

    std::unordered_set<std::string> defunct_wells;

    if ( wells )
    {
        cpgrid::WellConnections well_connections(*wells,
                                                 cpgdim,
                                                 cartesian_to_compressed);

        auto wells_on_proc =
            cpgrid::postProcessPartitioningForWells(cell_part,
                                                    *wells,
                                                    well_connections,
                                                    cc.size());
        defunct_wells = cpgrid::computeDefunctWellNames(wells_on_proc,
                                                        *wells,
                                                        cc,
                                                        0);
    }
#endif

    MPI_Comm new_comm = MPI_COMM_NULL;

    if(num_parts < cc.size())
    {
        std::vector<int> ranks(num_parts);
        for(int i=0; i<num_parts; ++i)
            ranks[i]=i;
        MPI_Group new_group;
        MPI_Group old_group;
        MPI_Comm_group(cc, &old_group);
        MPI_Group_incl(old_group, num_parts, &(ranks[0]), &new_group);

        // Not all procs take part in the parallel computation
        MPI_Comm_create(cc, new_group, &new_comm);
        cc=CollectiveCommunication(new_comm);
    }else{
        new_comm = cc;
    }
    if(my_num<cc.size())
    {
        distributed_data_.reset(new cpgrid::CpGridData(new_comm));
        distributed_data_->distributeGlobalGrid(*this,*this->current_view_data_, cell_part,
                                                overlapLayers);
        std::cout << "After loadbalancing process " << my_num << " has " <<
            distributed_data_->cell_to_face_.size() << " cells." << std::endl;

        // add an interface for gathering/scattering data with communication
        // forward direction will be scatter and backward gather
        cell_scatter_gather_interfaces_.reset(new InterfaceMap);

        auto rank = distributed_data_->ccobj_.rank();

        if ( rank == 0)
        {
            std::map<int, std::size_t> proc_to_no_cells;
            for(auto cell_owner = cell_part.begin(); cell_owner != cell_part.end();
                ++cell_owner)
            {
                ++proc_to_no_cells[*cell_owner];
            }

            for(const auto& proc_no_cells : proc_to_no_cells)
            {
                (*cell_scatter_gather_interfaces_)[proc_no_cells.first]
                    .first.reserve(proc_no_cells.second);
            }

            std::size_t cell_index = 0;

            for(auto cell_owner = cell_part.begin(); cell_owner != cell_part.end();
                ++cell_owner, ++cell_index)
            {
                auto& indices = (*cell_scatter_gather_interfaces_)[*cell_owner];
                indices.first.add(cell_index);
            }

        }

        (*cell_scatter_gather_interfaces_)[0].second
            .reserve(distributed_data_->cell_indexset_.size());

        for( auto& index: distributed_data_->cell_indexset_)
        {
            typedef typename cpgrid::CpGridData::AttributeSet AttributeSet;
            if ( index.local().attribute() == AttributeSet::owner)
            {
                auto& indices = (*cell_scatter_gather_interfaces_)[0];
                indices.second.add(index.local());
            }
        }
    }
    current_view_data_ = distributed_data_.get();
    return std::make_pair(true, defunct_wells);

#else // #if HAVE_MPI
    std::cerr << "CpGrid::scatterGrid() is non-trivial only with "
              << "MPI support and if the target Dune platform is "
              << "sufficiently recent.\n";
    return std::make_pair(false, std::unordered_set<std::string>());
#endif
}


    void CpGrid::createCartesian(const std::array<int, 3>& dims,
                                 const std::array<double, 3>& cellsize)
    {
        // Make the grdecl format arrays.
        // Pillar coords.
        std::vector<double> coord;
        coord.reserve(6*(dims[0] + 1)*(dims[1] + 1));
        double bot = 0.0;
        double top = dims[2]*cellsize[2];
        // i runs fastest for the pillars.
        for (int j = 0; j < dims[1] + 1; ++j) {
            double y = j*cellsize[1];
            for (int i = 0; i < dims[0] + 1; ++i) {
                double x = i*cellsize[0];
                double pillar[6] = { x, y, bot, x, y, top };
                coord.insert(coord.end(), pillar, pillar + 6);
            }
        }
        std::vector<double> zcorn(8*dims[0]*dims[1]*dims[2]);
        const int num_per_layer = 4*dims[0]*dims[1];
        double* offset = &zcorn[0];
        for (int k = 0; k < dims[2]; ++k) {
            double zlow = k*cellsize[2];
            std::fill(offset, offset + num_per_layer, zlow);
            offset += num_per_layer;
            double zhigh = (k+1)*cellsize[2];
            std::fill(offset, offset + num_per_layer, zhigh);
            offset += num_per_layer;
        }
        std::vector<int> actnum(dims[0]*dims[1]*dims[2], 1);

        // Process them.
        grdecl g;
        g.dims[0] = dims[0];
        g.dims[1] = dims[1];
        g.dims[2] = dims[2];
        g.coord = &coord[0];
        g.zcorn = &zcorn[0];
        g.actnum = &actnum[0];
        current_view_data_->processEclipseFormat(g, {}, 0.0, false, false);
    }

    void CpGrid::readSintefLegacyFormat(const std::string& grid_prefix)
    {
        current_view_data_->readSintefLegacyFormat(grid_prefix);
    }
    void CpGrid::writeSintefLegacyFormat(const std::string& grid_prefix) const
    {
        current_view_data_->writeSintefLegacyFormat(grid_prefix);
    }


#if HAVE_ECL_INPUT
    void CpGrid::processEclipseFormat(const Opm::EclipseGrid& ecl_grid,
                                      bool periodic_extension,
                                      bool turn_normals, bool clip_z,
                                      const std::vector<double>& poreVolume,
                                      const Opm::NNC& nncs)
    {
        current_view_data_->processEclipseFormat(ecl_grid, periodic_extension,
                                                 turn_normals, clip_z,
                                                 poreVolume, nncs);
    }
#endif

    void CpGrid::processEclipseFormat(const grdecl& input_data, double z_tolerance,
                                      bool remove_ij_boundary, bool turn_normals)
    {
        current_view_data_->processEclipseFormat(input_data, {}, z_tolerance, remove_ij_boundary, turn_normals);
    }

} // namespace Dune
