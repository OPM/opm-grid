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
#include <iomanip>

namespace
{

#if HAVE_MPI

using AttributeSet = Dune::cpgrid::CpGridData::AttributeSet;

template<typename Tuple, bool first>
void reserveInterface(const std::vector<Tuple>& list, Dune::CpGrid::InterfaceMap& interface,
                      const std::integral_constant<bool, first>&)
{
    std::map<int, std::size_t> proc_to_no_cells;
    for(const auto& entry: list)
    {
        ++proc_to_no_cells[std::get<1>(entry)];
    }
    for(const auto& proc: proc_to_no_cells)
    {
        auto& entry = interface[proc.first];
        if ( first )
            entry.first.reserve(proc.second);
        else
            entry.second.reserve(proc.second);
    }
}

void setupSendInterface(const std::vector<std::tuple<int, int, char> >& list, Dune::CpGrid::InterfaceMap& interface)
{
    reserveInterface(list, interface, std::integral_constant<bool, true>());
    int cellIndex=-1;
    int oldIndex = std::numeric_limits<int>::max();
    for(const auto& entry: list)
    {
        auto index = std::get<0>(entry);
        assert(oldIndex == std::numeric_limits<int>::max() || index >= oldIndex);

        if (index != oldIndex )
        {
            oldIndex = index;
            ++cellIndex;
        }
        interface[std::get<1>(entry)].first.add(cellIndex);
    }
}

void setupRecvInterface(const std::vector<std::tuple<int, int, char, int> >& list, Dune::CpGrid::InterfaceMap& interface)
{
    reserveInterface(list, interface, std::integral_constant<bool, false>());
    for(const auto& entry: list)
    {
        auto index = std::get<3>(entry);
        interface[std::get<1>(entry)].second.add(index);
    }
}
#endif // HAVE_MPI
}

namespace Dune
{

    CpGrid::CpGrid()
        : data_( new cpgrid::CpGridData(*this)),
          current_view_data_(data_.get()),
          distributed_data_(),
          cell_scatter_gather_interfaces_(new InterfaceMap),
          point_scatter_gather_interfaces_(new InterfaceMap)
    {}


    CpGrid::CpGrid(MPIHelper::MPICommunicator  comm)
        : data_( new cpgrid::CpGridData(comm)),
          current_view_data_(data_.get()),
          distributed_data_(),
          cell_scatter_gather_interfaces_(new InterfaceMap),
          point_scatter_gather_interfaces_(new InterfaceMap)
    {}



std::pair<bool, std::unordered_set<std::string> >
CpGrid::scatterGrid(EdgeWeightMethod method, bool ownersFirst, const std::vector<cpgrid::OpmWellType> * wells,
                    const double* transmissibilities, bool addCornerCells, int overlapLayers)
{
    // Silence any unused argument warnings that could occur with various configurations.
    static_cast<void>(wells);
    static_cast<void>(transmissibilities);
    static_cast<void>(overlapLayers);
    static_cast<void>(method);
    if(distributed_data_)
    {
        std::cerr<<"There is already a distributed version of the grid."
                 << " Maybe scatterGrid was called before?"<<std::endl;
        return std::make_pair(false, std::unordered_set<std::string>());
    }

#if HAVE_MPI
    auto& cc = data_->ccobj_;

    if (cc.size() > 1)
    {
#ifdef HAVE_ZOLTAN
        auto part_and_wells =
            cpgrid::zoltanGraphPartitionGridOnRoot(*this, wells, transmissibilities, cc, method, 0);
        using std::get;
        auto cell_part = std::get<0>(part_and_wells);
        auto defunct_wells = std::get<1>(part_and_wells);
        auto exportList = std::get<2>(part_and_wells);
        auto importList = std::get<3>(part_and_wells);
#else
        OPM_THROW(std::runtime_error, "Parallel runs depend on ZOLTAN. Please install!");
        // std::vector<int> cell_part(current_view_data_->global_cell_.size());
        // int  num_parts=-1;
        // std::array<int, 3> initial_split;
        // initial_split[1]=initial_split[2]=std::pow(cc.size(), 1.0/3.0);
        // initial_split[0]=cc.size()/(initial_split[1]*initial_split[2]);
        // partition(*this, initial_split, num_parts, cell_part, false, false);
        // const auto& cpgdim =  logicalCartesianSize();
        // std::vector<int> cartesian_to_compressed(cpgdim[0]*cpgdim[1]*cpgdim[2], -1);
        // for( int i=0; i < numCells(); ++i )
        // {
        //     cartesian_to_compressed[globalCell()[i]] = i;
        // }

        // std::unordered_set<std::string> defunct_wells;

        // if ( wells )
        // {
        //     cpgrid::WellConnections well_connections(*wells,
        //                                              cpgdim,
        //                                              cartesian_to_compressed);

        //     auto wells_on_proc =
        //         cpgrid::postProcessPartitioningForWells(cell_part,
        //                                                 *wells,
        //                                                 well_connections,
        //                                                 cc.size());
        //     defunct_wells = cpgrid::computeDefunctWellNames(wells_on_proc,
        //                                                     *wells,
        //                                                     cc,
        //                                                     0);
        // }
#endif

        // first create the overlap
        // map from process to global cell indices in overlap
        std::map<int,std::set<int> > overlap;
        auto noImportedOwner = addOverlapLayer(*this, cell_part, exportList, importList, cc, addCornerCells,
                                               transmissibilities);
        // importList contains all the indices that will be here.
        auto compareImport = [](const std::tuple<int,int,char,int>& t1,
                                const std::tuple<int,int,char,int>&t2)
                             {
                                 return std::get<0>(t1) < std::get<0>(t2);
                             };

        if ( ! ownersFirst )
        {
            // merge owner and overlap sorted by global index
            std::inplace_merge(importList.begin(), importList.begin()+noImportedOwner,
                               importList.end(), compareImport);
        }
        // assign local indices
        int localIndex = 0;
        for(auto&& entry: importList)
            std::get<3>(entry) = localIndex++;

        if ( ownersFirst )
        {
            // merge owner and overlap sorted by global index
            std::inplace_merge(importList.begin(), importList.begin()+noImportedOwner,
                               importList.end(), compareImport);
        }

        distributed_data_.reset(new cpgrid::CpGridData(cc));
        distributed_data_->setUniqueBoundaryIds(data_->uniqueBoundaryIds());
        // Just to be sure we assume that only master knows
        cc.broadcast(&distributed_data_->use_unique_boundary_ids_, 1, 0);

        // Create indexset
        distributed_data_->cell_indexset_.beginResize();
        for(const auto& entry: importList)
        {
            distributed_data_->cell_indexset_.add(std::get<0>(entry), ParallelIndexSet::LocalIndex(std::get<3>(entry), AttributeSet(std::get<2>(entry)), true));
        }
        distributed_data_->cell_indexset_.endResize();
        // add an interface for gathering/scattering data with communication
        // forward direction will be scatter and backward gather
        // Interface will communicate from owner to all
        setupSendInterface(exportList, *cell_scatter_gather_interfaces_);
        setupRecvInterface(importList, *cell_scatter_gather_interfaces_);

        distributed_data_->distributeGlobalGrid(*this,*this->current_view_data_, cell_part);


        // Compute the partition type for cell
        int owned_cells = 0;
        int overlap_cells = 0;
        for (const auto i : distributed_data_->cell_indexset_) {
            if (i.local().attribute() == AttributeSet::owner) {
                ++owned_cells;
            } else {
                ++overlap_cells;
            }
        }

        // owned cells
        std::vector<int> cc_owned_cells(cc.size(), 0);
        cc_owned_cells[cc.rank()] = owned_cells;
        cc.sum(cc_owned_cells.data(), cc.size());

        // overlap cells
        std::vector<int> cc_overlap_cells(cc.size(), 0);
        cc_overlap_cells[cc.rank()] = overlap_cells;
        cc.sum(cc_overlap_cells.data(), cc.size());

        // total cells
        const int total_cells = distributed_data_->cell_to_face_.size();
        std::vector<int> cc_total_cells(cc.size(), 0);
        cc_total_cells[cc.rank()] = total_cells;
        cc.sum(cc_total_cells.data(), cc.size());

        if (cc.rank() == 0) {
            std::ostringstream ostr;
            ostr << "\nLoad balancing distributes " << data_->size(0)
                << " active cells on " << cc.size() << " processes as follows:\n"; 
            ostr << "  rank   owned cells   overlap cells   total cells\n";
            ostr << "--------------------------------------------------\n";
            for (int i = 0; i < cc.size(); ++i) {
                ostr << std::setw(6) << i
                    << std::setw(14) << cc_owned_cells[i]
                    << std::setw(16) << cc_overlap_cells[i]
                    << std::setw(14) << cc_total_cells[i] << "\n";
            }
            ostr << "--------------------------------------------------\n";
            ostr << "   sum";
            auto sum_owned = std::accumulate(cc_owned_cells.begin(), cc_owned_cells.end(), 0);
            ostr << std::setw(14) << sum_owned;
            auto sum_overlap = std::accumulate(cc_overlap_cells.begin(), cc_overlap_cells.end(), 0);
            ostr << std::setw(16) << sum_overlap;
            auto sum_total = std::accumulate(cc_total_cells.begin(), cc_total_cells.end(), 0);
            ostr << std::setw(14) << sum_total << "\n";
            Opm::OpmLog::info(ostr.str());
        }
        for (const auto& nc : cc_owned_cells) {
            if (nc == 0) {
                throw std::runtime_error("At least one process has zero cells. Aborting.");
            }
        }

        current_view_data_ = distributed_data_.get();
        return std::make_pair(true, defunct_wells);
    }
    else
    {
        std::cerr << "CpGrid::scatterGrid() only makes sense in a parallel run. "
                  << "This run only uses one process.\n";
        return std::make_pair(false, std::unordered_set<std::string>());
    }
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
        if ( current_view_data_->ccobj_.rank() != 0 )
        {
            grdecl g;
            g.dims[0] = g.dims[1] = g.dims[2] = 0;
            current_view_data_->processEclipseFormat(g, {}, 0.0, false, false);
            // global grid only on rank 0
            return;
        }

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
        current_view_data_->ccobj_.broadcast(current_view_data_->logical_cartesian_size_.data(),
                                             current_view_data_->logical_cartesian_size_.size(),
                                             0);
    }
    void CpGrid::writeSintefLegacyFormat(const std::string& grid_prefix) const
    {
        current_view_data_->writeSintefLegacyFormat(grid_prefix);
        current_view_data_->ccobj_.broadcast(current_view_data_->logical_cartesian_size_.data(),
                                             current_view_data_->logical_cartesian_size_.size(),
                                             0);
    }


#if HAVE_ECL_INPUT
    void CpGrid::processEclipseFormat(const Opm::EclipseGrid* ecl_grid,
                                      bool periodic_extension,
                                      bool turn_normals, bool clip_z,
                                      const std::vector<double>& poreVolume,
                                      const Opm::NNC& nncs)
    {
        current_view_data_->processEclipseFormat(ecl_grid, periodic_extension,
                                                 turn_normals, clip_z,
                                                 poreVolume, nncs);
        current_view_data_->ccobj_.broadcast(current_view_data_->logical_cartesian_size_.data(),
                                             current_view_data_->logical_cartesian_size_.size(),
                                             0);
    }
#endif

    void CpGrid::processEclipseFormat(const grdecl& input_data, double z_tolerance,
                                      bool remove_ij_boundary, bool turn_normals)
    {
        current_view_data_->processEclipseFormat(input_data, {}, z_tolerance, remove_ij_boundary, turn_normals);
        current_view_data_->ccobj_.broadcast(current_view_data_->logical_cartesian_size_.data(),
                                             current_view_data_->logical_cartesian_size_.size(),
                                             0);
    }

} // namespace Dune
