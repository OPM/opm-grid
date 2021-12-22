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
#include <opm/grid/common/ZoltanGraphFunctions.hpp>
#include <opm/grid/common/GridPartitioning.hpp>
#include <opm/grid/common/WellConnections.hpp>

#include <opm/grid/common/CommunicationUtils.hpp>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <tuple>

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
        : data_( new cpgrid::CpGridData()),
          current_view_data_(data_.get()),
          distributed_data_(),
          cell_scatter_gather_interfaces_(new InterfaceMap),
          point_scatter_gather_interfaces_(new InterfaceMap),
          global_id_set_(*current_view_data_)
    {}


    CpGrid::CpGrid(MPIHelper::MPICommunicator comm)
        : data_( new cpgrid::CpGridData(comm)),
          current_view_data_(data_.get()),
          distributed_data_(),
          cell_scatter_gather_interfaces_(new InterfaceMap),
          point_scatter_gather_interfaces_(new InterfaceMap),
          global_id_set_(*current_view_data_)
    {}


std::pair<bool, std::vector<std::pair<std::string,bool> > >
CpGrid::scatterGrid(EdgeWeightMethod method,
                    [[maybe_unused]] bool ownersFirst,
                    const std::vector<cpgrid::OpmWellType> * wells,
                    [[maybe_unused]] bool serialPartitioning,
                    const double* transmissibilities,
                    [[maybe_unused]] bool addCornerCells,
                    int overlapLayers,
                    [[maybe_unused]] bool useZoltan,
                    double zoltanImbalanceTol,
                    [[maybe_unused]] bool allowDistributedWells,
                    [[maybe_unused]] const std::vector<int>& input_cell_part)
{
    // Silence any unused argument warnings that could occur with various configurations.
    static_cast<void>(wells);
    static_cast<void>(transmissibilities);
    static_cast<void>(overlapLayers);
    static_cast<void>(method);
    static_cast<void>(zoltanImbalanceTol);

    if(distributed_data_)
    {
        std::cerr<<"There is already a distributed version of the grid."
                 << " Maybe scatterGrid was called before?"<<std::endl;
        return std::make_pair(false, std::vector<std::pair<std::string,bool> >());
    }

#if HAVE_MPI
    auto& cc = data_->ccobj_;

    if (cc.size() > 1)
    {
        std::vector<int> computedCellPart;
        std::vector<std::pair<std::string,bool>> wells_on_proc;
        std::vector<std::tuple<int,int,char>> exportList;
        std::vector<std::tuple<int,int,char,int>> importList;
        cpgrid::WellConnections wellConnections;

        auto inputNumParts = input_cell_part.size();
        inputNumParts = this->comm().max(inputNumParts);

        if ( inputNumParts > 0 )
        {
            std::vector<int> errors;
            std::vector<std::string> errorMessages =
                { "More parts than MPI Communicator can handle",
                  "Indices of parts need to zero starting",
                  "Indices of parts need to be consecutive",
                  "Only rank 0 should provide partitioning information for each cell"};

            std::set<int> existingParts;

            if (comm().rank() == 0)
            {
                for(const auto& part: input_cell_part)
                {
                    existingParts.insert(part);
                }
                if (*input_cell_part.rbegin() >= comm().size())
                {
                    errors.push_back(0);
                }

                int i = 0;
                if (*existingParts.begin() != i)
                {
                    errors.push_back(1);
                }
                for (const auto& part: existingParts)
                {
                    if (part != i++)
                    {
                        errors.push_back(2);
                        break;
                    }
                }
                if (std::size_t(size(0)) != input_cell_part.size())
                {
                    errors.push_back(3);
                }
            }
            auto size = errors.size();
            comm().broadcast(&size, 1, 0);
            errors.resize(size);

            if (!errors.empty())
            {
                comm().broadcast(errors.data(), size, 0);
                std::string message("Loadbalance: ");
                for ( const auto& e: errors)
                {
                    message.append(errorMessages[e]).append(". ");
                }
                if (comm().rank() == 0)
                {
                    OPM_THROW(std::logic_error, message);
                }
                else
                {
                    OPM_THROW_NOLOG(std::logic_error, message);
                }
            }


            // Partitioning given externally
            std::tie(computedCellPart, wells_on_proc, exportList, importList, wellConnections) =
                cpgrid::createZoltanListsFromParts(*this, wells, nullptr, input_cell_part,
                                                   true);
        }
        else
        {
            if (useZoltan)
            {
#ifdef HAVE_ZOLTAN
                std::tie(computedCellPart, wells_on_proc, exportList, importList, wellConnections)
                    = serialPartitioning
                    ? cpgrid::zoltanSerialGraphPartitionGridOnRoot(*this, wells, transmissibilities, cc, method, 0, zoltanImbalanceTol, allowDistributedWells)
                    : cpgrid::zoltanGraphPartitionGridOnRoot(*this, wells, transmissibilities, cc, method, 0, zoltanImbalanceTol, allowDistributedWells);
#else
                OPM_THROW(std::runtime_error, "Parallel runs depend on ZOLTAN if useZoltan is true. Please install!");
#endif // HAVE_ZOLTAN
            }
            else
            {
                std::tie(computedCellPart, wells_on_proc, exportList, importList, wellConnections) =
                    cpgrid::vanillaPartitionGridOnRoot(*this, wells, transmissibilities, allowDistributedWells);
            }
        }
        comm().barrier();

        // first create the overlap
        // map from process to global cell indices in overlap
        std::map<int,std::set<int> > overlap;
        auto noImportedOwner = addOverlapLayer(*this, computedCellPart, exportList, importList, cc, addCornerCells,
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

        int procsWithZeroCells{};

        if (cc.rank()==0)
        {
            // Print some statistics without communication
            std::vector<int> ownedCells(cc.size(), 0);
            std::vector<int> overlapCells(cc.size(), 0);
            for (const auto& entry: exportList)
            {
                if(std::get<2>(entry) == AttributeSet::owner)
                {
                    ++ownedCells[std::get<1>(entry)];
                }
                else
                {
                    ++overlapCells[std::get<1>(entry)];
                }
            }

            for(const auto& cellsOnProc: ownedCells)
            {
                procsWithZeroCells += (cellsOnProc == 0);
            }
            std::ostringstream ostr;
            ostr << "\nLoad balancing distributes " << data_->size(0)
                 << " active cells on " << cc.size() << " processes as follows:\n";
            ostr << "  rank   owned cells   overlap cells   total cells\n";
            ostr << "--------------------------------------------------\n";
            for (int i = 0; i < cc.size(); ++i) {
                ostr << std::setw(6) << i
                     << std::setw(14) << ownedCells[i]
                     << std::setw(16) << overlapCells[i]
                     << std::setw(14) << ownedCells[i] + overlapCells[i] << "\n";
            }
            ostr << "--------------------------------------------------\n";
            ostr << "   sum";
            auto sumOwned = std::accumulate(ownedCells.begin(), ownedCells.end(), 0);
            ostr << std::setw(14) << sumOwned;
            auto sumOverlap = std::accumulate(overlapCells.begin(), overlapCells.end(), 0);
            ostr << std::setw(16) << sumOverlap;
            ostr << std::setw(14) << (sumOwned + sumOverlap) << "\n";
            Opm::OpmLog::info(ostr.str());
        }

        // Print well distribution
        std::vector<std::pair<int,int> > procWellPairs;

        // range filters would be nice here. so C++20.
        procWellPairs.reserve(std::count_if(std::begin(wells_on_proc),
                                            std::end(wells_on_proc),
                                            [](const std::pair<std::string, bool>& p){ return p.second; }));
        int wellIndex = 0;
        for ( const auto& well: wells_on_proc)
        {
            if ( well.second )
            {
                procWellPairs.emplace_back(cc.rank(), wellIndex);
            }
            ++wellIndex;
        }

        std::tie(procWellPairs, std::ignore) = Opm::gatherv(procWellPairs, cc, 0);

        if (cc.rank() == 0)
        {
            std::sort(std::begin(procWellPairs), std::end(procWellPairs),
                      [](const std::pair<int,int>& p1, const std::pair<int,int>& p2)
                      { return p1.second < p2.second;});
            std::ostringstream ostr;
            ostr << "\nLoad balancing distributed the wells as follows:\n"
                 << "  well name            ranks with perforated cells\n"
                 << "---------------------------------------------------\n";
            auto procWellPair = std::begin(procWellPairs);
            auto endProcWellPair = std::end(procWellPairs);
            int wellIdx = 0;
            for ( const auto& well: wells_on_proc)
            {
                ostr << std::setw(16) << well.first;
                while (procWellPair != endProcWellPair && procWellPair->second < wellIdx)
                {
                    ++procWellPair;
                }
                for ( ; procWellPair != endProcWellPair && procWellPair->second == wellIdx;
                      ++procWellPair)
                {
                    ostr << " "<< std::setw(7) << procWellPair->first;
                }
                ostr << "\n";
                ++wellIdx;
            }
            Opm::OpmLog::info(ostr.str());
        }

        procsWithZeroCells = cc.sum(procsWithZeroCells);

        if (procsWithZeroCells) {
            std::string msg = "At least one process has zero cells. Aborting. \n"
                     " Try decreasing the imbalance tolerance for zoltan with \n"
                     " --zoltan-imbalance-tolerance. The current value is "
                     + std::to_string(zoltanImbalanceTol);
            if (cc.rank()==0)
            {
                OPM_THROW(std::runtime_error, msg );
            }
            else
            {
                OPM_THROW_NOLOG(std::runtime_error, msg);
            }
        }

        distributed_data_.reset(new cpgrid::CpGridData(cc));
        distributed_data_->setUniqueBoundaryIds(data_->uniqueBoundaryIds());
        // Just to be sure we assume that only master knows
        cc.broadcast(&distributed_data_->use_unique_boundary_ids_, 1, 0);

        // Create indexset
        distributed_data_->cellIndexSet().beginResize();
        for(const auto& entry: importList)
        {
            distributed_data_->cellIndexSet().add(std::get<0>(entry), ParallelIndexSet::LocalIndex(std::get<3>(entry), AttributeSet(std::get<2>(entry)), true));
        }
        distributed_data_->cellIndexSet().endResize();
        // add an interface for gathering/scattering data with communication
        // forward direction will be scatter and backward gather
        // Interface will communicate from owner to all
        setupSendInterface(exportList, *cell_scatter_gather_interfaces_);
        setupRecvInterface(importList, *cell_scatter_gather_interfaces_);

        distributed_data_->distributeGlobalGrid(*this,*this->current_view_data_, computedCellPart);
        global_id_set_.insertIdSet(*distributed_data_);


        current_view_data_ = distributed_data_.get();
        return std::make_pair(true, wells_on_proc);
    }
    else
    {
        std::cerr << "CpGrid::scatterGrid() only makes sense in a parallel run. "
                  << "This run only uses one process.\n";
        return std::make_pair(false, std::vector<std::pair<std::string,bool>>());
    }
#else // #if HAVE_MPI
    std::cerr << "CpGrid::scatterGrid() is non-trivial only with "
              << "MPI support and if the target Dune platform is "
              << "sufficiently recent.\n";
    return std::make_pair(false, std::vector<std::pair<std::string,bool>>());
#endif
}


    void CpGrid::createCartesian(const std::array<int, 3>& dims,
                                 const std::array<double, 3>& cellsize)
    {
        if ( current_view_data_->ccobj_.rank() != 0 )
        {
            // global grid only on rank 0
            current_view_data_->ccobj_.broadcast(current_view_data_->logical_cartesian_size_.data(),
                                                 current_view_data_->logical_cartesian_size_.size(),
                                                 0);
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
        using NNCMap = std::set<std::pair<int, int>>;
        using NNCMaps = std::array<NNCMap, 2>;
        NNCMaps nnc;
        current_view_data_->processEclipseFormat(g, nullptr, nnc, false, false, false);
        // global grid only on rank 0
        current_view_data_->ccobj_.broadcast(current_view_data_->logical_cartesian_size_.data(),
                                             current_view_data_->logical_cartesian_size_.size(),
                                             0);
    }

    void CpGrid::readSintefLegacyFormat(const std::string& grid_prefix)
    {
        if ( current_view_data_->ccobj_.rank() == 0 )
        {
            current_view_data_->readSintefLegacyFormat(grid_prefix);
        }
        current_view_data_->ccobj_.broadcast(current_view_data_->logical_cartesian_size_.data(),
                                             current_view_data_->logical_cartesian_size_.size(),
                                             0);
    }
    void CpGrid::writeSintefLegacyFormat(const std::string& grid_prefix) const
    {
        // Only rank 0 has the full data. Use that for writing.
        if ( current_view_data_->ccobj_.rank() == 0 )
        {
            data_->writeSintefLegacyFormat(grid_prefix);
        }
    }


#if HAVE_ECL_INPUT
    std::vector<std::size_t> CpGrid::processEclipseFormat(const Opm::EclipseGrid* ecl_grid,
                                                          Opm::EclipseState* ecl_state,
                                                          bool periodic_extension,
                                                          bool turn_normals, bool clip_z,
                                                          bool pinchActive)
    {
        auto removed_cells = current_view_data_->processEclipseFormat(ecl_grid, ecl_state, periodic_extension,
                                                                      turn_normals, clip_z, pinchActive);
        current_view_data_->ccobj_.broadcast(current_view_data_->logical_cartesian_size_.data(),
                                             current_view_data_->logical_cartesian_size_.size(),
                                             0);
        return removed_cells;
    }

    std::vector<std::size_t> CpGrid::processEclipseFormat(const Opm::EclipseGrid* ecl_grid_ptr,
                                                              Opm::EclipseState* ecl_state,
                                                              bool periodic_extension, bool turn_normals, bool clip_z)
    {
        return processEclipseFormat(ecl_grid_ptr, ecl_state, periodic_extension, turn_normals, clip_z,
                             !ecl_grid_ptr || ecl_grid_ptr->isPinchActive());
    }

#endif

    void CpGrid::processEclipseFormat(const grdecl& input_data,
                                      bool remove_ij_boundary, bool turn_normals)
    {
        using NNCMap = std::set<std::pair<int, int>>;
        using NNCMaps = std::array<NNCMap, 2>;
        NNCMaps nnc;
        current_view_data_->processEclipseFormat(input_data, nullptr, nnc,
                                                 remove_ij_boundary, turn_normals, false);
        current_view_data_->ccobj_.broadcast(current_view_data_->logical_cartesian_size_.data(),
                                             current_view_data_->logical_cartesian_size_.size(),
                                             0);
    }

} // namespace Dune
