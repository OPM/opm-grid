//===========================================================================
//
// File: CpGrid.cpp
//
// Created: Thu Jun  4 12:55:28 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            B�rd Skaflestad     <bard.skaflestad@sintef.no>
//            Antonella Ritorto   <antonella.ritorto@opm-op.com>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2023 Equinor ASA.
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
#include <opm/grid/utility/platform_dependent/reenable_warnings.h>
#endif

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include "../CpGrid.hpp"
#include "LgrHelpers.hpp"
#include "ParentToChildrenCellGlobalIdHandle.hpp"
#include "NestedRefinementUtilities.hpp"
#include <opm/grid/common/MetisPartition.hpp>
#include <opm/grid/common/ZoltanPartition.hpp>
#include <opm/grid/GraphOfGridWrappers.hpp>
//#include <opm/grid/common/ZoltanGraphFunctions.hpp>
#include <opm/grid/common/GridPartitioning.hpp>
#include <opm/grid/cpgrid/refinement/conforming/RefinedDistribution.hpp>
//#include <opm/grid/common/WellConnections.hpp>
#include <opm/grid/common/CommunicationUtils.hpp>

//#include <fstream>
//#include <iostream>
#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <numeric>
#include <string>
#include <unordered_map>
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

/// Release memory resources from CpGrid::InterfaceMap.  Used as custom
/// deleter for std::shared_ptr<InterfaceMap>.
struct FreeInterfaces
{
#if !HAVE_MPI

    /// Release memory resources for InterfaceMap object
    ///
    /// \param[in] interfaces Object for which to release memory resources.
    void operator()([[maybe_unused]] Dune::CpGrid::InterfaceMap* interfaces) const
    {
        // Nothing to do in the sequential case as the CpGrid::InterfaceMap
        // handles interface deletion in its destructor in this case.
    }

#else // HAVE_MPI

    /// Release memory resources for InterfaceMap object
    ///
    /// \param[in] interfaces Object for which to release memory resources.
    void operator()(Dune::CpGrid::InterfaceMap* interfaces) const
    {
        if (interfaces == nullptr) {
            return;
        }

        for (auto& interface : *interfaces) {
            auto& [scatter, gather] = interface.second;
            scatter.free();
            gather.free();
        }
    }

#endif // HAVE_MPI
};
}

namespace Dune
{

CpGrid::CpGrid()
    : distributed_data_(),
      cell_scatter_gather_interfaces_(new InterfaceMap, FreeInterfaces{}),
      point_scatter_gather_interfaces_(new InterfaceMap, FreeInterfaces{}),
      global_id_set_ptr_()
{
    data_.push_back(std::make_shared<cpgrid::CpGridData>(data_));
    current_data_ = &data_;
    global_id_set_ptr_ = std::make_shared<cpgrid::GlobalIdSet>(*(current_data_->back()));
}

CpGrid::CpGrid(MPIHelper::MPICommunicator comm)
    : distributed_data_(),
      cell_scatter_gather_interfaces_(new InterfaceMap, FreeInterfaces{}),
      point_scatter_gather_interfaces_(new InterfaceMap, FreeInterfaces{}),
      global_id_set_ptr_()
{
    data_.push_back(std::make_shared<cpgrid::CpGridData>(comm, data_));
    current_data_ = &data_;
    global_id_set_ptr_ = std::make_shared<cpgrid::GlobalIdSet>(*(current_data_->back()));
}

std::vector<int>
CpGrid::zoltanPartitionWithoutScatter([[maybe_unused]] const std::vector<cpgrid::OpmWellType>* wells,
                                      [[maybe_unused]] const std::unordered_map<std::string, std::set<int>>& possibleFutureConnections,
                                      [[maybe_unused]] const double* transmissibilities,
                                      [[maybe_unused]] const int numParts,
                                      [[maybe_unused]] const double zoltanImbalanceTol) const
{
#if HAVE_MPI && HAVE_ZOLTAN
    const auto met = EdgeWeightMethod(1);

    return cpgrid::zoltanGraphPartitionGridForJac(*this, wells, possibleFutureConnections, transmissibilities,
                                                  this->data_[0]->ccobj_, met,
                                                  0, numParts, zoltanImbalanceTol);
#else
    return std::vector<int>(this->numCells(), 0);
#endif
}

std::pair<bool, std::vector<std::pair<std::string,bool> > >
CpGrid::scatterGrid(EdgeWeightMethod method,
                    [[maybe_unused]] bool ownersFirst,
                    const std::vector<cpgrid::OpmWellType> * wells,
                    [[maybe_unused]] const std::unordered_map<std::string, std::set<int>>& possibleFutureConnections,
                    [[maybe_unused]] bool serialPartitioning,
                    const double* transmissibilities,
                    [[maybe_unused]] bool addCornerCells,
                    int overlapLayers,
                    [[maybe_unused]] int partitionMethod,
                    double imbalanceTol,
                    [[maybe_unused]] bool allowDistributedWells,
                    [[maybe_unused]] const std::vector<int>& input_cell_part,
                    int level,
                    [[maybe_unused]] bool useTransToFilterOverlap)
{
    // Silence any unused argument warnings that could occur with various configurations.
    static_cast<void>(wells);
    static_cast<void>(transmissibilities);
    static_cast<void>(overlapLayers);
    static_cast<void>(method);
    static_cast<void>(imbalanceTol);
    static_cast<void>(level);

    if(!distributed_data_.empty())
    {
        std::cerr<<"There is already a distributed version of the grid."
                 << " Maybe scatterGrid was called before?"<<std::endl;
        return std::make_pair(false, std::vector<std::pair<std::string,bool> >());
    }

#if HAVE_MPI
    bool validLevel = (level>-1) && (level <= maxLevel());
    // If level == -1, leaf grid view should be distributed (with/without LGRs).
    // - without LGRs: leaf grid view coincides with level zero grid. Supported.
    // - with LGRs: not supported yet. Throw in that case.
    int selectedLevel = validLevel? level : 0;
    // Conforming refine-before-redistribute: partition level zero, then distribute the refined leaf.
    const bool refinedLeaf = (lgr_backend_ == Opm::Refinement::Backend::Conforming) && (maxLevel() > 0) && (level == -1);
    if (refinedLeaf) {
        level = 0;
    }
    if (validLevel && (level>0)) {
        if (comm().rank() == 0) {
            OPM_THROW(std::logic_error, "Loadbalancing a refined level grid is not supported, yet.");
        }
        else {
            OPM_THROW_NOLOG(std::logic_error, "Loadbalancing a refined level grid is not supported, yet.");
        }
    }

    if ( (maxLevel()>0) && (level==-1) ) {
        if (comm().rank() == 0) {
            OPM_THROW(std::logic_error, "Loadbalancing a leaf grid view with local refinement is not supported, yet.");
        }
        else {
            OPM_THROW_NOLOG(std::logic_error, "Loadbalancing a leaf grid view with local refinement is not supported, yet.");
        }
    }

    if ((maxLevel()>0) && (partitionMethod!= Dune::PartitionMethod::zoltanGoG) && !refinedLeaf) {
        if (comm().rank() == 0) {
            OPM_THROW(std::logic_error, "Loadbalancing level zero grid of a grid with local refinement is supported for ZOLTANGOG.");
        }
        else {
            OPM_THROW_NOLOG(std::logic_error, "Loadbalancing level zero grid of a grid with local refinement is supported for ZOLTANGOG.");
        }
    }

    auto& cc = data_[selectedLevel]->ccobj_;

    if (cc.size() > 1)
    {
        // The grid of a dual-continuum run is partitioned with its matrix cells and fracture
        // twins joined, as the cells of a well, and with Zoltan on the graph of the grid.
        // The graph of the wells gets the pairs, also if there are no wells.
        // The offset is known on the root only, where the global grid is.
        const int continuumOffset = cc.max(static_cast<int>(data_[selectedLevel]->continuum_cart_offset_));
        const bool dualContinuum = continuumOffset > 0;
        const std::vector<cpgrid::OpmWellType> noWells;
        if (dualContinuum) {
            if (wells == nullptr) {
                wells = &noWells;
            }
            if (partitionMethod == Dune::PartitionMethod::zoltanGoG) {
                partitionMethod = Dune::PartitionMethod::zoltan;
            }
        }

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
                if (std::ranges::any_of(existingParts,
                                        [&i](const auto& part)
                                        { return part != i++; }))
                {
                    errors.push_back(2);
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
                cpgrid::createListsFromParts(*this, wells, possibleFutureConnections, /* transmissibilities = */ nullptr, input_cell_part,
                                              /* allowDistributedWells = */ true, /* gridAndWells = */ nullptr, level);
        }
        else
        {
            if (partitionMethod == Dune::PartitionMethod::zoltan)
            {
#ifdef HAVE_ZOLTAN
                std::tie(computedCellPart, wells_on_proc, exportList, importList, wellConnections)
                    = serialPartitioning
                    ? cpgrid::zoltanSerialGraphPartitionGridOnRoot(*this, wells, possibleFutureConnections, transmissibilities, cc, method, 0, imbalanceTol, allowDistributedWells, partitioningParams)
                    : cpgrid::zoltanGraphPartitionGridOnRoot(*this, wells, possibleFutureConnections, transmissibilities, cc, method, 0, imbalanceTol, allowDistributedWells, partitioningParams);
#else
                OPM_THROW(std::runtime_error, "Parallel runs depend on ZOLTAN if useZoltan is true. Please install!");
#endif // HAVE_ZOLTAN
            }
            else if (partitionMethod == Dune::PartitionMethod::metis)
            {
#ifdef HAVE_METIS
                if (!serialPartitioning)
                    OPM_MESSAGE("Warning: Serial partitioning is set to false and METIS was selected to partition the grid, but METIS is a serial partitioner. Continuing with serial partitioning...");
                std::tie(computedCellPart, wells_on_proc, exportList, importList, wellConnections) = cpgrid::metisSerialGraphPartitionGridOnRoot(*this, wells, possibleFutureConnections, transmissibilities, cc, method, 0, imbalanceTol, allowDistributedWells, partitioningParams);
#else
                OPM_THROW(std::runtime_error, "Parallel runs depend on METIS if useMetis is true. Please install!");
#endif // HAVE_METIS

            }
            else if (partitionMethod == Dune::PartitionMethod::zoltanGoG)
            {
#ifdef HAVE_ZOLTAN
                std::tie(computedCellPart, wells_on_proc, exportList, importList, wellConnections)
                    = serialPartitioning
                    ? Opm::zoltanSerialPartitioningWithGraphOfGrid(*this, wells, possibleFutureConnections, transmissibilities, cc, method, 0, imbalanceTol, allowDistributedWells, partitioningParams)
                    : Opm::zoltanPartitioningWithGraphOfGrid(*this, wells, possibleFutureConnections, transmissibilities, cc, method, 0, imbalanceTol, allowDistributedWells, partitioningParams, level);
#else
                OPM_THROW(std::runtime_error, "Parallel runs depend on ZOLTAN if useZoltan is true. Please install!");
#endif // HAVE_ZOLTAN
            }
            else
            {
                std::tie(computedCellPart, wells_on_proc, exportList, importList, wellConnections) =
                    cpgrid::vanillaPartitionGridOnRoot(*this, wells, possibleFutureConnections, transmissibilities, allowDistributedWells);
            }
        }
        // A dual-continuum grid keeps each matrix cell on the process of its fracture twin: the
        // coupling between them is a connection like any other, and the cells have to be
        // neighbours of each other for the transmissibilities to be computed.
        if (dualContinuum) {
            std::vector<int> twinParts;
            if (cc.rank() == 0) {
                const auto offset = data_[selectedLevel]->continuum_cart_offset_;
                const auto& globalCell = data_[selectedLevel]->global_cell_;
                const auto localIndex = [&globalCell](const std::size_t cartesian)
                {
                    const auto it = std::lower_bound(globalCell.begin(), globalCell.end(),
                                                     static_cast<int>(cartesian));
                    return (it != globalCell.end() && *it == static_cast<int>(cartesian))
                        ? static_cast<int>(it - globalCell.begin())
                        : -1;
                };

                // Cells that have to be on one process form a group: a matrix cell and its
                // twin, and, unless wells may be distributed, the cells of a well. A group
                // that the partitioner has split, as a twin pair of two different wells can
                // be, is moved as a whole, so that no well is split by the pairing. It goes
                // to the process of the first cell of a well in it, or else of its first
                // matrix cell.
                std::vector<int> group(globalCell.size());
                std::iota(group.begin(), group.end(), 0);
                const auto find = [&group](int c) {
                    while (group[c] != c) {
                        group[c] = group[group[c]];
                        c = group[c];
                    }
                    return c;
                };
                const auto join = [&group, &find](const int c1, const int c2) {
                    const auto r1 = find(c1);
                    const auto r2 = find(c2);
                    // The smaller index is the root, so that it is the first cell of the group.
                    group[std::max(r1, r2)] = std::min(r1, r2);
                };

                std::vector<char> wellCell(globalCell.size(), 0);
                std::unordered_map<std::string, int> firstCellOfWell;
                const auto addToWell = [&](const std::string& name, const std::size_t cartesian) {
                    const auto c = localIndex(cartesian);
                    if (c < 0) {
                        return;
                    }
                    wellCell[c] = 1;
                    if (allowDistributedWells) {
                        return;
                    }
                    if (const auto it = firstCellOfWell.find(name); it != firstCellOfWell.end()) {
                        join(it->second, c);
                    } else {
                        firstCellOfWell.emplace(name, c);
                    }
                };
                if (wells != nullptr) {
                    for (const auto& well : *wells) {
                        for (const auto& connection : well.getConnections()) {
                            addToWell(well.name(), connection.global_index());
                        }
                    }
                }
                for (const auto& future : possibleFutureConnections) {
                    for (const auto cartesian : future.second) {
                        addToWell(future.first, cartesian);
                    }
                }

                for (std::size_t c = 0; c < globalCell.size(); ++c) {
                    if (static_cast<std::size_t>(globalCell[c]) >= offset) {
                        break;      // the matrix cells come first
                    }
                    if (const auto twin = localIndex(globalCell[c] + offset); twin >= 0) {
                        join(c, twin);
                    }
                }

                // The part of a group: that of its first well cell, otherwise of its first cell.
                std::vector<int> groupPart(globalCell.size(), -1);
                for (std::size_t c = 0; c < globalCell.size(); ++c) {
                    if (wellCell[c] && groupPart[find(c)] < 0) {
                        groupPart[find(c)] = computedCellPart[c];
                    }
                }
                twinParts = computedCellPart;
                for (std::size_t c = 0; c < globalCell.size(); ++c) {
                    const auto root = find(c);
                    twinParts[c] = (groupPart[root] >= 0) ? groupPart[root] : computedCellPart[root];
                }
            }

            std::tie(computedCellPart, wells_on_proc, exportList, importList, wellConnections) =
                cpgrid::createListsFromParts(*this, wells, possibleFutureConnections, transmissibilities,
                                             twinParts, (inputNumParts > 0) || allowDistributedWells,
                                             /* gridAndWells = */ nullptr, level);
        }

        comm().barrier();

        if (refinedLeaf) {
            const auto leafPart = Opm::Refinement::leafPartitionFromLevelZero(*this, computedCellPart);
            Opm::Refinement::prepareLeafForScatter(*data_.back(), leafPart, cc);
            std::tie(computedCellPart, wells_on_proc, exportList, importList, wellConnections) =
                cpgrid::createListsFromParts(*this, wells, possibleFutureConnections, /* transmissibilities = */ nullptr,
                                             leafPart, allowDistributedWells, /* gridAndWells = */ nullptr, /* level = */ -1);
            computedCellPart = leafPart;
            selectedLevel = static_cast<int>(data_.size()) - 1;
            level = -1;
        }

        // first create the overlap
        auto noImportedOwner = addOverlapLayer(*this,
                                               computedCellPart,
                                               exportList,
                                               importList,
                                               cc,
                                               addCornerCells,
                                               transmissibilities,
                                               useTransToFilterOverlap,
                                               1 /*layers*/,
                                               level);
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

            procsWithZeroCells =
                    std::accumulate(ownedCells.begin(), ownedCells.end(), 0,
                                    [](const auto acc, const auto cellsOnProc)
                                    { return acc + (cellsOnProc == 0); });
            std::ostringstream ostr;
            ostr << "\nLoad balancing distributes level " << selectedLevel << " with " << data_[selectedLevel]->size(0)
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
        procWellPairs.reserve(std::ranges::count_if(wells_on_proc,
                                                    [](const auto& p)
                                                    { return p.second; }));
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
            std::ranges::sort(procWellPairs,
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
                " Try decreasing the imbalance tolerance with the argument \n"
                " --imbalance-tol. The current value is "
                + std::to_string(imbalanceTol);
            if (cc.rank()==0)
            {
                OPM_THROW(std::runtime_error, msg );
            }
            else
            {
                OPM_THROW_NOLOG(std::runtime_error, msg);
            }
        }

        // distributed_data should be empty at this point.
        distributed_data_.push_back(std::make_shared<cpgrid::CpGridData>(cc, distributed_data_));
        distributed_data_[0]->setUniqueBoundaryIds(data_[selectedLevel]->uniqueBoundaryIds());

        // Just to be sure we assume that only master knows
        cc.broadcast(&distributed_data_[0]->use_unique_boundary_ids_, 1, 0);

        // Create indexset
        distributed_data_[0]->cellIndexSet().beginResize();
        for(const auto& entry: importList)
        {
            distributed_data_[0]->cellIndexSet()
                .add(std::get<0>(entry),ParallelIndexSet::LocalIndex(std::get<3>(entry),AttributeSet(std::get<2>(entry)), true));
        }
        distributed_data_[0]->cellIndexSet().endResize();
        // add an interface for gathering/scattering data with communication
        // forward direction will be scatter and backward gather
        // Interface will communicate from owner to all
        setupSendInterface(exportList, *cell_scatter_gather_interfaces_);
        setupRecvInterface(importList, *cell_scatter_gather_interfaces_);

        distributed_data_[0]->distributeGlobalGrid(*this,*this->data_[selectedLevel], computedCellPart);
        distributed_data_[0]->continuum_cart_offset_ = continuumOffset;
        (*global_id_set_ptr_).insertIdSet(*distributed_data_[0]);
        distributed_data_[0]-> index_set_.reset(new cpgrid::IndexSet(distributed_data_[0]->cell_to_face_.size(),
                                                                     distributed_data_[0]-> geomVector<3>().size()));

        current_data_ = &distributed_data_;
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
                             const std::array<double, 3>& cellsize,
                             const std::array<int, 3>& shift)
{
    if ( current_data_->back()->ccobj_.rank() != 0 )
    {
        // global grid only on rank 0
        current_data_->back()->ccobj_.broadcast(current_data_->back()->logical_cartesian_size_.data(),
                                                current_data_->back()->logical_cartesian_size_.size(),
                                                0);
        return;
    }

    // Make the grdecl format arrays.
    // Pillar coords.
    std::vector<double> coord;
    coord.reserve(6*(dims[0] + 1)*(dims[1] + 1));
    double bot = 0.0+shift[2]*cellsize[2];
    double top = (dims[2]+shift[2])*cellsize[2];
    // i runs fastest for the pillars.
    for (int j = 0; j < dims[1] + 1; ++j) {
        double y = (j+shift[1])*cellsize[1];
        for (int i = 0; i < dims[0] + 1; ++i) {
            double x = (i+shift[0])*cellsize[0];
            double pillar[6] = { x, y, bot, x, y, top };
            coord.insert(coord.end(), pillar, pillar + 6);
        }
    }
    std::vector<double> zcorn(8*dims[0]*dims[1]*dims[2]);
    const int num_per_layer = 4*dims[0]*dims[1];
    double* offset = &zcorn[0];
    for (int k = 0; k < dims[2]; ++k) {
        double zlow = (k+shift[2])*cellsize[2];
        std::fill_n(offset, num_per_layer, zlow);
        offset += num_per_layer;
        double zhigh = (k+1+shift[2])*cellsize[2];
        std::fill_n(offset, num_per_layer, zhigh);
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

    // Note: This is a Cartesian, matching grid which is edge-conforming
    // regardless of the edge_conformal flag.
    current_data_->back()->processEclipseFormat(g,
                                                /* ecl_state = */ nullptr,
                                                nnc,
                                                /* remove_ij_boundary = */ false,
                                                /* turn_normals = */ false,
                                                /* pinchActive = */ false,
                                                /* tolerance_unique_ponts = */ 0.0,
                                                /* edge_conformal = */ false);

    // global grid only on rank 0
    current_data_->back()->ccobj_.broadcast(current_data_->back()->logical_cartesian_size_.data(),
                                            current_data_->back()->logical_cartesian_size_.size(),
                                            0);
}

const std::array<int, 3>& CpGrid::logicalCartesianSize() const
{
    // Temporary. For a grid with LGRs, we set the logical cartesian size of the LeafGridView as the one for level 0.
    //            Goal: CartesianIndexMapper well-defined for CpGrid LeafView with LGRs.
    return current_data_->front() -> logical_cartesian_size_;
}

const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& CpGrid::currentData() const
{
    return *current_data_;
}

std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& CpGrid::currentData()
{
    return *current_data_;
}
const Dune::cpgrid::CpGridData& CpGrid::currentLeafData() const
{
    return *current_data_->back();
}

Dune::cpgrid::CpGridData& CpGrid::currentLeafData()
{
    return *current_data_->back();
}

const std::vector<int>& CpGrid::globalCell() const
{
    // Temporary. For a grid with LGRs, we set the globalCell() of the as the one for level 0.
    //            Goal: CartesianIndexMapper well-defined for CpGrid LeafView with LGRs.
    return currentLeafData().global_cell_;
}

void CpGrid::getIJK(const int c, std::array<int,3>& ijk) const
{
    current_data_->back()->getIJK(c, ijk);
}

bool CpGrid::uniqueBoundaryIds() const
{
    return current_data_->back()->uniqueBoundaryIds();
}

void CpGrid::setUniqueBoundaryIds(bool uids)
{
    current_data_->back()->setUniqueBoundaryIds(uids);
}

std::string CpGrid::name() const
{
    return "CpGrid";
}

int CpGrid::maxLevel() const
{
    if (currentData().size() == 1){
        return 0; // "GLOBAL" grid is the unique one
    }
    else {  // There are multiple LGRs
        return this -> currentData().size() - 2; // last entry is leafView, and it starts in level 0 = GLOBAL grid.
    }
}

template<int codim>
typename CpGridTraits::template Codim<codim>::LevelIterator CpGrid::lbegin (int level) const
{
    if (level<0 || level>maxLevel())
        DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
    return  cpgrid::Iterator<codim, All_Partition>( *(*current_data_)[level], 0, true);
}
template typename CpGridTraits::template Codim<0>::LevelIterator CpGrid::lbegin<0>(int) const;
template typename CpGridTraits::template Codim<1>::LevelIterator CpGrid::lbegin<1>(int) const;
template typename CpGridTraits::template Codim<3>::LevelIterator CpGrid::lbegin<3>(int) const;

template<int codim>
typename CpGridTraits::template Codim<codim>::LevelIterator CpGrid::lend (int level) const
{
    if (level<0 || level>maxLevel())
        DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
    return  cpgrid::Iterator<codim, All_Partition>( *(*current_data_)[level], size(level, codim), true);
}
template typename CpGridTraits::template Codim<0>::LevelIterator CpGrid::lend<0>(int) const;
template typename CpGridTraits::template Codim<1>::LevelIterator CpGrid::lend<1>(int) const;
template typename CpGridTraits::template Codim<3>::LevelIterator CpGrid::lend<3>(int) const;

template<int codim>
typename CpGridTraits::template Codim<codim>::LeafIterator CpGrid::leafbegin() const
{
    return cpgrid::Iterator<codim, All_Partition>(*(current_data_->back()), 0, true);
}
template typename CpGridTraits::template Codim<0>::LeafIterator CpGrid::leafbegin<0>() const;
template typename CpGridTraits::template Codim<1>::LeafIterator CpGrid::leafbegin<1>() const;
template typename CpGridTraits::template Codim<3>::LeafIterator CpGrid::leafbegin<3>() const;

template<int codim>
typename CpGridTraits::template Codim<codim>::LeafIterator CpGrid::leafend() const
{
    return cpgrid::Iterator<codim, All_Partition>(*(current_data_->back()), size(codim), true);
}
template typename CpGridTraits::template Codim<0>::LeafIterator CpGrid::leafend<0>() const;
template typename CpGridTraits::template Codim<1>::LeafIterator CpGrid::leafend<1>() const;
template typename CpGridTraits::template Codim<3>::LeafIterator CpGrid::leafend<3>() const;

template<int codim, PartitionIteratorType PiType>
typename CpGridTraits::template Codim<codim>::template Partition<PiType>::LevelIterator CpGrid::lbegin (int level) const
{
    if (level<0 || level>maxLevel())
        DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
    return  cpgrid::Iterator<codim, PiType>( *(*current_data_)[level], 0, true);
}
template typename CpGridTraits::template Codim<0>::template Partition<Dune::Interior_Partition>::LevelIterator
CpGrid::lbegin<0,Dune::Interior_Partition>(int) const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::InteriorBorder_Partition>::LevelIterator
CpGrid::lbegin<0,Dune::InteriorBorder_Partition>(int) const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::Overlap_Partition>::LevelIterator
CpGrid::lbegin<0,Dune::Overlap_Partition>(int) const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::OverlapFront_Partition>::LevelIterator
CpGrid::lbegin<0,Dune::OverlapFront_Partition>(int) const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::All_Partition>::LevelIterator
CpGrid::lbegin<0,Dune::All_Partition>(int) const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::Ghost_Partition>::LevelIterator
CpGrid::lbegin<0,Dune::Ghost_Partition>(int) const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::Interior_Partition>::LevelIterator
CpGrid::lbegin<1,Dune::Interior_Partition>(int) const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::InteriorBorder_Partition>::LevelIterator
CpGrid::lbegin<1,Dune::InteriorBorder_Partition>(int) const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::Overlap_Partition>::LevelIterator
CpGrid::lbegin<1,Dune::Overlap_Partition>(int) const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::OverlapFront_Partition>::LevelIterator
CpGrid::lbegin<1,Dune::OverlapFront_Partition>(int) const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::All_Partition>::LevelIterator
CpGrid::lbegin<1,Dune::All_Partition>(int) const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::Ghost_Partition>::LevelIterator
CpGrid::lbegin<1,Dune::Ghost_Partition>(int) const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::Interior_Partition>::LevelIterator
CpGrid::lbegin<3,Dune::Interior_Partition>(int) const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::InteriorBorder_Partition>::LevelIterator
CpGrid::lbegin<3,Dune::InteriorBorder_Partition>(int) const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::Overlap_Partition>::LevelIterator
CpGrid::lbegin<3,Dune::Overlap_Partition>(int) const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::OverlapFront_Partition>::LevelIterator
CpGrid::lbegin<3,Dune::OverlapFront_Partition>(int) const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::All_Partition>::LevelIterator
CpGrid::lbegin<3,Dune::All_Partition>(int) const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::Ghost_Partition>::LevelIterator
CpGrid::lbegin<3,Dune::Ghost_Partition>(int) const;

template<int codim, PartitionIteratorType PiType>
typename CpGridTraits::template Codim<codim>::template Partition<PiType>::LevelIterator CpGrid::lend (int level) const
{
    if (level<0 || level>maxLevel())
        DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
    return  cpgrid::Iterator<codim, PiType>( *(*current_data_)[level], size(level, codim), true);
}
template typename CpGridTraits::template Codim<0>::template Partition<Dune::Interior_Partition>::LevelIterator
CpGrid::lend<0,Dune::Interior_Partition>(int) const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::InteriorBorder_Partition>::LevelIterator
CpGrid::lend<0,Dune::InteriorBorder_Partition>(int) const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::Overlap_Partition>::LevelIterator
CpGrid::lend<0,Dune::Overlap_Partition>(int) const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::OverlapFront_Partition>::LevelIterator
CpGrid::lend<0,Dune::OverlapFront_Partition>(int) const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::All_Partition>::LevelIterator
CpGrid::lend<0,Dune::All_Partition>(int) const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::Ghost_Partition>::LevelIterator
CpGrid::lend<0,Dune::Ghost_Partition>(int) const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::Interior_Partition>::LevelIterator
CpGrid::lend<1,Dune::Interior_Partition>(int) const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::InteriorBorder_Partition>::LevelIterator
CpGrid::lend<1,Dune::InteriorBorder_Partition>(int) const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::Overlap_Partition>::LevelIterator
CpGrid::lend<1,Dune::Overlap_Partition>(int) const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::OverlapFront_Partition>::LevelIterator
CpGrid::lend<1,Dune::OverlapFront_Partition>(int) const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::All_Partition>::LevelIterator
CpGrid::lend<1,Dune::All_Partition>(int) const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::Ghost_Partition>::LevelIterator
CpGrid::lend<1,Dune::Ghost_Partition>(int) const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::Interior_Partition>::LevelIterator
CpGrid::lend<3,Dune::Interior_Partition>(int) const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::InteriorBorder_Partition>::LevelIterator
CpGrid::lend<3,Dune::InteriorBorder_Partition>(int) const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::Overlap_Partition>::LevelIterator
CpGrid::lend<3,Dune::Overlap_Partition>(int) const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::OverlapFront_Partition>::LevelIterator
CpGrid::lend<3,Dune::OverlapFront_Partition>(int) const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::All_Partition>::LevelIterator
CpGrid::lend<3,Dune::All_Partition>(int) const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::Ghost_Partition>::LevelIterator
CpGrid::lend<3,Dune::Ghost_Partition>(int) const;

template<int codim, PartitionIteratorType PiType>
typename CpGridFamily::Traits::template Codim<codim>::template Partition<PiType>::LeafIterator CpGrid::leafbegin() const
{
    return cpgrid::Iterator<codim, PiType>(*(current_data_->back()), 0, true);
}
template typename CpGridTraits::template Codim<0>::template Partition<Dune::Interior_Partition>::LeafIterator
CpGrid::leafbegin<0,Dune::Interior_Partition>() const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::InteriorBorder_Partition>::LeafIterator
CpGrid::leafbegin<0,Dune::InteriorBorder_Partition>() const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::Overlap_Partition>::LeafIterator
CpGrid::leafbegin<0,Dune::Overlap_Partition>() const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::OverlapFront_Partition>::LeafIterator
CpGrid::leafbegin<0,Dune::OverlapFront_Partition>() const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::All_Partition>::LeafIterator
CpGrid::leafbegin<0,Dune::All_Partition>() const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::Ghost_Partition>::LeafIterator
CpGrid::leafbegin<0,Dune::Ghost_Partition>() const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::Interior_Partition>::LeafIterator
CpGrid::leafbegin<1,Dune::Interior_Partition>() const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::InteriorBorder_Partition>::LeafIterator
CpGrid::leafbegin<1,Dune::InteriorBorder_Partition>() const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::Overlap_Partition>::LeafIterator
CpGrid::leafbegin<1,Dune::Overlap_Partition>() const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::OverlapFront_Partition>::LeafIterator
CpGrid::leafbegin<1,Dune::OverlapFront_Partition>() const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::All_Partition>::LeafIterator
CpGrid::leafbegin<1,Dune::All_Partition>() const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::Ghost_Partition>::LeafIterator
CpGrid::leafbegin<1,Dune::Ghost_Partition>() const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::Interior_Partition>::LeafIterator
CpGrid::leafbegin<3,Dune::Interior_Partition>() const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::InteriorBorder_Partition>::LeafIterator
CpGrid::leafbegin<3,Dune::InteriorBorder_Partition>() const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::Overlap_Partition>::LeafIterator
CpGrid::leafbegin<3,Dune::Overlap_Partition>() const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::OverlapFront_Partition>::LeafIterator
CpGrid::leafbegin<3,Dune::OverlapFront_Partition>() const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::All_Partition>::LeafIterator
CpGrid::leafbegin<3,Dune::All_Partition>() const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::Ghost_Partition>::LeafIterator
CpGrid::leafbegin<3,Dune::Ghost_Partition>() const;

template<int codim, PartitionIteratorType PiType>
typename CpGridFamily::Traits::template Codim<codim>::template Partition<PiType>::LeafIterator CpGrid::leafend() const
{
    return cpgrid::Iterator<codim, PiType>(*(current_data_->back()), size(codim), true);
}
template typename CpGridTraits::template Codim<0>::template Partition<Dune::Interior_Partition>::LeafIterator
CpGrid::leafend<0,Dune::Interior_Partition>() const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::InteriorBorder_Partition>::LeafIterator
CpGrid::leafend<0,Dune::InteriorBorder_Partition>() const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::Overlap_Partition>::LeafIterator
CpGrid::leafend<0,Dune::Overlap_Partition>() const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::OverlapFront_Partition>::LeafIterator
CpGrid::leafend<0,Dune::OverlapFront_Partition>() const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::All_Partition>::LeafIterator
CpGrid::leafend<0,Dune::All_Partition>() const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::Ghost_Partition>::LeafIterator
CpGrid::leafend<0,Dune::Ghost_Partition>() const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::Interior_Partition>::LeafIterator
CpGrid::leafend<1,Dune::Interior_Partition>() const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::InteriorBorder_Partition>::LeafIterator
CpGrid::leafend<1,Dune::InteriorBorder_Partition>() const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::Overlap_Partition>::LeafIterator
CpGrid::leafend<1,Dune::Overlap_Partition>() const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::OverlapFront_Partition>::LeafIterator
CpGrid::leafend<1,Dune::OverlapFront_Partition>() const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::All_Partition>::LeafIterator
CpGrid::leafend<1,Dune::All_Partition>() const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::Ghost_Partition>::LeafIterator
CpGrid::leafend<1,Dune::Ghost_Partition>() const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::Interior_Partition>::LeafIterator
CpGrid::leafend<3,Dune::Interior_Partition>() const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::InteriorBorder_Partition>::LeafIterator
CpGrid::leafend<3,Dune::InteriorBorder_Partition>() const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::Overlap_Partition>::LeafIterator
CpGrid::leafend<3,Dune::Overlap_Partition>() const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::OverlapFront_Partition>::LeafIterator
CpGrid::leafend<3,Dune::OverlapFront_Partition>() const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::All_Partition>::LeafIterator
CpGrid::leafend<3,Dune::All_Partition>() const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::Ghost_Partition>::LeafIterator
CpGrid::leafend<3,Dune::Ghost_Partition>() const;

int CpGrid::size (int level, int codim) const
{
    if (level<0 || level>maxLevel())
        DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
    return currentData()[level]->size(codim);
}

int CpGrid::size (int codim) const
{
    return currentLeafData().size(codim);
}

int CpGrid::size (int level, GeometryType type) const
{
    if (level<0 || level>maxLevel())
        DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
    return currentData()[level]->size(type);
}

int CpGrid::size (GeometryType type) const
{
    return currentLeafData().size(type);
}

const CpGridFamily::Traits::GlobalIdSet& CpGrid::globalIdSet() const
{
    return  *global_id_set_ptr_;
}

const CpGridFamily::Traits::LocalIdSet& CpGrid::localIdSet() const
{
    return *global_id_set_ptr_;
}

const CpGridFamily::Traits::LevelIndexSet& CpGrid::levelIndexSet(int level) const
{
    if (level<0 || level>maxLevel())
        DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
    return *currentData()[level]->index_set_;
}

const CpGridFamily::Traits::LeafIndexSet& CpGrid::leafIndexSet() const
{
    return *current_data_->back()->index_set_;
}

const std::vector< Dune :: GeometryType >& CpGrid::geomTypes( const int codim ) const
{
    return leafIndexSet().geomTypes( codim );
}

template <int codim>
cpgrid::Entity<codim> CpGrid::entity( const cpgrid::Entity< codim >& seed ) const
{
    return cpgrid::Entity<codim>( *(this->current_data_->back()), seed );
}

template cpgrid::Entity<0> CpGrid::entity<0>( const cpgrid::Entity<0>&) const;
template cpgrid::Entity<3> CpGrid::entity<3>( const cpgrid::Entity<3>&) const;

/// \brief Size of the overlap on the leaf level
unsigned int CpGrid::overlapSize(int) const {
    return 1;
}

/// \brief Size of the ghost cell layer on the leaf level
unsigned int CpGrid::ghostSize(int) const {
    return 0;
}

/// \brief Size of the overlap on a given level
unsigned int CpGrid::overlapSize(int, int) const {
    return 1;
}

/// \brief Size of the ghost cell layer on a given level
unsigned int CpGrid::ghostSize(int, int) const {
    return 0;
}

unsigned int CpGrid::numBoundarySegments() const
{
    if( uniqueBoundaryIds() )
    {
        return current_data_->back()->unique_boundary_ids_.size();
    }
    else
    {
        unsigned int numBndSegs = 0;
        const int num_faces = numFaces();
        for (int i = 0; i < num_faces; ++i) {
            cpgrid::EntityRep<1> face(i, true);
            if (current_data_->back()->face_to_cell_[face].size() == 1) {
                ++numBndSegs;
            }
        }
        return numBndSegs;
    }
}

void CpGrid::setPartitioningParams(const std::map<std::string,std::string>& params)
{
    partitioningParams = params;
}

const typename CpGridTraits::Communication& Dune::CpGrid::comm () const
{
    return current_data_->back()->ccobj_;
}

const typename CpGridTraits::Communication& Dune::CpGrid::getLevelZeroComm() const
{
    return current_data_->front()->ccobj_;
}

//

const std::vector<double>& CpGrid::zcornData() const {
    return current_data_->back()->zcornData();
}

int CpGrid::numCells(int level) const
{
    bool validLevel = (level>-1) && (level<= maxLevel());
    return validLevel? data_[level]->cell_to_face_.size() : current_data_->back()->cell_to_face_.size();
}
/// \brief Get the number of faces.
int CpGrid::numFaces(int level) const
{
    bool validLevel = (level>-1) && (level<= maxLevel());
    return validLevel? data_[level]->face_to_cell_.size() : current_data_->back()->face_to_cell_.size();
}
/// \brief Get The number of vertices.
int CpGrid::numVertices() const
{
    return current_data_->back()->geomVector<3>().size();
}

int CpGrid::numCellFaces(int cell, int level) const
{
    bool validLevel = (level>-1) && (level<= maxLevel());
    return validLevel? data_[level]->cell_to_face_[cpgrid::EntityRep<0>(cell, true)].size()
        : current_data_->back()->cell_to_face_[cpgrid::EntityRep<0>(cell, true)].size();
}

int CpGrid::cellFace(int cell, int local_index, int level) const
{
    bool validLevel = (level>-1) && (level<= maxLevel());
    return validLevel? data_[level]-> cell_to_face_[cpgrid::EntityRep<0>(cell, true)][local_index].index()
        : current_data_->back()->cell_to_face_[cpgrid::EntityRep<0>(cell, true)][local_index].index();
}

const cpgrid::OrientedEntityTable<0,1>::row_type CpGrid::cellFaceRow(int cell) const
{
    return current_data_->back()->cell_to_face_[cpgrid::EntityRep<0>(cell, true)];
}

int CpGrid::faceCell(int face, int local_index, int level) const
{
    // Get the correct neigbour relation (row from face->cell table).
    const bool validLevel = (level > -1) && (level <= maxLevel());
    const cpgrid::OrientedEntityTable<1,0>::row_type r = validLevel
        ? data_[level]->face_to_cell_[cpgrid::EntityRep<1>(face, true)]
        : current_data_->back()->face_to_cell_[cpgrid::EntityRep<1>(face, true)];

    // Deal with special case of single element in row.
    if (r.size() == 1) {
        const bool use_first = (r[0].orientation() == (local_index == 0));
        return use_first ? r[0].index() : -1;
    }

    // In the parallel case we store non-existent cells for faces along the front region.
    // These are marked with EntityRep<0>::InvalidIndex; replace them
    // with -1 for this interface's sake.
    auto getLocalIndexOrMinusOne = [](const cpgrid::EntityRep<0>& cell) {
        return cell.isValid() ? cell.index() : -1;
    };

    // We will use a helper array to store the ordered face->cell relation requested in
    // this interface, i.e. the face (normal) is oriented from nb[0] to nb[1].
    std::array<int, 2> nbs{ getLocalIndexOrMinusOne(r[0]), getLocalIndexOrMinusOne(r[1]) };

    // The entities in r are in arbitrary order, and we must use the orientations to
    // get the correct order. If the first entity has orientation true, we are fine,
    // otherwise we must flip the order.
    //
    // Note that this requires that also the "other process" neighbours have proper
    // orientations!
    if (!r[0].orientation()) {
        std::swap(nbs[0], nbs[1]);
    }

    return nbs[local_index];
}

int CpGrid::numCellFaces() const
{
    return current_data_->back()->cell_to_face_.dataSize();
}

int CpGrid::numFaceVertices(int face) const
{
    return current_data_->back()->face_to_point_[face].size();
}

int CpGrid::faceVertex(int face, int local_index) const
{
    return current_data_->back()->face_to_point_[face][local_index];
}

double CpGrid::cellCenterDepth(int cell_index) const
{
    // Here cell center depth is computed as a raw average of cell corner depths.
    // This generally gives slightly different results than using the cell centroid.
    double zz = 0.0;
    const int nv = current_data_->back()->cell_to_point_[cell_index].size();
    const int nd = 3;
    for (int i=0; i<nv; ++i) {
        zz += vertexPosition(current_data_->back()->cell_to_point_[cell_index][i])[nd-1];
    }
    return zz/nv;
}

const Dune::FieldVector<double,3> CpGrid::faceCenterEcl(int cell_index, int face, const Dune::cpgrid::Intersection& intersection) const
{
    // This method is an alternative to the method faceCentroid(...).
    // The face center is computed as a raw average of cell corners.
    // For faulted cells this gives different results then average of face nodes
    // that seems to agree more with eclipse.
    // This assumes the cell nodes are ordered
    // 6---7
    // | T |
    // 4---5
    //   2---3
    //   | B |
    //   0---1

    // this follows the DUNE reference cube
    static const int faceVxMap[ 6 ][ 4 ] = { {0, 2, 4, 6}, // face 0 - I_FACE false
                                             {1, 3, 5, 7}, // face 1 - I_FACE true
                                             {0, 1, 4, 5}, // face 2 - J_FACE false
                                             {2, 3, 6, 7}, // face 3 - J_FACE true
                                             {0, 1, 2, 3}, // face 4 - K_FACE false
                                             {4, 5, 6, 7}  // face 5 - K_FACE true
    };

    assert (current_data_->back()->cell_to_point_[cell_index].size() == 8);
    if (lgr_backend_ == Opm::Refinement::Backend::Conforming && maxLevel() > 0) {
        return faceCenterEclConforming_(cell_index, face, intersection);
    }
    Dune::FieldVector<double,3> center(0.0);

    bool isCoarseCellInside = (intersection.inside().level() == 0);
    bool isCoarseCellOutside = false;
    if (intersection.neighbor()){
        isCoarseCellOutside = (intersection.outside().level() == 0);
    }
    bool twoCoarseNeighboringCells = isCoarseCellInside && isCoarseCellOutside;
    bool isOnGridBoundary_coarseNeighboringCell = intersection.boundary() && isCoarseCellInside && (!intersection.neighbor());

    // For CpGrid with LGRs, a refined face with a coarse neighboring cell and a refined neighboring cell
    // (that is when the face belongs to the boundary of an LGR and is located in the interior of the grid),
    // unfortunately leads us to a different order of the faces, in cell_to_face_, depending on if the
    // neighboring cell, here with cell_index index, is the coarse one or the refined one. Preceisely,
    // cell_to_face_[cell_index - coarse neighboring cell] = { left, right, front, back, bottom, top} = {0,1,2,3,4,5} with
    // the notation above, and
    // cell_to_face_[cell_index - refined neighboring cell] = {bottom, front, left, right, back, top} = {2,3,1,4,0,5} with
    // the notation used in faceVxMap.

    for( int i=0; i<4; ++i ) {
        if ((maxLevel() == 0) || twoCoarseNeighboringCells || isOnGridBoundary_coarseNeighboringCell) {
            center += vertexPosition(current_data_->back()->cell_to_point_[cell_index][ faceVxMap[ face ][ i ] ]);
        }
        else { //  (refined) intersection with one coarse neighboring cell and one refined neighboring cell
            center += vertexPosition(current_data_->back()->face_to_point_[intersection.id()][i]);
        }
    }

    for (int i=0; i<3; ++i) {
        center[i] /= 4;
    }
    return center;

}

const Dune::FieldVector<double,3> CpGrid::faceAreaNormalEcl(int face) const
{
    return faceAreaNormalEcl(face, static_cast<int>(current_data_->size()) - 1);
}

const Dune::FieldVector<double,3> CpGrid::faceAreaNormalEcl(int face, int level) const
{
    assert((level >= 0) && (level < static_cast<int>(current_data_->size())));
    const auto& data = *(*current_data_)[level];
    assert((face >= 0) && (face < static_cast<int>(data.face_to_point_.size())));
    const auto pos = [&data](int vertex) -> const Dune::FieldVector<double,3>& {
        return data.geomVector<3>()[cpgrid::EntityRep<3>(vertex, true)].center();
    };
    // same implementation as ResInsight
    const int nd = Dune::FieldVector<double,3>::dimension;
    const int nv = data.face_to_point_[face].size();
    switch (nv)
    {
    case 0:
    case 1:
    case 2:
        {
            return Dune::FieldVector<double,3>(0.0);
        }
        break;
    case 3:
        {
            Dune::FieldVector<double,3> a = pos(data.face_to_point_[face][0])
                - pos(data.face_to_point_[face][2]);
            Dune::FieldVector<double,3> b = pos(data.face_to_point_[face][1])
                - pos(data.face_to_point_[face][2]);
            Dune::FieldVector<double,3> areaNormal = cross(a,b);
            for (int i=0; i<nd; ++i) {
                areaNormal[i] /= 2;
            }
            return areaNormal;
        }
        break;
    case 4:
        {
            Dune::FieldVector<double,3> a = pos(data.face_to_point_[face][0])
                - pos(data.face_to_point_[face][2]);
            Dune::FieldVector<double,3> b = pos(data.face_to_point_[face][1])
                - pos(data.face_to_point_[face][3]);
            Dune::FieldVector<double,3> areaNormal = cross(a,b);
            areaNormal *= 0.5;
            return areaNormal;
        }
        break;
    default:
        {
            int h = (nv - 1)/2;
            int k = (nv % 2) ? 0 : nv - 1;

            Dune::FieldVector<double,3> areaNormal(0.0);
            // First quads
            for (int i = 1; i < h; ++i)
            {
                Dune::FieldVector<double,3> a = pos(data.face_to_point_[face][2*i])
                    - pos(data.face_to_point_[face][0]);
                Dune::FieldVector<double,3> b = pos(data.face_to_point_[face][2*i+1])
                    - pos(data.face_to_point_[face][2*i-1]);
                areaNormal += cross(a,b);
            }

            // Last triangle or quad
            Dune::FieldVector<double,3> a = pos(data.face_to_point_[face][2*h])
                - pos(data.face_to_point_[face][0]);
            Dune::FieldVector<double,3> b = pos(data.face_to_point_[face][k])
                - pos(data.face_to_point_[face][2*h-1]);
            areaNormal += cross(a,b);

            areaNormal *= 0.5;

            return areaNormal;
        }

    }
}

const Dune::FieldVector<double,3>& CpGrid::vertexPosition(int vertex) const
{
    return current_data_->back()->geomVector<3>()[cpgrid::EntityRep<3>(vertex, true)].center();
}

double CpGrid::faceArea(int face) const
{
    return current_data_->back()->geomVector<1>()[cpgrid::EntityRep<1>(face, true)].volume();
}

const Dune::FieldVector<double,3>& CpGrid::faceCentroid(int face) const
{
    return current_data_->back()->geomVector<1>()[cpgrid::EntityRep<1>(face, true)].center();
}

const Dune::FieldVector<double,3>& CpGrid::faceNormal(int face) const
{
    return current_data_->back()->face_normals_.get(face);
}

double CpGrid::cellVolume(int cell) const
{
    return current_data_->back()->geomVector<0>()[cpgrid::EntityRep<0>(cell, true)].volume();
}

const Dune::FieldVector<double,3>& CpGrid::cellCentroid(int cell) const
{
    return current_data_->back()->geomVector<0>()[cpgrid::EntityRep<0>(cell, true)].center();
}

CpGrid::CentroidIterator<0> CpGrid::beginCellCentroids() const
{
    return CentroidIterator<0>(current_data_->back()->geomVector<0>().begin());
}

CpGrid::CentroidIterator<1> CpGrid::beginFaceCentroids() const
{
    return CentroidIterator<1>(current_data_->back()->geomVector<1>().begin());
}

const std::vector<int>& CpGrid::sortedNumAquiferCells() const{
    return current_data_->back()->sortedNumAquiferCells();
}

int CpGrid::boundaryId(int face) const
{
    // Note that this relies on the following implementation detail:
    // The grid is always construct such that the faces where
    // orientation() returns true are oriented along the positive IJK
    // direction. Oriented means that the first cell attached to face
    // has the lower index.
    int ret = 0;
    cpgrid::EntityRep<1> f(face, true);
    if (current_data_->back()->face_to_cell_[f].size() == 1) {
        if (current_data_->back()->uniqueBoundaryIds()) {
            // Use the unique boundary ids.
            ret = current_data_->back()->unique_boundary_ids_[f];
        } else {
            // Use the face tag based ids, i.e. 1-6 for i-, i+, j-, j+, k-, k+.
            const bool normal_is_in =
                !(current_data_->back()->face_to_cell_[f][0].orientation());
            enum face_tag tag = current_data_->back()->face_tag_[f];
            switch (tag) {
            case I_FACE:
                //                   LEFT : RIGHT
                ret = normal_is_in ? 1    : 2; // min(I) : max(I)
                break;
            case J_FACE:
                //                   BACK : FRONT
                ret = normal_is_in ? 3    : 4; // min(J) : max(J)
                break;
            case K_FACE:
                // Note: TOP at min(K) as 'z' measures *depth*.
                //                   TOP  : BOTTOM
                ret = normal_is_in ? 5    : 6; // min(K) : max(K)
                break;
            case NNC_FACE:
                // This should not be possible, as NNC "faces" always
                // have two cell neighbours.
                OPM_THROW(std::logic_error, "NNC face at boundary. This should never happen!");
            }
        }
    }
    return ret;
}

const CpGrid::InterfaceMap& CpGrid::cellScatterGatherInterface() const
{
    return *cell_scatter_gather_interfaces_;
}

const CpGrid::InterfaceMap& CpGrid::pointScatterGatherInterface() const
{
    return *point_scatter_gather_interfaces_;
}

void CpGrid::switchToGlobalView()
{
    current_data_ = &data_;
}

void CpGrid::switchToDistributedView()
{
    if (distributed_data_.empty()) {
        OPM_THROW(std::logic_error, "No distributed view available in grid");
    } else {
        current_data_ = &distributed_data_;
    }
}

#if HAVE_MPI

const cpgrid::CpGridDataTraits::CommunicationType& CpGrid::cellCommunication() const
{
    return current_data_->back()->cellCommunication();
}

cpgrid::CpGridDataTraits::ParallelIndexSet& CpGrid::getCellIndexSet()
{
    return current_data_->back()->cellIndexSet();
}

cpgrid::CpGridDataTraits::RemoteIndices& CpGrid::getCellRemoteIndices()
{
    return current_data_->back()->cellRemoteIndices();
}

const cpgrid::CpGridDataTraits::ParallelIndexSet& CpGrid::getCellIndexSet() const
{
    return current_data_->back()->cellIndexSet();
}

const cpgrid::CpGridDataTraits::RemoteIndices& CpGrid::getCellRemoteIndices() const
{
    return current_data_->back()->cellRemoteIndices();
}

#endif

std::vector<std::size_t>
CpGrid::processEclipseFormat(const Opm::EclipseGrid* ecl_grid,
                             Opm::EclipseState* ecl_state,
                             const bool periodic_extension,
                             const bool turn_normals,
                             const bool clip_z,
                             const bool pinchActive,
                             const bool edge_conformal)
{
    if (ecl_grid != nullptr && ecl_grid->dualPorosity()) {
        return this->processDualContinuumEclipseFormat(*ecl_grid, ecl_state, periodic_extension,
                                                       turn_normals, clip_z, pinchActive,
                                                       edge_conformal);
    }

    auto removed_cells = current_data_->back()->
        processEclipseFormat(ecl_grid, ecl_state,
                             periodic_extension,
                             turn_normals,
                             clip_z,
                             pinchActive,
                             edge_conformal);

    current_data_->back()->ccobj_.broadcast(current_data_->back()->logical_cartesian_size_.data(),
                                            current_data_->back()->logical_cartesian_size_.size(),
                                            0);
    return removed_cells;
}

std::vector<std::size_t>
CpGrid::processDualContinuumEclipseFormat(const Opm::EclipseGrid& ecl_grid,
                                          Opm::EclipseState* ecl_state,
                                          const bool periodic_extension,
                                          const bool turn_normals,
                                          const bool clip_z,
                                          const bool pinchActive,
                                          const bool edge_conformal)
{
    const auto continuumSize = ecl_grid.getCartesianSize() / 2;

    // The corner-point preprocessing cannot make a grid without cells.
    const auto hasActiveCells = [&ecl_grid, continuumSize](const std::size_t first)
    {
        for (std::size_t g = first; g < first + continuumSize; ++g) {
            if (ecl_grid.cellActive(g)) {
                return true;
            }
        }
        return false;
    };
    if (!hasActiveCells(0) || !hasActiveCells(continuumSize)) {
        OPM_THROW(std::invalid_argument,
                  "A dual-continuum grid needs active cells in both the matrix and the fracture continuum.");
    }

    // One grid per continuum, each with its own cells, activity, MINPV and
    // geometry, joined into one grid with the matrix cells first. The copies of
    // the input and the grid of the fracture continuum are released as soon as
    // they are used, before the NNC faces are added.
    //
    // The pore volumes, permeabilities and multipliers used for the MINPV and PINCH
    // processing are those of the dual-continuum grid, indexed by its cells.
    cpgrid::ContinuumSlice matrixSlice;
    matrixSlice.broadcastSuccess = false;
    cpgrid::ContinuumSlice fractureSlice;
    fractureSlice.cartOffset = static_cast<int>(continuumSize);

    auto removed_cells = [&]
    {
        const auto matrix = ecl_grid.continuumGrid(Opm::EclipseGrid::Continuum::Matrix);
        return current_data_->back()->
            processEclipseFormat(&matrix, ecl_state, periodic_extension, turn_normals,
                                 clip_z, pinchActive, edge_conformal, &matrixSlice);
    }();

    {
        CpGrid fracture_grid(current_data_->back()->ccobj_);
        {
            const auto fracture = ecl_grid.continuumGrid(Opm::EclipseGrid::Continuum::Fracture);
            const auto removed_fracture_cells = fracture_grid.current_data_->back()->
                processEclipseFormat(&fracture, ecl_state, periodic_extension, turn_normals,
                                     clip_z, pinchActive, edge_conformal, &fractureSlice);

            // The cells removed from the fracture continuum are numbered in that continuum.
            for (const auto cell : removed_fracture_cells) {
                removed_cells.push_back(cell + continuumSize);
            }
        }

        this->appendDisjointGrid(fracture_grid);
    }

    current_data_->back()->continuum_cart_offset_ = continuumSize;

    if (ecl_state != nullptr) {
        // The pinch-out connections of both continua, and the pruning of the global
        // property arrays that follows, are what processing a single grid does.
        std::vector<Opm::NNCdata> pinchNnc = std::move(matrixSlice.pinchNnc);
        pinchNnc.insert(pinchNnc.end(), fractureSlice.pinchNnc.begin(), fractureSlice.pinchNnc.end());
        ecl_state->setPinchNNC(std::move(pinchNnc));
        ecl_state->prune_global_for_schedule_run();

        // The explicit NNCs include the coupling between each matrix cell and its
        // fracture twin. They connect cells of different continua, so they could not be
        // added when processing the continua one by one.
        std::vector<std::pair<int, int>> nnc;
        for (const auto& entry : ecl_state->getInputNNC().input()) {
            nnc.emplace_back(entry.cell1, entry.cell2);
        }
        this->addNnc(nnc);
    }

    current_data_->back()->ccobj_.broadcast(current_data_->back()->logical_cartesian_size_.data(),
                                            current_data_->back()->logical_cartesian_size_.size(),
                                            0);
    return removed_cells;
}

std::vector<std::size_t>
CpGrid::processEclipseFormat(const Opm::EclipseGrid* ecl_grid_ptr,
                             Opm::EclipseState* ecl_state,
                             const bool periodic_extension,
                             const bool turn_normals,
                             const bool clip_z,
                             const bool edge_conformal)
{
    return processEclipseFormat(ecl_grid_ptr, ecl_state, periodic_extension, turn_normals, clip_z,
                                !ecl_grid_ptr || ecl_grid_ptr->isPinchActive(), edge_conformal);
}

void CpGrid::processEclipseFormat(const grdecl& input_data,
                                  const bool remove_ij_boundary,
                                  const bool turn_normals,
                                  const bool edge_conformal)
{
    using NNCMap = std::set<std::pair<int, int>>;
    using NNCMaps = std::array<NNCMap, 2>;
    NNCMaps nnc;
    current_data_->back()->processEclipseFormat(input_data,
                                                nullptr,
                                                nnc,
                                                remove_ij_boundary,
                                                turn_normals,
                                                /* pinchActive = */ false,
                                                /* tolerance_unique_ponts = */ 0.0,
                                                edge_conformal);

    current_data_->back()->ccobj_.broadcast(current_data_->back()->logical_cartesian_size_.data(),
                                            current_data_->back()->logical_cartesian_size_.size(),
                                            0);
}

void CpGrid::appendDisjointGrid(const CpGrid& other)
{
    const auto isPlainSerialGrid = [](const CpGrid& grid)
    {
        return grid.data_.size() == 1 && grid.distributed_data_.empty()
            && grid.current_data_ == &grid.data_;
    };

    if (!isPlainSerialGrid(*this) || !isPlainSerialGrid(other)) {
        OPM_THROW(std::invalid_argument,
                  "appendDisjointGrid() requires two grids that are not distributed, without local grid refinement. "
                  "Extend the grid before load balancing.");
    }

    data_.back()->appendDisjointGrid(*other.data_.back());
}

void CpGrid::addNnc(const std::vector<std::pair<int, int>>& nnc)
{
    if (data_.size() != 1 || !distributed_data_.empty() || current_data_ != &data_) {
        OPM_THROW(std::invalid_argument,
                  "addNnc() requires a grid that is not distributed, without local grid refinement.");
    }

    auto& data = *data_.back();
    const auto& global_cell = data.global_cell_;
    // The compressed index of the cell with the Cartesian index, or -1 if that cell
    // is inactive.
    const auto compressedIndex = [&global_cell](const int cartesian)
    {
        // The active cells are numbered in increasing Cartesian order.
        const auto it = std::lower_bound(global_cell.begin(), global_cell.end(), cartesian);
        return (it != global_cell.end() && *it == cartesian)
            ? static_cast<int>(it - global_cell.begin())
            : -1;
    };
    const auto areNeighbours = [this](const int c1, const int c2)
    {
        const auto faces = this->cellFaceRow(c1);
        for (int i = 0; i < faces.size(); ++i) {
            const int face = faces[i].index();
            const int other = (this->faceCell(face, 0) == c1) ? this->faceCell(face, 1)
                                                               : this->faceCell(face, 0);
            if (other == c2) {
                return true;
            }
        }
        return false;
    };

    // An index outside the logical grid is an error, not an inactive cell.
    const auto& logicalSize = data.logical_cartesian_size_;
    const auto numCartesian = static_cast<std::int64_t>(logicalSize[0]) * logicalSize[1] * logicalSize[2];
    const auto checkCartesian = [numCartesian](const int cartesian)
    {
        if (cartesian < 0 || cartesian >= numCartesian) {
            OPM_THROW(std::invalid_argument,
                      "addNnc(): the Cartesian index " + std::to_string(cartesian)
                      + " is outside the logical grid of " + std::to_string(numCartesian) + " cells.");
        }
    };

    // Ordered and unique by Cartesian index, as the NNC faces made when processing the
    // Eclipse format.
    std::set<std::pair<int, int>> unique_nnc;
    for (const auto& [a, b] : nnc) {
        checkCartesian(a);
        checkCartesian(b);
        unique_nnc.emplace(std::min(a, b), std::max(a, b));
    }

    std::vector<std::pair<int, int>> cell_pairs;
    std::size_t numInactive = 0;
    for (const auto& [a, b] : unique_nnc) {
        const int c1 = compressedIndex(a);
        const int c2 = compressedIndex(b);
        if (c1 < 0 || c2 < 0) {
            // A connection to an inactive cell has nothing to connect, as for a grid
            // processed in one piece.
            ++numInactive;
            continue;
        }
        if (areNeighbours(c1, c2)) {
            continue;
        }
        cell_pairs.emplace_back(c1, c2);
    }

    if (numInactive > 0) {
        Opm::OpmLog::info(std::to_string(numInactive) +
                          " non-neighbouring connection(s) skipped because a cell is inactive.");
    }

    data.appendNncFaces(cell_pairs);
}

void CpGrid::processEclipseFormat(const grdecl& input_data,
                                  const std::set<std::pair<int,int>>& pinchNnc,
                                  const bool pinchActive,
                                  const bool edge_conformal)
{
    std::array<std::set<std::pair<int, int>>, 2> nnc;
    nnc[0] = pinchNnc;
    current_data_->back()->processEclipseFormat(input_data,
                                                nullptr,
                                                nnc,
                                                /* remove_ij_boundary = */ false,
                                                /* turn_normals = */ false,
                                                pinchActive,
                                                /* tolerance_unique_ponts = */ 0.0,
                                                edge_conformal);

    current_data_->back()->ccobj_.broadcast(current_data_->back()->logical_cartesian_size_.data(),
                                            current_data_->back()->logical_cartesian_size_.size(),
                                            0);
}

template<int dim>
cpgrid::Entity<dim> createEntity(const CpGrid& grid,int index,bool orientation)
{
    return cpgrid::Entity<dim>(*(grid.current_data_->back()), index, orientation);
}
template cpgrid::Entity<0> createEntity(const CpGrid&, int, bool);
template cpgrid::Entity<3> createEntity(const CpGrid&, int, bool);
template cpgrid::Entity<1> createEntity(const CpGrid&, int, bool); // needed in distribution_test.cpp

std::array<double,3> CpGrid::getEclCentroid(const int& elemIdx) const
{
    return this-> current_data_->back() -> computeEclCentroid(elemIdx);
}

std::array<double,3> CpGrid::getEclCentroid(const cpgrid::Entity<0>& elem) const
{
    return this-> getEclCentroid(elem.index());
}

} // namespace Dune
