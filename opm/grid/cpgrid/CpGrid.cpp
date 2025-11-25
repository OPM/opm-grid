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
#include "mpi.h"
#include <opm/grid/utility/platform_dependent/reenable_warnings.h>
#endif

#if HAVE_ECL_INPUT
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#endif

#include "../CpGrid.hpp"
#include "LgrHelpers.hpp"
#include "ParentToChildrenCellGlobalIdHandle.hpp"
#include "NestedRefinementUtilities.hpp"
#include <opm/grid/common/MetisPartition.hpp>
#include <opm/grid/common/ZoltanPartition.hpp>
#include <opm/grid/GraphOfGridWrappers.hpp>
//#include <opm/grid/common/ZoltanGraphFunctions.hpp>
#include <opm/grid/common/GridPartitioning.hpp>
//#include <opm/grid/common/WellConnections.hpp>
#include <opm/grid/common/CommunicationUtils.hpp>

//#include <fstream>
//#include <iostream>
#include <algorithm>
#include <iomanip>
#include <numeric>
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
    : current_view_data_(),
      distributed_data_(),
      cell_scatter_gather_interfaces_(new InterfaceMap, FreeInterfaces{}),
      point_scatter_gather_interfaces_(new InterfaceMap, FreeInterfaces{}),
      global_id_set_ptr_()
{
    data_.push_back(std::make_shared<cpgrid::CpGridData>(data_));
    current_view_data_ = data_[0].get();
    current_data_ = &data_;
    global_id_set_ptr_ = std::make_shared<cpgrid::GlobalIdSet>(*current_view_data_);
}

CpGrid::CpGrid(MPIHelper::MPICommunicator comm)
    : current_view_data_(),
      distributed_data_(),
      cell_scatter_gather_interfaces_(new InterfaceMap, FreeInterfaces{}),
      point_scatter_gather_interfaces_(new InterfaceMap, FreeInterfaces{}),
      global_id_set_ptr_()
{
    data_.push_back(std::make_shared<cpgrid::CpGridData>(comm, data_));
    current_view_data_ = data_[0].get();
    current_data_ = &data_;
    global_id_set_ptr_ = std::make_shared<cpgrid::GlobalIdSet>(*current_view_data_);
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
                    int level)
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

    if ((maxLevel()>0) && (partitionMethod!= Dune::PartitionMethod::zoltanGoG)) {
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
                if (std::any_of(existingParts.begin(), existingParts.end(),
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
        comm().barrier();

        // first create the overlap
        auto noImportedOwner = addOverlapLayer(*this,
                                               computedCellPart,
                                               exportList,
                                               importList,
                                               cc,
                                               addCornerCells,
                                               transmissibilities,
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
                " Try decreasing the imbalance tolerance with the argument \n"
                " --imbalance-tolerance. The current value is "
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
        (*global_id_set_ptr_).insertIdSet(*distributed_data_[0]);
        distributed_data_[0]-> index_set_.reset(new cpgrid::IndexSet(distributed_data_[0]->cell_to_face_.size(),
                                                                     distributed_data_[0]-> geomVector<3>().size()));



        current_view_data_ = distributed_data_[0].get();
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
        std::fill(offset, offset + num_per_layer, zlow);
        offset += num_per_layer;
        double zhigh = (k+1+shift[2])*cellsize[2];
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

    // Note: This is a Cartesian, matching grid which is edge-conforming
    // regardless of the edge_conformal flag.
    current_view_data_->processEclipseFormat(g,
#if HAVE_ECL_INPUT
                                             /* ecl_state = */ nullptr,
#endif
                                             nnc,
                                             /* remove_ij_boundary = */ false,
                                             /* turn_normals = */ false,
                                             /* pinchActive = */ false,
                                             /* tolerance_unique_ponts = */ 0.0,
                                             /* edge_conformal = */ false);

    // global grid only on rank 0
    current_view_data_->ccobj_.broadcast(current_view_data_->logical_cartesian_size_.data(),
                                         current_view_data_->logical_cartesian_size_.size(),
                                         0);
}

const std::array<int, 3>& CpGrid::logicalCartesianSize() const
{
    // Temporary. For a grid with LGRs, we set the logical cartesian size of the LeafGridView as the one for level 0.
    //            Goal: CartesianIndexMapper well-defined for CpGrid LeafView with LGRs.
    return current_view_data_ -> logical_cartesian_size_;
}

const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& CpGrid::currentData() const
{
    return *current_data_;
}

std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& CpGrid::currentData()
{
    return *current_data_;
}

const std::vector<int>& CpGrid::globalCell() const
{
    // Temporary. For a grid with LGRs, we set the globalCell() of the as the one for level 0.
    //            Goal: CartesianIndexMapper well-defined for CpGrid LeafView with LGRs.
    return currentData().back() -> global_cell_;
}

void CpGrid::computeGlobalCellLgr(const int& level, const std::array<int,3>& startIJK, std::vector<int>& global_cell_lgr)
{
    assert(level);
    for (const auto& element : elements(levelGridView(level))) {
        // Element belogns to an LGR, therefore has a father. Get IJK of the father in the level grid the father was born.
        // For CARFIN, parent cells belong to level 0.
        std::array<int,3> parentIJK = {0,0,0};
        currentData()[element.father().level()]->getIJK(element.father().index(), parentIJK);
        // Each parent cell has been refined in cells_per_dim[0]*cells_per_dim[1]*cells_per_dim[2] child cells.
        // element has certain 'position' inside its parent cell that can be described with 'IJK' indices, let's denote them by ijk,
        // where 0<= i < cells_per_dim[0], 0<= j < cells_per_dim[1], 0<= k < cells_per_dim[2].
        const auto& cells_per_dim = currentData()[level]->cells_per_dim_;
        //
        // Refined cell (here 'element') has "index in parent cell": k*cells_per_dim[0]*cells_per_dim[1] + j*cells_per_dim[0] + i
        // and it's stored in  cell_to_idxInParentCell_.
        auto idx_in_parent_cell =  currentData()[level]-> cell_to_idxInParentCell_[element.index()];
        // Find ijk.
        std::array<int,3> childIJK = Opm::Lgr::getIJK(idx_in_parent_cell, cells_per_dim);
        // The corresponding lgrIJK can be computed as follows:
        const std::array<int,3>& lgrIJK = { ( (parentIJK[0] - startIJK[0])*cells_per_dim[0] ) + childIJK[0],  // Shift parent index according to the startIJK of the LGR.
                                            ( (parentIJK[1] - startIJK[1])*cells_per_dim[1] ) + childIJK[1],
                                            ( (parentIJK[2] - startIJK[2])*cells_per_dim[2] ) + childIJK[2] };
        // Dimensions of the "patch of cells" formed when providing startIJK and endIJK for an LGR
        const auto& lgr_logical_cartesian_size = currentData()[level]->logical_cartesian_size_;
        global_cell_lgr[element.index()] = (lgrIJK[2]*lgr_logical_cartesian_size[0]*lgr_logical_cartesian_size[1]) + (lgrIJK[1]*lgr_logical_cartesian_size[0]) + lgrIJK[0];
    }
}

void CpGrid::computeGlobalCellLeafGridViewWithLgrs(std::vector<int>& global_cell_leaf)
{
    for (const auto& element: elements(leafGridView())) {
        // In the context of allowed nested refinement, we lookup for the oldest ancestor, belonging to level-zero-grid.
        auto ancestor = element.getOrigin();
        int origin_in_level_zero = ancestor.index();
        assert(origin_in_level_zero < currentData().front()->size(0));
        global_cell_leaf[element.index()] = currentData().front()-> global_cell_[origin_in_level_zero];
    }
}

std::vector<std::unordered_map<std::size_t, std::size_t>> CpGrid::mapLocalCartesianIndexSetsToLeafIndexSet() const
{
    std::vector<std::unordered_map<std::size_t, std::size_t>> localCartesianIdxSets_to_leafIdx(maxLevel()+1); // Plus level 0
    for (const auto& element : elements(leafGridView())) {
        const auto& global_cell_level = currentData()[element.level()]->globalCell()[element.getLevelElem().index()];
        localCartesianIdxSets_to_leafIdx[element.level()][global_cell_level] = element.index();
    }
    return localCartesianIdxSets_to_leafIdx;
}

std::vector<std::array<int,2>> CpGrid::mapLeafIndexSetToLocalCartesianIndexSets() const
{
    std::vector<std::array<int,2>> leafIdx_to_localCartesianIdxSets(currentData().back()->size(0));
    for (const auto& element : elements(leafGridView())) {
        const auto& global_cell_level = currentData()[element.level()]->globalCell()[element.getLevelElem().index()];
        leafIdx_to_localCartesianIdxSets[element.index()] = {element.level(), global_cell_level};
    }
    return leafIdx_to_localCartesianIdxSets;
}

void CpGrid::getIJK(const int c, std::array<int,3>& ijk) const
{
    current_view_data_->getIJK(c, ijk);
}

bool CpGrid::uniqueBoundaryIds() const
{
    return current_view_data_->uniqueBoundaryIds();
}

void CpGrid::setUniqueBoundaryIds(bool uids)
{
    current_view_data_->setUniqueBoundaryIds(uids);
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
    return cpgrid::Iterator<codim, All_Partition>(*current_view_data_, 0, true);
}
template typename CpGridTraits::template Codim<0>::LeafIterator CpGrid::leafbegin<0>() const;
template typename CpGridTraits::template Codim<1>::LeafIterator CpGrid::leafbegin<1>() const;
template typename CpGridTraits::template Codim<3>::LeafIterator CpGrid::leafbegin<3>() const;


template<int codim>
typename CpGridTraits::template Codim<codim>::LeafIterator CpGrid::leafend() const
{
    return cpgrid::Iterator<codim, All_Partition>(*current_view_data_, size(codim), true);
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
    return cpgrid::Iterator<codim, PiType>(*current_view_data_, 0, true);
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
    return cpgrid::Iterator<codim, PiType>(*current_view_data_, size(codim), true);
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
    return currentData().back()->size(codim);
}

int CpGrid::size (int level, GeometryType type) const
{
    if (level<0 || level>maxLevel())
        DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
    return currentData()[level]->size(type);
}

int CpGrid::size (GeometryType type) const
{
    return currentData().back()->size(type);
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
    return *current_view_data_->index_set_;
}

void CpGrid::globalRefine (int refCount)
{
    if (refCount < 0) {
        OPM_THROW(std::logic_error, "Invalid argument. Provide a nonnegative integer for global refinement.");
    }

    // Throw if strict local refinement is detected.
    if(this->maxLevel()) {
        // For strict local refinement, these sizes are identical. Global refinement
        // results in a difference in at least one direction (x, y, or z).
        bool sameNX = currentData().back()->logicalCartesianSize()[0] == currentData().front()->logicalCartesianSize()[0];
        bool sameNY = currentData().back()->logicalCartesianSize()[1] == currentData().front()->logicalCartesianSize()[1];
        bool sameNZ = currentData().back()->logicalCartesianSize()[2] == currentData().front()->logicalCartesianSize()[2];
        if (sameNX && sameNY && sameNZ) {
            OPM_THROW(std::logic_error, "Global refinement of a mixed grid with coarse and refined cells is not supported yet.");
        }
    }
    if (refCount>0) {
        for (int refinedLevel = 0; refinedLevel < refCount; ++refinedLevel) {

            std::vector<int> assignRefinedLevel(current_view_data_-> size(0));
            const auto& preAdaptMaxLevel = this ->maxLevel();
            std::vector<std::string> lgr_name_vec = { "GR" + std::to_string(preAdaptMaxLevel +1) };
            const std::array<int,3>& endIJK = currentData().back()->logicalCartesianSize();

            for(const auto& element: elements(this-> leafGridView())) {
                // Mark all the elements of the current leaf grid view for refinement
                mark(1, element);
                assignRefinedLevel[element.index()] = preAdaptMaxLevel +1;
            }

            preAdapt();
            refineAndUpdateGrid(/* cells_per_dim_vec = */ {{2,2,2}}, assignRefinedLevel, lgr_name_vec, {{0,0,0}}, {endIJK});
            postAdapt();
        }
    }
}

const std::vector< Dune :: GeometryType >& CpGrid::geomTypes( const int codim ) const
{
    return leafIndexSet().geomTypes( codim );
}

template <int codim>
cpgrid::Entity<codim> CpGrid::entity( const cpgrid::Entity< codim >& seed ) const
{
    return cpgrid::Entity<codim>( *(this->current_view_data_), seed );
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
        return current_view_data_->unique_boundary_ids_.size();
    }
    else
    {
        unsigned int numBndSegs = 0;
        const int num_faces = numFaces();
        for (int i = 0; i < num_faces; ++i) {
            cpgrid::EntityRep<1> face(i, true);
            if (current_view_data_->face_to_cell_[face].size() == 1) {
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
    return current_view_data_->ccobj_;
}

//

const std::vector<double>& CpGrid::zcornData() const {
    return current_view_data_->zcornData();
}

int CpGrid::numCells(int level) const
{
    bool validLevel = (level>-1) && (level<= maxLevel());
    return validLevel? data_[level]->cell_to_face_.size() : current_view_data_->cell_to_face_.size();
}
/// \brief Get the number of faces.
int CpGrid::numFaces(int level) const
{
    bool validLevel = (level>-1) && (level<= maxLevel());
    return validLevel? data_[level]->face_to_cell_.size() : current_view_data_->face_to_cell_.size();
}
/// \brief Get The number of vertices.
int CpGrid::numVertices() const
{
    return current_view_data_->geomVector<3>().size();
}

int CpGrid::numCellFaces(int cell, int level) const
{
    bool validLevel = (level>-1) && (level<= maxLevel());
    return validLevel? data_[level]->cell_to_face_[cpgrid::EntityRep<0>(cell, true)].size()
        : current_view_data_->cell_to_face_[cpgrid::EntityRep<0>(cell, true)].size();
}

int CpGrid::cellFace(int cell, int local_index, int level) const
{
    bool validLevel = (level>-1) && (level<= maxLevel());
    return validLevel? data_[level]-> cell_to_face_[cpgrid::EntityRep<0>(cell, true)][local_index].index()
        : current_view_data_->cell_to_face_[cpgrid::EntityRep<0>(cell, true)][local_index].index();
}

const cpgrid::OrientedEntityTable<0,1>::row_type CpGrid::cellFaceRow(int cell) const
{
    return current_view_data_->cell_to_face_[cpgrid::EntityRep<0>(cell, true)];
}

int CpGrid::faceCell(int face, int local_index, int level) const
{
    // Get the correct neigbour relation (row from face->cell table).
    const bool validLevel = (level > -1) && (level <= maxLevel());
    const cpgrid::OrientedEntityTable<1,0>::row_type r = validLevel
        ? data_[level]->face_to_cell_[cpgrid::EntityRep<1>(face, true)]
        : current_view_data_->face_to_cell_[cpgrid::EntityRep<1>(face, true)];

    // Deal with special case of single element in row.
    if (r.size() == 1) {
        const bool use_first = (r[0].orientation() == (local_index == 0));
        return use_first ? r[0].index() : -1;
    }

    // In the parallel case we store non-existent cells for faces along the front region.
    // These are marked with index std::numeric_limits<int>::max(), we will replace them
    // with -1 for this interface's sake.
    auto getIndexConvertMaxToNegOne = [](const cpgrid::EntityRep<0>& cell) {
        return cell.index() == std::numeric_limits<int>::max() ? -1 : cell.index();
    };

    // We will use a helper array to store the ordered face->cell relation requested in
    // this interface, i.e. the face (normal) is oriented from nb[0] to nb[1].
    std::array<int, 2> nbs{ getIndexConvertMaxToNegOne(r[0]), getIndexConvertMaxToNegOne(r[1]) };

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
    return current_view_data_->cell_to_face_.dataSize();
}

int CpGrid::numFaceVertices(int face) const
{
    return current_view_data_->face_to_point_[face].size();
}

int CpGrid::faceVertex(int face, int local_index) const
{
    return current_view_data_->face_to_point_[face][local_index];
}

Dune::cpgrid::Intersection CpGrid::getParentIntersectionFromLgrBoundaryFace(const Dune::cpgrid::Intersection& intersection) const
{
    if ( intersection.neighbor()) {
        if ((intersection.inside().level() != intersection.outside().level())) {
            // one coarse and one refined neighboring cell
            /** Now, it could also be two refined cells. In that case, any of them will fit to search for the parent face */
            const auto& cellIn = intersection.inside();
            const auto& cellOut = intersection.outside();

            // Identify the coarse and the refined neighboring cell
            const auto coarseCell =  (cellIn.level() == 0) ? cellIn : cellOut;
            const auto refinedCell =  (coarseCell == cellIn) ? cellOut : cellIn;
            assert(coarseCell.level() != refinedCell.level());

            // Get parent cell (on level zero) of the refined cell
            const auto& parentCell = refinedCell.father();
            assert(refinedCell.father().level() == 0);

            // Get the index inside and orientation from the leaf grid (refined) face
            const auto& intersectionIdxInInside = intersection.indexInInside();

            for(const auto& parentIntersection : intersections(this->levelGridView(0), parentCell)){
                // Get the inInsideIdx and orientation from the parent intersection
                const auto& parentIdxInInside = parentIntersection.indexInInside();
                if (parentIdxInInside == intersectionIdxInInside) {
                    return parentIntersection;
                }
            }
        }
        OPM_THROW(std::invalid_argument, "Parent intersection not found for face with index: " + std::to_string(intersection.id()) +
                  " and index in inside: " + std::to_string(intersection.indexInInside()));
    }
    OPM_THROW(std::invalid_argument, "Face is on the boundary of the grid");
}

void CpGrid::markElemAssignLevelDetectActiveLgrs(const std::vector<std::array<int,3>>& startIJK_vec,
                                                 const std::vector<std::array<int,3>>& endIJK_vec,
                                                 std::vector<int>& assignRefinedLevel,
                                                 std::vector<int>& lgr_with_at_least_one_active_cell)
{
    auto assignAndDetect = [this, &assignRefinedLevel, &lgr_with_at_least_one_active_cell](const cpgrid::Entity<0>& element, int level)
    {
        mark(1, element);
        assignRefinedLevel[element.index()] = level+1;
        // shifted since starting grid is level 0, and refined grids levels are >= 1.
        lgr_with_at_least_one_active_cell[level] = 1;
    };
    Opm::Lgr::computeOnLgrParents(*this, startIJK_vec, endIJK_vec, assignAndDetect);
}

void CpGrid::populateCellIndexSetRefinedGrid([[maybe_unused]] int level)
{
#if HAVE_MPI
    const auto& level_data = currentData()[level];
    const auto& level_global_id_set =  level_data->global_id_set_;
    auto& level_index_set =  level_data->cellIndexSet();
    // Clean up level cell index set - needed e.g. for synchronization of cell ids.
    level_index_set = ParallelIndexSet();

    level_index_set.beginResize();
    // ParallelIndexSet::LocalIndex( elemIdx, attribute /* owner or copy */, true/false)
    // The bool indicates whether the global index might exist on other processes with other attributes.
    // RemoteIndices::rebuild has a Boolean as a template parameter telling the method whether to ignore this
    // Boolean on the indices or not when building.
    //
    // For refined level grids, we check if the parent cell is fully interior. Then, its children won't be seen
    // by any other process. Therefore, the boolean is set to false.
    for(const auto& element : elements(levelGridView(level))) {
        if ( element.partitionType() == InteriorEntity) {

            // Check if it has an overlap neighbor
            bool parentFullyInterior = true;
            const auto& parent_cell = element.father();

            const auto& intersections = Dune::intersections(levelGridView(parent_cell.level()), parent_cell);
            for (const auto& intersection : intersections) {
                if ( intersection.neighbor() ) {
                    const auto& neighborPartitionType = intersection.outside().partitionType();
                    // To help detection of fully interior cells, i.e., without overlap neighbors
                    if (neighborPartitionType == OverlapEntity )  {
                        parentFullyInterior = false;
                        // Interior cells with overlap neighbor may appear in other processess -> true
                        level_index_set.add( level_global_id_set->id(element),
                                             ParallelIndexSet::LocalIndex(element.index(), AttributeSet(AttributeSet::owner), true));
                        // Store it only once
                        break;
                    }
                }
            }
            if(parentFullyInterior) { // Fully interior cells do not appear in any other process -> false
                level_index_set.add( level_global_id_set->id(element),
                                     ParallelIndexSet::LocalIndex(element.index(), AttributeSet(AttributeSet::owner), false));
            }
        }
        else { // overlap cell
            assert(element.partitionType() == OverlapEntity);
            level_index_set.add( level_global_id_set->id(element),
                                 ParallelIndexSet::LocalIndex(element.index(), AttributeSet(AttributeSet::copy), true));
        }
    }
    level_index_set.endResize();

    currentData()[level]->cellRemoteIndices().template rebuild<true>(); // true->ignore bool in ParallelIndexSet!
    // Check the size
    assert(static_cast<std::size_t>(level_data->cellIndexSet().size()) == static_cast<std::size_t>(level_data->size(0)) );
#endif
}

void CpGrid::populateCellIndexSetLeafGridView()
{
#if HAVE_MPI
    auto& leaf_index_set =  (*current_data_).back()->cellIndexSet();
    // Clean up leaf cell index set - needed e.g. for synchronization of cell ids.
    leaf_index_set = ParallelIndexSet();

    leaf_index_set.beginResize();

    // ParallelIndexSet::LocalIndex( elemIdx, attribute /* owner or copy */, true/false)
    // The bool indicates whether the global index might exist on other processes with other attributes.
    // RemoteIndices::rebuild has a Boolean as a template parameter telling the method whether to ignore this
    // Boolean on the indices or not when building.
    for(const auto& element : elements(leafGridView())) {
        const auto& elemPartitionType = element.getLevelElem().partitionType();
        if ( elemPartitionType == InteriorEntity) {
            // Check if it has an overlap neighbor
            bool isFullyInterior = true;
            for (const auto& intersection : intersections(leafGridView(), element)) {
                if ( intersection.neighbor() ) {
                    const auto& neighborPartitionType = intersection.outside().getLevelElem().partitionType();
                    // To help detection of fully interior cells, i.e., without overlap neighbors
                    if (neighborPartitionType == OverlapEntity )  {
                        isFullyInterior = false;
                        // Interior cells with overlap neighbor may appear in other processess -> false
                        leaf_index_set.add(globalIdSet().id(element),
                                           ParallelIndexSet::LocalIndex(element.index(), AttributeSet(AttributeSet::owner), true));
                        // Store it only once
                        break;
                    }
                }
            }
            if(isFullyInterior) { // Fully interior cells do not appear in any other process -> false
                leaf_index_set.add(globalIdSet().id(element),
                                   ParallelIndexSet::LocalIndex(element.index(), AttributeSet(AttributeSet::owner), false));
            }
        }
        else { // overlap cell
            assert(elemPartitionType == OverlapEntity);
            leaf_index_set.add(globalIdSet().id(element),
                               ParallelIndexSet::LocalIndex(element.index(), AttributeSet(AttributeSet::copy), true));
        }
    }
    leaf_index_set.endResize();

    (*current_data_).back()->cellRemoteIndices().template rebuild<true>(); // true->ignore bool in ParallelIndex!

#endif
}

void CpGrid::populateLeafGlobalIdSet()
{
    // Global id for the cells in leaf grid view
    std::vector<int> leafCellIds(current_data_->back()->size(0));
    for(const auto& element: elements(leafGridView())){
        // Notice that for level zero cells the global_id_set_ is given, for refined level grids was defined
        // under the assumption of each lgr being fully contained in the interior of a process.
        // Therefore, it is not needed here to distingish between owned and overlap cells.
        auto equivElem = element.getLevelElem();
        leafCellIds[element.index()] = (*current_data_)[element.level()]->global_id_set_->id(equivElem);
    }

    // Global id for the faces in leaf grid view. Empty vector (Entity<1> not supported for CpGrid).
    std::vector<int> leafFaceIds{};

    // Global id for the points in leaf grid view
    std::vector<int> leafPointIds(current_data_->back()->size(3));
    for(const auto& point : vertices(leafGridView())){
        const auto& level_pointLevelIdx = current_data_->back()->corner_history_[point.index()];
        assert(level_pointLevelIdx[0] != -1);
        assert(level_pointLevelIdx[1] != -1);
        const auto& pointLevelEntity =  cpgrid::Entity<3>(*( (*current_data_)[level_pointLevelIdx[0]]), level_pointLevelIdx[1], true);
        leafPointIds[point.index()] = (*current_data_)[level_pointLevelIdx[0]]->global_id_set_->id(pointLevelEntity);
    }

    current_data_->back()->global_id_set_->swap(leafCellIds, leafFaceIds, leafPointIds);
}

double CpGrid::cellCenterDepth(int cell_index) const
{
    // Here cell center depth is computed as a raw average of cell corner depths.
    // This generally gives slightly different results than using the cell centroid.
    double zz = 0.0;
    const int nv = current_view_data_->cell_to_point_[cell_index].size();
    const int nd = 3;
    for (int i=0; i<nv; ++i) {
        zz += vertexPosition(current_view_data_->cell_to_point_[cell_index][i])[nd-1];
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


    assert (current_view_data_->cell_to_point_[cell_index].size() == 8);
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
            center += vertexPosition(current_view_data_->cell_to_point_[cell_index][ faceVxMap[ face ][ i ] ]);
        }
        else { //  (refined) intersection with one coarse neighboring cell and one refined neighboring cell
            center += vertexPosition(current_view_data_->face_to_point_[intersection.id()][i]);
        }
    }

    for (int i=0; i<3; ++i) {
        center[i] /= 4;
    }
    return center;

}

const Dune::FieldVector<double,3> CpGrid::faceAreaNormalEcl(int face) const
{
    // same implementation as ResInsight
    const int nd = Dune::FieldVector<double,3>::dimension;
    const int nv =  numFaceVertices(face);
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
            Dune::FieldVector<double,3> a = vertexPosition(current_view_data_->face_to_point_[face][0])
                - vertexPosition(current_view_data_->face_to_point_[face][2]);
            Dune::FieldVector<double,3> b = vertexPosition(current_view_data_->face_to_point_[face][1])
                - vertexPosition(current_view_data_->face_to_point_[face][2]);
            Dune::FieldVector<double,3> areaNormal = cross(a,b);
            for (int i=0; i<nd; ++i) {
                areaNormal[i] /= 2;
            }
            return areaNormal;
        }
        break;
    case 4:
        {
            Dune::FieldVector<double,3> a = vertexPosition(current_view_data_->face_to_point_[face][0])
                - vertexPosition(current_view_data_->face_to_point_[face][2]);
            Dune::FieldVector<double,3> b = vertexPosition(current_view_data_->face_to_point_[face][1])
                - vertexPosition(current_view_data_->face_to_point_[face][3]);
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
                Dune::FieldVector<double,3> a = vertexPosition(current_view_data_->face_to_point_[face][2*i])
                    - vertexPosition(current_view_data_->face_to_point_[face][0]);
                Dune::FieldVector<double,3> b = vertexPosition(current_view_data_->face_to_point_[face][2*i+1])
                    - vertexPosition(current_view_data_->face_to_point_[face][2*i-1]);
                areaNormal += cross(a,b);
            }

            // Last triangle or quad
            Dune::FieldVector<double,3> a = vertexPosition(current_view_data_->face_to_point_[face][2*h])
                - vertexPosition(current_view_data_->face_to_point_[face][0]);
            Dune::FieldVector<double,3> b = vertexPosition(current_view_data_->face_to_point_[face][k])
                - vertexPosition(current_view_data_->face_to_point_[face][2*h-1]);
            areaNormal += cross(a,b);

            areaNormal *= 0.5;

            return areaNormal;
        }

    }
}

const Dune::FieldVector<double,3>& CpGrid::vertexPosition(int vertex) const
{
    return current_view_data_->geomVector<3>()[cpgrid::EntityRep<3>(vertex, true)].center();
}

double CpGrid::faceArea(int face) const
{
    return current_view_data_->geomVector<1>()[cpgrid::EntityRep<1>(face, true)].volume();
}

const Dune::FieldVector<double,3>& CpGrid::faceCentroid(int face) const
{
    return current_view_data_->geomVector<1>()[cpgrid::EntityRep<1>(face, true)].center();
}

const Dune::FieldVector<double,3>& CpGrid::faceNormal(int face) const
{
    return current_view_data_->face_normals_.get(face);
}

double CpGrid::cellVolume(int cell) const
{
    return current_view_data_->geomVector<0>()[cpgrid::EntityRep<0>(cell, true)].volume();
}

const Dune::FieldVector<double,3>& CpGrid::cellCentroid(int cell) const
{
    return current_view_data_->geomVector<0>()[cpgrid::EntityRep<0>(cell, true)].center();
}

CpGrid::CentroidIterator<0> CpGrid::beginCellCentroids() const
{
    return CentroidIterator<0>(current_view_data_->geomVector<0>().begin());
}

CpGrid::CentroidIterator<1> CpGrid::beginFaceCentroids() const
{
    return CentroidIterator<1>(current_view_data_->geomVector<1>().begin());
}

const std::vector<int>& CpGrid::sortedNumAquiferCells() const{
           return current_view_data_->sortedNumAquiferCells();
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
    if (current_view_data_->face_to_cell_[f].size() == 1) {
        if (current_view_data_->uniqueBoundaryIds()) {
            // Use the unique boundary ids.
            ret = current_view_data_->unique_boundary_ids_[f];
        } else {
            // Use the face tag based ids, i.e. 1-6 for i-, i+, j-, j+, k-, k+.
            const bool normal_is_in =
                !(current_view_data_->face_to_cell_[f][0].orientation());
            enum face_tag tag = current_view_data_->face_tag_[f];
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
    current_view_data_ = data_.back().get();
    current_data_ = &data_;
}

void CpGrid::switchToDistributedView()
{
    if (distributed_data_.empty()) {
        OPM_THROW(std::logic_error, "No distributed view available in grid");
    } else {
        current_view_data_ = distributed_data_.back().get();
        current_data_ = &distributed_data_;
    }
}

#if HAVE_MPI

const cpgrid::CpGridDataTraits::CommunicationType& CpGrid::cellCommunication() const
{
    return current_view_data_->cellCommunication();
}

cpgrid::CpGridDataTraits::ParallelIndexSet& CpGrid::getCellIndexSet()
{
    return current_view_data_->cellIndexSet();
}

cpgrid::CpGridDataTraits::RemoteIndices& CpGrid::getCellRemoteIndices()
{
    return current_view_data_->cellRemoteIndices();
}

const cpgrid::CpGridDataTraits::ParallelIndexSet& CpGrid::getCellIndexSet() const
{
    return current_view_data_->cellIndexSet();
}

const cpgrid::CpGridDataTraits::RemoteIndices& CpGrid::getCellRemoteIndices() const
{
    return current_view_data_->cellRemoteIndices();
}

#endif

#if HAVE_ECL_INPUT
std::vector<std::size_t>
CpGrid::processEclipseFormat(const Opm::EclipseGrid* ecl_grid,
                             Opm::EclipseState* ecl_state,
                             const bool periodic_extension,
                             const bool turn_normals,
                             const bool clip_z,
                             const bool pinchActive,
                             const bool edge_conformal)
{
    auto removed_cells = current_view_data_->
        processEclipseFormat(ecl_grid, ecl_state,
                             periodic_extension,
                             turn_normals,
                             clip_z,
                             pinchActive,
                             edge_conformal);

    current_view_data_->ccobj_.broadcast(current_view_data_->logical_cartesian_size_.data(),
                                         current_view_data_->logical_cartesian_size_.size(),
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

#endif

void CpGrid::processEclipseFormat(const grdecl& input_data,
                                  const bool remove_ij_boundary,
                                  const bool turn_normals,
                                  const bool edge_conformal)
{
    using NNCMap = std::set<std::pair<int, int>>;
    using NNCMaps = std::array<NNCMap, 2>;
    NNCMaps nnc;
    current_view_data_->processEclipseFormat(input_data,
#if HAVE_ECL_INPUT
                                             nullptr,
#endif
                                             nnc,
                                             remove_ij_boundary,
                                             turn_normals,
                                             /* pinchActive = */ false,
                                             /* tolerance_unique_ponts = */ 0.0,
                                             edge_conformal);

    current_view_data_->ccobj_.broadcast(current_view_data_->logical_cartesian_size_.data(),
                                         current_view_data_->logical_cartesian_size_.size(),
                                         0);
}

template<int dim>
cpgrid::Entity<dim> createEntity(const CpGrid& grid,int index,bool orientation)
{
    return cpgrid::Entity<dim>(*grid.current_view_data_, index, orientation);
}
template cpgrid::Entity<0> createEntity(const CpGrid&, int, bool);
template cpgrid::Entity<3> createEntity(const CpGrid&, int, bool);
template cpgrid::Entity<1> createEntity(const CpGrid&, int, bool); // needed in distribution_test.cpp

bool CpGrid::mark(int refCount, const cpgrid::Entity<0>& element)
{
    // Throw if element has a neighboring cell from a different level.
    // E.g., a coarse cell touching the boundary of an LGR, or
    // a refined cell with a coarser/finner neighboring cell.
    for (const auto& intersection : Dune::intersections(leafGridView(), element)){
        if (intersection.neighbor() && (intersection.outside().level() != element.level()) && (element.level()==0))
            OPM_THROW(std::logic_error, "Refinement of cells at LGR boundaries is not supported, yet.");
    }
    // For serial run, mark elements also in the level they were born.
    if(currentData().size()>1) {
        // Mark element in its level
        currentData()[element.level()] -> mark(refCount, element.getLevelElem());
    }
    // Mark element (also in the serial run case) in current_view_data_. Note that if scatterGrid has been invoked, then
    // current_view_data_ == distributed_data_[0].
    return current_view_data_-> mark(refCount, element);
}

int CpGrid::getMark(const cpgrid::Entity<0>& element) const
{
    return current_view_data_->getMark(element);
}

bool CpGrid::preAdapt()
{
    // Check if elements in pre-adapt existing grids have been marked for refinment.
    // Serial run: currentData() = data_. Parallel run: currentData() = distributed_data_.
    bool isPreAdapted = false; // 0
    for (const auto& preAdaptGrid : currentData()) {
        isPreAdapted = std::max(isPreAdapted, preAdaptGrid -> preAdapt()); // could be 0 or 1
    }
    // If at least one process has marked elements, return true.
    return this->comm().max(isPreAdapted);
}

bool CpGrid::adapt()
{
    if(!preAdapt()) { // marked cells set can be empty
        return false; // the grid does not change at all.
    }

    const std::vector<std::array<int,3>>& cells_per_dim_vec = {{2,2,2}}; // Arbitrary chosen values.
    std::vector<int> assignRefinedLevel(current_view_data_-> size(0));
    const auto& preAdaptMaxLevel = this ->maxLevel();

    int local_marked_elem_count = 0;
    for (int elemIdx = 0; elemIdx < current_view_data_->size(0); ++elemIdx) {
        const auto& element = cpgrid::Entity<0>(*current_view_data_, elemIdx, true);
        assignRefinedLevel[elemIdx] = (this->getMark(element) == 1) ? (preAdaptMaxLevel +1) : 0;
        if (this->getMark(element) == 1) {
            ++local_marked_elem_count;
        }
    }

    std::vector<std::string> lgr_name_vec = { "LGR" + std::to_string(preAdaptMaxLevel +1) };

    auto global_marked_elem_count = comm().sum(local_marked_elem_count);
    auto global_cell_count_before_adapt = comm().sum(current_view_data_-> size(0)); // Recall overlap cells are also marked
    // Check if its a global refinement
    bool is_global_refine = (global_marked_elem_count == global_cell_count_before_adapt);
    if (is_global_refine) { // parallel or sequential
        // Rewrite the lgr name (GR stands for GLOBAL REFINEMET)
        lgr_name_vec = { "GR" + std::to_string(preAdaptMaxLevel +1) };
        const std::array<int,3>& endIJK = currentData().back()->logicalCartesianSize();
        return this->refineAndUpdateGrid(cells_per_dim_vec, assignRefinedLevel, lgr_name_vec, {{0,0,0}}, {endIJK});
    }
    return this-> refineAndUpdateGrid(cells_per_dim_vec, assignRefinedLevel, lgr_name_vec);
}

bool CpGrid::refineAndUpdateGrid(const std::vector<std::array<int,3>>& cells_per_dim_vec,
                                 const std::vector<int>& assignRefinedLevel,
                                 const std::vector<std::string>& lgr_name_vec,
                                 const std::vector<std::array<int,3>>& startIJK_vec,
                                 const std::vector<std::array<int,3>>& endIJK_vec)
{
    // To do: support coarsening.
    assert( static_cast<int>(assignRefinedLevel.size()) == currentData().back()->size(0));
    assert(cells_per_dim_vec.size() == lgr_name_vec.size());

    auto& data = currentData(); // data pointed by current_view_data_ (data_ or distributed_data_[if loadBalance() has been invoked before adapt()]).
    // Logical Cartesian Size before adapting the grid - to be used in case the entire grid will be refined.
    const auto& lcs =  data.back()->logical_cartesian_size_;

    bool isCARFIN = !startIJK_vec.empty();
    bool isGlobalRefine = isCARFIN; // One way of global-refine the grid is via startIJK and endIJK.
    for (int c = 0; c<3; ++c) {
        isGlobalRefine = isGlobalRefine && (startIJK_vec[0][c] == 0 ) && (endIJK_vec[0][c] == lcs[c]);
        if (!isGlobalRefine)
            break;
    }

    // Each marked element has its assigned level where its refined entities belong.
    const int& levels = cells_per_dim_vec.size();
    // Notice that "levels" represents also the total amount of new (after calling adapt) refined level grids.
    const int& preAdaptMaxLevel = this->maxLevel();
    // Copy corner history - needed to compute later ids, empty vector if the grid to be adapted is level 0 grid, or the grid has been distributed.
    const auto& preAdaptGrid_corner_history = (preAdaptMaxLevel>0) ? current_view_data_->corner_history_ : std::vector<std::array<int,2>>();

    if (!global_id_set_ptr_) {
        global_id_set_ptr_ = std::make_shared<cpgrid::GlobalIdSet>(*data.back());

        for (int preAdaptLevel = 0; preAdaptLevel <= preAdaptMaxLevel; ++preAdaptLevel) {
            global_id_set_ptr_->insertIdSet(*data[preAdaptLevel]);
        }
    }

    // To determine if an LGR is not empty in a given process, we set
    // lgr_with_at_least_one_active_cell[in that level] to 1 if it contains
    // at least one active cell, and to 0 otherwise.
    std::vector<int> lgr_with_at_least_one_active_cell(levels);
    Opm::Lgr::detectActiveLgrs(*this, startIJK_vec, endIJK_vec, lgr_with_at_least_one_active_cell);

    // To store/build refined level grids.
    std::vector<std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>> refined_data_vec(levels, data);
    std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>> refined_grid_ptr_vec(levels);

    std::vector<Dune::cpgrid::DefaultGeometryPolicy> refined_geometries_vec(levels);
    std::vector<std::vector<std::array<int,8>>> refined_cell_to_point_vec(levels);
    std::vector<cpgrid::OrientedEntityTable<0,1>> refined_cell_to_face_vec(levels);
    std::vector<Opm::SparseTable<int>> refined_face_to_point_vec(levels);
    std::vector<cpgrid::OrientedEntityTable<1,0>> refined_face_to_cell_vec(levels);

    // Mutable containers for refined corners, faces, cells, face tags, and face normals.
    std::vector<Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<0,3>>> refined_corners_vec(levels);
    std::vector<Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<2,3>>> refined_faces_vec(levels);
    std::vector<Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<3,3>>> refined_cells_vec(levels);
    std::vector<Dune::cpgrid::EntityVariableBase<enum face_tag>> mutable_refined_face_tags_vec(levels);
    typedef Dune::FieldVector<double,3> PointType;
    std::vector<Dune::cpgrid::EntityVariableBase<PointType>> mutable_refined_face_normals_vec(levels);

    std::vector<std::vector<int>> refined_global_cell_vec(levels);


    // To store adapted grid
    std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& adaptedData = data;
#if HAVE_MPI
    auto adaptedGrid_ptr =
        std::make_shared<Dune::cpgrid::CpGridData>((*(data[0])).ccobj_, adaptedData);
#else
    // DUNE 2.7 is missing convertion to NO_COMM
    auto adaptedGrid_ptr = std::make_shared<Dune::cpgrid::CpGridData>(adaptedData);
#endif
    auto& adaptedGrid = *adaptedGrid_ptr;
    Dune::cpgrid::DefaultGeometryPolicy&                         adapted_geometries = adaptedGrid.geometry_;
    std::vector<std::array<int,8>>&                              adapted_cell_to_point = adaptedGrid.cell_to_point_;
    cpgrid::OrientedEntityTable<0,1>&                            adapted_cell_to_face = adaptedGrid.cell_to_face_;
    Opm::SparseTable<int>&                                       adapted_face_to_point = adaptedGrid.face_to_point_;
    cpgrid::OrientedEntityTable<1,0>&                            adapted_face_to_cell = adaptedGrid.face_to_cell_;
    cpgrid::EntityVariable<enum face_tag,1>&                     adapted_face_tags = adaptedGrid.face_tag_;
    cpgrid::SignedEntityVariable<Dune::FieldVector<double,3>,1>& adapted_face_normals = adaptedGrid.face_normals_;
    // Mutable containers for adapted corners, faces, cells, face tags, and face normals.
    Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<0,3>>& adapted_corners =
        *(adapted_geometries.geomVector(std::integral_constant<int,3>()));
    Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<2,3>>& adapted_faces =
        *(adapted_geometries.geomVector(std::integral_constant<int,1>()));
    Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<3,3>>& adapted_cells =
        *(adapted_geometries.geomVector(std::integral_constant<int,0>()));
    Dune::cpgrid::EntityVariableBase<enum face_tag>& mutable_face_tags = adapted_face_tags;
    typedef Dune::FieldVector<double,3> PointType;
    Dune::cpgrid::EntityVariableBase<PointType>& mutable_face_normals = adapted_face_normals;



    // Refine marked cells and provide marked-corner/face/cell - refined-corner/faces/cells relations.
    //
    // ------------------------ Marked elements parameters
    // -- markedElem_to_itsLgr :
    // Each marked element gets refined and we store this "auxiliary markedElementLGR", to later
    // build a unique level containing all the refined entities from all the marked elements.
    std::vector<std::shared_ptr<Dune::cpgrid::CpGridData> > markedElem_to_itsLgr;
    markedElem_to_itsLgr.resize(current_view_data_->size(0));
    // -- markedElem_count: Total amount of marked elements to be refined. It will be used to print grid info.
    int markedElem_count = 0;
    // -- cornerInMarkedElemWithEquivRefinedCorner :
    // For each corner from level zero, we store the marked elements where the corner appears and its equivalent
    // refined corner in  each auxiliary marked-element-lgr. Example: corner with index 5 appears in marked
    // elements 0 and 1, with refined equivalent corner indices 8 and 2 respectively. Then,
    // cornerInMarkedElemWithEquivRefinedCorner[5] = {{0, 8}, {1, 2}}.
    // For corners not appearing in any marked element, empty vector.
    std::vector<std::vector<std::array<int,2>>> cornerInMarkedElemWithEquivRefinedCorner;
    cornerInMarkedElemWithEquivRefinedCorner.resize(current_view_data_->size(3) );
    // -- markedElemAndEquivRefinedCorner_to_corner :
    // To correctly build the level-refined and adapted-grid topology features, we need to keep track of the
    // corners that got replaced by equivalent refined corners, in each marked element where the corner appeared,
    // not only in its last appearance. The last appearance will be used to avoid repetition when storing.
    // Following the example above,
    // markedElemAndEquivRefinedCorner_to_corner[{0, 8}] = 5;
    // markedElemAndEquivRefinedCorner_to_corner[{1, 2}] = 5;
    std::map< std::array<int,2>, int > markedElemAndEquivRefinedCorn_to_corner;
    // -- faceInMarkedElemAndRefinedFaces :
    // For each face from level zero, we store the marked elements where the face appears (maximum 2 cells)
    // and its new-born refined faces from each auxiliary marked-element-lgr. Example: face with index 9
    // appears in marked elements 0 and 1. Then,
    // faceInMarkedElemAndRefinedFaces[9] = {{0, {refinedFace0_0, ..., refinedFaceN_0}},
    //                                       {1, {refinedFace0_1, ..., refinedFaceM_1}}}.
    // For faces not appearing in any marked element, empty vector.
    std::vector<std::vector<std::pair<int, std::vector<int>>>> faceInMarkedElemAndRefinedFaces;
    faceInMarkedElemAndRefinedFaces.resize(current_view_data_->face_to_cell_.size());
    // ------------------------ Refined cells parameters
    // --- Refined cells and PreAdapt cells relations ---
    std::map<std::array<int,2>,std::array<int,2>> elemLgrAndElemLgrCell_to_refinedLevelAndRefinedCell;
    std::map<std::array<int,2>,std::array<int,2>> refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell;
    // Integer to count only REFINED cells (new-born refined cells from ANY marked element).
    std::vector<int> refined_cell_count_vec(levels, 0);
    // -- Parent-child relations --
    // Relation between the grids before adapt() and new refined cells (on the refined grid level - not in each individual lgr).
    std::vector<std::vector<std::tuple<int,std::vector<int>>>> preAdapt_parent_to_children_cells_vec(preAdaptMaxLevel +1);
    // ------------------------ Adapted cells parameters
    // --- Adapted cells and PreAdapt cells relations ---
    std::map<std::array<int,2>,int>           elemLgrAndElemLgrCell_to_adaptedCell;
    std::unordered_map<int,std::array<int,2>> adaptedCell_to_elemLgrAndElemLgrCell;
    // Integer to count adapted cells (mixed between cells from level0 (not involved in LGRs), and (new-born) refined cells).
    int cell_count = 0;
    // -- Some extra indices relations between preAdapt-grid and adapted-grid --
    // Relation between the grids before adapt() and leafview cell indices.
    std::vector<std::vector<int>> preAdapt_level_to_leaf_cells_vec(preAdaptMaxLevel +1);
    for (int preAdaptLevel = 0; preAdaptLevel < preAdaptMaxLevel +1; ++preAdaptLevel) {
        // Resize with the corresponding amount of cells of the preAdapt level. Deafualt {-1, empty vector} when the cell has no children.
        if ( (*data[preAdaptLevel]).parent_to_children_cells_.empty()){
            preAdapt_parent_to_children_cells_vec[preAdaptLevel].resize(data[preAdaptLevel]->size(0), std::make_pair(-1, std::vector<int>{}));
        }
        else {
            preAdapt_parent_to_children_cells_vec[preAdaptLevel] =  (*data[preAdaptLevel]).parent_to_children_cells_;
        }
        // Resize with the corresponding amount of cell of the preAdapt level. Dafualt -1 when the cell vanished and does not appear on the leaf grid view.
        // In entry 'level cell index', we store 'leafview cell index', or -1 when the cell vanished.
        preAdapt_level_to_leaf_cells_vec[preAdaptLevel].resize(data[preAdaptLevel]->size(0), -1);
    }
    //
    Opm::Lgr::refineAndProvideMarkedRefinedRelations( *this,
                                                      /* Marked elements parameters */
                                                      markedElem_to_itsLgr,
                                                      markedElem_count,
                                                      cornerInMarkedElemWithEquivRefinedCorner,
                                                      markedElemAndEquivRefinedCorn_to_corner,
                                                      faceInMarkedElemAndRefinedFaces,
                                                      /* Refined cells parameters */
                                                      elemLgrAndElemLgrCell_to_refinedLevelAndRefinedCell,
                                                      refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                                      refined_cell_count_vec,
                                                      assignRefinedLevel,
                                                      preAdapt_parent_to_children_cells_vec,
                                                      /* Adapted cells parameters */
                                                      elemLgrAndElemLgrCell_to_adaptedCell,
                                                      adaptedCell_to_elemLgrAndElemLgrCell,
                                                      cell_count,
                                                      preAdapt_level_to_leaf_cells_vec,
                                                      /* Additional parameters */
                                                      cells_per_dim_vec);

#if HAVE_MPI
    auto global_markedElem_count = comm().sum(markedElem_count);
    if ( global_markedElem_count == 0 ) {
        return false;
    }
#else
    if ( markedElem_count == 0 ) {
        return false;
    }
#endif

    // Update/define parent_to_children_cells_ and level_to_leaf_cells_ for all the existing level grids (level 0, 1, ..., preAdaptMaxLevel), before this call of adapt.
    for (int preAdaptLevel = 0; preAdaptLevel < preAdaptMaxLevel +1; ++preAdaptLevel) {
        (*data[preAdaptLevel]).parent_to_children_cells_ = preAdapt_parent_to_children_cells_vec[preAdaptLevel];
        (*data[preAdaptLevel]).level_to_leaf_cells_ =  preAdapt_level_to_leaf_cells_vec[preAdaptLevel];
    }

    // -- Child-parent relations --
    // refined_child_to_parent_cells_vec:   Refined child cells and their parents. Entry is {-1,-1} when cell has no father.
    //                                      Otherwise, {level parent cell, parent cell index}. Each entry represents a refined level.
    // refined_cell_to_idxInParentCell_vec: Each refined child cell has a unique index in its parent cell, to be used to build geometryInFather().
    //                                      Each entry represents a refined level.
    // adapted_child_to_parent_cells:       Adapted child cells and their parents. Entry is {-1,-1} when cell has no father. Otherwise, {level parent cell, parent cell index}
    // adapted_cell_to_idxInParentCell:     Each refined child cell has a unique index in its parent cell, to be used to build geometryInFather().
    //                                      When the cell has not been refined, -1.
    const auto& [refined_child_to_parent_cells_vec,
                 refined_cell_to_idxInParentCell_vec,
                 adapted_child_to_parent_cells,
                 adapted_cell_to_idxInParentCell] = Opm::Lgr::defineChildToParentAndIdxInParentCell(*currentData().back(),
                                                                                                    preAdaptMaxLevel,
                                                                                                    refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                                                                                    refined_cell_count_vec,adaptedCell_to_elemLgrAndElemLgrCell,
                                                                                                    cell_count);

    // -- Refined to Adapted cells and Adapted-cells to {level where the cell was born, cell index on that level} --
    // refined_level_to_leaf_cells_vec:  Relation between the refined grid and leafview cell indices.
    // leaf_to_level_cells:              Relation between an adapted cell and its equivalent cell coming either from current_view_data_ or from the refined grid (level)
    const auto& [refined_level_to_leaf_cells_vec,
                 leaf_to_level_cells] = Opm::Lgr::defineLevelToLeafAndLeafToLevelCells(*currentData().back(),
                                                                                       preAdaptMaxLevel,
                                                                                       elemLgrAndElemLgrCell_to_refinedLevelAndRefinedCell,
                                                                                       refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                                                                       refined_cell_count_vec,
                                                                                       elemLgrAndElemLgrCell_to_adaptedCell,
                                                                                       adaptedCell_to_elemLgrAndElemLgrCell,
                                                                                       cell_count);

    // CORNERS
    // Stablish relationships between PreAdapt corners and refined or adapted ones ---
    //
    // --- Refined corners and PreAdapt corners relations ---
    std::map<std::array<int,2>,std::array<int,2>> elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner;
    std::map<std::array<int,2>,std::array<int,2>> refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner;
    std::map<std::array<int,2>,std::array<int,2>> vanishedRefinedCorner_to_itsLastAppearance;
    // Integer to count only refined corners.
    std::vector<int> refined_corner_count_vec(levels, 0);
    Opm::Lgr::identifyRefinedCornersPerLevel(*currentData().back(),
                                             preAdaptMaxLevel,
                                             /* Refined grid parameters */
                                             elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                                             refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner,
                                             refined_corner_count_vec,
                                             vanishedRefinedCorner_to_itsLastAppearance,
                                             /* Additional parameters */
                                             markedElem_to_itsLgr,
                                             assignRefinedLevel,
                                             cornerInMarkedElemWithEquivRefinedCorner,
                                             faceInMarkedElemAndRefinedFaces,
                                             cells_per_dim_vec);

    // --- Adapted corners and PreAdapt corners relations ---
    std::map<std::array<int,2>,int>           elemLgrAndElemLgrCorner_to_adaptedCorner;
    std::unordered_map<int,std::array<int,2>> adaptedCorner_to_elemLgrAndElemLgrCorner;
    // Integer to count adapted corners (mixed between corners from current_view_data_ (not involved in LGRs), and (new-born) refined corners).
    int corner_count = 0;
    Opm::Lgr::identifyLeafGridCorners(*currentData().back(),
                                      preAdaptMaxLevel,
                                      /* Adapted grid parameters */
                                      elemLgrAndElemLgrCorner_to_adaptedCorner,
                                      adaptedCorner_to_elemLgrAndElemLgrCorner,
                                      corner_count,
                                      /* Additional parameters */
                                      markedElem_to_itsLgr,
                                      assignRefinedLevel,
                                      cornerInMarkedElemWithEquivRefinedCorner,
                                      vanishedRefinedCorner_to_itsLastAppearance,
                                      faceInMarkedElemAndRefinedFaces,
                                      cells_per_dim_vec);

    // FACES
    // Stablish relationships between PreAdapt faces and refined or adapted ones ---
    // --- Refined faces and PreAdapt faces relations ---
    std::map<std::array<int,2>,std::array<int,2>> elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace;
    std::map<std::array<int,2>,std::array<int,2>> refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace;
    // Integer to count adapted faces (mixed between faces from level0 (not involved in LGRs), and (new-born) refined faces).
    std::vector<int> refined_face_count_vec(levels, 0);
    Opm::Lgr::identifyRefinedFacesPerLevel(*currentData().back(),
                                           preAdaptMaxLevel,
                                           elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace,
                                           refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace,
                                           refined_face_count_vec,
                                           markedElem_to_itsLgr,
                                           assignRefinedLevel,
                                           faceInMarkedElemAndRefinedFaces,
                                           cells_per_dim_vec);

    // --- Adapted faces and PreAdapt faces relations ---
    std::map< std::array<int,2>, int >           elemLgrAndElemLgrFace_to_adaptedFace;
    std::unordered_map< int, std::array<int,2> > adaptedFace_to_elemLgrAndElemLgrFace;
    // Integer to count adapted faces (mixed between faces from current_view_data_ (not involved in LGRs), and (new-born) refined faces).
    int face_count = 0;
    Opm::Lgr::identifyLeafGridFaces(*currentData().back(),
                                    preAdaptMaxLevel,
                                    elemLgrAndElemLgrFace_to_adaptedFace,
                                    adaptedFace_to_elemLgrAndElemLgrFace,
                                    face_count,
                                    markedElem_to_itsLgr,
                                    assignRefinedLevel,
                                    faceInMarkedElemAndRefinedFaces,
                                    cells_per_dim_vec);

    // Set refined level grids geometries
    // --- Refined corners  ---
    Opm::Lgr::populateRefinedCorners(refined_corners_vec,
                                     refined_corner_count_vec,
                                     markedElem_to_itsLgr,
                                     preAdaptMaxLevel,
                                     refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner);
    // --- Refined faces  ---
    Opm::Lgr::populateRefinedFaces(refined_faces_vec,
                                   mutable_refined_face_tags_vec,
                                   mutable_refined_face_normals_vec,
                                   refined_face_to_point_vec,
                                   refined_face_count_vec,
                                   refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace,
                                   elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                                   vanishedRefinedCorner_to_itsLastAppearance,
                                   markedElem_to_itsLgr,
                                   preAdaptMaxLevel,
                                   cornerInMarkedElemWithEquivRefinedCorner,
                                   markedElemAndEquivRefinedCorn_to_corner);
    // --- Refined cells  ---
    Opm::Lgr::populateRefinedCells(*currentData().back(),
                                   refined_cells_vec,
                                   refined_cell_to_point_vec,
                                   refined_global_cell_vec,
                                   refined_cell_count_vec,
                                   refined_cell_to_face_vec,
                                   refined_face_to_cell_vec,
                                   refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                   elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace,
                                   faceInMarkedElemAndRefinedFaces,
                                   refined_geometries_vec,
                                   elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                                   vanishedRefinedCorner_to_itsLastAppearance,
                                   markedElem_to_itsLgr,
                                   assignRefinedLevel,
                                   preAdaptMaxLevel,
                                   markedElemAndEquivRefinedCorn_to_corner,
                                   cornerInMarkedElemWithEquivRefinedCorner,
                                   cells_per_dim_vec);

    // Update leaf grid geometries
    // --- Adapted corners ---
    Opm::Lgr::populateLeafGridCorners(*currentData().back(),
                                      adapted_corners,
                                      corner_count,
                                      markedElem_to_itsLgr,
                                      adaptedCorner_to_elemLgrAndElemLgrCorner);
    // --- Adapted faces ---
    Opm::Lgr::populateLeafGridFaces(*currentData().back(),
                                    adapted_faces,
                                    mutable_face_tags,
                                    mutable_face_normals,
                                    adapted_face_to_point,
                                    face_count,
                                    adaptedFace_to_elemLgrAndElemLgrFace,
                                    elemLgrAndElemLgrCorner_to_adaptedCorner,
                                    vanishedRefinedCorner_to_itsLastAppearance,
                                    markedElem_to_itsLgr,
                                    assignRefinedLevel,
                                    markedElemAndEquivRefinedCorn_to_corner,
                                    cornerInMarkedElemWithEquivRefinedCorner,
                                    cells_per_dim_vec,
                                    preAdaptMaxLevel);
    // --- Adapted cells ---
    Opm::Lgr::populateLeafGridCells(*currentData().back(),
                                    adapted_cells,
                                    adapted_cell_to_point,
                                    cell_count,
                                    adapted_cell_to_face,
                                    adapted_face_to_cell,
                                    adaptedCell_to_elemLgrAndElemLgrCell,
                                    elemLgrAndElemLgrFace_to_adaptedFace,
                                    faceInMarkedElemAndRefinedFaces,
                                    adapted_geometries,
                                    elemLgrAndElemLgrCorner_to_adaptedCorner,
                                    vanishedRefinedCorner_to_itsLastAppearance,
                                    markedElem_to_itsLgr,
                                    assignRefinedLevel,
                                    markedElemAndEquivRefinedCorn_to_corner,
                                    cornerInMarkedElemWithEquivRefinedCorner,
                                    cells_per_dim_vec,
                                    preAdaptMaxLevel);

    for (int level = 0; level < levels; ++level) {
        const int refinedLevelGridIdx = level + preAdaptMaxLevel +1;
#if HAVE_MPI
        refined_grid_ptr_vec[level] = std::make_shared<Dune::cpgrid::CpGridData>((*(data[0])).ccobj_, refined_data_vec[level]);
#else
        // DUNE 2.7 is missing convertion to NO_COMM
        refined_grid_ptr_vec[level] = std::make_shared<Dune::cpgrid::CpGridData>(refined_data_vec[level]);
#endif
        // Store refined grid
        if ((level == 0) && (preAdaptMaxLevel>0)) { // Overwrite the leaf-grid-view with the first new-refined-level-grid
            data[preAdaptMaxLevel+1] = refined_grid_ptr_vec[level];
        }
        else {
            data.push_back(refined_grid_ptr_vec[level]);
        }

        Dune::cpgrid::DefaultGeometryPolicy&  refinedLevel_geometries = (*data[refinedLevelGridIdx]).geometry_;
        // Mutable containers for adapted corners, faces, cells, face tags, and face normals.
        Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<0,3>>& level_corners =
            *(refinedLevel_geometries.geomVector(std::integral_constant<int,3>()));
        Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<2,3>>& level_faces =
            *(refinedLevel_geometries.geomVector(std::integral_constant<int,1>()));
        Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<3,3>>& level_cells =
            *(refinedLevel_geometries.geomVector(std::integral_constant<int,0>()));

        level_corners.swap(refined_corners_vec[level]);
        level_faces.swap(refined_faces_vec[level]);
        level_cells.swap(refined_cells_vec[level]);

        (*data[refinedLevelGridIdx]).cell_to_point_.swap(refined_cell_to_point_vec[level]);
        (*data[refinedLevelGridIdx]).cell_to_face_.swap(refined_cell_to_face_vec[level]);

        (*data[refinedLevelGridIdx]).face_to_point_.swap(refined_face_to_point_vec[level]);
        (*data[refinedLevelGridIdx]).face_to_cell_.swap(refined_face_to_cell_vec[level]);

        cpgrid::EntityVariable<enum face_tag,1>& level_face_tags =   (*data[refinedLevelGridIdx]).face_tag_;
        Dune::cpgrid::EntityVariableBase<enum face_tag>& level_mutable_face_tags = level_face_tags;
        level_mutable_face_tags.swap(mutable_refined_face_tags_vec[level]);

        cpgrid::SignedEntityVariable<Dune::FieldVector<double,3>,1>&  level_face_normals =   (*data[refinedLevelGridIdx]).face_normals_;
        Dune::cpgrid::EntityVariableBase<PointType>& level_mutable_face_normals = level_face_normals;
        level_mutable_face_normals.swap(mutable_refined_face_normals_vec[level]);

        // Further Refined grid Attributes
        //
        // Populate some attributes of the level LGR
        (*data[refinedLevelGridIdx]).level_data_ptr_ = &(this -> currentData());
        (*data[refinedLevelGridIdx]).level_ = refinedLevelGridIdx;
        this -> lgr_names_[lgr_name_vec[level]] = refinedLevelGridIdx; // {"name_lgr", level}
        (*data[refinedLevelGridIdx]).child_to_parent_cells_ = refined_child_to_parent_cells_vec[level];
        (*data[refinedLevelGridIdx]).cell_to_idxInParentCell_ = refined_cell_to_idxInParentCell_vec[level];
        (*data[refinedLevelGridIdx]).level_to_leaf_cells_ =  refined_level_to_leaf_cells_vec[level];
        (*data[refinedLevelGridIdx]).index_set_ = std::make_unique<cpgrid::IndexSet>(data[refinedLevelGridIdx]->size(0),
                                                                                     data[refinedLevelGridIdx]->size(3));
        (*data[refinedLevelGridIdx]).refinement_max_level_ = levels + preAdaptMaxLevel;
        // Determine the amount of cells per direction, per parent cell, of the corresponding LGR.
        (*data[refinedLevelGridIdx]).cells_per_dim_ = cells_per_dim_vec[level];
        // TO DO: This new code for refinement do not assume Cartesian Shape. How does logical_cartesian_size_ should be defined then?
        // When the refined level grid has been originated from a block of cells, then its logical Cartesian size
        // corresponds to the inner product between cells_per_dim_vec[level] and the dimension of the block (amount of cells in each direction).
        // In the case of a block of cells, e.g., when CARFIN keyword is used, we need the following:
        if (isCARFIN) {
            const auto& blockDim = Opm::Lgr::getPatchDim(startIJK_vec[level], endIJK_vec[level]);
            (*data[refinedLevelGridIdx]).logical_cartesian_size_ = { cells_per_dim_vec[level][0]*blockDim[0],
                                                                     cells_per_dim_vec[level][1]*blockDim[1],
                                                                     cells_per_dim_vec[level][2]*blockDim[2] };
        }
        else {
            (*data[refinedLevelGridIdx]).logical_cartesian_size_ = (*data[0]).logical_cartesian_size_;
            (*data[refinedLevelGridIdx]).global_cell_.swap(refined_global_cell_vec[level]);
        }
        // One alternative definition for logical_cartesian_size_ in the case where the marked elements for refinement do not form a block of cells,
        // therefore, are not associated with the keyword CARFIN, is to imagine that we put all the marked elements one next to the other, along
        // the x-axis. Then, the "imaginary" logical Cartesian size of the refined level grid would be
        // { (# marked elemnts)x cells_per_dim_vec[level][0], cells_per_dim_vec[level][1], cells_per_dim_vec[level][2]}.
        /** To do: how the definition of refined level grids logical_cartesian_size_ affects LookUpData class (and LookUpCartesianData)*/
    }

    // Store adapted grid
    data.push_back(adaptedGrid_ptr);

    // Further Adapted  grid Attributes
    (*data[levels + preAdaptMaxLevel +1]).child_to_parent_cells_ = adapted_child_to_parent_cells;
    (*data[levels + preAdaptMaxLevel +1]).cell_to_idxInParentCell_ = adapted_cell_to_idxInParentCell;
    (*data[levels + preAdaptMaxLevel +1]).leaf_to_level_cells_ =  leaf_to_level_cells;
    (*data[levels + preAdaptMaxLevel +1]).index_set_ = std::make_unique<cpgrid::IndexSet>(data[levels + preAdaptMaxLevel +1]->size(0),
                                                                                          data[levels + preAdaptMaxLevel +1]->size(3));
    (*data[levels + preAdaptMaxLevel +1]).refinement_max_level_ = levels + preAdaptMaxLevel;

    if (isGlobalRefine) {
        assert(cells_per_dim_vec.size() == 1);
        (*data[levels + preAdaptMaxLevel +1]).logical_cartesian_size_ =  { lcs[0]*cells_per_dim_vec[0][0],
                                                                           lcs[1]*cells_per_dim_vec[0][1],
                                                                           lcs[2]*cells_per_dim_vec[0][2] };
    }
    else {
        (*data[levels + preAdaptMaxLevel +1]).logical_cartesian_size_ =  (*data[0]).logical_cartesian_size_;
    }

    // Update the leaf grid view
    current_view_data_ = data.back().get();

    // When the refinement is determined by startIJK and endIJK values, the LGR has a (local) Cartesian size.
    // Therefore, each refined cell belonging to the LGR can be associated with a (local) IJK and its (local) Cartesian index.
    // If the LGR has NXxNYxNZ dimension, then the Cartesian indices take values
    // k*NN*NY + j*NX + i, where i<NX, j<Ny, k<NZ.
    // This index is stored in <refined-level-grid>.global_cell_[ refined cell index (~element.index()) ] =  k*NN*NY + j*NX + i.
    if (isCARFIN) {
        for (int level = 0; level < levels; ++level) {
            const int refinedLevelGridIdx = level + preAdaptMaxLevel +1;
            std::vector<int> global_cell_lgr(data[refinedLevelGridIdx]->size(0));
            computeGlobalCellLgr(refinedLevelGridIdx, startIJK_vec[level], global_cell_lgr);
            (*data[refinedLevelGridIdx]).global_cell_.swap(global_cell_lgr);
        }
    }

    std::vector<int> global_cell_leaf( data[levels + preAdaptMaxLevel +1]->size(0));
    computeGlobalCellLeafGridViewWithLgrs(global_cell_leaf);
    (*data[levels + preAdaptMaxLevel +1]).global_cell_.swap(global_cell_leaf);

    updateCornerHistoryLevels(cornerInMarkedElemWithEquivRefinedCorner,
                              elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                              adaptedCorner_to_elemLgrAndElemLgrCorner,
                              corner_count,
                              preAdaptGrid_corner_history,
                              preAdaptMaxLevel,
                              levels);

    // Insert the new id sets into the grid global_id_set_ptr_
    for (int level = 0; level < levels; ++level) {
        const int refinedLevelGridIdx = level + preAdaptMaxLevel +1;
        this->global_id_set_ptr_->insertIdSet(*data[refinedLevelGridIdx]);
    }
    this->global_id_set_ptr_->insertIdSet(*data.back());

    // Only for parallel runs
    // - Define global ids for refined level grids (level 1, 2, ..., maxLevel)
    // - Define GlobalIdMapping (cellMapping, faceMapping, pointMapping required per level)
    // - Define ParallelIndex for overlap cells and their neighbors
    if(comm().size()>1) {
        globalIdsPartitionTypesLgrAndLeafGrids(assignRefinedLevel,
                                               cells_per_dim_vec,
                                               lgr_with_at_least_one_active_cell);
    }

    // Print total amount of cells on the adapted grid
    Opm::OpmLog::info(std::to_string(markedElem_count) + " elements have been marked (in " + std::to_string(comm().rank()) + " rank).\n");
    Opm::OpmLog::info(std::to_string(levels)  + " (new) refined level grid(s) (in " + std::to_string(comm().rank()) + " rank).\n");
    Opm::OpmLog::info(std::to_string(cell_count)  + " total cells on the leaf grid view (in " + std::to_string(comm().rank()) + " rank).\n");

    return (preAdaptMaxLevel < this->maxLevel()); // true if at least one entity was refined
}

void CpGrid::postAdapt()
{
    // - Resize with the new amount of cells on the leaf grid view
    // - Set marks equal to zero (representing 'doing nothing')
    current_view_data_ -> postAdapt();
}

void CpGrid::globalIdsPartitionTypesLgrAndLeafGrids([[maybe_unused]] const std::vector<int>& assignRefinedLevel,
                                                    [[maybe_unused]] const std::vector<std::array<int,3>>& cells_per_dim_vec,
                                                    [[maybe_unused]] const std::vector<int>& lgr_with_at_least_one_active_cell)
{
#if HAVE_MPI
    // Prediction min cell and point global ids per process
    //
    // Predict how many new cells/points (born in refined level grids) need new globalIds, so we can assign unique
    // new ids ( and anticipate the maximum).
    // The grid is already refined according to the LGR specification
    // At this point, neither cell_index_set_ nor partition_type_indicator_ are populated.
    // Refined level grid cells:
    //    1. Inherit their partition type from their parent cell (i.e., element.father().partitionType()).
    //    2. Assign global ids only for interior cells.
    //    3. Communicate the already assigned cell global ids from interior to overlap refined cells.
    // Refined level grid points/vertices:
    // There are 4 partition types: interior, border, front, overlap. This classification requires that both
    // cell and face partition types are already defined, not available yet for refined level grids.
    //    1. Assign for all partition type points a 'candidate of global id' (unique in each process).
    //       Except the points that coincide with a point from level zero.
    //    2. Re-write the values for points that are corners of overlap refined cells, via communication.
    // Under the assumption of LGRs fully-interior, no communication is needed. In the general case, communication will be used
    // to populate overlap cell/point global ids on the refined level grids.
    /** Warning: due to the overlap layer size (equal to 1) cells that share corners or edges (not faces) with interior cells
        are not included/seen by the process. This, in some cases, ends up in multiple ids for the same point. */

    int min_globalId_cell_in_proc = 0;
    int min_globalId_point_in_proc = 0;
    Opm::Lgr::predictMinCellAndPointGlobalIdPerProcess(*this,
                                                       assignRefinedLevel,
                                                       cells_per_dim_vec,
                                                       lgr_with_at_least_one_active_cell,
                                                       min_globalId_cell_in_proc,
                                                       min_globalId_point_in_proc);

    // Only for level 1,2,.., maxLevel grids.
    // For each level, define the local-to-global maps for cells and points (for faces: empty).
    // 1) Assignment of new global ids is done only for owned cells and non-overlap points.
    // 2) For overlap cells and points: communicate. Not needed under the assumption of fully interior LGRs.
    std::vector<std::vector<int>> localToGlobal_cells_per_level(cells_per_dim_vec.size());
    std::vector<std::vector<int>> localToGlobal_points_per_level(cells_per_dim_vec.size());
    // Ignore faces - empty vectors.
    std::vector<std::vector<int>> localToGlobal_faces_per_level(cells_per_dim_vec.size());

    Opm::Lgr::assignCellIdsAndCandidatePointIds(*this,
                                                localToGlobal_cells_per_level,
                                                localToGlobal_points_per_level,
                                                min_globalId_cell_in_proc,
                                                min_globalId_point_in_proc,
                                                cells_per_dim_vec);


    const auto& parent_to_children = current_data_->front()->parent_to_children_cells_;
    ParentToChildrenCellGlobalIdHandle parentToChildrenGlobalId_handle(parent_to_children, localToGlobal_cells_per_level);
    currentData().front()->communicate(parentToChildrenGlobalId_handle,
                                       Dune::InteriorBorder_All_Interface,
                                       Dune::ForwardCommunication );

    // After assigning global IDs to points in refined-level grids, a single point may have
    // a "unique" global ID in each local leaf grid view for every process to which it belongs.
    // To ensure true uniqueness, since global IDs must be distinct across the global leaf view
    // and consistent across each refined-level grid, we will rewrite the entries in
    // localToGlobal_points_per_level.
    //
    // This correction is done using cell_to_point_ across all refined cells through
    // communication: gathering the 8 corner points of each interior cell and scattering the
    // 8 corner points of overlapping cells, for all child cells of a parent cell in level zero grid.
    //
    // Design decision: Why we communicate via level zero grid instead of in each refined level grid.
    // The reason is that how children are stored (the ordering) in parent_to_children_cells_
    // is always the same, accross all processes.
    // Even though the ordering of the corners in cell_to_point_ is the same accross all processes,
    // this may not be enough to correctly overwrite the "winner" point global ids for refined cells.
    //
    /** Current approach avoids duplicated point ids when
     // 1. the LGR is distributed in P_{i_0}, ..., P_{i_n}, with n+1 < comm().size(),
     // AND
     // 2. there is no coarse cell seen by a process P with P != P_{i_j}, j = 0, ..., n.
     // Otherwise, there will be duplicated point ids.
     //
     // Reason: neighboring cells that only share corners (not faces) are NOT considered in the
     // overlap layer of the process.*/
    Opm::Lgr::selectWinnerPointIds(*this,
                                   localToGlobal_points_per_level,
                                   parent_to_children,
                                   cells_per_dim_vec);

    for (std::size_t level = 1; level < cells_per_dim_vec.size()+1; ++level) {
        // Global id set for each (refined) level grid.
        if(lgr_with_at_least_one_active_cell[level-1]>0) { // Check if LGR is active in currect process.
            (*current_data_)[level]->global_id_set_->swap(localToGlobal_cells_per_level[level-1],
                                                          localToGlobal_faces_per_level[level-1],
                                                          localToGlobal_points_per_level[level-1]);
        }
    }

    for (std::size_t level = 1; level < cells_per_dim_vec.size()+1; ++level) {

        populateCellIndexSetRefinedGrid(level);
        // Compute the partition type for cell
        currentData()[level]->computeCellPartitionType();
        // Compute the partition type for point
        currentData()[level]->computePointPartitionType();
        // Now we can compute the communication interface.
        currentData()[level]->computeCommunicationInterfaces(currentData()[level]->size(3));
    }

    populateLeafGlobalIdSet();

    // Insert the new id sets into the grid global_id_set_ptr_
    for (std::size_t level = 0; level < cells_per_dim_vec.size()+1; ++level) {
        this->global_id_set_ptr_->insertIdSet(*(*current_data_)[level]);
    }
    this->global_id_set_ptr_->insertIdSet(*currentData().back());


    populateCellIndexSetLeafGridView();

    // Compute the partition type for cell
    (*current_data_).back()->computeCellPartitionType();

    // Compute the partition type for point
    (*current_data_).back()->computePointPartitionType();

    // Now we can compute the communication interface.
    current_data_->back()->computeCommunicationInterfaces(current_data_->back()->size(3));
    assert(static_cast<std::size_t>(current_data_->back()->cellIndexSet().size()) == static_cast<std::size_t>(current_data_->back()->size(0)) );
#endif
}

void CpGrid::getFirstChildGlobalIds([[maybe_unused]] std::vector<int>& parentToFirstChildGlobalIds)
{
#if HAVE_MPI
    switchToGlobalView();
    Opm::Lgr::getFirstChildGlobalIds(*this, parentToFirstChildGlobalIds);
#endif
}

void CpGrid::syncDistributedGlobalCellIds()
{
#if HAVE_MPI
    std::vector<int> parentToFirstChildGlobalIds;
    getFirstChildGlobalIds(parentToFirstChildGlobalIds);

    auto size = comm().max(parentToFirstChildGlobalIds.size());

    switchToDistributedView();

    parentToFirstChildGlobalIds.resize(size);
    comm().broadcast(parentToFirstChildGlobalIds.data(), parentToFirstChildGlobalIds.size(), 0);

    const int maxLevel = this->maxLevel();

    // Preallocate syncCellIds (and vertexIds, which will NOT be synchronized)
    std::vector<std::vector<int>> syncCellIds(maxLevel);
    std::vector<std::vector<int>> vertexIds(maxLevel);
    for (int level = 1; level <= maxLevel; ++level) {
        syncCellIds[level-1].resize(currentData()[level]->size(0));
        vertexIds[level-1].resize(currentData()[level]->size(3));
    }

    const auto& globalIdSet = this->globalIdSet();

    // Populate syncCellIds and vertexIds
    for (int level = 1; level <= maxLevel; ++level) {
        const auto& elements = Dune::elements(levelGridView(level));
        for (const auto& element : elements) {
            const int parent_globalId = globalIdSet.id(element.father());
            const int idx_in_parent = element.getIdxInParentCell();
            const int first_child_id = parentToFirstChildGlobalIds[parent_globalId];
            const int new_elem_globalId = first_child_id + idx_in_parent;

            syncCellIds[element.level()-1][element.index()] = new_elem_globalId;
        }

        for (const auto& vertex : Dune::vertices(levelGridView(level))){
            vertexIds[level-1][vertex.index()] = globalIdSet.id(vertex);
        }
    }

    // Re-assign new cell global ids for all refined level grids
    std::vector<int> faceIds; // empty for all
    for (int level = 1; level <= maxLevel; ++level) {
        if(currentData()[level]->size(0)) { // Check if LGR is active in currect process.
            currentData()[level]->global_id_set_->swap(syncCellIds[level-1],
                                                       faceIds,
                                                       vertexIds[level-1]);

            populateCellIndexSetRefinedGrid(level);
            // Insert the new id sets into the grid global_id_set_ptr_
            this->global_id_set_ptr_->insertIdSet(*(*current_data_)[level]);
        }
    }

    populateLeafGlobalIdSet();
    this->global_id_set_ptr_->insertIdSet(*currentData().back());

    assert(static_cast<std::size_t>(current_data_->back()->cellIndexSet().size()) == static_cast<std::size_t>(current_data_->back()->size(0)) );
#endif
}

void CpGrid::addLgrsUpdateLeafView(const std::vector<std::array<int,3>>& cells_per_dim_vec,
                                   const std::vector<std::array<int,3>>& startIJK_vec,
                                   const std::vector<std::array<int,3>>& endIJK_vec,
                                   const std::vector<std::string>& lgr_name_vec,
                                   const std::vector<std::string>& lgr_parent_grid_name_vec)
{
    // For parallel run, level zero grid is stored in distributed_data_[0]. If CpGrid::scatterGrid has been invoked,
    // then current_view_data_ == distributed_data_[0].
    // For serial run, level zero grid is stored in data_[0]. In this case, current_view_data_ == data_[0].
    // Note: currentData() returns data_ (if grid is not distributed) or distributed_data_ otherwise.

    Opm::Lgr::validStartEndIJKs(startIJK_vec, endIJK_vec);

    // If no parent grid name vector has been provided, then default "GLOBAL" for all (new) level grids.
    std::vector<std::string> parent_grid_names = lgr_parent_grid_name_vec;
    if (parent_grid_names.size() == 0){ // No parent grid name given->default "GLOBAL" parent grid
        parent_grid_names.resize(cells_per_dim_vec.size(), "GLOBAL");
    }

    // Sizes of provided vectors (number of subivisions per cells and lgrs name) should coincide.
    bool matchingSizeHasFailed = false;
    if ( (cells_per_dim_vec.size() != startIJK_vec.size()) ||
         (lgr_name_vec.size() != startIJK_vec.size()) ||
         (parent_grid_names.size() != startIJK_vec.size())) {
        matchingSizeHasFailed = true;
    }
    matchingSizeHasFailed = comm().max(matchingSizeHasFailed);
    if (matchingSizeHasFailed) {
        OPM_THROW(std::invalid_argument, "Sizes of provided vectors with subdivisions per cell and LGR names need to match.");
    }

    // Discard LGRs whose subdivisions do not trigger actual refinement, i.e., cells_per_dim_ = {1,1,1}
    const auto [filtered_cells_per_dim_vec,
                filtered_startIJK_vec,
                filtered_endIJK_vec,
                filtered_lgr_name_vec,
                filtered_lgr_parent_grid_name_vec] = Opm::Lgr::excludeFakeSubdivisions(cells_per_dim_vec,
                                                                                       startIJK_vec,
                                                                                       endIJK_vec,
                                                                                       lgr_name_vec,
                                                                                       parent_grid_names);
    if (filtered_cells_per_dim_vec.size() == 0) { // if all LGRs expect 1 child per direction, then no refinement will be done.
        return;
    }

    if (!Opm::areParentGridsAvailableBeforeTheirLgrs(getLgrNameToLevel(),
                                                     filtered_lgr_name_vec,
                                                     filtered_lgr_parent_grid_name_vec)) {
        OPM_THROW(std::invalid_argument, "Parent grid (name) must exist before its LGRs.");
    }

    // Refinement proceeds in steps. In each step, for every parent grid:
    //   1. Gather the data of its child LGRs.
    //   2. Create the corresponding LGRs.
    //   3. Update the leaf grid view.
    //
    // To achieve this, we first collect the set of unique parent grid names
    // (avoiding duplicates).
    std::set<std::string> non_repeated_parent_grid_names(filtered_lgr_parent_grid_name_vec.begin(),
                                                         filtered_lgr_parent_grid_name_vec.end());
    int tmp_maxLevel = this->maxLevel();

    for (const auto& parent_grid_name : non_repeated_parent_grid_names) {
        //   1. Gather the data of its child LGRs.
        auto [cells_per_dim_vec_parent_grid,
              startIJK_vec_parent_grid,
              endIJK_vec_parent_grid,
              lgr_name_vec_parent_grid] =  Opm::filterLgrDataPerParentGridName(filtered_cells_per_dim_vec,
                                                                               filtered_startIJK_vec,
                                                                               filtered_endIJK_vec,
                                                                               filtered_lgr_name_vec,
                                                                               filtered_lgr_parent_grid_name_vec,
                                                                               parent_grid_name);


        int parent_grid_index = getLgrNameToLevel().at(parent_grid_name);
        // Determine the assigned level for the refinement of each marked cell
        std::vector<int> assignRefinedLevel(currentData().back()->size(0));

        // Compatibility of numbers of subdivisions of neighboring LGRs".
        // The method compatibleSubdivision returns a bool. We convert it into an int since MPI within DUNE does not support bool directly.
        int compatibleSubdivisions = Opm::Lgr::compatibleSubdivisions(filtered_cells_per_dim_vec,
                                                                      filtered_startIJK_vec,
                                                                      filtered_endIJK_vec,
                                                                      currentData()[parent_grid_index]->logicalCartesianSize());
        compatibleSubdivisions = comm().min(compatibleSubdivisions); // 0 when at least one process returns false (un-compatible subdivisions).
        if(!compatibleSubdivisions) {
            if (comm().rank()==0){
                OPM_THROW(std::logic_error, "Subdivisions of neighboring LGRs sharing at least one face do not coincide. Not suppported yet.");
            }
            else{
                OPM_THROW_NOLOG(std::logic_error, "Subdivisions of neighboring LGRs sharing at least one face do not coincide. Not suppported yet.");
            }
        }

        // To determine if an LGR is not empty in a given process, for each
        // parent grid, we set at_least_one_active_parent[in that level] to 1
        // if it contains at least one active cell, and to 0 otherwise.
        std::vector<int> at_least_one_active_parent(startIJK_vec_parent_grid.size());

        // Find out which (ACTIVE) elements belong to the block cells defined by startIJK and endIJK values.
        for(const auto& element: elements(this->leafGridView())) {
            std::array<int,3> ijk;
            // Skip the element if its level is not equal to parent_grid_index
            // Note: Inconsistency with DUNE Grid interface. element.level()
            //       returns the index to access the level grid where the entity was born.
            //       This means that, in general,
            //       elment.level() != element.father().level() + 1/
            //       It can happen that | element.level() - element.father().level()| >1.
            if (parent_grid_index != element.level())
                continue;
            currentData()[element.level()]->getIJK(element.getLevelElem().index(), ijk);

            for (std::size_t level = 0; level < startIJK_vec_parent_grid.size(); ++level) {
                bool belongsToLevel = true;
                for (int c = 0; c < 3; ++c) {
                    belongsToLevel = belongsToLevel && ( (ijk[c] >= startIJK_vec_parent_grid[level][c]) && (ijk[c] < endIJK_vec_parent_grid[level][c]) );
                    if (!belongsToLevel)
                        break;
                }
                if(belongsToLevel) {
                    this->mark(1, element);
                    at_least_one_active_parent[level] = 1;
                    assignRefinedLevel[element.index()] = tmp_maxLevel + level +1;
                }
            }
        }
        tmp_maxLevel = tmp_maxLevel + startIJK_vec_parent_grid.size(); // Update the maxLevel

        //   2. Create the corresponding LGRs. and  3. Update the leaf grid view.
        refineAndUpdateGrid(cells_per_dim_vec_parent_grid,
                            assignRefinedLevel,
                            lgr_name_vec_parent_grid,
                            startIJK_vec_parent_grid,
                            endIJK_vec_parent_grid);

        int non_empty_lgrs = 0;
        for (std::size_t level = 0; level < startIJK_vec_parent_grid.size(); ++level){
            // Do not throw if all cells of an LGR are inactive in a parallel run (The process might not 'see' those cells.)
            if (at_least_one_active_parent[level]) {
                ++non_empty_lgrs;
            }
            if ((comm().max(at_least_one_active_parent[level]) == 0) && (comm().rank() == 0)) {
                Opm::OpmLog::warning(lgr_name_vec_parent_grid[level]+ " contains only inactive cells (in rank " + std::to_string(comm().rank()) +").\n");
            }
        }

        // Notice that in a parallel run, non_empty_lgrs represents the local active lgrs, i.e. the lgrs containing active cells which also belong
        // to the current process. Considered per parent grid.
        auto globalActiveLgrs_currentParentGrid = comm().sum(non_empty_lgrs);
        if((globalActiveLgrs_currentParentGrid == 0) && (comm().rank() == 0)) {
            Opm::OpmLog::warning("All the LGRs with parent grid " + parent_grid_name + " contain only inactive cells.\n");
        }
    }
}

void CpGrid::autoRefine(const std::array<int,3>& nxnynz)
{
    // Refinement factors must be odd and positive.
    for (const auto& nd : nxnynz) {
        if (nd<=0 || (nd%2==0)) {
            OPM_THROW(std::invalid_argument, "Refinement factor must be odd and positive.\n");
        }
    }
    const auto endIJK = this->logicalCartesianSize();
    addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {nxnynz},
                          /* startIJK_vec = */ {{0,0,0}},
                          /* endIJK_vec = */ {endIJK},
                          /* lgr_name_vec = */ {"GLOBAL_REFINED"});
}

const std::map<std::string,int>& CpGrid::getLgrNameToLevel() const{
    return lgr_names_;
}

std::array<double,3> CpGrid::getEclCentroid(const int& elemIdx) const
{
    return this-> current_view_data_ -> computeEclCentroid(elemIdx);
}

std::array<double,3> CpGrid::getEclCentroid(const cpgrid::Entity<0>& elem) const
{
    return this-> getEclCentroid(elem.index());
}

void CpGrid::updateCornerHistoryLevels(const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                       const std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                                       const std::unordered_map<int,std::array<int,2>>& adaptedCorner_to_elemLgrAndElemLgrCorner,
                                       const int& corner_count,
                                       const std::vector<std::array<int,2>>& preAdaptGrid_corner_history,
                                       const int& preAdaptMaxLevel,
                                       const int& newLevels)
{
    for (int level = preAdaptMaxLevel+1; level < preAdaptMaxLevel + newLevels+1; ++level) {
        currentData()[level]->corner_history_.resize( currentData()[level] ->size(3), std::array<int,2>({-1,-1}));
    }
    // corner_history_ for levels 0, level 1, ..., preAdapt-maxLevel (maximum level before calling (again) adapt) should be already populated
    // corner_history_[ new corner ] = {-1,-1}
    // corner_history_[ corner equivalent to a corner in a previous level ] = { level where the corner was born, its index in that level grid}.
    for (std::size_t corner = 0; corner < cornerInMarkedElemWithEquivRefinedCorner.size(); ++corner) {
        if (!cornerInMarkedElemWithEquivRefinedCorner[corner].empty()) {
            const auto& [refinedLevel, refinedCorner] = elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner.at(cornerInMarkedElemWithEquivRefinedCorner[corner].back());
            currentData()[refinedLevel]->corner_history_[refinedCorner] = preAdaptGrid_corner_history.empty() ? std::array<int,2>{{0, static_cast<int>(corner)}} :  preAdaptGrid_corner_history[corner];
        }
    }

    // corner_history_ leaf grid view
    for ( int leafCorner = 0; leafCorner < corner_count; ++leafCorner){
        currentData().back()->corner_history_.resize(corner_count);
        const auto& [elemLgr, elemLgrCorner] = adaptedCorner_to_elemLgrAndElemLgrCorner.at(leafCorner);
        if (elemLgr != -1) {
            const auto& [refinedLevel, refinedCorner] = elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner.at({elemLgr, elemLgrCorner});
            currentData().back()->corner_history_[leafCorner] = { refinedLevel, refinedCorner};
        }
        else {
            currentData().back()->corner_history_[leafCorner] =  preAdaptGrid_corner_history.empty() ? std::array<int,2>{{0, elemLgrCorner}} : preAdaptGrid_corner_history[elemLgrCorner];
        }
    }
}

} // namespace Dune
