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
#include <opm/grid/common/ZoltanPartition.hpp>
//#include <opm/grid/common/ZoltanGraphFunctions.hpp>
#include <opm/grid/common/GridPartitioning.hpp>
//#include <opm/grid/common/WellConnections.hpp>
#include <opm/grid/common/CommunicationUtils.hpp>

//#include <fstream>
//#include <iostream>
#include <iomanip>
//#include <tuple>

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
    current_view_data_= data_[0].get();
    global_id_set_ptr_ = std::make_shared<cpgrid::GlobalIdSet>(*current_view_data_);
    
}

std::vector<int>
CpGrid::zoltanPartitionWithoutScatter([[maybe_unused]] const std::vector<cpgrid::OpmWellType>* wells,
                                      [[maybe_unused]] const double* transmissibilities,
                                      [[maybe_unused]] const int numParts,
                                      [[maybe_unused]] const double zoltanImbalanceTol) const
{
#if HAVE_MPI && HAVE_ZOLTAN
    const auto met = EdgeWeightMethod(1);

    return cpgrid::zoltanGraphPartitionGridForJac(*this, wells, transmissibilities,
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

    if(!distributed_data_.empty())
    {
        std::cerr<<"There is already a distributed version of the grid."
                 << " Maybe scatterGrid was called before?"<<std::endl;
        return std::make_pair(false, std::vector<std::pair<std::string,bool> >());
    }

    if (data_.size() > 1)
    {
        if (comm().rank() == 0)
        {
            OPM_THROW(std::logic_error, "Loadbalancing a grid with local grid refinement is not supported, yet.");
        }
        else
        {
            OPM_THROW_NOLOG(std::logic_error, "Loadbalancing a grid with local grid refinement is not supported, yet.");
        }
    }

#if HAVE_MPI
    auto& cc = data_[0]->ccobj_;

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
                    ? cpgrid::zoltanSerialGraphPartitionGridOnRoot(*this, wells, transmissibilities, cc, method, 0, zoltanImbalanceTol, allowDistributedWells, zoltanParams)
                    : cpgrid::zoltanGraphPartitionGridOnRoot(*this, wells, transmissibilities, cc, method, 0, zoltanImbalanceTol, allowDistributedWells, zoltanParams);
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
            ostr << "\nLoad balancing distributes " << data_[0]->size(0)
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


        // distributed_data should be empty at this point.
        distributed_data_.push_back(std::make_shared<cpgrid::CpGridData>(cc, distributed_data_)); 
        distributed_data_[0]->setUniqueBoundaryIds(data_[0]->uniqueBoundaryIds());
       
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

        distributed_data_[0]->distributeGlobalGrid(*this,*this->current_view_data_, computedCellPart);
        // global_id_set_.insertIdSet(*distributed_data_[0]);
        (*global_id_set_ptr_).insertIdSet(*distributed_data_[0]);
        distributed_data_[0]-> index_set_.reset(new cpgrid::IndexSet(distributed_data_[0]->cell_to_face_.size(),
                                                                     distributed_data_[0]-> geomVector<3>().size()));
       


        current_view_data_ = distributed_data_[0].get();
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
    current_view_data_->processEclipseFormat(g,
#if HAVE_ECL_INPUT
                                             nullptr,
#endif
                                             nnc, false, false, false, 0.0);
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

const std::vector<int>& CpGrid::globalCell() const
{
    // Temporary. For a grid with LGRs, we set the globalCell() of the as the one for level 0.
    //            Goal: CartesianIndexMapper well-defined for CpGrid LeafView with LGRs.
    if (current_view_data_ == this-> data_.back().get()){
        return current_view_data_ -> global_cell_;
    }
    else{
        return this -> distributed_data_[0] ->global_cell_;
    }
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
    if (!distributed_data_.empty()){
        return 0;
    }
    if (data_.size() == 1){
        return 0; // "GLOBAL" grid is the unique one
    }
    else {  // There are multiple LGRs
        return double(this -> data_.size() - 2); // last entry is leafView, and it starts in level 0 = GLOBAL grid.
    }
}

template<int codim>
typename CpGridTraits::template Codim<codim>::LevelIterator CpGrid::lbegin (int level) const{
    if (level<0 || level>maxLevel())
        DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
    if (!distributed_data_.empty()){
        return cpgrid::Iterator<codim, All_Partition>(*current_view_data_, 0, true);
    }
    else{
        return cpgrid::Iterator<codim, All_Partition>(*data_[level], 0, true);
    }
}
template typename CpGridTraits::template Codim<0>::LevelIterator CpGrid::lbegin<0>(int) const;
template typename CpGridTraits::template Codim<1>::LevelIterator CpGrid::lbegin<1>(int) const;
template typename CpGridTraits::template Codim<3>::LevelIterator CpGrid::lbegin<3>(int) const;

template<int codim>
typename CpGridTraits::template Codim<codim>::LevelIterator CpGrid::lend (int level) const
{
    if (level<0 || level>maxLevel())
        DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
    if (!distributed_data_.empty()){
        return cpgrid::Iterator<codim, All_Partition>(*current_view_data_, size(codim), true);
    }
    else{
        return cpgrid::Iterator<codim,All_Partition>(*data_[level], size(level, codim), true );
    }
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
    if (!distributed_data_.empty()){
        return cpgrid::Iterator<codim,PiType>(*current_view_data_, 0, true);
    }
    else{
        return cpgrid::Iterator<codim,PiType>(*data_[level], 0, true);
    }
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
    if (!distributed_data_.empty()){
        return cpgrid::Iterator<codim,PiType>(*current_view_data_, size(codim), true);
    }
    else{
        return cpgrid::Iterator<codim,PiType>(*data_[level], size(level, codim), true);
    }

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
    return data_[level]-> size(codim);
}

int CpGrid::size (int codim) const
{
    return current_view_data_->size(codim);
}

int CpGrid::size (int level, GeometryType type) const
{
    if (level<0 || level>maxLevel())
        DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
    return data_[level] -> size(type);
}

int CpGrid::size (GeometryType type) const
{
    return current_view_data_->size(type);
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
    if (current_view_data_ == data_.back().get())
    {
        return *data_[level] -> index_set_;
    }
    else
    {
        return *distributed_data_[level] -> index_set_;
    }
}

const CpGridFamily::Traits::LeafIndexSet& CpGrid::leafIndexSet() const
{
    return *current_view_data_->index_set_;
}

void CpGrid::globalRefine (int)
{
    std::cout << "Warning: Global refinement not implemented, yet." << std::endl;
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

void CpGrid::setZoltanParams(const std::map<std::string,std::string>& params)
{
    zoltanParams = params;
}

const typename CpGridTraits::Communication& Dune::CpGrid::comm () const
{
    return current_view_data_->ccobj_;
}

//

const std::vector<double>& CpGrid::zcornData() const {
    return current_view_data_->zcornData();
}

int CpGrid::numCells() const
{
    return current_view_data_->cell_to_face_.size();
}
/// \brief Get the number of faces.
int CpGrid::numFaces() const
{
    return current_view_data_->face_to_cell_.size();
}
/// \brief Get The number of vertices.
int CpGrid::numVertices() const
{
    return current_view_data_->geomVector<3>().size();
}

int CpGrid::numCellFaces(int cell) const
{
    return current_view_data_->cell_to_face_[cpgrid::EntityRep<0>(cell, true)].size();
}

int CpGrid::cellFace(int cell, int local_index) const
{
    return current_view_data_->cell_to_face_[cpgrid::EntityRep<0>(cell, true)][local_index].index();
}

const cpgrid::OrientedEntityTable<0,1>::row_type CpGrid::cellFaceRow(int cell) const
{
    return current_view_data_->cell_to_face_[cpgrid::EntityRep<0>(cell, true)];
}

int CpGrid::faceCell(int face, int local_index) const
{
    // In the parallel case we store non-existent cells for faces along
    // the front region. Theses marked with index std::numeric_limits<int>::max(),
    // orientation might be arbitrary, though.
    cpgrid::OrientedEntityTable<1,0>::row_type r
        = current_view_data_->face_to_cell_[cpgrid::EntityRep<1>(face, true)];
    bool a = (local_index == 0);
    bool b = r[0].orientation();
    bool use_first = a ? b : !b;
    // The number of valid cells.
    int r_size = r.size();
    // In the case of only one valid cell, this is the index of it.
    int index = 0;
    if(r[0].index()==std::numeric_limits<int>::max()){
        assert(r_size==2);
        --r_size;
        index=1;
    }
    if(r.size()>1 && r[1].index()==std::numeric_limits<int>::max())
    {
        assert(r_size==2);
        --r_size;
    }
    if (r_size == 2) {
        return use_first ? r[0].index() : r[1].index();
    } else {
        return use_first ? r[index].index() : -1;
    }
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

const Dune::FieldVector<double,3> CpGrid::faceCenterEcl(int cell_index, int face) const
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
    static const int faceVxMap[ 6 ][ 4 ] = { {0, 2, 4, 6}, // face 0
                                             {1, 3, 5, 7}, // face 1
                                             {0, 1, 4, 5}, // face 2
                                             {2, 3, 6, 7}, // face 3
                                             {0, 1, 2, 3}, // face 4
                                             {4, 5, 6, 7}  // face 5
    };


    assert (current_view_data_->cell_to_point_[cell_index].size() == 8);
    Dune::FieldVector<double,3> center(0.0);
    for( int i=0; i<4; ++i )
    {
        center += vertexPosition(current_view_data_->cell_to_point_[cell_index][ faceVxMap[ face ][ i ] ]);
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
    current_view_data_=data_[0].get();
}

void CpGrid::switchToDistributedView()
{
    if (distributed_data_.empty())
        OPM_THROW(std::logic_error, "No distributed view available in grid");
    current_view_data_=distributed_data_[0].get();
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

//
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
        data_[0]->writeSintefLegacyFormat(grid_prefix);
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
    current_view_data_->processEclipseFormat(input_data,
#if HAVE_ECL_INPUT
                                             nullptr,
#endif
                                             nnc,
                                             remove_ij_boundary, turn_normals, false, 0.0);
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


void CpGrid::addLgrUpdateLeafView(const std::array<int,3>& cells_per_dim, const std::array<int,3>& startIJK,
                                  const std::array<int,3>& endIJK, const std::string& lgr_name)
{
    this -> addLgrsUpdateLeafView({cells_per_dim}, {startIJK}, {endIJK}, {lgr_name});
}


void CpGrid::addLgrsUpdateLeafView(const std::vector<std::array<int,3>>& cells_per_dim_vec,
                                   const std::vector<std::array<int,3>>& startIJK_vec,
                                   const std::vector<std::array<int,3>>& endIJK_vec,
                                   const std::vector<std::string>& lgr_name_vec)
{
    // Check startIJK_vec and endIJK_vec have same size, and "startIJK[patch][coordinate] < endIJK[patch][coordinate]"
    (*data_[0]).validStartEndIJKs(startIJK_vec, endIJK_vec);
    if (!distributed_data_.empty()){
        if (comm().rank()==0){
            OPM_THROW(std::logic_error, "Adding LGRs to a distributed grid is not supported, yet.");
        }
        else{
            OPM_THROW_NOLOG(std::logic_error, "Adding LGRs to a distributed grid is not supported, yet.");
        }
    }
    // Check LGRs are disjoint
    if (startIJK_vec.size() > 0 && !(*data_[0]).disjointPatches(startIJK_vec, endIJK_vec)){
        if (comm().rank()==0){
            OPM_THROW(std::logic_error, "LGRs are not disjoint.");
        }
        else{
            OPM_THROW_NOLOG(std::logic_error, "LGRs are not disjoint.");
        }
    }
    // Check grid is Cartesian
    const std::array<int,3>& coarseGrid_dim =  (*data_[0]).logical_cartesian_size_;
    long unsigned int coarseGridXYZ = coarseGrid_dim[0]*coarseGrid_dim[1]*coarseGrid_dim[2];
    if ((*data_[0]).global_cell_.size() != coarseGridXYZ){
        if (comm().rank()==0){
            OPM_THROW(std::logic_error, "Grid is not Cartesian. This type of refinement is not supported yet.");
        }
        else {
            // No cells on rank > 0
            return;
        }
    }
    // Check all the cells to be refined have no NNC (no neighbouring connections).
    std::vector<int> all_patch_cells = (*data_[0]).getPatchesCells(startIJK_vec, endIJK_vec);
    if ((*data_[0]).hasNNCs(all_patch_cells)){
        OPM_THROW(std::logic_error, "NNC face on a cell containing LGR is not supported yet.");
    }
    // Total amount of patches:
    const int& num_patches = startIJK_vec.size();
    assert(cells_per_dim_vec.size() == startIJK_vec.size());
    assert(cells_per_dim_vec.size() == endIJK_vec.size());
    assert(cells_per_dim_vec.size() == lgr_name_vec.size());
    // Map to relate boundary patches corners with their equivalent refined/new-born ones. {0,oldCornerIdx} -> {level,newCornerIdx}
    std::map<std::array<int,2>, std::array<int,2>> old_to_new_boundaryPatchCorners;
    // Map to relate boundary patch faces with their children refined/new-born ones. {0,oldFaceIdx} -> {level,{newFaceIdx0, ...}}
    std::map<std::array<int,2>, std::tuple<int, std::vector<int>>> old_to_new_boundaryPatchFaces;
    //
    // For level0, attach children to each parent cell. For no parents, entry {-1, {}} representing {no level, {no children}}
    auto& l0_parent_to_children_cells = (*data_[0]).parent_to_children_cells_;
    l0_parent_to_children_cells.resize(data_[0]-> size(0), std::make_tuple(-1, std::vector<int>()));
    //
    // Get patches corner and face indices. We instantiate them with level1 info and then insert other levels info.
    std::vector<int> all_patch_corners = (*data_[0]).getPatchCorners(startIJK_vec[0], endIJK_vec[0]);
    std::vector<int> all_patch_faces = (*data_[0]).getPatchFaces(startIJK_vec[0], endIJK_vec[0]);
    for (int patch = 0; patch < num_patches; ++patch){
        if (patch+1 < num_patches) { // Populate last pacht information at the end of the for-loop
            const auto& next_patch_corners = (*data_[0]).getPatchCorners(startIJK_vec[patch+1], endIJK_vec[patch+1]);
            const auto& next_patch_faces = (*data_[0]).getPatchFaces(startIJK_vec[patch+1], endIJK_vec[patch+1]);
            all_patch_corners.insert(all_patch_corners.end(), next_patch_corners.begin(), next_patch_corners.end());
            all_patch_faces.insert(all_patch_faces.end(), next_patch_faces.begin(), next_patch_faces.end());
        }
        // Build each LGR from the selected patche of cells from level0 (level0 = this->data_[0]).
        const auto& [level_ptr, boundary_old_to_new_corners, boundary_old_to_new_faces, parent_to_children_faces,
                     parent_to_children_cells, child_to_parent_faces, child_to_parent_cells]
            = (*(this-> data_[0])).refinePatch(cells_per_dim_vec[patch], startIJK_vec[patch], endIJK_vec[patch]);
        //
        // Add each LGR to data_ in entry [patch +1] (shifted +1 since level0 is coarse grid. Levels are 1,2,..., num_patches).
        (this-> data_).push_back(level_ptr);
        //
        // Populate some attributes of the LGR
        //          level_data_ptr_
        (*data_[patch +1]).level_data_ptr_ = &(this -> data_);
        //          level_
        (*data_[patch +1]).level_ = patch +1;
        // Add the name of each LGR in this->lgr_names_
        this -> lgr_names_[lgr_name_vec[patch]] = patch +1; // {"name_lgr", level}
        std::vector<int> l_global_cell(data_[patch+1]->size(0), 0); 
        std::iota(l_global_cell.begin()+1, l_global_cell.end(), 1); // from entry[1], adds +1 per entry: {0,1,2,3,...}
        (*data_[patch+1]).global_cell_ = l_global_cell;
        //          index_set_
        (*data_[patch+1]).index_set_ = std::make_unique<cpgrid::IndexSet>(data_[patch+1]->size(0), data_[patch+1]->size(3));
        //          local_id_set_
        (*data_[patch+1]).local_id_set_ = std::make_shared<const cpgrid::IdSet>(*data_[patch+1]);
        //          cells_per_dim_ Determine the amount of cells per direction, per parent cell, of the corresponding LGR. 
        (*data_[patch +1]).cells_per_dim_ = cells_per_dim_vec[patch];
        //          logical_cartesian_size_ Assuming Cartesian Grid Shape (GLOBAL grid is required to be Cartesian)
        (*data_[patch+1]).logical_cartesian_size_ = {cells_per_dim_vec[patch][0]*(endIJK_vec[patch][0]-startIJK_vec[patch][0]),
            cells_per_dim_vec[patch][1]*(endIJK_vec[patch][1]-startIJK_vec[patch][1]),
            cells_per_dim_vec[patch][2]*(endIJK_vec[patch][2]-startIJK_vec[patch][2])};
        //
        // POPULATING (*data_[0]).parent_to_children_cells_
        // POPULATING (*data_[patch +1]).child_to_parent_cells_
        //    True parent cells have entries: {level of the LGR the parent cell has its children, {child0, child1, ...}}
        //    False parent cells (NO PARENT CELLS) have {-1,{}} entries.
        //
        //    True child cells have entries: {0, parent cell index} (0 represents the "GLOBAL" coarse grid)
        //    False child cells (with no parent) have {-1,-1} entries.
        //    child_to_parent_cells entries look like {child index in the LGR, parent cell index}
        //
        // Re-write entries of actual parent/Create the ones for child cells in each level (not default value {-1,-1} needed for LGRs).
        assert(!parent_to_children_cells.empty());
        (*data_[patch +1]).child_to_parent_cells_.resize(child_to_parent_cells.size());
        (*data_[patch +1]).global_cell_.resize(child_to_parent_cells.size());
        for (const auto& [trueParent, children_list] : parent_to_children_cells){
            l0_parent_to_children_cells[trueParent] = std::make_tuple(patch +1, children_list); // {level/LGR, {child0, child1, ...}}
            assert(!children_list.empty());
            for (const auto& child : children_list){
                (*data_[patch +1]).child_to_parent_cells_[child] = {0, trueParent}; //{level parent-cell, parent-cell-index}
            }
        }
        // Populate old_to_new_boundaryPatchCorners
        for (const auto& [oldCorner, newCorner] : boundary_old_to_new_corners) {
            old_to_new_boundaryPatchCorners[{0, oldCorner}] = {patch +1, newCorner};
            // (shifted) [patch +1] since coarse grid is level 0, levels/LGRs are 1,2, ..., num_patches.
        }
        // Populate old_to_new_boundaryPatchFaces
        for (const auto& [face, children_list] : boundary_old_to_new_faces) {
            old_to_new_boundaryPatchFaces[{0,face}] = {patch+1, children_list};
        }
    } // end-patch-forloop
    // Last patch
    const auto& last_patch_corners = (*data_[0]).getPatchCorners(startIJK_vec[num_patches-1], endIJK_vec[num_patches -1]);
    const auto& last_patch_faces = (*data_[0]).getPatchFaces(startIJK_vec[num_patches -1], endIJK_vec[num_patches -1]);
    const auto& last_patch_cells = (*data_[0]).getPatchCells(startIJK_vec[num_patches -1], endIJK_vec[num_patches -1]);
    all_patch_corners.insert(all_patch_corners.end(), last_patch_corners.begin(), last_patch_corners.end());
    all_patch_faces.insert(all_patch_faces.end(), last_patch_faces.begin(), last_patch_faces.end());
    all_patch_cells.insert(all_patch_cells.end(), last_patch_cells.begin(), last_patch_cells.end()); 
    // Relation between level and leafview cell indices.
    std::vector<int>& l0_to_leaf_cells = (*data_[0]).level_to_leaf_cells_;
    l0_to_leaf_cells.resize(data_[0]->size(0));
    // To store the leaf view (mixed grid: with (non parents) coarse and (children) refined entities).
    typedef Dune::FieldVector<double,3> PointType;
    std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& leaf_data = this -> data_;
#if HAVE_MPI
    auto leaf_view_ptr =
        std::make_shared<Dune::cpgrid::CpGridData>((*(this-> data_[0])).ccobj_, leaf_data);
#else
    // DUNE 2.7 is missing convertion to NO_COMM
    auto leaf_view_ptr = std::make_shared<Dune::cpgrid::CpGridData>(leaf_data);
#endif
    auto& leaf_view = *leaf_view_ptr;
    Dune::cpgrid::DefaultGeometryPolicy& leaf_geometries = leaf_view.geometry_;
    std::vector<std::array<int,8>>& leaf_cell_to_point = leaf_view.cell_to_point_;
    cpgrid::OrientedEntityTable<0,1>& leaf_cell_to_face = leaf_view.cell_to_face_;
    Opm::SparseTable<int>& leaf_face_to_point = leaf_view.face_to_point_;
    cpgrid::OrientedEntityTable<1,0>& leaf_face_to_cell = leaf_view.face_to_cell_;
    cpgrid::EntityVariable<enum face_tag,1>& leaf_face_tags = leaf_view.face_tag_;
    cpgrid::SignedEntityVariable<Dune::FieldVector<double,3>,1>& leaf_face_normals = leaf_view.face_normals_;
    //
    std::vector<std::array<int,2>>& leaf_to_level_cells = leaf_view.leaf_to_level_cells_; // {level, cell index in that level}
    // leaf_child_to_parent_cells[ cell index ] must be {-1,-1} when the cell has no father.
    std::vector<std::array<int,2>>& leaf_child_to_parent_cells = leaf_view.child_to_parent_cells_;
    // All active cells. leaf_global_cell[ leaf cell idx] = origin cell idx in level 0 (parent/equiv cell in level 0)
    std::vector<int>& leaf_global_cell = leaf_view.global_cell_;
    //
    // Mutable containers for leaf view corners, faces, cells, face tags, and face normals.
    Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<0,3>>& leaf_corners =
        *(leaf_geometries.geomVector(std::integral_constant<int,3>()));
    Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<2,3>>& leaf_faces =
        *(leaf_geometries.geomVector(std::integral_constant<int,1>()));
    Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<3,3>>& leaf_cells =
        *(leaf_geometries.geomVector(std::integral_constant<int,0>()));
    Dune::cpgrid::EntityVariableBase<enum face_tag>& mutable_face_tags = leaf_face_tags;
    Dune::cpgrid::EntityVariableBase<PointType>& mutable_face_normals = leaf_face_normals;
    //
    // Integer to count leaf view corners (mixed between corners from level0 not involved in LGRs, and new-born-corners).
    int corner_count = 0;
    // Relation between {level0/level1, old-corner-index/new-born-corner-index} and its corresponding leafview-corner-index.
    std::vector<std::vector<int>> level_to_leaf_corners; // level_to_leaf_corners[level][corner] = leaf_corner_index
    level_to_leaf_corners.resize(num_patches +1); // level 0 + LGRs
    level_to_leaf_corners[0].resize(this -> data_[0]-> size(3),-1); // Entry '-1' for corners not appearing in the leafview
    // Corners coming from the level0, excluding patch_corners, i.e., the old-corners involved in the LGR.
    for (int corner = 0; corner < this-> data_[0]->size(3); ++corner) {
        // Auxiliary bool to discard ANY CORNER COMING FROM ANY PATCH
        bool is_there_allPatchCorn = false;
        for(const auto& patchCorn : all_patch_corners) {
            is_there_allPatchCorn = is_there_allPatchCorn || (corner == patchCorn); //true->corn coincides with one patch corner
            if (is_there_allPatchCorn)
                break;
        }
        if(!is_there_allPatchCorn) { // corner is not involved in any refinement, so we store it.
            level_to_leaf_corners[0][corner] = corner_count; // Only write entries of corners that appear in the LeafView
            corner_count +=1;
        }
    } // end corner-forloop
    // Corners coming from LGRs, i.e. refined (new-born) corners.
    for (int patch = 0; patch < num_patches; ++patch){
        level_to_leaf_corners[patch +1].resize(this -> data_[patch +1] -> size(3));
        for (int corner = 0; corner < this -> data_[patch+1]->size(3); ++corner) {
            level_to_leaf_corners[patch + 1][corner] = corner_count;
            corner_count +=1;
        }
    }
    // Resize the container of the leaf view corners.
    leaf_corners.resize(corner_count);
    for (long unsigned int corner = 0; corner < level_to_leaf_corners[0].size(); ++corner) {
        if (level_to_leaf_corners[0][corner] != -1){ // ONLY NEEDED FOR LEVEL 0
            leaf_corners[level_to_leaf_corners[0][corner]]
                = (*(this->data_[0])).geometry_.geomVector(std::integral_constant<int,3>()) -> get(corner);
        }
    }
    for (int l = 1; l < num_patches +1; ++l) {
        for (long unsigned int corner = 0; corner < level_to_leaf_corners[l].size(); ++corner) {
            const auto& level_data = *(this->data_[l]);
            leaf_corners[level_to_leaf_corners[l][corner]]
                = level_data.geometry_.geomVector(std::integral_constant<int,3>()) -> get(corner);
        }
    }
    // Integer to count leaf view faces (mixed between faces from level0 not involved in LGR, and new-born-faces).
    int face_count = 0;
    // Relation between {level0/level1/..., old-face-index/new-born-face-index} and its corresponding leafview-face-index.
    std::vector<std::vector<int>> level_to_leaf_faces; // level_to_leaf_faces[level][face] = leaf_face_index
    level_to_leaf_faces.resize(num_patches +1); // Level 0 + LGRs
    level_to_leaf_faces[0].resize(this -> data_[0]->face_to_cell_.size(), -1); // Entry -1 for faces that do not appear in leafView
    // Faces coming from the level0, that do not belong to any patch.
    for (int face = 0; face < this->data_[0]->face_to_cell_.size(); ++face) {
        // Auxiliary bool to discard patch faces OF ANY PATCH
        bool is_there_allPatchFace = false;
        for(const auto& patchFace : all_patch_faces) {
            is_there_allPatchFace = is_there_allPatchFace || (face == patchFace); //true->face coincides with one patch face
            if (is_there_allPatchFace)
                break;
        }
        if(!is_there_allPatchFace) { // false-> face was not involved in any LGRs, so we store it.
            level_to_leaf_faces[0][face] = face_count;
            face_count +=1;
        }
    } //end-face-forloop
    // Faces coming from LGRs, i.e. refined faces.
    for (int patch = 0; patch < num_patches; ++patch) {
        level_to_leaf_faces[patch +1].resize(this -> data_[patch +1]->face_to_cell_.size());
        for (int face = 0; face < this->data_[patch +1]-> face_to_cell_.size(); ++face) {
            level_to_leaf_faces[patch +1][face] = face_count;
            face_count +=1;
        }
    }
    // Resize leaf_faces, mutable_face_tags, and mutable_face_normals.
    leaf_faces.resize(face_count);
    mutable_face_tags.resize(face_count);
    mutable_face_normals.resize(face_count);
    // Auxiliary integer to count all the points in leaf_face_to_point.
    int num_points = 0;
    // Auxiliary vector to store face_to_point with non consecutive indices.
    std::vector<std::vector<int>> aux_face_to_point;
    aux_face_to_point.resize(face_count);
    // FACES COMING FROM LEVEL 0
    for (int face = 0; face < static_cast<int>(this -> data_[0]->face_to_cell_.size()); ++face){
        if (level_to_leaf_faces[0][face] != -1){ // ONLY NEEDED FOR LEVEL 0
            const auto& leafFaceIdx =  level_to_leaf_faces[0][face];
            // Get the level data.
            const auto& level_data =  *(this->data_[0]);
            // Get the (face) entity (from level data).
            const auto& entity =  Dune::cpgrid::EntityRep<1>(face, true);
            // Get the face geometry.
            leaf_faces[leafFaceIdx] = (*(level_data.geometry_.geomVector(std::integral_constant<int,1>())))[entity];
            // Get the face tag.
            mutable_face_tags[leafFaceIdx] = level_data.face_tag_[entity];
            // Get the face normal.
            mutable_face_normals[leafFaceIdx] = level_data.face_normals_[entity];
            // Get old_face_to_point.
            auto old_face_to_point = level_data.face_to_point_[face];
            aux_face_to_point[leafFaceIdx].reserve(old_face_to_point.size());
            // Add the amount of points to the count num_points.
            num_points += old_face_to_point.size();
            for (int corn = 0; corn < 4; ++corn) {
                // Auxiliary bool to identify boundary patch corners.
                bool is_there_allPatchBoundCorn = false;
                for(const auto& [l0_oldCorner, level_newCorner] : old_to_new_boundaryPatchCorners){
                    is_there_allPatchBoundCorn = is_there_allPatchBoundCorn || (corn == l0_oldCorner[1]);
                    //true-> boundary patch corner
                    if (is_there_allPatchBoundCorn) {
                        aux_face_to_point[leafFaceIdx].push_back(level_to_leaf_corners[level_newCorner[0]][level_newCorner[1]]);
                        break; // Go to the next corner (on the boundary of a patch)
                    }
                }
                if (!is_there_allPatchBoundCorn) {// Corner not involded in any LGR.
                    aux_face_to_point[leafFaceIdx].push_back(level_to_leaf_corners[0][old_face_to_point[corn]]);
                }
            } // end-corn-forloop
        } // end-if
    } // end-face-for-loop
    // FACES COMING FROM LGRS
    for (int level = 1; level < num_patches +1; ++level) {
        for (int face = 0; face < static_cast<int>(this -> data_[level]->face_to_cell_.size()); ++face){
            const auto& leafFaceIdx =  level_to_leaf_faces[level][face];
            // Get the level data.
            const auto& level_data =  *(this->data_[level]);
            // Get the (face) entity (from level data).
            const auto& entity =  Dune::cpgrid::EntityRep<1>(face, true);
            // Get the face geometry.
            leaf_faces[leafFaceIdx] = (*(level_data.geometry_.geomVector(std::integral_constant<int,1>())))[entity];
            // Get the face tag.
            mutable_face_tags[leafFaceIdx] = level_data.face_tag_[entity];
            // Get the face normal.
            mutable_face_normals[leafFaceIdx] = level_data.face_normals_[entity];
            // Get old_face_to_point.
            auto old_face_to_point = level_data.face_to_point_[face];
            aux_face_to_point[leafFaceIdx].reserve(old_face_to_point.size());
            // Add the amount of points to the count num_points.
            num_points += old_face_to_point.size();
            // Face comes from level1/leavel2/....
            for (int corn = 0; corn < static_cast<int>(old_face_to_point.size()); ++corn) {
                aux_face_to_point[leafFaceIdx].push_back(level_to_leaf_corners[level][old_face_to_point[corn]]);
            }
        } // end-face-for-loop
    } // end-level-forloop
    // Leaf view face_to_point.
    leaf_face_to_point.reserve(face_count, num_points);
    for (int face = 0; face < face_count; ++face) {
        leaf_face_to_point.appendRow(aux_face_to_point[face].begin(), aux_face_to_point[face].end());
    }
    // Integer to count leaf view cells (mixed between cells from level0 not involved in LGR, and new-born-cells).
    int cell_count = 0;
    // Map between leafCellIndices and {level0/level1/.., old-cell-index/new-born-cell-index}.
    leaf_to_level_cells.reserve((this ->data_[0] ->size(0)) + all_patch_cells.size()); // MORE ENTRIES THAT ACTUALLY NEEDED
    // Cells coming from the level0, that do not belong to the patch.
    for (int cell = 0; cell < this->data_[0]-> size(0); ++cell) {
        // Auxiliary bool to identify cells of the patch1.
        bool is_there_allPatchCell = false;
        for(const auto& patch_cell : all_patch_cells) {
            is_there_allPatchCell = is_there_allPatchCell || (cell == patch_cell); //true-> coincides with one patch-cell
            if (is_there_allPatchCell)
                break;
        }
        if(!is_there_allPatchCell) {// Cell does not belong to any patch, so we store it.
            l0_to_leaf_cells[cell] = cell_count;
            leaf_to_level_cells.push_back({0,cell});
            cell_count +=1;
        }
    }
    // Cells coming from LGRs, i.e. refined cells.
    for (int patch = 0; patch < num_patches; ++patch){
        (*data_[patch +1]).level_to_leaf_cells_.resize(this->data_[patch+1]-> size(0));
        // shifted +1 since "GLOBAL" is level0 (coarse grid). Levels are 1,2,..., num_patches.
        for (int cell = 0; cell < this->data_[patch+1]-> size(0); ++cell) {
            (*data_[patch +1]).level_to_leaf_cells_[cell] = cell_count;
            leaf_to_level_cells.push_back({patch +1, cell});
            cell_count +=1;
        }
    }
    leaf_cells.resize(cell_count);
    leaf_cell_to_point.resize(cell_count);
    // For cells that do not have a parent, we set {-1,-1} by defualt and rewrite later for actual children
    leaf_child_to_parent_cells.resize(cell_count, std::array<int,2>({-1,-1}));
    leaf_global_cell.resize(cell_count);
    for (int leafCellIdx = 0; leafCellIdx < cell_count; ++leafCellIdx){
        const auto& level_cellIdx = leaf_to_level_cells[leafCellIdx]; // {level, cellIdx}
        const auto& level_data =  *(this->data_[level_cellIdx[0]]);
        const auto& entity =  Dune::cpgrid::EntityRep<0>(level_cellIdx[1], true);
        // Get the cell geometry.
        leaf_cells[leafCellIdx] = (*(level_data.geometry_.geomVector(std::integral_constant<int,0>())))[entity];
        // Get old corners of the cell that will be replaced with leaf view ones.
        auto old_cell_to_point = level_data.cell_to_point_[level_cellIdx[1]];
        // Get old faces of the cell that will be replaced with leaf view ones.
        auto old_cell_to_face = level_data.cell_to_face_[entity];
        // Auxiliary cell_to_face
        std::vector<cpgrid::EntityRep<1>> aux_cell_to_face;
        if (level_cellIdx[0] == 0) { // Cell comes from level0
            leaf_global_cell[leafCellIdx] = (*data_[0]).global_cell_[level_cellIdx[1]]; // global_cell_[origin cell] in level0
            // Cell to point.
            for (int corn = 0; corn < 8; ++corn) {
                // Auxiliary bool to identity boundary patch corners
                bool is_there_allPatchBoundCorn = false;
                for(const auto& [l0_oldCorner, level_newCorner] : old_to_new_boundaryPatchCorners) {
                    is_there_allPatchBoundCorn = is_there_allPatchBoundCorn || (old_cell_to_point[corn] == l0_oldCorner[1]);
                    if (is_there_allPatchBoundCorn) { //true-> coincides with one boundary patch corner
                        leaf_cell_to_point[leafCellIdx][corn] = level_to_leaf_corners[level_newCorner[0]][level_newCorner[1]];
                        break; // Go to the next corner (of the bondary of a patch)
                    }
                }
                if(!is_there_allPatchBoundCorn) { // Corner does not belong to any patch boundary.
                    leaf_cell_to_point[leafCellIdx][corn] = level_to_leaf_corners[0][old_cell_to_point[corn]];
                }
            } // end-corn-forloop
            // Cell to face.
            for (const auto& face : old_cell_to_face) {   // Auxiliary bool to identity boundary patch faces
                bool is_there_allPatchBoundFace = false;
                for (const auto& [l0_boundFace, level_childrenList] : old_to_new_boundaryPatchFaces) {
                    // l0_boundFace = {0,face},  level_childrenList = {patch +1, children_list}
                    is_there_allPatchBoundFace = is_there_allPatchBoundFace || (face.index() == l0_boundFace[1]);
                    //true-> coincides with one boundary patch face
                    if (is_there_allPatchBoundFace) { // Face belongs to one of the patch boundaries.
                        for (const auto& new_face : std::get<1>(level_childrenList)) {
                            aux_cell_to_face.push_back({new_face, face.orientation()});
                        }
                        is_there_allPatchBoundFace = true;
                        break; // Go to the next corner (on the boundary of a patch)
                    }
                }
                if (!is_there_allPatchBoundFace) { // Face does not belong to any of the patch boundaries.
                    aux_cell_to_face.push_back({level_to_leaf_faces[0][face.index()], face.orientation()});
                }
            } // end-old_cell_to_face-forloop
        }
        else { // Refined cells. (Cell comes from LGRs)
            // Get level where cell was created and its local index, to later deduce its parent.
            auto& [level, levelIdx]  = leaf_to_level_cells[leafCellIdx]; // {level, cell index in that level} (level != 0)
            leaf_child_to_parent_cells[leafCellIdx] = (*data_[level]).child_to_parent_cells_[levelIdx]; //{0, parent cell index}
            leaf_global_cell[leafCellIdx] = (*data_[0]).global_cell_[leaf_child_to_parent_cells[leafCellIdx][1]];
            // global_cell_[parentCell in l0]
            // Cell to point.
            for (int corn = 0; corn < 8; ++corn) {
                leaf_cell_to_point[leafCellIdx][corn] = level_to_leaf_corners[level][old_cell_to_point[corn]];
            }
            // Cell to face.
            for (auto& face : old_cell_to_face) {
                aux_cell_to_face.push_back({level_to_leaf_faces[level][face.index()], face.orientation()});
            }
        }
        // Leaf view cell to face.
        leaf_cell_to_face.appendRow(aux_cell_to_face.begin(), aux_cell_to_face.end());
    }
    // Leaf view face to cell.
    leaf_cell_to_face.makeInverseRelation(leaf_face_to_cell);
    //  Add Leaf View to data_.
    (this-> data_).push_back(leaf_view_ptr);
    current_view_data_ = data_[num_patches +1].get();
    // Leaf  index_set_
    (*data_[num_patches +1]).index_set_ = std::make_unique<cpgrid::IndexSet>(data_[num_patches+1]->size(0), data_[num_patches+1]->size(3));
    // Leaf local_id_set_
    (*data_[num_patches +1]).local_id_set_ = std::make_shared<const cpgrid::IdSet>(*data_[num_patches+1]);
}


const std::map<std::string,int>& CpGrid::getLgrNameToLevel() const{
    return lgr_names_;
}

Dune::cpgrid::CpGridData* CpGrid::getCurrentViewData() const
{
    return  current_view_data_;
}

std::array<double,3> CpGrid::getEclCentroid(const int& elemIdx) const
{
    return this-> current_view_data_ -> computeEclCentroid(elemIdx);
}

std::array<double,3> CpGrid::getEclCentroid(const cpgrid::Entity<0>& elem) const
{
    return this-> getEclCentroid(elem.index());
}

} // namespace Dune
