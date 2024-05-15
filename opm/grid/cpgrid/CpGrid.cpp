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

const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& CpGrid::chooseData() const
{
    if (current_view_data_ == this-> data_.back().get()){
        return data_;
    }
    else{
        return distributed_data_;
    }
}

const std::vector<int>& CpGrid::globalCell() const
{
    // Temporary. For a grid with LGRs, we set the globalCell() of the as the one for level 0.
    //            Goal: CartesianIndexMapper well-defined for CpGrid LeafView with LGRs.
    return chooseData().back() -> global_cell_;
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
    return *chooseData()[level] -> index_set_;
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

Dune::cpgrid::Intersection CpGrid::getParentIntersectionFromLgrBoundaryFace(const Dune::cpgrid::Intersection& intersection) const
{
    if ( intersection.neighbor()) {
        if ((intersection.inside().level() != intersection.outside().level())) {
            // one coarse and one refined neighboring cell
            const auto& cellIn = intersection.inside();
            const auto& cellOut = intersection.outside();

            // Identify the coarse and the refined neighboring cell
            const auto coarseCell =  (cellIn.level() == 0) ? cellIn : cellOut;
            const auto refinedCell =  (coarseCell == cellIn) ? cellOut : cellIn;
            assert(coarseCell.level() == 0);
            assert(refinedCell.level() > 0);

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

bool CpGrid::mark(int refCount, const cpgrid::Entity<0>& element)
{
    // chooseData() is equal to 'data_' when the grid has not been distributed,
    //                          'distributed_data_' otherwise.
    return current_view_data_-> mark(refCount, element);
}

int CpGrid::getMark(const cpgrid::Entity<0>& element) const
{
    return current_view_data_->getMark(element);
}

bool CpGrid::preAdapt()
{
    // Set the flags mighVanish for elements that have been marked for refinement/coarsening.
    return current_view_data_-> preAdapt();
}

bool CpGrid::adapt(const std::array<int,3>& cells_per_dim)
{
    // To do: support coarsening.

    /** Change std::size_t by int since -1 represents current_view_data_ and might lead to narrowing errors */

    
    if (!distributed_data_.empty()){
        if (comm().rank()==0){
            OPM_THROW(std::logic_error, "Adding LGRs to a distributed grid is not supported, yet.");
        }
        else{
            OPM_THROW_NOLOG(std::logic_error, "Adding LGRs to a distributed grid is not supported, yet.");
        }
    }

    // To store/build one unique level with the refinement of all the marked elements.
    std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>> refined_data = this -> data_;
#if HAVE_MPI
    auto refined_grid_ptr = std::make_shared<Dune::cpgrid::CpGridData>((*(this-> data_[0])).ccobj_, refined_data);
#else
    // DUNE 2.7 is missing convertion to NO_COMM
    auto refined_grid_ptr = std::make_shared<Dune::cpgrid::CpGridData>(refined_data);
#endif
    auto& refined_grid = *refined_grid_ptr;
    Dune::cpgrid::DefaultGeometryPolicy&                         refined_geometries = refined_grid.geometry_;
    std::vector<std::array<int,8>>&                              refined_cell_to_point = refined_grid.cell_to_point_;
    cpgrid::OrientedEntityTable<0,1>&                            refined_cell_to_face = refined_grid.cell_to_face_;
    Opm::SparseTable<int>&                                       refined_face_to_point = refined_grid.face_to_point_;
    cpgrid::OrientedEntityTable<1,0>&                            refined_face_to_cell = refined_grid.face_to_cell_;
    cpgrid::EntityVariable<enum face_tag,1>&                     refined_face_tags = refined_grid.face_tag_;
    cpgrid::SignedEntityVariable<Dune::FieldVector<double,3>,1>& refined_face_normals = refined_grid.face_normals_;
    // Mutable containers for refined corners, faces, cells, face tags, and face normals.
    Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<0,3>>& refined_corners =
        *(refined_geometries.geomVector(std::integral_constant<int,3>()));
    Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<2,3>>& refined_faces =
        *(refined_geometries.geomVector(std::integral_constant<int,1>()));
    Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<3,3>>& refined_cells =
        *(refined_geometries.geomVector(std::integral_constant<int,0>()));
    Dune::cpgrid::EntityVariableBase<enum face_tag>& mutable_refined_face_tags = refined_face_tags;
    typedef Dune::FieldVector<double,3> PointType;
    Dune::cpgrid::EntityVariableBase<PointType>& mutable_refined_face_normals = refined_face_normals;

    // To store adapted grid
    std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& adaptedData = this -> data_;
#if HAVE_MPI
    auto adaptedGrid_ptr =
        std::make_shared<Dune::cpgrid::CpGridData>((*(this-> data_[0])).ccobj_, adaptedData);
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

    // -- markedElem_to_itsLgr :
    // Each marked element gets refined and we store this "auxiliary markedElementLGR", to later
    // build a unique level containing all the refined entities from all the marked elements.
    std::vector<std::shared_ptr<Dune::cpgrid::CpGridData> > markedElem_to_itsLgr;
    markedElem_to_itsLgr.resize(current_view_data_->size(0));
    int markedElem_count = 0;

    // -- cornerInMarkedElemWithEquivRefinedCorner :
    // For each corner from level zero, we store the marked elements where the corner appears and its equivalent
    // refined corner in  each auxiliary marked-element-lgr. Example: corner with index 5 appears in marked
    // elements 0 and 1, with refined equivalent corner indices 8 and 2 respectively. Then,
    // cornerInMarkedElemWithEquivRefinedCorner[5] = {{0, 8}, {1, 2}}.
    // For corners not appearing in any marked element, empty vector.
    std::vector<std::vector<std::array<int,2>>> cornerInMarkedElemWithEquivRefinedCorner;
    cornerInMarkedElemWithEquivRefinedCorner.resize(current_view_data_->size(3));
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

    // --- Adapted cells and PreAdapt cells relations ---
    // Integer to count adapted cells (mixed between cells from level0 (not involved in LGRs), and (new-born) refined cells).
    std::size_t cell_count = 0;
    std::map<std::array<int,2>,int>           elemLgrAndElemLgrCell_to_adaptedCell;
    std::unordered_map<int,std::array<int,2>> adaptedCell_to_elemLgrAndElemLgrCell;
    // --- Refined cells and PreAdapt cells relations ---
    // Integer to count only REFINED cells (new-born refined cells from ANY marked element).
    std::size_t refined_cell_count = 0;
    std::vector<std::size_t> refined_cell_count_vec = {refined_cell_count};
    std::map<std::array<int,2>,std::array<int,2>> elemLgrAndElemLgrCell_to_refinedLevelAndRefinedCell;
    std::map<std::array<int,2>,std::array<int,2>> refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell;

    // Auxiliary integer to get the entry where the level with all refined entities will be stored.
    int refinedLevel = this->data_.size();
    // The unique refined grid created with all refined entities will be stored in data_[refinedLevel].
    // The adapted grid containing coarse and refined cells will be stored in data_[refinedLevel +1].
    // Notice that current_view_data_ == data_.back() == data_[refinedLevel-1] before adapt().
    
    // There is only one refined level grid. The entries of non-marked cells will get modified. 
    std::vector<int> assignRefinedLevel(current_view_data_->size(0), refinedLevel);
    // Not necessary
    const std::vector<std::array<int,3>>& cells_per_dim_vec = {cells_per_dim};
    
    // -- Parent-child relations --
    // Relation between the grid before adapt() ("current_view_data_") and refined cells (on the refined grid level - not in each individual lgr).
    std::vector<std::tuple<int,std::vector<int>>>& parent_to_refinedChildCells = current_view_data_ -> parent_to_children_cells_;
    parent_to_refinedChildCells.resize(current_view_data_->size(0), std::make_pair(-1, std::vector<int>{}));
    // Refined child cells and their parents. Entry is {-1,-1} when cell has no father. // {level parent cell, parent cell index}
    std::vector<std::array<int,2>> refinedChild_to_parentCell;
    // Each entry represents a refined level. Here, there is only one.
    std::vector<std::vector<std::array<int,2>>> refinedChild_to_parentCell_vec = {refinedChild_to_parentCell};
    // Adapted child cells and their parents. Entry is {-1,-1} when cell has no father. // {level parent cell, parent cell index}
    std::vector<std::array<int,2>> adaptedChild_to_parentCell;
    std::vector<int> adaptedChild_to_idxInParentCell;
    std::vector<int> refinedChild_to_idxInParentCell;
    // Each entry represents a refined level. Here, there is only one.
    std::vector<std::vector<int>> refinedChild_to_idxInParentCell_vec = {refinedChild_to_idxInParentCell};

    // -- Some extra indices relations between preAdapt-grid, refined-grid, and adapted-grid --
    // Relation between the grid before adapt() ("current_view_data_") and leafview cell indices.
    std::vector<int>& preAdaptCells_to_adaptedCells = current_view_data_ -> level_to_leaf_cells_;
    preAdaptCells_to_adaptedCells.resize(current_view_data_->size(0),-1);
    // Relation between the refined grid ("current_view_data_") and leafview cell indices.
    std::vector<int> refinedCells_to_adaptedCells; // = data_[level] -> level_to_leaf_cells_;
    std::vector<std::vector<int>> refinedCells_to_adaptedCells_vec = {refinedCells_to_adaptedCells};
    // Relation between an adapted cell and its equivalent cell coming either from current_view_data_ or from the refined grid (level)
    std::vector<std::array<int,2>> adaptedCell_to_levelAndLevelCell; // data_.back() -> leaf_to_level_cells_


    // Refine marked cells and provide marked-corner/face/cell - refined-corner/faces/cells relations.
    refineAndProvideMarkedRefinedRelations( markedElem_to_itsLgr,
                                            assignRefinedLevel,
                                            elemLgrAndElemLgrCell_to_adaptedCell,
                                            adaptedCell_to_elemLgrAndElemLgrCell,
                                            cell_count,
                                            elemLgrAndElemLgrCell_to_refinedLevelAndRefinedCell,
                                            refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                            refined_cell_count_vec,
                                            markedElem_count,
                                            cornerInMarkedElemWithEquivRefinedCorner,
                                            markedElemAndEquivRefinedCorn_to_corner,
                                            faceInMarkedElemAndRefinedFaces,
                                            cells_per_dim_vec,
                                            parent_to_refinedChildCells,
                                            preAdaptCells_to_adaptedCells,
                                            refinedLevel); // Notice that current_view_data_ == data_.back() == data_[refinedLevel-1] before adapt().

    // -- Child-parent relations --
    defineChildToParentRelation(adaptedCell_to_elemLgrAndElemLgrCell,
                                cell_count,
                                refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                refined_cell_count_vec,
                                adaptedChild_to_parentCell,
                                refinedChild_to_parentCell_vec,
                                adaptedChild_to_idxInParentCell,
                                refinedChild_to_idxInParentCell_vec,
                                refinedLevel);

    // -- Refined to Adapted cells and Adapted-cells to {level where the cell was born, cell index on that level} --
    defineRefinedAdaptedCellsRelation(  elemLgrAndElemLgrCell_to_adaptedCell,
                                        adaptedCell_to_elemLgrAndElemLgrCell,
                                        cell_count,
                                        elemLgrAndElemLgrCell_to_refinedLevelAndRefinedCell, 
                                        refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                        refined_cell_count_vec,
                                        refinedCells_to_adaptedCells_vec,
                                        adaptedCell_to_levelAndLevelCell,
                                        refinedLevel);

    // --- Adapted corners and PreAdapt corners relations ---
    // Integer to count adapted corners (mixed between corners from current_view_data_ (not involved in LGRs), and (new-born) refined corners).
    std::size_t corner_count = 0;
    std::map<std::array<int,2>,int>           elemLgrAndElemLgrCorner_to_adaptedCorner;
    std::unordered_map<int,std::array<int,2>> adaptedCorner_to_elemLgrAndElemLgrCorner;
    // --- Refined corners and PreAdapt corners relations ---
    // Integer to count only refined corners.
    std::size_t refined_corner_count = 0;
    std::vector<std::size_t> refined_corner_count_vec = {refined_corner_count};
    std::map<std::array<int,2>,std::array<int,2>> elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner;
    std::map<std::array<int,2>,std::array<int,2>> refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner;
    std::map<std::array<int,2>,std::array<int,2>> vanishedRefinedCorner_to_itsLastAppearance;

    definePreAdaptRefinedAndAdaptedCornerRelations(elemLgrAndElemLgrCorner_to_adaptedCorner,
                                                   adaptedCorner_to_elemLgrAndElemLgrCorner,
                                                   corner_count,
                                                   elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                                                   refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner,
                                                   refined_corner_count_vec,
                                                   vanishedRefinedCorner_to_itsLastAppearance,
                                                   markedElem_to_itsLgr,
                                                   assignRefinedLevel,
                                                   cornerInMarkedElemWithEquivRefinedCorner,
                                                   faceInMarkedElemAndRefinedFaces,
                                                   cells_per_dim_vec);

    // --- Adapted corners ---
    populateAdaptedCorners(adapted_corners,
                           corner_count,
                           markedElem_to_itsLgr,
                           adaptedCorner_to_elemLgrAndElemLgrCorner);
    // --- Refined corners  ---
    std::vector<Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<0,3>>> refined_corners_vec = {refined_corners};
    populateRefinedCorners(refined_corners_vec,
                           refined_corner_count_vec,
                           markedElem_to_itsLgr,
                           refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner);

    // --- Adapted faces and PreAdapt faces relations ---
    // Integer to count adapted faces (mixed between faces from current_view_data_ (not involved in LGRs), and (new-born) refined faces).
    std::size_t face_count = 0;
    std::map< std::array<int,2>, int >           elemLgrAndElemLgrFace_to_adaptedFace;
    std::unordered_map< int, std::array<int,2> > adaptedFace_to_elemLgrAndElemLgrFace;
    // --- Refined faces and PreAdapt faces relations ---
    // Integer to count adapted faces (mixed between faces from level0 (not involved in LGRs), and (new-born) refined faces).
    std::size_t refined_face_count = 0;
    std::vector<std::size_t> refined_face_count_vec = {refined_face_count};
    std::map<std::array<int,2>,std::array<int,2>> elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace;
    std::map<std::array<int,2>,std::array<int,2>> refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace;

    definePreAdaptRefinedAndAdaptedFaceRelations( elemLgrAndElemLgrFace_to_adaptedFace,
                                                  adaptedFace_to_elemLgrAndElemLgrFace,
                                                  face_count,
                                                  elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace,
                                                  refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace,
                                                  refined_face_count_vec,
                                                  markedElem_to_itsLgr,
                                                  assignRefinedLevel,
                                                  faceInMarkedElemAndRefinedFaces,
                                                  cells_per_dim_vec);

    // --- Adapted faces ---
    populateAdaptedFaces(adapted_faces,
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
                         cells_per_dim_vec);

    
  
    
    // --- Refined faces  ---
    std::vector<Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<2,3>>> refined_faces_vec = {refined_faces};
    std::vector<Dune::cpgrid::EntityVariableBase<enum face_tag>> mutable_refined_face_tags_vec = {mutable_refined_face_tags}; 
    std::vector<Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>> mutable_refined_face_normals_vec = {mutable_refined_face_normals};
    std::vector<Opm::SparseTable<int>> refined_face_to_point_vec = {refined_face_to_point};
    populateRefinedFaces(refined_faces_vec,
                         mutable_refined_face_tags_vec,
                         mutable_refined_face_normals_vec,
                         refined_face_to_point_vec,
                         refined_face_count_vec,
                         refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace,
                         elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                         vanishedRefinedCorner_to_itsLastAppearance,
                         markedElem_to_itsLgr,
                         cornerInMarkedElemWithEquivRefinedCorner,
                         markedElemAndEquivRefinedCorn_to_corner);

    std::vector<int> adapted_global_cell(cell_count, 0);
    // --- Adapted cells ---
    populateAdaptedCells(adapted_cells,
                         adapted_cell_to_point,
                         adapted_global_cell,
                         cell_count,
                         adapted_cell_to_face,
                         adapted_face_to_cell,
                         adaptedCell_to_elemLgrAndElemLgrCell,
                         elemLgrAndElemLgrFace_to_adaptedFace,
                         faceInMarkedElemAndRefinedFaces,
                         elemLgrAndElemLgrCorner_to_adaptedCorner,
                         vanishedRefinedCorner_to_itsLastAppearance,
                         markedElem_to_itsLgr,
                         assignRefinedLevel,
                         markedElemAndEquivRefinedCorn_to_corner,
                         cornerInMarkedElemWithEquivRefinedCorner,
                         cells_per_dim_vec);

   
    // --- Refined cells  ---
    std::vector<Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<3,3>>> refined_cells_vec = {refined_cells};
    std::vector<std::vector<std::array<int,8>>> refined_cell_to_point_vec = {refined_cell_to_point}; 
    std::vector<int> refined_global_cell(refined_cell_count, 0);
    std::vector<std::vector<int>> refined_global_cell_vec = {refined_global_cell};
    std::vector<cpgrid::OrientedEntityTable<0,1>> refined_cell_to_face_vec = { refined_cell_to_face};
    std::vector<cpgrid::OrientedEntityTable<1,0>> refined_face_to_cell_vec = { refined_face_to_cell};
    populateRefinedCells(refined_cells_vec,
                         refined_cell_to_point_vec,
                         refined_global_cell_vec,
                         refined_cell_count_vec,
                         refined_cell_to_face_vec,
                         refined_face_to_cell_vec,
                         refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                         elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace,
                         faceInMarkedElemAndRefinedFaces,
                         elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                         vanishedRefinedCorner_to_itsLastAppearance,
                         markedElem_to_itsLgr,
                         assignRefinedLevel,
                         markedElemAndEquivRefinedCorn_to_corner,
                         cornerInMarkedElemWithEquivRefinedCorner,
                         cells_per_dim_vec);


    (this-> data_).push_back(refined_grid_ptr);
    // Populate some attributes of the level LGR
    (*data_[refinedLevel]).level_data_ptr_ = &(this -> data_);
    (*data_[refinedLevel]).level_ = refinedLevel;
    this -> lgr_names_["LGR"] = refinedLevel; // {"name_lgr", level}
    (*data_[refinedLevel]).child_to_parent_cells_ = refinedChild_to_parentCell_vec[0];
    (*data_[refinedLevel]).cell_to_idxInParentCell_ = refinedChild_to_idxInParentCell_vec[0];
    (*data_[refinedLevel]).level_to_leaf_cells_ =  refinedCells_to_adaptedCells_vec[0];
    (*data_[refinedLevel]).global_cell_ = refined_global_cell_vec[0];
    (*data_[refinedLevel]).index_set_ = std::make_unique<cpgrid::IndexSet>(data_[refinedLevel]->size(0), data_[refinedLevel]->size(3));
    (*data_[refinedLevel]).local_id_set_ = std::make_shared<const cpgrid::IdSet>(*data_[refinedLevel]);
    // Determine the amount of cells per direction, per parent cell, of the corresponding LGR.
    (*data_[refinedLevel]).cells_per_dim_ = cells_per_dim;
    // TO DO: This new code for refinement do not assume Cartesian Shape. How does logical_cartesian_size_ should be defined then?
     (*data_[refinedLevel]).logical_cartesian_size_ =  (*data_[0]).logical_cartesian_size_;

    (this-> data_).push_back(adaptedGrid_ptr);
    current_view_data_ = data_.back().get();
    (*data_[refinedLevel+1]).child_to_parent_cells_ = adaptedChild_to_parentCell;
    (*data_[refinedLevel+1]).cell_to_idxInParentCell_ = adaptedChild_to_idxInParentCell;
    (*data_[refinedLevel+1]).leaf_to_level_cells_ =  adaptedCell_to_levelAndLevelCell;
    (*data_[refinedLevel+1]).global_cell_ = adapted_global_cell;
    (*data_[refinedLevel+1]).index_set_ = std::make_unique<cpgrid::IndexSet>(data_[refinedLevel+1]->size(0), data_[refinedLevel+1]->size(3));
    (*data_[refinedLevel+1]).local_id_set_ = std::make_shared<const cpgrid::IdSet>(*data_[refinedLevel+1]);
     // Determine the amount of cells per direction, per parent cell, of the corresponding LGR.
    (*data_[refinedLevel+1]).cells_per_dim_ = cells_per_dim;
    // TO DO: How to modified logical_cartesian_size_
    (*data_[refinedLevel +1]).logical_cartesian_size_ =  (*data_[0]).logical_cartesian_size_;

    // Print total amount of cells on the leaf grid view
    Opm::OpmLog::info(std::to_string(markedElem_count) + " marked elements have been refined.\n");
    Opm::OpmLog::info(std::to_string(cell_count)  + " total cells on the adapted grid.\n");

    return preAdapt();
}

void CpGrid::postAdapt()
{
    // - Resize with the new amount of cells on the leaf grid view
    // - Set marks equal to zero (representing 'doing nothing')
    current_view_data_ -> postAdapt();
}

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
    /* Begin new approach 
    // Check all the cells to be refined have no NNC (no neighbouring connections).
    std::vector<int> markedCells = (*data_[0]).getPatchesCells(startIJK_vec, endIJK_vec);
    for (const auto& elemIdx : markedCells) {
        const auto& elem =  Dune::cpgrid::Entity<0>(*(coarse_grid.chooseData()[0]), elemIdx, true);
        coarse_grid.mark(1, elem);
    }


    End new approach */
    
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
    // Check LGRs are disjoint (sharing corners allowed, sharing faces not allowed)
    if (startIJK_vec.size() > 0 && (*data_[0]).patchesShareFace(startIJK_vec, endIJK_vec)){
        if (comm().rank()==0){
            OPM_THROW(std::logic_error, "LGRs share at least one face.");
        }
        else{
            OPM_THROW_NOLOG(std::logic_error, "LGRs share at least one face.");
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
    //
    // Total amount of patches:
    const int& num_patches = startIJK_vec.size();
    assert(cells_per_dim_vec.size() == startIJK_vec.size());
    assert(cells_per_dim_vec.size() == endIJK_vec.size());
    assert(cells_per_dim_vec.size() == lgr_name_vec.size());
    // Leaf grid view corners are unique. They are stored with consecutives starting from zero.
    // To update the leaf grid view corners after refinement, we store:
    // 1) the level zero corners that do not belong to any LGR.
    // 2) the refined corners from each LGR.
    //
    // Attention: LGRs are allowed to share corners (NOT faces!), that means that at least one corner from
    //            level zero can appear on more than one LGR boundary.
    //            With this in mind, we need to keep track of the relationships:
    //            {level zero, corner index in level zero } <---> { level l, (refined) equivalent corner index in level l}
    //            [for corners laying on boundary of LGRs that coincide with corners on level zero].
    //
    //            A corner from level zero is equivalent to a refined corner when their geometry().center() coincide.
    //            (they have the same values in the x-,y-,and z-coordinates).
    //
    // We introduce a few containers to store level zero corners and their equivalent corners in different LGRs.
    //
    // Map to relate level zero corners laying on LGRs boundaries. Each entry looks like
    // {0, corner index in level zero} -> {level, equivalent (refined) corner index in that level}
    //    Notice that when the same corner from level zero appears in more than one LGR, the entry for
    //    that corner, {0, corner index in level zero}, will be overwritten. In the end, this level zero
    //    corner will be associated with the last LGR where it apperas, together with the corresponding
    //    equivalent refined corner index in that last level.
    //    Example: if the grid has 3 LGRs in total, and the level zero corner with index 67 appears
    //    in LGR1 and LGR3 boundaries, with lgr1_corner_index 34 and lgr3_corner_index 23, then
    //    {0, 67} -> {3, 23}.
    //    This map will be key to count the correct amount of leaf grid view corners since it does not allow
    //    any repetition. Recall that the corners on the leaf grid view are unique and store with consecutive
    //    indices starting from 0.
    // Attention: these containers involve only refined corners on the boundary of LGRs that coincide
    //            with a corner in level zero. Namely, there might be refined corners that coincide
    //            and lay on boundary of LGRs, but do not have an equivalent corner in level zero.
    //            These are not taken into account (yet). A situation where this can happen is when
    //            two LGRs share an edge, have similar widths, lengths, and heights, and the same
    //            amount of children in the axis-direction that involve the shared edge.
    std::map<std::array<int,2>, std::array<int,2>> levelZeroToLGRsBoundaryCorners_oneToOne;
    // We need also the inverse relationship, namely, for each refined equivalent corner (stored
    // ONLY ONCE), we associte the corresponding level zero corner, that is:
    // {last level where that corner appears, its index in that level} -> {0, equivalent corner in level zero}.
    // Example: if the grid has 3 LGRs in total, and the level zero corner with index 67 appears
    //    in LGR1 and LGR3 boundaries, with lgr1_corner_index 34 and lgr3_corner_index 23, then
    //    {3, 23}-> {0, 67}.
    std::map<std::array<int,2>, std::array<int,2>> lgrsToLevelZeroBoundaryCorners_oneToOne;
    // To be able to track the origin of a refined corner laying on the boundary an LGR,
    // we store for each of those corners, the corresponding equivalent corner on level zero.
    // For the same value {0, corner index in level zero}, there might be more than one key
    // {level, refined equivalent corner index in that level}.
    // Example: if the grid has 3 LGRs in total, and the level zero corner with index 67 appears
    // in LGR1 and LGR3 boundaries, with lgr1_corner_index 34, and lgr3_corner_index 23, then:
    // {1, 34} -> {0, 67}
    // {3, 23} -> {0, 67}.
    std::map<std::array<int,2>, std::array<int,2>> lgrsToLevelZeroBoundaryCorners_oneToMoreThanOne;
    // The last container of this sequence to distinguish and relate corners laying on boundary of LGRs,
    // separate such information per level. Namely, for each level, we store the indeces of the refined
    // corners laying on the boundary of the corresponding LGR level, and associate it with the equivalent
    // corner in level zero.
    // Example: if the grid has 3 LGRs in total, and the level zero corner with index 67 appears
    // in LGR1 and LGR3 boundaries, with lgr1_corner_index 34, and lgr3_corner_index 23, then:
    // lgrBoundaryCornerWithEquivalentInLevelZero[0] contains 34
    // lgrBoundaryCornerWithEquivalentInLevelZero[2] contains 23.
    // Attetnion: LGR1 information is stored in entry 0, LGR2, in entry 1, LGR3 in 2, and so on.
    std::vector<std::vector<int>> lgrBoundaryCornerWithEquivalentInLevelZero;
    lgrBoundaryCornerWithEquivalentInLevelZero.resize(num_patches);

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
        /** Logical Cartesian size is equal to cells_per_dim_ in the adapt-method */
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
        // Populate levelZeroToLGRsBoundaryCorners_oneToOne
        // Populate lgrsToLevelZeroBoundaryCorners_oneToMoreThanOne
        // Populate lgrBoundaryCornerWithEquivalentInLevelZero
        // For a description of these containers, see above.
        for (const auto& [cornerIdxLevelZero, cornerIdxLgr] : boundary_old_to_new_corners) {
            levelZeroToLGRsBoundaryCorners_oneToOne[{0, cornerIdxLevelZero}] = {patch +1, cornerIdxLgr};
            // (shifted) [patch +1] since coarse grid is level 0, levels/LGRs are 1,2, ..., num_patches.
            lgrsToLevelZeroBoundaryCorners_oneToMoreThanOne[{patch +1, cornerIdxLgr}] = {0, cornerIdxLevelZero};
            lgrBoundaryCornerWithEquivalentInLevelZero[patch].push_back(cornerIdxLgr);
        }
        // Populate old_to_new_boundaryPatchFaces
        for (const auto& [face, children_list] : boundary_old_to_new_faces) {
            old_to_new_boundaryPatchFaces[{0,face}] = {patch+1, children_list};
        }
    } // end-patch-forloop

    // Populate lgrsToLevelZeroBoundaryCorners_oneToOne
    // For a description of these containers, see above.
    for (const auto [zeroLevelInfo, newLevelInfo] : levelZeroToLGRsBoundaryCorners_oneToOne) {
        // zeroLevelInfo = { 0, corner index in level zero }
        // newLevelInfo = { last level where the equivalent refined corner appears, its index in that level}
        lgrsToLevelZeroBoundaryCorners_oneToOne[newLevelInfo] = zeroLevelInfo;
    }

    // Last patch cornes and faces.
    const auto& last_patch_corners = (*data_[0]).getPatchCorners(startIJK_vec[num_patches-1], endIJK_vec[num_patches -1]);
    const auto& last_patch_faces = (*data_[0]).getPatchFaces(startIJK_vec[num_patches -1], endIJK_vec[num_patches -1]);
    all_patch_corners.insert(all_patch_corners.end(), last_patch_corners.begin(), last_patch_corners.end());
    all_patch_faces.insert(all_patch_faces.end(), last_patch_faces.begin(), last_patch_faces.end());
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
            if (is_there_allPatchCorn) {
                break;
            }
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
            // Auxiliary bool to detect if a refined corner (born in an LGR) lays on the boundary
            // of the LGR and coincides with a corner from level zero.
            bool isEquivCornerLevelZero_laysOnLgrBoundary = false;
            for (const auto& boundaryCorner : lgrBoundaryCornerWithEquivalentInLevelZero[patch]) {
                isEquivCornerLevelZero_laysOnLgrBoundary = isEquivCornerLevelZero_laysOnLgrBoundary || (corner == boundaryCorner);
            }
            // We store the refined corner
            // 1. lays on the LGR boundary, but does not coincide with any corner from level zero, or
            // 2. lays on the LGR boundary, it does coincide with a corner from level zero, BUT its
            //    equivalent corner in level zero appeared in a later LGR, or
            // 3. is in the interior of the LGR.
            //
            // Notice that condition 2. can be check via the lgrsToLevelZeroBoundaryConerns_oneToOne:
            // for each refined equivalent corner (stored ONLY ONCE), we associte the corresponding level zero corner:
            // {last level where that corner appears, its index in that level} -> {0, equivalent corner in level zero}.
            // This means that
            // - lgrsToLevelZeroBoundaryCorners_oneToOne.count({patch+1, corner}) takes the value 0
            // when the refined corner appears in a latter level ( and lays on the lgr boundary, being equivalent to a
            // corner in level zero).
            // - lgrToLevelZeroBoundaryCorners_oneToOne.count({patch+1, corner}) takes the value 1
            // when "patch+1" is the last level where the equivalent corner from level zero appears on an LGR boundary.
            // In this way, we avoid storing repeated corners on the leaf grid view, except from the ones mentioned
            // above levelZeroToLGRsBoundaryCorners_oneToOne declaration. See "Attention" comment.
            if (!isEquivCornerLevelZero_laysOnLgrBoundary || (lgrsToLevelZeroBoundaryCorners_oneToOne.count({patch+1, corner}) == 1) ){
                level_to_leaf_corners[patch + 1][corner] = corner_count;
                corner_count +=1;
            }
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
                if (level_to_leaf_corners[0][old_face_to_point[corn]] == -1) {
                    // In this case, the corner from level zero got replaced by a refined one.
                    // Detect the corresponding LGR (if the corner appears in more than one LGR, then it's the last one)
                    // and the refined corner index in that LGR, which is equivalen (meaning that both - corner
                    // from level zero and refined corner - have exactly the same values in x-,y-, and z-coordinates.
                    const auto [lgr, lgrCornIdx] = levelZeroToLGRsBoundaryCorners_oneToOne[{0, old_face_to_point[corn]}];
                    aux_face_to_point[leafFaceIdx].push_back(level_to_leaf_corners[lgr][lgrCornIdx]);
                }
                else{// Corner not involded in any LGR.
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
                // If the corner is on the LGR boundary and coincides with a corner from level zero,
                // the corresponding equivalent leaf grid view corner could be associated to another level.
                // Namely, the last LGR where this special corner appeared.
                //
                // We can distinguish the following three cases:
                //
                // 1. lays on the LGR boundary, it does coincide with a corner from level zero, BUT its
                //    equivalent corner in level zero appeared in a later LGR, or
                // 2. lays on the LGR boundary, but does not coincide with any corner from level zero, or
                // 3. is in the interior of the LGR.
                //
                // Notice that part of conditions 1. and 2. can be check via the lgrsToLevelZeroBoundaryConerns_oneToOne:
                // for each refined equivalent corner (stored ONLY ONCE), we associte the corresponding level zero corner:
                // {last level where that corner appears, its index in that level} -> {0, equivalent corner in level zero}.
                // This means that
                // - lgrsToLevelZeroBoundaryCorners_oneToOne.count({level, old_face_to_point[corn]}) takes the value 0
                // when the refined corner appears in a latter level ( and lays on the lgr boundary, being equivalent to a
                // corner in level zero).
                // - lgrToLevelZeroBoundaryCorners_oneToOne.count({level, old_face_to_point[corn]}) takes the value 1
                // when "level" is the last level where the equivalent corner from level zero appears on an LGR boundary.
                // In this way, we avoid storing repeated corners on the leaf grid view, except from the ones mentioned
                // above levelZeroToLGRsBoundaryCorners_oneToOne declaration. See "Attention" comment.
                //
                // Auxiliary bool to detect if a refined corner (born in an LGR) lays on the boundary
                // of the LGR and coincides with a corner from level zero.
                bool isEquivCornerLevelZero_laysOnLgrBoundary = false;
                for (const auto& boundaryCorner : lgrBoundaryCornerWithEquivalentInLevelZero[level-1]) {
                    isEquivCornerLevelZero_laysOnLgrBoundary = isEquivCornerLevelZero_laysOnLgrBoundary
                        || (old_face_to_point[corn] == boundaryCorner);
                    // Case 1. described above.
                    if (isEquivCornerLevelZero_laysOnLgrBoundary &&
                        (lgrsToLevelZeroBoundaryCorners_oneToOne.count({level, old_face_to_point[corn]}) == 0)) {
                        // --- Find corner in the later level ---
                        // Get the equivalent corner in level zero (its index in level zero)
                        const auto& levelZeroInfo = lgrsToLevelZeroBoundaryCorners_oneToMoreThanOne[{level, old_face_to_point[corn]}];
                        // Use the info from the equivalent corner in level zero to find the (last) level where
                        // the corner appeared for last time, and its index in that level.
                        const auto& [updateLevel, newCornerIdx ] = levelZeroToLGRsBoundaryCorners_oneToOne[levelZeroInfo];
                        // Get the leaf grid view corner, with the "updatedLevel" and the corresponding corner index in that level.
                        aux_face_to_point[leafFaceIdx].push_back(level_to_leaf_corners[updateLevel][newCornerIdx]);
                        break;
                    }
                }
                // Cases 2. and 3. described above.
                if (!isEquivCornerLevelZero_laysOnLgrBoundary
                    || (lgrsToLevelZeroBoundaryCorners_oneToOne.count({level, old_face_to_point[corn]}) == 1) ){
                    aux_face_to_point[leafFaceIdx].push_back(level_to_leaf_corners[level][old_face_to_point[corn]]);
                }
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
                for(const auto& [l0_oldCorner, level_newCorner] : levelZeroToLGRsBoundaryCorners_oneToOne) {
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
                        const auto levelTouched = std::get<0>(level_childrenList);
                        for (const auto& new_face : std::get<1>(level_childrenList)) {
                            //  level_to_leaf_faces[level][face] = leaf_face_index
                            aux_cell_to_face.push_back({level_to_leaf_faces[levelTouched][new_face], face.orientation()});
                        }
                        is_there_allPatchBoundFace = true;
                        break; // Go to the next face (on the boundary of a patch)
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
                // If the corner is on the LGR boundary and coincides with a corner from level zero,
                // the corresponding equivalent leaf grid view corner could be associated to another level.
                // Namely, the last LGR where this special corner appeared.
                //
                // We can distinguish the following three cases:
                //
                // 1. lays on the LGR boundary, it does coincide with a corner from level zero, BUT its
                //    equivalent corner in level zero appeared in a later LGR, or
                // 2. lays on the LGR boundary, but does not coincide with any corner from level zero, or
                // 3. is in the interior of the LGR.
                //
                // Notice that part of conditions 1. and 2. can be check via the lgrsToLevelZeroBoundaryConerns_oneToOne:
                // for each refined equivalent corner (stored ONLY ONCE), we associte the corresponding level zero corner:
                // {last level where that corner appears, its index in that level} -> {0, equivalent corner in level zero}.
                // This means that
                // - lgrsToLevelZeroBoundaryCorners_oneToOne.count({level, old_face_to_point[corn]}) takes the value 0
                // when the refined corner appears in a latter level ( and lays on the lgr boundary, being equivalent to a
                // corner in level zero).
                // - lgrToLevelZeroBoundaryCorners_oneToOne.count({level, old_face_to_point[corn]}) takes the value 1
                // when "level" is the last level where the equivalent corner from level zero appears on an LGR boundary.
                // In this way, we avoid storing repeated corners on the leaf grid view, except from the ones mentioned
                // above levelZeroToLGRsBoundaryCorners_oneToOne declaration. See "Attention" comment.
                //
                // Auxiliary bool to detect if a refined corner (born in an LGR) lays on the boundary
                // of the LGR and coincides with a corner from level zero.
                bool isEquivCornerLevelZero_laysOnLgrBoundary = false;
                for (const auto& boundaryCorner : lgrBoundaryCornerWithEquivalentInLevelZero[level-1]) {
                    isEquivCornerLevelZero_laysOnLgrBoundary = isEquivCornerLevelZero_laysOnLgrBoundary
                        || (old_cell_to_point[corn] == boundaryCorner);
                    // Case 1. described above.
                    if (isEquivCornerLevelZero_laysOnLgrBoundary &&
                        (lgrsToLevelZeroBoundaryCorners_oneToOne.count({level, old_cell_to_point[corn]}) == 0)) {
                        // --- Find corner in the later level ---
                        // Get the equivalent corner in level zero (its index in level zero)
                        const auto& levelZeroInfo = lgrsToLevelZeroBoundaryCorners_oneToMoreThanOne[{level, old_cell_to_point[corn]}];
                        // Use the info from the equivalent corner in level zero to find the (last) level where
                        // the corner appeared for last time, and its index in that level.
                        const auto& [updateLevel, newCornerIdx ] = levelZeroToLGRsBoundaryCorners_oneToOne[levelZeroInfo];
                        // Get the leaf grid view corner, with the "updatedLevel" and the corresponding corner index in that level.
                        leaf_cell_to_point[leafCellIdx][corn] = level_to_leaf_corners[updateLevel][newCornerIdx];
                        break;
                    }
                }
                // Case 2. and 3. described above.
                if (!isEquivCornerLevelZero_laysOnLgrBoundary ||
                    (lgrsToLevelZeroBoundaryCorners_oneToOne.count({level, old_cell_to_point[corn]}) == 1) ){
                    leaf_cell_to_point[leafCellIdx][corn] = level_to_leaf_corners[level][old_cell_to_point[corn]];
                }
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
    // Leaf logical_cartesian_size_
    (*data_[num_patches +1]).logical_cartesian_size_ =  (*data_[0]).logical_cartesian_size_;

    // Print total amount of cells on the leaf grid view
    Opm::OpmLog::info(std::to_string(num_patches) + " LGRs applied to global grid.\n");
    Opm::OpmLog::info(std::to_string(cell_count)  + " total cells on the refined grid.\n");
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

void CpGrid::refineAndProvideMarkedRefinedRelations(  std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                                      std::vector<int>& assignRefinedLevel,
                                                      std::map<std::array<int,2>,int>& elemLgrAndElemLgrCell_to_adaptedCell,
                                                      std::unordered_map<int,std::array<int,2>>& adaptedCell_to_elemLgrAndElemLgrCell,
                                                      std::size_t& cell_count,
                                                      std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrCell_to_refinedLevelAndRefinedCell,
                                                      std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                                      std::vector<std::size_t>& refined_cell_count_vec,
                                                      int& markedElem_count,
                                                      std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                                      std::map<std::array<int,2>,int>& markedElemAndEquivRefinedCorn_to_corner,
                                                      std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                                      const std::vector<std::array<int,3>>& cells_per_dim_vec,
                                                      std::vector<std::tuple<int,std::vector<int>>>& parent_to_refinedChildCells,
                                                      std::vector<int>& preAdaptCells_to_adaptedCells,
                                                      const int refinedLevel)
{   // Each marked element for refinement (mark equal to 1), will be refined individuality, creating its own Lgr. The element index will
    // be also used to identify its lgr. Even though, in the end, all the refined entities will belong to a unique level grid.
    // For this reason, we associate "-1" with those elements that are not involved in any refinement and will appear
    // as "coarse" cells in the leaf-grid-view (adapted-grid).
    for (int elemIdx = 0; elemIdx <  current_view_data_ -> size(0); ++elemIdx) {
        const auto element = Dune::cpgrid::Entity<0>(*current_view_data_, elemIdx, true);
        const auto elemMark = getMark(element);
        // When the element is marked with 0 ("doing nothing"), it will appear in the adapted grid with same geometrical features (center, volume).
        if (elemMark ==  0) {
            elemLgrAndElemLgrCell_to_adaptedCell[{-1, elemIdx}] = cell_count;
            adaptedCell_to_elemLgrAndElemLgrCell[cell_count] = {-1, elemIdx};
            cell_count +=1;
            preAdaptCells_to_adaptedCells[elemIdx] = cell_count;
        }
        // When the element is marked for refinement, we also mark its corners and faces
        // since they will get replaced by refined ones.
        if (elemMark ==  1) {
            markedElem_count +=1;
            const auto& markedElemLevel = assignRefinedLevel[elemIdx];
            assert(markedElemLevel >= refinedLevel);
            // Shift the markedElemRefinedLevel to access data containers
            const auto& shiftedLevel = markedElemLevel - refinedLevel;
            // Build auxiliary LGR for the refinement of this element
            const auto& [elemLgr_ptr,
                         parentCorners_to_equivalentRefinedCorners,
                         parentFace_to_itsRefinedFaces,
                         parentCell_to_itsRefinedCells,
                         refinedFace_to_itsParentFace,
                         refinedCell_to_itsParentCell]
                = current_view_data_-> refineSingleCell(cells_per_dim_vec[shiftedLevel], elemIdx);
            markedElem_to_itsLgr[ elemIdx ] = elemLgr_ptr;

            const auto& childrenCount = cells_per_dim_vec[shiftedLevel][0]*cells_per_dim_vec[shiftedLevel][1]*cells_per_dim_vec[shiftedLevel][2];
            std::vector<int> refinedChildrenList(childrenCount);

            for (int refinedCell = 0; refinedCell < childrenCount; ++refinedCell) {

                elemLgrAndElemLgrCell_to_adaptedCell[{elemIdx, refinedCell}] = cell_count;
                adaptedCell_to_elemLgrAndElemLgrCell[cell_count] = {elemIdx, refinedCell};
                cell_count +=1;

                elemLgrAndElemLgrCell_to_refinedLevelAndRefinedCell[{elemIdx, refinedCell}] = { markedElemLevel, refined_cell_count_vec[shiftedLevel]};
                refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell[{markedElemLevel, refined_cell_count_vec[shiftedLevel]}] = {elemIdx, refinedCell};
                refinedChildrenList[refinedCell] = refined_cell_count_vec[shiftedLevel];
                refined_cell_count_vec[shiftedLevel] +=1;

            }
            parent_to_refinedChildCells[elemIdx] = std::make_pair( markedElemLevel, refinedChildrenList);
            for (const auto& [markedCorner, lgrEquivCorner] : parentCorners_to_equivalentRefinedCorners) {
                cornerInMarkedElemWithEquivRefinedCorner[markedCorner].push_back({elemIdx, lgrEquivCorner});
                markedElemAndEquivRefinedCorn_to_corner[ {elemIdx, lgrEquivCorner}] = markedCorner;
            }
            for (const auto& [markedFace, itsRefinedFaces] : parentFace_to_itsRefinedFaces) {
                faceInMarkedElemAndRefinedFaces[markedFace].push_back({elemIdx, itsRefinedFaces});
            }
        } // end-if-elemMark==1
    } // end-elem-for-loop
}

void CpGrid::defineChildToParentRelation(std::unordered_map<int,std::array<int,2>> adaptedCell_to_elemLgrAndElemLgrCell,
                                         const std::size_t cell_count,
                                         std::map<std::array<int,2>,std::array<int,2>> refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                         const std::vector<std::size_t> refined_cell_count_vec,
                                         std::vector<std::array<int,2>>& adaptedChild_to_parentCell, // {level parent cell, parent cell index}
                                         std::vector<std::vector<std::array<int,2>>>& refinedChild_to_parentCell_vec,
                                         std::vector<int>& adaptedChild_to_idxInParentCell, // {level parent cell, parent cell index}
                                         std::vector<std::vector<int>>& refinedChild_to_idxInParentCell_vec,
                                         const int refinedLevel)
{
    adaptedChild_to_parentCell.resize(cell_count, std::array<int,2>{-1,-1});
    adaptedChild_to_idxInParentCell.resize(cell_count, -1);
    // Rewrite only the entries of adapted cells that have a parent cell
    for (std::size_t cell = 0; cell < cell_count; ++cell) {
        const auto& elemLgr = adaptedCell_to_elemLgrAndElemLgrCell[cell][0];
        // -1 represents current_view_data_ which has level "refinedLevel-1"
        if (elemLgr != -1) {
            adaptedChild_to_parentCell[cell] = {refinedLevel-1, elemLgr};
            adaptedChild_to_idxInParentCell[cell] = adaptedCell_to_elemLgrAndElemLgrCell[cell][1];
        }
    }
    
    for (std::size_t shiftedLevel = 0; shiftedLevel < refined_cell_count_vec.size(); ++shiftedLevel) {
        refinedChild_to_parentCell_vec[shiftedLevel].resize(refined_cell_count_vec[shiftedLevel]);
        refinedChild_to_idxInParentCell_vec[shiftedLevel].resize(refined_cell_count_vec[shiftedLevel]);
        const auto& level =  shiftedLevel + refinedLevel; /** refinedLevel should be called firstRefinedLevel */
        // Every refined cell has a parent cell
        for (std::size_t cell = 0; cell < refined_cell_count_vec[shiftedLevel]; ++cell) {
            const auto& elemLgr = refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell[{level, cell}][0];
            // -1 represents current_view_data_ which has level "refinedLevel-1"
            assert(elemLgr != -1);
            refinedChild_to_parentCell_vec[shiftedLevel][cell] = {refinedLevel-1, elemLgr};
            refinedChild_to_idxInParentCell_vec[shiftedLevel][cell] =  refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell[{level,cell}][1];
        }
    }
}

void CpGrid::defineRefinedAdaptedCellsRelation(  std::map<std::array<int,2>,int> elemLgrAndElemLgrCell_to_adaptedCell,
                                                 std::unordered_map<int,std::array<int,2>> adaptedCell_to_elemLgrAndElemLgrCell,
                                                 const std::size_t cell_count,
                                                 std::map<std::array<int,2>,std::array<int,2>> elemLgrAndElemLgrCell_to_refinedLevelAndRefinedCell,
                                                 std::map<std::array<int,2>,std::array<int,2>> refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                                 const std::vector<std::size_t> refined_cell_count_vec,
                                                 std::vector<std::vector<int>>& refinedCells_to_adaptedCells_vec,
                                                 std::vector<std::array<int,2>>& adaptedCell_to_levelAndLevelCell,
                                                 const int refinedLevel)
{
    assert(refinedLevel>0);
    // -- Adapted to {level, cell index in that level}  --
    adaptedCell_to_levelAndLevelCell.resize(cell_count);
    for (std::size_t cell = 0; cell < cell_count; ++cell) {
        const auto& elemLgr = adaptedCell_to_elemLgrAndElemLgrCell[cell][0];
        const auto& elemLgrCell = adaptedCell_to_elemLgrAndElemLgrCell[cell][1];
        // -1 represents current_view_data_ which has level "refinedLevel-1"
        const auto& levelWhereCellWasBorn = (elemLgr == -1) ? (refinedLevel-1) : elemLgrAndElemLgrCell_to_refinedLevelAndRefinedCell[{elemLgr, elemLgrCell}][0];
        const auto& levelCellIdx =  (elemLgr == -1) ? elemLgrCell : elemLgrAndElemLgrCell_to_refinedLevelAndRefinedCell[{elemLgr, elemLgrCell}][1];
        adaptedCell_to_levelAndLevelCell[cell] = { levelWhereCellWasBorn, levelCellIdx};
    }
    // -- Refined to adapted cells --
    for (std::size_t shiftedLevel = 0; shiftedLevel < refined_cell_count_vec.size(); ++shiftedLevel){
        refinedCells_to_adaptedCells_vec[shiftedLevel].resize(refined_cell_count_vec[shiftedLevel]);
        for (std::size_t cell = 0; cell < refined_cell_count_vec[shiftedLevel]; ++cell) {
            const auto& elemLgr = refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell[{shiftedLevel + refinedLevel, cell}][0];
            const auto& elemLgrCell = refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell[{shiftedLevel + refinedLevel, cell}][1];
            const auto& adaptedCellIdx = elemLgrAndElemLgrCell_to_adaptedCell[{elemLgr, elemLgrCell}];
            refinedCells_to_adaptedCells_vec[shiftedLevel][cell] = adaptedCellIdx;
        }
    }
}

void CpGrid::definePreAdaptRefinedAndAdaptedCornerRelations(std::map<std::array<int,2>,int>& elemLgrAndElemLgrCorner_to_adaptedCorner,
                                                            std::unordered_map<int,std::array<int,2>>& adaptedCorner_to_elemLgrAndElemLgrCorner,
                                                            std::size_t& corner_count,
                                                            std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                                                            std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner,
                                                            std::vector<std::size_t>& refined_corner_count_vec,
                                                            std::map<std::array<int,2>, std::array<int,2>>& vanishedRefinedCorner_to_itsLastAppearance,
                                                            const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                                            const std::vector<int>& assignRefinedLevel,
                                                            const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                                            const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                                            const std::vector<std::array<int,3>>& cells_per_dim_vec)
{
    // Step 1. Select/store the corners from level 0 not involved in any LGR.
    //         Replace the corners from level zero involved in LGR by the equivalent ones, born in LGRs.
    //         In this case, we avoid repetition considering the last appearance of the level zero corner
    //         in the LGRs.
    for (int corner = 0; corner < current_view_data_->size(3); ++corner) {
        if (cornerInMarkedElemWithEquivRefinedCorner[corner].empty()) { // save it
            // Note: Since we are associating each LGR with its parent cell index, and this index can take
            //       the value 0, we will represent the grid before adapt is called with the value -1.
            elemLgrAndElemLgrCorner_to_adaptedCorner[{-1, corner}] = corner_count;
            adaptedCorner_to_elemLgrAndElemLgrCorner[corner_count] = {-1, corner};
            corner_count +=1;
        }
        else { // corner involved in refinement, so we search it in one LGR (the last one where it appears)
            assert(!cornerInMarkedElemWithEquivRefinedCorner[corner].empty());
            // Get the lgr corner that replaces the marked corner from level zero.
            // Note: Recall that lgr coincides with the marked element index from level 0 that got refined.
            //       Since the container is a map, the lgr and the lgr corner index correspond to the last
            //       appearance of the marked corner (from level 0).
            const auto& lastAppearanceLgr = cornerInMarkedElemWithEquivRefinedCorner[corner].back()[0];
            const auto& lastAppearanceLgrCorner = cornerInMarkedElemWithEquivRefinedCorner[corner].back()[1];

            // Build the relationships between adapted corner and level corner, for future search due topology aspects.
            elemLgrAndElemLgrCorner_to_adaptedCorner[{lastAppearanceLgr, lastAppearanceLgrCorner}] = corner_count;
            adaptedCorner_to_elemLgrAndElemLgrCorner[corner_count] = {lastAppearanceLgr, lastAppearanceLgrCorner};
            corner_count +=1;

            const auto& level = assignRefinedLevel[lastAppearanceLgr];
            assert(level>= this->data_.size());
            const auto& shiftedLevel = level - (this->data_.size());

            elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner[{lastAppearanceLgr, lastAppearanceLgrCorner}] = {level, refined_corner_count_vec[shiftedLevel]};
            refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner[{level, refined_corner_count_vec[shiftedLevel]}] = {lastAppearanceLgr, lastAppearanceLgrCorner};
            refined_corner_count_vec[shiftedLevel] +=1;
        }
    } // end corner-forloop

    for (int elemIdx = 0; elemIdx < current_view_data_->size(0); ++elemIdx) {
        if (markedElem_to_itsLgr[elemIdx]!= nullptr) {
            const auto& level = assignRefinedLevel[elemIdx];
            assert(level>= this->data_.size());
            const auto& shiftedLevel = level - (this->data_.size());
            for (int corner = 0; corner < markedElem_to_itsLgr[elemIdx] ->size(3); ++corner) {
                // Discard marked corners. Store (new born) refined corners

                // INTERIOR
                if (isRefinedCornerInInteriorLgr(cells_per_dim_vec[shiftedLevel], corner)) { // It's a refined interior corner, so we store it.
                    // In this case, the corner is a new born refined corner that does not
                    // coincide with any corner from the GLOBAL grid (level 0). Therefore,
                    // it has to be stored.
                    elemLgrAndElemLgrCorner_to_adaptedCorner[{elemIdx, corner}] = corner_count;
                    adaptedCorner_to_elemLgrAndElemLgrCorner[corner_count] = {elemIdx, corner};
                    corner_count += 1;

                    elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner[{elemIdx, corner}] = {level, refined_corner_count_vec[shiftedLevel]};
                    refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner[{level, refined_corner_count_vec[shiftedLevel]}] = {elemIdx, corner};
                    refined_corner_count_vec[shiftedLevel] +=1;
                }

                // LAYING ON EDGES
                //
                // Refined corners laying on edges - Refined edge has a 'coarse' parent edge (line between 2 corners of the parent cell)
                // To avoid repetition, we distinguish the case where the refined corner lays on an edge of its parent cell.
                // We detect the two coarse faces involved (Notice that the extremes of the parent cell have been stored previously).
                // When the marked faces appeares only once, we store the corner now. Otherwise, we store the refined corner on its
                // last apparance assocaited to one of these parent faces, taking also into account the elemLgr. For example, when
                // the refined corners lays on an edge connecting I_FACE false and K_FACE true of the parent cell, let's say iFaceIdx,
                // kFaceIdx, with each of those faces appearing twice (maximum) :
                // iFaceIdx appearing in current "elem" and elemLgr1
                // kFaceIdx appearing in current "elem" and elemLgr2
                // Then, we take the max(elemLgr1, elemLgr2) and store the refined corner only if this maximum equals elem.
                if (newRefinedCornerLaysOnEdge(cells_per_dim_vec[shiftedLevel], corner)) {
                    const auto& markedFacesTouchingEdge = getParentFacesAssocWithNewRefinedCornLayingOnEdge(cells_per_dim_vec[shiftedLevel], corner, elemIdx);
                    const auto& markedFace1 = markedFacesTouchingEdge[0];
                    const auto& markedFace2 = markedFacesTouchingEdge[1];

                    int lastAppearanceMarkedFace1 = faceInMarkedElemAndRefinedFaces[markedFace1].back().first; // elemLgr1
                    int lastAppearanceMarkedFace2 = faceInMarkedElemAndRefinedFaces[markedFace2].back().first; // elemLgr2

                    int maxLastAppearance = std::max(lastAppearanceMarkedFace1, lastAppearanceMarkedFace2);
                    int faceAtMaxLastAppearance = (maxLastAppearance == lastAppearanceMarkedFace1) ? markedFace1 : markedFace2;

                    // Save the relationship between the vanished refined corner and its last appearance
                    bool atLeastOneFaceAppearsTwice = (faceInMarkedElemAndRefinedFaces[markedFace1].size()>1) ||
                        (faceInMarkedElemAndRefinedFaces[markedFace2].size()>1);
                    if (atLeastOneFaceAppearsTwice && (maxLastAppearance != elemIdx)) {
                        // modify the method below to take into account face!
                        /** To be done: if the neighboring cells sharing the face do not belong to the same refined level, then
                            their cells_per_dim_[ levelA ]  and cells_per_dim_[ levelB ] might differ. Modify the code to take this into account.
                            Namely, when they differ, store both. */
                        const auto& neighboringLgrCornerIdx = replaceLgr1CornerIdxByLgr2CornerIdx(cells_per_dim_vec[shiftedLevel], corner, elemIdx, faceAtMaxLastAppearance);
                        vanishedRefinedCorner_to_itsLastAppearance[{elemIdx, corner}] = {maxLastAppearance, neighboringLgrCornerIdx};
                        // Notice that, when we use these container to locate vanished corners, we might need a while-loop,
                        // since {elem, corner} leads to {lastMaxAppearance, neighboringLgrCornerIdx}, which can also vanish.
                        // So we need something like:
                        // if (elemLgrAndElemLgrCorner_to_adapted/refinedCorner.count({elem, corner}) == 0)
                        //    int updateElemLgr =  vanishedRefinedCorner_to_itsLastAppearance[{elem, corner}][0];
                        //    int updateElemLgrCorner =  vanishedRefinedCorner_to_itsLastAppearance[{elem, corner}][1];
                        //     while (elemLgrAndElemLgrCorner_to_adapted/refinedCorner.count({updateElemLgr, updateElemLgCorner}) == 0)
                        //        int tempElemLgr =  updateElemLgr;
                        //        int tempElemLgrCorner =  updateElemLgrCorner;
                        //        updateElemLgr =  vanishedRefinedCorner_to_itsLastAppearance[{ tempElemLgr ,  tempElemLgrCorner}][0];
                        //        updateElemLgrCorner =  vanishedRefinedCorner_to_itsLastAppearance[{ tempElemLgr ,  tempElemLgrCorner}][1];
                        // Then, use the lastest update to search for the corner in teh refined/adapted grid (which would be the one that
                        // gives elemLgrAndElemLgrCorner_to_adapted/refinedCorner.count({updateElemLgr, updateElemLgCorner}) == 1).
                    }
                    if (maxLastAppearance == elemIdx) {
                        // Store the refined corner in its last appearence - to avoid repetition.
                        elemLgrAndElemLgrCorner_to_adaptedCorner[{elemIdx, corner}] = corner_count;
                        adaptedCorner_to_elemLgrAndElemLgrCorner[corner_count] = {elemIdx, corner};
                        corner_count += 1;

                        elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner[{elemIdx, corner}] = {level, refined_corner_count_vec[shiftedLevel]};
                        refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner[{level,refined_corner_count_vec[shiftedLevel]}] = {elemIdx, corner};
                        refined_corner_count_vec[shiftedLevel] +=1;
                    }
                }

                // LAYING ON BOUNDARY LGR - NOT ON AN EDGE - NOT COINCIDING WITH A MARKED CORNER
                //
                // If the refined corner lays on the boundary of the LGR, e.i., it was born on one of the faces
                // of the marked element that got refined, then, we have two cases:
                // - the marked face appears only in one marked element -> then, we store this corner now.
                // - the marked face appears twice (maximum times) in two marked elements -> we store it later.
                if ( isRefinedNewBornCornerOnLgrBoundary(cells_per_dim_vec[shiftedLevel], corner) && !newRefinedCornerLaysOnEdge(cells_per_dim_vec[shiftedLevel], corner)) {
                    // Get the index of the marked face where the refined corner was born.
                    const auto& markedFace = getParentFaceWhereNewRefinedCornerLaysOn(cells_per_dim_vec[shiftedLevel], corner, elemIdx);
                    // check how many times marked face appearn
                    // Get the last LGR (marked element) where the marked face appeared.
                    int lastLgrWhereMarkedFaceAppeared = faceInMarkedElemAndRefinedFaces[markedFace].back().first;
                    // Save the relationship between the vanished refined corner and its last appearance
                    if ((faceInMarkedElemAndRefinedFaces[markedFace].size()>1) && (lastLgrWhereMarkedFaceAppeared != elemIdx)) {
                        /** To be done: if the neighboring cells sharing the face do not belong to the same refined level, then
                            their cells_per_dim_[ levelA ]  and cells_per_dim_[ levelB ] might differ. Modify the code to take this into account.
                            Namely, when they differ, store both. */
                        const auto& neighboringLgrCornerIdx = replaceLgr1CornerIdxByLgr2CornerIdx(cells_per_dim_vec[shiftedLevel], corner);
                        vanishedRefinedCorner_to_itsLastAppearance[{elemIdx, corner}] = {lastLgrWhereMarkedFaceAppeared, neighboringLgrCornerIdx};
                    }
                    if (lastLgrWhereMarkedFaceAppeared == elemIdx) {
                        // Store the refined corner in its last appearence - to avoid repetition.
                        elemLgrAndElemLgrCorner_to_adaptedCorner[{elemIdx, corner}] = corner_count;
                        adaptedCorner_to_elemLgrAndElemLgrCorner[corner_count] = {elemIdx, corner};
                        corner_count += 1;

                        elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner[{elemIdx, corner}] = {level, refined_corner_count_vec[shiftedLevel]};
                        refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner[{level, refined_corner_count_vec[shiftedLevel]}] = {elemIdx, corner};
                        refined_corner_count_vec[shiftedLevel] +=1;
                    }
                }
            } // end-corner-for-loop
        } // end-if-nullptr
    } // end-elem-for-loop
}

void  CpGrid::definePreAdaptRefinedAndAdaptedFaceRelations(  std::map<std::array<int,2>,int>& elemLgrAndElemLgrFace_to_adaptedFace,
                                                             std::unordered_map<int,std::array<int,2>>& adaptedFace_to_elemLgrAndElemLgrFace,
                                                             std::size_t& face_count,
                                                             std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace,
                                                             std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace,
                                                             std::vector<std::size_t>& refined_face_count_vec,
                                                             const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                                             const std::vector<int>& assignRefinedLevel,
                                                             const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                                             const std::vector<std::array<int,3>>& cells_per_dim_vec)
{
    // Step 1. Add the LGR faces, for each LGR
    for (int elem = 0; elem < current_view_data_->size(0); ++elem) {
        if (markedElem_to_itsLgr[elem]!=nullptr)  {
            const auto& level = assignRefinedLevel[elem];
            assert(level>= this->data_.size());
            const auto& shiftedLevel = level - (this->data_.size());
            for (int face = 0; face < markedElem_to_itsLgr[elem] ->face_to_cell_.size(); ++face) {
                // Discard marked faces. Store (new born) refined faces
                bool isNewRefinedFaceOnLgrBoundary = isRefinedFaceOnLgrBoundary(cells_per_dim_vec[shiftedLevel], face, markedElem_to_itsLgr[elem]);
                if (!isNewRefinedFaceOnLgrBoundary) { // It's a refined interior face, so we store it
                    //  if (isRefinedFaceInInteriorLgr(cells_per_dim, face, markedElem_to_itsLgr[elem])) {
                    // In this case, the face is a new born refined face that does not
                    // have any "parent face" from the GLOBAL grid (level 0).
                    elemLgrAndElemLgrFace_to_adaptedFace[{elem, face}] = face_count;
                    adaptedFace_to_elemLgrAndElemLgrFace[face_count] = {elem, face};
                    face_count += 1;

                    elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace[{elem, face}] = {level, refined_face_count_vec[shiftedLevel]};
                    refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace[{level, refined_face_count_vec[shiftedLevel]}] = {elem, face};
                    refined_face_count_vec[shiftedLevel] +=1;
                }
                // If the refined face lays on the boundary of the LGR, e.i., it was born on one of the faces
                // of the marked element that got refined, then, we have two cases:
                // - the marked face appears only in one marked element -> then, we store this face now.
                // - the marked face appears twice (maximum times) in two marked elements -> we store it later.
                else {
                    // Get the index of the marked face where the refined corner was born.
                    int markedFace = getParentFaceWhereNewRefinedFaceLaysOn(cells_per_dim_vec[shiftedLevel], face, markedElem_to_itsLgr[elem], elem);
                    assert(!faceInMarkedElemAndRefinedFaces[markedFace].empty());
                    // Get the last LGR (marked element) where the marked face appeared.
                    int lastLgrWhereMarkedFaceAppeared = faceInMarkedElemAndRefinedFaces[markedFace].back().first;
                    if (lastLgrWhereMarkedFaceAppeared == elem) {
                        // Store the refined face in its last appearence - to avoid repetition.
                        elemLgrAndElemLgrFace_to_adaptedFace[{elem, face}] = face_count;
                        adaptedFace_to_elemLgrAndElemLgrFace[face_count] = {elem, face};
                        face_count += 1;

                        elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace[{elem, face}] = {level, refined_face_count_vec[shiftedLevel]};
                        refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace[{level, refined_face_count_vec[shiftedLevel]}] = {elem, face};
                        refined_face_count_vec[shiftedLevel] +=1;
                    }
                }
            } // end-face-for-loop
        } // end-if-nullptr
    } // end-elem-for-loop
    // Step 2. Select/store the faces from level 0 not involved in any LGR.
    //         Replace the faces from level zero involved in LGR by the equivalent ones, born in LGRs.
    //         In this case, we avoid repetition considering the last appearance of the level zero corner
    //         in the LGRs.
    for (int face = 0; face < current_view_data_->face_to_cell_.size(); ++face) {
        if (faceInMarkedElemAndRefinedFaces[face].empty()) { // save it
            // Note: Since we are associating each LGR with its parent cell index, and this index can take
            //       the value 0, we will represent the current_view_data_ with the value -1
            elemLgrAndElemLgrFace_to_adaptedFace[{-1, face}] = face_count;
            adaptedFace_to_elemLgrAndElemLgrFace[face_count] = {-1, face};
            face_count +=1;

        }
    } // end face-forloop
}

void CpGrid::populateAdaptedCorners(Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<0,3>>& adapted_corners,
                                    const std::size_t& corner_count,
                                    const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                    std::unordered_map<int,std::array<int,2>> adaptedCorner_to_elemLgrAndElemLgrCorner)
{
    adapted_corners.resize(corner_count);
    for (std::size_t corner = 0; corner < corner_count; ++corner) {
        const auto& elemLgr = adaptedCorner_to_elemLgrAndElemLgrCorner[corner][0];
        const auto& elemLgrCorner = adaptedCorner_to_elemLgrAndElemLgrCorner[corner][1];
        // Note: Since we are associating each LGR with its parent cell index, and this index can take
        //       the value 0, we will represent the level 0 (current_view_data_) with the value -1
        if (elemLgr == -1){ // Corner from current_view_data_
            adapted_corners[corner] = current_view_data_->geometry_.geomVector(std::integral_constant<int,3>()) -> get(elemLgrCorner);
        }
        else {
            const auto& elemLgrData = markedElem_to_itsLgr[elemLgr];
            adapted_corners[corner] = elemLgrData -> geometry_.geomVector(std::integral_constant<int,3>()) -> get(elemLgrCorner);
        }
    }
}

void CpGrid::populateRefinedCorners(std::vector<Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<0,3>>>& refined_corners_vec,
                                    const std::vector<std::size_t>& refined_corner_count_vec,
                                    const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                    std::map<std::array<int,2>,std::array<int,2>> refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner)
{
    for (std::size_t shiftedLevel = 0; shiftedLevel < refined_corner_count_vec.size(); ++shiftedLevel) {
        refined_corners_vec[shiftedLevel].resize(refined_corner_count_vec[shiftedLevel]);
        const auto& level = shiftedLevel + (this->data_.size());
        for (std::size_t corner = 0; corner < refined_corner_count_vec[shiftedLevel]; ++corner) {
            const auto& elemLgr = refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner[{level,corner}][0];
            const auto& elemLgrCorner = refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner[{level,corner}][1];
            const auto& elemLgrData = markedElem_to_itsLgr[elemLgr];
            refined_corners_vec[shiftedLevel][corner] = elemLgrData -> geometry_.geomVector(std::integral_constant<int,3>()) -> get(elemLgrCorner);
        }
    }
}

void CpGrid::populateAdaptedFaces(Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<2,3>>& adapted_faces,
                                  Dune::cpgrid::EntityVariableBase<enum face_tag>& mutable_face_tags,
                                  Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>& mutable_face_normals,
                                  Opm::SparseTable<int>& adapted_face_to_point,
                                  const std::size_t& face_count,
                                  std::unordered_map<int,std::array<int,2>> adaptedFace_to_elemLgrAndElemLgrFace,
                                  std::map<std::array<int,2>,int> elemLgrAndElemLgrCorner_to_adaptedCorner,
                                  std::map<std::array<int,2>, std::array<int,2>> vanishedRefinedCorner_to_itsLastAppearance,
                                  const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                  const std::vector<int>& assignRefinedLevel,
                                  std::map<std::array<int,2>,int> markedElemAndEquivRefinedCorn_to_corner,
                                  const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                  const std::vector<std::array<int,3>>& cells_per_dim_vec)
{
    adapted_faces.resize(face_count);
    mutable_face_tags.resize(face_count);
    mutable_face_normals.resize(face_count);
    //
    // Auxiliary integer to count all the points in leaf_face_to_point.
    std::size_t num_points = 0;
    // Auxiliary vector to store face_to_point with non consecutive indices.
    std::vector<std::vector<int>> aux_face_to_point;
    aux_face_to_point.resize(face_count);
    for (std::size_t face = 0; face < face_count; ++face) {

        const auto& elemLgr = adaptedFace_to_elemLgrAndElemLgrFace[face][0];
        // Note: Since we are associating each LGR with its parent cell index, and this index can take
        //       the value 0, we will represent the level 0 (GLOBAL grid) with the value -1.
        const auto& elemLgrFace = adaptedFace_to_elemLgrAndElemLgrFace[face][1];
        const auto& elemLgrFaceEntity =  Dune::cpgrid::EntityRep<1>(elemLgrFace, true);
        Opm::SparseTable<int>::mutable_row_type preAdapt_face_to_point;

        if (elemLgr == -1) { // The value -1 represents the current_grid_view_
            // Get the face geometry.
            adapted_faces[face] = (*(current_view_data_->geometry_.geomVector(std::integral_constant<int,1>())))[elemLgrFaceEntity];
            // Get the face tag.
            mutable_face_tags[face] = current_view_data_->face_tag_[elemLgrFaceEntity];
            // Get the face normal.
            mutable_face_normals[face] = current_view_data_->face_normals_[elemLgrFaceEntity];
            // Get face_to_point_ before adapting - we need to replace the level corners by the adapted ones.
            preAdapt_face_to_point = current_view_data_->face_to_point_[elemLgrFace];
            // Add the amount of points to the count num_points.
            num_points += preAdapt_face_to_point.size();
        }
        else {
            const auto& elemLgrData = markedElem_to_itsLgr[elemLgr]; // Recall elemLgr == marked elemIdx
            // Get the face geometry.
            adapted_faces[face] = (*(elemLgrData->geometry_.geomVector(std::integral_constant<int,1>())))[elemLgrFaceEntity];
            // Get the face tag.
            mutable_face_tags[face] = elemLgrData->face_tag_[elemLgrFaceEntity];
            // Get the face normal.
            mutable_face_normals[face] = elemLgrData->face_normals_[elemLgrFaceEntity];
            // Get face_to_point_ before adapting - we need to replace the level corners by the adapted ones.
            preAdapt_face_to_point = elemLgrData->face_to_point_[elemLgrFace];
            // Add the amount of points to the count num_points.
            num_points += preAdapt_face_to_point.size();
        }
        for (int corn = 0; corn < static_cast<int>(preAdapt_face_to_point.size()); ++corn) {
            int adaptedCorn;
            const auto& elemLgrCorn = preAdapt_face_to_point[corn];
            // Corner is stored in adapted_corners
            if ( elemLgrAndElemLgrCorner_to_adaptedCorner.count({elemLgr, elemLgrCorn}) == 1)  {
                adaptedCorn =  elemLgrAndElemLgrCorner_to_adaptedCorner[{elemLgr, elemLgrCorn}];
            }
            else{
                // Corner might have vanished - Search its equivalent lgr-corner in that case -
                // last lgr where the corner appears -
                int lastAppearanceLgr = 0; // It'll get rewritten
                int lastAppearanceLgrEquivCorner = 0; // It'll get rewritten
                if ((elemLgr ==-1) && (!cornerInMarkedElemWithEquivRefinedCorner[elemLgrCorn].empty())) {
                    lastAppearanceLgr = cornerInMarkedElemWithEquivRefinedCorner[elemLgrCorn].back()[0];
                    lastAppearanceLgrEquivCorner = cornerInMarkedElemWithEquivRefinedCorner[elemLgrCorn].back()[1];
                }
                if (elemLgr>-1) { // It represents the refinement of the element with index "elemIdx == elemLgr"
                    if (markedElemAndEquivRefinedCorn_to_corner.count({elemLgr, elemLgrCorn}) == 1) {
                        const auto& markedCorner = markedElemAndEquivRefinedCorn_to_corner[{elemLgr, elemLgrCorn}];
                        lastAppearanceLgr = cornerInMarkedElemWithEquivRefinedCorner[markedCorner].back()[0];
                        lastAppearanceLgrEquivCorner = cornerInMarkedElemWithEquivRefinedCorner[markedCorner].back()[1];
                    }
                    else {
                        const auto& shiftedLevel = assignRefinedLevel[elemLgr] - (this->data_.size());
                        bool isNewRefinedCornInInteriorLgr = isRefinedCornerInInteriorLgr(cells_per_dim_vec[shiftedLevel], elemLgrCorn);
                        assert(!isNewRefinedCornInInteriorLgr); // and  elemLgrAndElemLgrCorner_to_adaptedCorner.count({elemLgr, elemLgrCorn}) == 0.,
                        // To locate vanished corners, we need a while-loop, since {elemLgr, elemLgrcorner} leads to
                        // {neighboringElemLgr, neighboringElemLgrCornerIdx}, which might have also vanished.
                        // Then, use the lastest appearance of the current corner, meaning, the first (and unique one - by construction) that
                        // gives elemLgrAndElemLgrCorner_to_adaptedCorner.count({lastAppearanceElemLgr, lastAppearanceElemLgCorner}) == 1).
                        lastAppearanceLgr = vanishedRefinedCorner_to_itsLastAppearance[{elemLgr, elemLgrCorn}][0];
                        // this corner lays on the area occupied by a coarse face that got refined and belonged to two marked elements.
                        // Get the index of this corner with respect to the greatest marked element index.
                        lastAppearanceLgrEquivCorner = vanishedRefinedCorner_to_itsLastAppearance[{elemLgr, elemLgrCorn}][1];
                        while (elemLgrAndElemLgrCorner_to_adaptedCorner.count({lastAppearanceLgr, lastAppearanceLgrEquivCorner}) == 0) {
                            int tempElemLgr =  lastAppearanceLgr;
                            int tempElemLgrCorner =  lastAppearanceLgrEquivCorner;
                            lastAppearanceLgr =  vanishedRefinedCorner_to_itsLastAppearance[{ tempElemLgr ,  tempElemLgrCorner}][0];
                            lastAppearanceLgrEquivCorner =  vanishedRefinedCorner_to_itsLastAppearance[{ tempElemLgr ,  tempElemLgrCorner}][1];
                        }
                    }
                } // end-if(elemLgr>-1)
                adaptedCorn =   elemLgrAndElemLgrCorner_to_adaptedCorner[{lastAppearanceLgr, lastAppearanceLgrEquivCorner}];
            }
            aux_face_to_point[face].push_back(adaptedCorn);
        }
    } // end-adapted-faces
    // Adapted/Leaf-Grid-View face_to_point.
    adapted_face_to_point.reserve(face_count, num_points);
    for (std::size_t face = 0; face < face_count; ++face) {
        adapted_face_to_point.appendRow(aux_face_to_point[face].begin(), aux_face_to_point[face].end());
    }
}

   
void CpGrid::populateRefinedFaces(std::vector<Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<2,3>>>& refined_faces_vec,
                                  std::vector<Dune::cpgrid::EntityVariableBase<enum face_tag>>& mutable_refined_face_tags_vec,
                                  std::vector<Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>>& mutable_refined_face_normals_vec,
                                  std::vector<Opm::SparseTable<int>>& refined_face_to_point_vec,
                                  const std::vector<std::size_t>& refined_face_count_vec,
                                  std::map<std::array<int,2>,std::array<int,2>> refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace,
                                  std::map<std::array<int,2>,std::array<int,2>> elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                                  std::map<std::array<int,2>, std::array<int,2>> vanishedRefinedCorner_to_itsLastAppearance,
                                  const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                  const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                  std::map<std::array<int,2>,int> markedElemAndEquivRefinedCorn_to_corner)
{
    for (std::size_t shiftedLevel = 0; shiftedLevel < refined_face_count_vec.size(); ++shiftedLevel) {

        const auto& level = shiftedLevel + (this->data_.size());
        // Store the refined faces
        refined_faces_vec[shiftedLevel].resize(refined_face_count_vec[shiftedLevel]);
        mutable_refined_face_tags_vec[shiftedLevel].resize(refined_face_count_vec[shiftedLevel]);
        mutable_refined_face_normals_vec[shiftedLevel].resize(refined_face_count_vec[shiftedLevel]);
        //
        // Auxiliary integer to count all the points in refined_face_to_point.
        std::size_t refined_num_points = 0;
        // Auxiliary vector to store refined_face_to_point with non consecutive indices.
        std::vector<std::vector<int>> aux_refined_face_to_point;
        aux_refined_face_to_point.resize(refined_face_count_vec[shiftedLevel]);
        for (std::size_t face = 0; face < refined_face_count_vec[shiftedLevel]; ++face) {

            const auto& elemLgr = refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace[{level,face}][0];
            const auto& elemLgrFace = refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace[{level,face}][1];
            const auto& elemLgrFaceEntity =  Dune::cpgrid::EntityRep<1>(elemLgrFace, true);
            Opm::SparseTable<int>::mutable_row_type preAdapt_face_to_point;
            const auto& elemLgrData = markedElem_to_itsLgr[elemLgr];

            // Get the face geometry.
            refined_faces_vec[shiftedLevel][face] = (*(elemLgrData->geometry_.geomVector(std::integral_constant<int,1>())))[elemLgrFaceEntity];
            // Get the face tag.
            mutable_refined_face_tags_vec[shiftedLevel][face] = elemLgrData->face_tag_[elemLgrFaceEntity];
            // Get the face normal.
            mutable_refined_face_normals_vec[shiftedLevel][face] = elemLgrData->face_normals_[elemLgrFaceEntity];
            // Get face_to_point_ before adapting - we need to replace the level corners by the adapted ones.
            preAdapt_face_to_point = elemLgrData->face_to_point_[elemLgrFace];
            // Add the amount of points to the count num_points.
            refined_num_points += preAdapt_face_to_point.size();

            // Face_to_point
            for (int corn = 0; corn < static_cast<int>(preAdapt_face_to_point.size()); ++corn) {
                const auto& elemLgrCorn = preAdapt_face_to_point[corn];
                int refinedCorn;
                // Corner is stored in adapted_corners
                if (elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner.count({elemLgr, elemLgrCorn}) == 1)  {
                    refinedCorn = elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner[{elemLgr, elemLgrCorn}][1]; 
                }
                else{
                    // Corner might have vanished - Search its equivalent lgr-corner in that case -
                    // last lgr where the corner appears -
                    int lastAppearanceLgr = 0; // It'll get rewritten
                    int lastAppearanceLgrEquivCorner = 0; // It'll get rewritten
                    if(markedElemAndEquivRefinedCorn_to_corner.count({elemLgr, elemLgrCorn}) == 1) {
                        const auto& markedCorner = markedElemAndEquivRefinedCorn_to_corner[{elemLgr, elemLgrCorn}];
                        lastAppearanceLgr = cornerInMarkedElemWithEquivRefinedCorner[markedCorner].back()[0];
                        lastAppearanceLgrEquivCorner = cornerInMarkedElemWithEquivRefinedCorner[markedCorner].back()[1];
                    }
                    else {
                        // To locate vanished corners, we need a while-loop, since {elemLgr, elemLgrcorner} leads to
                        // {neighboringElemLgr, neighboringElemLgrCornerIdx}, which might have also vanished.
                        // Then, use the lastest appearance of the current corner, meaning, the first (and unique one - by construction) that
                        // gives elemLgrAndElemLgrCorner_to_refinedCorner.count({lastAppearanceElemLgr, lastAppearanceElemLgCorner}) == 1).
                        lastAppearanceLgr = vanishedRefinedCorner_to_itsLastAppearance[{elemLgr, elemLgrCorn}][0];
                        // this corner lays on the area occupied by a coarse face that got refined and belonged to two marked elements.
                        // Get the index of this corner with respect to the greatest marked element index.
                        lastAppearanceLgrEquivCorner = vanishedRefinedCorner_to_itsLastAppearance[{elemLgr, elemLgrCorn}][1];
                        while (elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner.count({lastAppearanceLgr, lastAppearanceLgrEquivCorner}) == 0) {
                            int tempElemLgr =  lastAppearanceLgr;
                            int tempElemLgrCorner =  lastAppearanceLgrEquivCorner;
                            lastAppearanceLgr =  vanishedRefinedCorner_to_itsLastAppearance[{ tempElemLgr ,  tempElemLgrCorner}][0];
                            lastAppearanceLgrEquivCorner =  vanishedRefinedCorner_to_itsLastAppearance[{ tempElemLgr ,  tempElemLgrCorner}][1];
                        }
                    }
                    refinedCorn =  elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner[{lastAppearanceLgr, lastAppearanceLgrEquivCorner}][1];
                }
                aux_refined_face_to_point[face].push_back(refinedCorn);
            }
        }
        // Refined face_to_point.
        refined_face_to_point_vec[shiftedLevel].reserve(refined_face_count_vec[shiftedLevel], refined_num_points);
        for (std::size_t face = 0; face < refined_face_count_vec[shiftedLevel]; ++face) {
            refined_face_to_point_vec[shiftedLevel].appendRow(aux_refined_face_to_point[face].begin(), aux_refined_face_to_point[face].end());
        }
    }
}

void CpGrid::populateAdaptedCells(Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<3,3>>& adapted_cells,
                                  std::vector<std::array<int,8>>& adapted_cell_to_point,
                                  std::vector<int>& adapted_global_cell,
                                  const std::size_t& cell_count,
                                  cpgrid::OrientedEntityTable<0,1>& adapted_cell_to_face,
                                  cpgrid::OrientedEntityTable<1,0>& adapted_face_to_cell,
                                  std::unordered_map<int,std::array<int,2>> adaptedCell_to_elemLgrAndElemLgrCell,
                                  std::map<std::array<int,2>,int> elemLgrAndElemLgrFace_to_adaptedFace,
                                  const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                  std::map<std::array<int,2>,int> elemLgrAndElemLgrCorner_to_adaptedCorner, 
                                  std::map<std::array<int,2>, std::array<int,2>> vanishedRefinedCorner_to_itsLastAppearance, 
                                  const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                  const std::vector<int>& assignRefinedLevel,
                                  std::map<std::array<int,2>,int> markedElemAndEquivRefinedCorn_to_corner,
                                  const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                  const std::vector<std::array<int,3>>& cells_per_dim_vec)
{

    // --- Adapted cells ---
    // Store the adapted cells. Main difficulty: to lookup correctly the indices of the corners and faces of each cell.
    adapted_cells.resize(cell_count);
    adapted_cell_to_point.resize(cell_count);
    for (std::size_t cell = 0; cell < cell_count; ++cell) {

        const auto& elemLgr = adaptedCell_to_elemLgrAndElemLgrCell[cell][0];
        const auto& elemLgrCell = adaptedCell_to_elemLgrAndElemLgrCell[cell][1];
        const auto& elemLgrCellEntity =  Dune::cpgrid::EntityRep<0>(elemLgrCell, true);

        adapted_global_cell[cell] = current_view_data_->global_cell_[(elemLgr == -1) ? elemLgrCell : elemLgr];

        std::array<int,8> preAdapt_cell_to_point;
        Dune::cpgrid::OrientedEntityRange<1> preAdapt_cell_to_face;
        // Auxiliary cell_to_face
        std::vector<cpgrid::EntityRep<1>> aux_cell_to_face;

        if (elemLgr == -1) { // The value -1 represents the current_view_data_
            // Get the cell geometry.
            adapted_cells[cell] =  (*(current_view_data_->geometry_.geomVector(std::integral_constant<int,0>())))[elemLgrCellEntity];
            // Get pre-adapt corners of the cell that will be replaced with leaf view ones.
            preAdapt_cell_to_point = current_view_data_->cell_to_point_[elemLgrCell];
            // Get pre-adapt faces of the cell that will be replaced with leaf view ones.
            preAdapt_cell_to_face = current_view_data_->cell_to_face_[elemLgrCellEntity];
        }
        else {
            const auto& elemLgrData = markedElem_to_itsLgr[elemLgr]; // Recall elemLgr == elemIdx (from an element marked for refinement)
            // Get the cell geometry.
            adapted_cells[cell] =  (*(elemLgrData->geometry_.geomVector(std::integral_constant<int,0>())))[elemLgrCellEntity];
            // Get pre-adapt corners of the cell that will be replaced with leaf view ones.
            preAdapt_cell_to_point = elemLgrData->cell_to_point_[elemLgrCell];
            // Get pre-adapt faces of the cell that will be replaced with leaf view ones.
            preAdapt_cell_to_face = elemLgrData->cell_to_face_[elemLgrCellEntity];
        }
        // Cell to point.
        for (int corn = 0; corn < 8; ++corn) {
            int adaptedCorn;
            const auto& preAdaptCorn = preAdapt_cell_to_point[corn];
            if (  elemLgrAndElemLgrCorner_to_adaptedCorner.count({elemLgr, preAdaptCorn}) == 0 ) {
                // Corner might have vanished - Search its equivalent lgr-corner in that case -
                // last lgr where the corner appears -
                int lastAppearanceLgr = 0; // It'll get rewritten
                int lastAppearanceLgrEquivCorner = 0; // It'll get rewritten
                if ((elemLgr ==-1) && (!cornerInMarkedElemWithEquivRefinedCorner[preAdaptCorn].empty())) {
                    lastAppearanceLgr = cornerInMarkedElemWithEquivRefinedCorner[preAdaptCorn].back()[0];
                    lastAppearanceLgrEquivCorner = cornerInMarkedElemWithEquivRefinedCorner[preAdaptCorn].back()[1];
                }
                if (elemLgr > -1) {
                    // Corner might have vanished - Search its equivalent lgr-corner in that case -
                    // last lgr where the corner appears -
                    if(markedElemAndEquivRefinedCorn_to_corner.count({elemLgr, preAdaptCorn}) == 1) {
                        const auto& markedCorner = markedElemAndEquivRefinedCorn_to_corner[{elemLgr, preAdaptCorn}];
                        lastAppearanceLgr = cornerInMarkedElemWithEquivRefinedCorner[markedCorner].back()[0];
                        lastAppearanceLgrEquivCorner = cornerInMarkedElemWithEquivRefinedCorner[markedCorner].back()[1];
                    }
                    else {
                        // To locate vanished corners, we need a while-loop, since {elemLgr, elemLgrcorner} leads to
                        // {neighboringElemLgr, neighboringElemLgrCornerIdx}, which might have also vanished.
                        // Then, use the lastest appearance of the current corner, meaning, the first (and unique one - by construction) that
                        // gives elemLgrAndElemLgrCorner_to_adaptedCorner.count({lastAppearanceElemLgr, lastAppearanceElemLgCorner}) == 1).
                        lastAppearanceLgr = vanishedRefinedCorner_to_itsLastAppearance[{elemLgr, preAdaptCorn}][0];
                        // this corner lays on the area occupied by a coarse face that got refined and belonged to two marked elements.
                        // Get the index of this corner with respect to the greatest marked element index.
                        lastAppearanceLgrEquivCorner = vanishedRefinedCorner_to_itsLastAppearance[{elemLgr, preAdaptCorn}][1];
                        while (elemLgrAndElemLgrCorner_to_adaptedCorner.count({lastAppearanceLgr, lastAppearanceLgrEquivCorner}) == 0) {
                            int tempElemLgr =  lastAppearanceLgr;
                            int tempElemLgrCorner =  lastAppearanceLgrEquivCorner;
                            lastAppearanceLgr =  vanishedRefinedCorner_to_itsLastAppearance[{ tempElemLgr ,  tempElemLgrCorner}][0];
                            lastAppearanceLgrEquivCorner =  vanishedRefinedCorner_to_itsLastAppearance[{ tempElemLgr ,  tempElemLgrCorner}][1];
                        }
                    }
                }
                adaptedCorn =  elemLgrAndElemLgrCorner_to_adaptedCorner[{lastAppearanceLgr, lastAppearanceLgrEquivCorner}];
            }
            // Corner is stored in adapted_corners
            else {
                assert(  elemLgrAndElemLgrCorner_to_adaptedCorner.count({elemLgr, preAdaptCorn}) == 1);
                adaptedCorn =  elemLgrAndElemLgrCorner_to_adaptedCorner[{elemLgr, preAdaptCorn}];
            }
            adapted_cell_to_point[cell][corn] = adaptedCorn;
        } // end-cell_to_point
        
        // Cell to face.
        for (const auto& face : preAdapt_cell_to_face) {
            const auto& preAdaptFace = face.index();
            int adaptedFace;
            // Face might have vanished - Search its refined lgr-children faces in that case -
            // last lgr where the face appears
            if (elemLgrAndElemLgrFace_to_adaptedFace.count({elemLgr, preAdaptFace}) == 0) {
                if (elemLgr ==-1) { // Coarse face got replaced by its children - from the last appearance of the marked face. 
                    assert(!faceInMarkedElemAndRefinedFaces[preAdaptFace].empty());
                    const auto& lastAppearanceLgr = faceInMarkedElemAndRefinedFaces[preAdaptFace].back().first;
                    const auto& lastAppearanceLgrFaces = faceInMarkedElemAndRefinedFaces[preAdaptFace].back().second;
                    for (const auto& refinedFace : lastAppearanceLgrFaces) {
                        adaptedFace = elemLgrAndElemLgrFace_to_adaptedFace[{lastAppearanceLgr, refinedFace}];
                        aux_cell_to_face.push_back({adaptedFace, face.orientation()});
                    }
                }
                if (elemLgr>-1) { // Refined face vanished and its equivalent refined face from a neighboring lgr got stored.
                    // Get shifted level
                    const auto& shiftedLevel = assignRefinedLevel[elemLgr] - (this->data_.size());
                    // Get the index of the marked face where the refined face was born.
                    const auto& markedFace = getParentFaceWhereNewRefinedFaceLaysOn(cells_per_dim_vec[shiftedLevel], preAdaptFace,
                                                                                    markedElem_to_itsLgr[elemLgr], elemLgr);
                    // Get the last LGR (marked element) where the marked face appeared.
                    const auto& lastLgrWhereMarkedFaceAppeared = faceInMarkedElemAndRefinedFaces[markedFace].back().first;
                    /** To do: modify the method to handle (and code) the case when  shifted level != shifted level last lgr. */
                    const auto& lastAppearanceLgrEquivFace = replaceLgr1FaceIdxByLgr2FaceIdx(cells_per_dim_vec[shiftedLevel], preAdaptFace,
                                                                                             markedElem_to_itsLgr[elemLgr]);
                    adaptedFace = elemLgrAndElemLgrFace_to_adaptedFace[{lastLgrWhereMarkedFaceAppeared, lastAppearanceLgrEquivFace}];
                    aux_cell_to_face.push_back({adaptedFace, face.orientation()});
                }
            }
            // Face is stored in adapted_faces
            else {
                assert( elemLgrAndElemLgrFace_to_adaptedFace.count({elemLgr, preAdaptFace}) == 1);
                adaptedFace = elemLgrAndElemLgrFace_to_adaptedFace[{elemLgr, preAdaptFace}];
                aux_cell_to_face.push_back({adaptedFace, face.orientation()});
            }
        } // end-cell_to_face
        // Adapted/Leaf-grid-view cell to face.
        adapted_cell_to_face.appendRow(aux_cell_to_face.begin(), aux_cell_to_face.end());
    } // adapted_cells

    // Adapted/Leaf-grid-view face to cell.
    adapted_cell_to_face.makeInverseRelation(adapted_face_to_cell);
}

   
void CpGrid::populateRefinedCells(std::vector<Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<3,3>>>& refined_cells_vec,
                                  std::vector<std::vector<std::array<int,8>>>& refined_cell_to_point_vec,
                                  std::vector<std::vector<int>>& refined_global_cell_vec,
                                  const std::vector<std::size_t>& refined_cell_count_vec,
                                  std::vector<cpgrid::OrientedEntityTable<0,1>>& refined_cell_to_face_vec,
                                  std::vector<cpgrid::OrientedEntityTable<1,0>>& refined_face_to_cell_vec,
                                  std::map<std::array<int,2>,std::array<int,2>> refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                  std::map<std::array<int,2>,std::array<int,2>> elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace,
                                  const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                  std::map<std::array<int,2>,std::array<int,2>> elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner, 
                                  std::map<std::array<int,2>, std::array<int,2>> vanishedRefinedCorner_to_itsLastAppearance, 
                                  const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                  const std::vector<int>& assignRefinedLevel,
                                  std::map<std::array<int,2>,int> markedElemAndEquivRefinedCorn_to_corner,
                                  const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                  const std::vector<std::array<int,3>>&  cells_per_dim_vec)
{
    // --- Refined cells ---
    for (std::size_t shiftedLevel = 0; shiftedLevel < refined_cell_count_vec.size(); ++shiftedLevel) {

        const auto& level = shiftedLevel + (this->data_.size());
        refined_cells_vec[shiftedLevel].resize(refined_cell_count_vec[shiftedLevel]);
        refined_cell_to_point_vec[shiftedLevel].resize(refined_cell_count_vec[shiftedLevel]);
          refined_global_cell_vec[shiftedLevel].resize(refined_cell_count_vec[shiftedLevel]);
    
        for (std::size_t cell = 0; cell < refined_cell_count_vec[shiftedLevel]; ++cell) {
        
            const auto& elemLgr = refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell[{level, cell}][0];
            assert(elemLgr >-1);
            const auto& elemLgrCell = refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell[{level, cell}][1];
            const auto& elemLgrCellEntity =  Dune::cpgrid::EntityRep<0>(elemLgrCell, true);
            const auto& elemLgrData = markedElem_to_itsLgr[elemLgr];
            assert( markedElem_to_itsLgr[elemLgr]!=nullptr);

            std::array<int,8> preAdapt_cell_to_point;
            Dune::cpgrid::OrientedEntityRange<1> preAdapt_cell_to_face;
            // Auxiliary cell_to_face
            std::vector<cpgrid::EntityRep<1>> aux_refined_cell_to_face;

            refined_global_cell_vec[shiftedLevel][cell] = cell; // current_view_data_ -> global_cell_[elemLgr]; instead?
            // Get the cell geometry.
            refined_cells_vec[shiftedLevel][cell] =  (*(elemLgrData->geometry_.geomVector(std::integral_constant<int,0>())))[elemLgrCellEntity];
            // Get pre-adapt corners of the cell that will be replaced with leaf view ones.
            preAdapt_cell_to_point = elemLgrData->cell_to_point_[elemLgrCell];
            // Get pre-adapt faces of the cell that will be replaced with leaf view ones.
            preAdapt_cell_to_face = elemLgrData->cell_to_face_[elemLgrCellEntity];

           

            // Cell to point.
            for (int corn = 0; corn < 8; ++corn) {
                int refinedCorn;
                const auto& preAdaptCorn = preAdapt_cell_to_point[corn];
                if ( elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner.count({elemLgr, preAdaptCorn}) == 0 ) { 
                    // Corner might have vanished - Search its equivalent lgr-corner in that case -
                    // last lgr where the corner appears -
                    int lastAppearanceLgr = 0; // It'll be rewritten
                    int lastAppearanceLgrEquivCorner = 0; // It'll be rewritten
                    if(markedElemAndEquivRefinedCorn_to_corner.count({elemLgr, preAdaptCorn}) == 1) {
                        const auto& markedCorner = markedElemAndEquivRefinedCorn_to_corner[{elemLgr, preAdaptCorn}];
                        lastAppearanceLgr = cornerInMarkedElemWithEquivRefinedCorner[markedCorner].back()[0];
                        lastAppearanceLgrEquivCorner = cornerInMarkedElemWithEquivRefinedCorner[markedCorner].back()[1];
                    }
                    else {
                        // To locate vanished corners, we need a while-loop, since {elemLgr, elemLgrcorner} leads to
                        // {neighboringElemLgr, neighboringElemLgrCornerIdx}, which might have also vanished.
                        // Then, use the lastest appearance of the current corner, meaning, the first (and unique one - by construction) that
                        // gives elemLgrAndElemLgrCorner_to_adaptedCorner.count({lastAppearanceElemLgr, lastAppearanceElemLgCorner}) == 1).
                        lastAppearanceLgr = vanishedRefinedCorner_to_itsLastAppearance[{elemLgr, preAdaptCorn}][0];
                        // this corner lays on the area occupied by a coarse face that got refined and belonged to two marked elements.
                        // Get the index of this corner with respect to the greatest marked element index.
                        lastAppearanceLgrEquivCorner = vanishedRefinedCorner_to_itsLastAppearance[{elemLgr, preAdaptCorn}][1];
                        while (elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner.count({lastAppearanceLgr, lastAppearanceLgrEquivCorner}) == 0) {
                            int tempElemLgr =  lastAppearanceLgr;
                            int tempElemLgrCorner =  lastAppearanceLgrEquivCorner;
                            lastAppearanceLgr =  vanishedRefinedCorner_to_itsLastAppearance[{ tempElemLgr ,  tempElemLgrCorner}][0];
                            lastAppearanceLgrEquivCorner =  vanishedRefinedCorner_to_itsLastAppearance[{ tempElemLgr ,  tempElemLgrCorner}][1];
                        }
                    }
                    /** Potential different level... */
                    refinedCorn = elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner[{lastAppearanceLgr, lastAppearanceLgrEquivCorner}][1];
                }
                // Corner is stored in adapted_corners
                else {
                    assert( elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner.count({elemLgr, preAdaptCorn}) == 1);
                    refinedCorn = elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner[{elemLgr, preAdaptCorn}][1];
                }
                refined_cell_to_point_vec[shiftedLevel][cell][corn] = refinedCorn;
            } // end-cell_to_point

            // Cell to face.
            for (const auto& face : preAdapt_cell_to_face) {
                const auto& preAdaptFace = face.index();
                int refinedFace;
                // Face might have vanished - Search its refined lgr-children faces in that case -
                // last lgr where the face appears
                if (elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace.count({elemLgr, preAdaptFace}) == 0) {
                    // Get the index of the marked face where the refined face was born.
                    const auto& markedFace = getParentFaceWhereNewRefinedFaceLaysOn(cells_per_dim_vec[shiftedLevel], preAdaptFace,
                                                                                    markedElem_to_itsLgr[elemLgr], elemLgr);
                    // Get the last LGR (marked element) where the marked face appeared.
                    int lastLgrWhereMarkedFaceAppeared = faceInMarkedElemAndRefinedFaces[markedFace].back().first;
                    /** To do: modify the method (and the code here) to take into account the case where the shifted level and
                        the shifted level of the lastAppearance do not coincide.
                        Alternatively, do not allow neighboring cells (sharing a face) to belong to different levels. */
                    const auto& lastAppearanceLgrEquivFace = replaceLgr1FaceIdxByLgr2FaceIdx(cells_per_dim_vec[shiftedLevel], preAdaptFace,
                                                                                             markedElem_to_itsLgr[elemLgr]);
                    /** To be modified! Refined face could belong to other level! Check similar aproach in other auxiliary methods */
                    refinedFace = elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace[{lastLgrWhereMarkedFaceAppeared, lastAppearanceLgrEquivFace}][1];
                    aux_refined_cell_to_face.push_back({refinedFace, face.orientation()});
                }
                // Face is stored in adapted_faces
                else {
                    assert( elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace.count({elemLgr, preAdaptFace}) == 1);
                    /** To be modified! Refined face could belong to other level! Check similar aproach in other auxiliary methods */
                    refinedFace =  elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace[{elemLgr, preAdaptFace}][1];
                    aux_refined_cell_to_face.push_back({refinedFace, face.orientation()});
                }
            } // end-cell_to_face
            // Refined cell to face.
            refined_cell_to_face_vec[shiftedLevel].appendRow(aux_refined_cell_to_face.begin(), aux_refined_cell_to_face.end());
        } // refined_cells
        // Refined face to cell.
        refined_cell_to_face_vec[shiftedLevel].makeInverseRelation(refined_face_to_cell_vec[shiftedLevel]);
    } // end-shiftedLevel-for-loop
}


std::array<int,3>  CpGrid::getRefinedCornerIJK(const std::array<int,3>& cells_per_dim, int cornerIdxInLgr) const
{
    // Order defined in Geometry::refine
    //  (j*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (i*(cells_per_dim[2]+1)) + k
    std::array<int,3> ijk;
    ijk[2] = cornerIdxInLgr % (cells_per_dim[2] +1);
    cornerIdxInLgr -= ijk[2];
    cornerIdxInLgr /= (cells_per_dim[2] +1);
    ijk[0] = cornerIdxInLgr % (cells_per_dim[0]+1);
    cornerIdxInLgr -=ijk[0];
    ijk[1] = cornerIdxInLgr / (cells_per_dim[0]+1);
    return ijk;
}

std::array<int,3>  CpGrid::getRefinedFaceIJK(const std::array<int,3>& cells_per_dim, int faceIdxInLgr,
                                             const std::shared_ptr<cpgrid::CpGridData>& elemLgr_ptr) const
{
    // Order defined in Geometry::refine
    // K_FACES  (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) + i
    // I_FACES  (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1))
    //           + (i*cells_per_dim[1]*cells_per_dim[2]) + (k*cells_per_dim[1]) + j
    // J_FACES  (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2] +1))
    //                    + ((cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2])
    //                    + (j*cells_per_dim[0]*cells_per_dim[2]) + (i*cells_per_dim[2]) + k
    const auto& faceEntity =  Dune::cpgrid::EntityRep<1>(faceIdxInLgr, true);
    const auto& faceTag = elemLgr_ptr ->face_tag_[faceEntity];
    std::array<int,3> ijk;
    switch (faceTag) {
    case I_FACE:
        faceIdxInLgr -= (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1));
        // faceIdxInLgr =  (i*cells_per_dim[1]*cells_per_dim[2]) + (k*cells_per_dim[1]) + j
        ijk[1] = faceIdxInLgr % cells_per_dim[1];
        faceIdxInLgr -= ijk[1]; // (i*cells_per_dim[1]*cells_per_dim[2]) + (k*cells_per_dim[1])
        faceIdxInLgr /= cells_per_dim[1]; // (i*cells_per_dim[2]) + k
        ijk[2] = faceIdxInLgr % cells_per_dim[2];
        faceIdxInLgr -=ijk[2]; // i*cells_per_dim[2]
        ijk[0] = faceIdxInLgr / cells_per_dim[2];   
        break;
    case J_FACE:
        faceIdxInLgr -=  (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2] +1))
            + ((cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2]);
        // faceIdxInLgr =  (j*cells_per_dim[0]*cells_per_dim[2]) + (i*cells_per_dim[2]) + k
        ijk[2] = faceIdxInLgr % cells_per_dim[2];
        faceIdxInLgr -= ijk[2]; // (j*cells_per_dim[0]*cells_per_dim[2]) + (i*cells_per_dim[2])
        faceIdxInLgr /= cells_per_dim[2]; // (j*cells_per_dim[0]) + i
        ijk[0] = faceIdxInLgr % cells_per_dim[0];
        faceIdxInLgr -=ijk[0]; // j*cells_per_dim[0]
        ijk[1] = faceIdxInLgr / cells_per_dim[0];  
        break;
    case K_FACE:
        //  (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) + i
        ijk[0] = faceIdxInLgr % cells_per_dim[0];
        faceIdxInLgr -= ijk[0]; // (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0])
        faceIdxInLgr /= cells_per_dim[0]; // (k*cells_per_dim[1]) + j
        ijk[1] = faceIdxInLgr % cells_per_dim[1];
        faceIdxInLgr -=ijk[1]; // k*cells_per_dim[1]
        ijk[2] = faceIdxInLgr / cells_per_dim[1];
        break;
    default:
        OPM_THROW(std::logic_error, "FaceTag is not I, J, or K!");
    }
    return ijk;
}

bool CpGrid::isRefinedCornerInInteriorLgr(const std::array<int,3>& cells_per_dim, int cornerIdxInLgr) const
{
    const auto& ijk = getRefinedCornerIJK(cells_per_dim, cornerIdxInLgr);
    return ((ijk[0]%cells_per_dim[0] > 0) &&  (ijk[1]%cells_per_dim[1]>0) && (ijk[2]%cells_per_dim[2]>0));
}

bool CpGrid::isRefinedFaceInInteriorLgr(const std::array<int,3>& cells_per_dim, int faceIdxInLgr, const std::shared_ptr<cpgrid::CpGridData>& elemLgr_ptr) const
{

    int refined_k_faces = cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1);
    int refined_i_faces = (cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2];
    
    bool isKface = (faceIdxInLgr < refined_k_faces);
    bool isIface = (faceIdxInLgr >= refined_k_faces) && (faceIdxInLgr < refined_k_faces + refined_i_faces);
    bool isJface = (faceIdxInLgr >= refined_k_faces + refined_i_faces);

    const auto& ijk = getRefinedFaceIJK(cells_per_dim, faceIdxInLgr, elemLgr_ptr);
    return ((ijk[0]%cells_per_dim[0] > 0 && isIface) ||  (ijk[1]%cells_per_dim[1]>0 && isJface) || (ijk[2]%cells_per_dim[2]>0 && isKface));
}

bool CpGrid::isRefinedNewBornCornerOnLgrBoundary(const std::array<int,3>& cells_per_dim, int cornerIdxInLgr) const
{
    const auto& ijk = getRefinedCornerIJK(cells_per_dim, cornerIdxInLgr);
    bool isOnParentCell_I_FACEfalse_and_newBornCorn = ( (ijk[0] == 0) && ((ijk[1] % cells_per_dim[1] != 0) || (ijk[2] % cells_per_dim[2] !=0) ));
    bool isOnParentCell_I_FACEtrue_and_newBornCorn = ( (ijk[0] == cells_per_dim[0]) && ((ijk[1] % cells_per_dim[1] != 0) || (ijk[2] % cells_per_dim[2] !=0) ));
    bool isOnParentCell_J_FACEfalse_and_newBornCorn = ( (ijk[1] == 0) && ((ijk[0] % cells_per_dim[0] != 0) || (ijk[2] % cells_per_dim[2] !=0) ));
    bool isOnParentCell_J_FACEtrue_and_newBornCorn = ( (ijk[1] == cells_per_dim[1]) && ((ijk[0] % cells_per_dim[0] != 0) || (ijk[2] % cells_per_dim[2] !=0) ));
    bool isOnParentCell_K_FACEfalse_and_newBornCorn = ( (ijk[2] == 0) && ((ijk[1] % cells_per_dim[1] != 0) || (ijk[0] % cells_per_dim[0] !=0) ));
    bool isOnParentCell_K_FACEtrue_and_newBornCorn = ( (ijk[2] == cells_per_dim[2]) && ((ijk[1] % cells_per_dim[1] != 0) || (ijk[0] % cells_per_dim[0] !=0) ));
    bool isOnParentCell_I_FACE = isOnParentCell_I_FACEfalse_and_newBornCorn || isOnParentCell_I_FACEtrue_and_newBornCorn;
    bool isOnParentCell_J_FACE = isOnParentCell_J_FACEfalse_and_newBornCorn || isOnParentCell_J_FACEtrue_and_newBornCorn;
    bool isOnParentCell_K_FACE = isOnParentCell_K_FACEfalse_and_newBornCorn || isOnParentCell_K_FACEtrue_and_newBornCorn;
    return (isOnParentCell_I_FACE || isOnParentCell_J_FACE || isOnParentCell_K_FACE);
}

bool CpGrid::newRefinedCornerLaysOnEdge(const std::array<int,3>& cells_per_dim, int cornerIdxInLgr) const
{
    const auto& ijk = getRefinedCornerIJK(cells_per_dim, cornerIdxInLgr);
    // Edges laying on bottom face
    bool isNewBornOnEdge01 = (ijk[0] % cells_per_dim[0] != 0) && (ijk[1] == 0) && (ijk[2] == 0);
    bool isNewBornOnEdge23 = (ijk[0] % cells_per_dim[0] != 0) && (ijk[1] == cells_per_dim[1]) && ( ijk[2] == 0);
    bool isNewBornOnEdge02 = (ijk[0] == 0) && (ijk[1] % cells_per_dim[1] != 0) && (ijk[2] == 0);
    bool isNewBornOnEdge13 = (ijk[0] == cells_per_dim[0]) && (ijk[1] % cells_per_dim[1] != 0) && (ijk[2] == 0);

    // Edges connecting bottom and top faces
    bool isNewBornOnEdge04 = (ijk[0] == 0) && (ijk[1] == 0) && (ijk[2] % cells_per_dim[2] != 0);
    bool isNewBornOnEdge26 = (ijk[0] == 0) && (ijk[1] == cells_per_dim[1]) && (ijk[2] % cells_per_dim[2] != 0);
    bool isNewBornOnEdge15 = (ijk[0] == cells_per_dim[0]) && (ijk[1] == 0) && (ijk[2] % cells_per_dim[2] != 0);
    bool isNewBornOnEdge37 = (ijk[0] == cells_per_dim[0]) && (ijk[1] == cells_per_dim[1]) && (ijk[2] % cells_per_dim[2] != 0);

    // Edges laying on top face
    bool isNewBornOnEdge45 = (ijk[0] % cells_per_dim[0] != 0) && (ijk[1] == 0) && (ijk[2] == cells_per_dim[2]);
    bool isNewBornOnEdge67 = (ijk[0] % cells_per_dim[0] != 0) && (ijk[1] == cells_per_dim[1]) && ( ijk[2] == cells_per_dim[2]);
    bool isNewBornOnEdge46 = (ijk[0] == 0) && (ijk[1] % cells_per_dim[1] != 0) && (ijk[2] == cells_per_dim[2]);
    bool isNewBornOnEdge57 = (ijk[0] == cells_per_dim[0]) && (ijk[1] % cells_per_dim[1] != 0) && (ijk[2] == cells_per_dim[2]);

    bool isOnEdge = isNewBornOnEdge01 || isNewBornOnEdge23 || isNewBornOnEdge02 || isNewBornOnEdge13 ||
        isNewBornOnEdge04 || isNewBornOnEdge26 || isNewBornOnEdge15 || isNewBornOnEdge37 ||
        isNewBornOnEdge45 || isNewBornOnEdge67 || isNewBornOnEdge46 || isNewBornOnEdge57;

    return isOnEdge;
}


bool CpGrid::isRefinedFaceOnLgrBoundary(const std::array<int,3>& cells_per_dim, int faceIdxInLgr,
                                        const std::shared_ptr<cpgrid::CpGridData>& elemLgr_ptr) const
{
    const auto& ijk = getRefinedFaceIJK(cells_per_dim, faceIdxInLgr, elemLgr_ptr);

    int refined_k_faces = cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1);
    int refined_i_faces = (cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2];

    bool isKface = (faceIdxInLgr < refined_k_faces);
    bool isIface = (faceIdxInLgr >= refined_k_faces) && (faceIdxInLgr < refined_k_faces + refined_i_faces);
    bool isJface = (faceIdxInLgr >= refined_k_faces + refined_i_faces);

    bool isOnParentCell_I_FACE = isIface && (ijk[0] % cells_per_dim[0] == 0) && (ijk[1]<cells_per_dim[1]) && (ijk[2]<cells_per_dim[2]);
    bool isOnParentCell_J_FACE = isJface && (ijk[1] % cells_per_dim[1] == 0) && (ijk[0]<cells_per_dim[0]) && (ijk[2]<cells_per_dim[2]);
    bool isOnParentCell_K_FACE = isKface && (ijk[2] % cells_per_dim[2] == 0) && (ijk[0]<cells_per_dim[0]) && (ijk[1]<cells_per_dim[1]);

    return (isOnParentCell_I_FACE || isOnParentCell_J_FACE || isOnParentCell_K_FACE);
}

std::array<int,2> CpGrid::getParentFacesAssocWithNewRefinedCornLayingOnEdge(const std::array<int,3>& cells_per_dim, int cornerIdxInLgr, int elemLgr) const
{
    assert(newRefinedCornerLaysOnEdge(cells_per_dim, cornerIdxInLgr));

    const auto& parentCell_to_face = current_view_data_->cell_to_face_[cpgrid::EntityRep<0>(elemLgr, true)];
    if(parentCell_to_face.size()>6){
        OPM_THROW(std::logic_error, "The associted parent cell has more than six faces. Refinment/Adaptivity not supported yet.");
    }
    // Corners Order defined in Geometry::refine  (j*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (i*(cells_per_dim[2]+1)) + k
    const auto& ijk = getRefinedCornerIJK(cells_per_dim, cornerIdxInLgr);
    // Edges laying on bottom face
    bool isNewBornOnEdge01 = (ijk[0] % cells_per_dim[0] != 0) && (ijk[1] == 0) && (ijk[2] == 0);
    bool isNewBornOnEdge23 = (ijk[0] % cells_per_dim[0] != 0) && (ijk[1] == cells_per_dim[1]) && ( ijk[2] == 0);
    bool isNewBornOnEdge02 = (ijk[0] == 0) && (ijk[1] % cells_per_dim[1] != 0) && (ijk[2] == 0);
    bool isNewBornOnEdge13 = (ijk[0] == cells_per_dim[0]) && (ijk[1] % cells_per_dim[1] != 0) && (ijk[2] == 0);

    // Edges connecting bottom and top faces
    bool isNewBornOnEdge04 = (ijk[0] == 0) && (ijk[1] == 0) && (ijk[2] % cells_per_dim[2] != 0);
    bool isNewBornOnEdge26 = (ijk[0] == 0) && (ijk[1] == cells_per_dim[1]) && (ijk[2] % cells_per_dim[2] != 0);
    bool isNewBornOnEdge15 = (ijk[0] == cells_per_dim[0]) && (ijk[1] == 0) && (ijk[2] % cells_per_dim[2] != 0);
    bool isNewBornOnEdge37 = (ijk[0] == cells_per_dim[0]) && (ijk[1] == cells_per_dim[1]) && (ijk[2] % cells_per_dim[2] != 0);

    // Edges laying on top face
    bool isNewBornOnEdge45 = (ijk[0] % cells_per_dim[0] != 0) && (ijk[1] == 0) && (ijk[2] == cells_per_dim[2]);
    bool isNewBornOnEdge67 = (ijk[0] % cells_per_dim[0] != 0) && (ijk[1] == cells_per_dim[1]) && ( ijk[2] == cells_per_dim[2]);
    bool isNewBornOnEdge46 = (ijk[0] == 0) && (ijk[1] % cells_per_dim[1] != 0) && (ijk[2] == cells_per_dim[2]);
    bool isNewBornOnEdge57 = (ijk[0] == cells_per_dim[0]) && (ijk[1] % cells_per_dim[1] != 0) && (ijk[2] == cells_per_dim[2]);

    std::vector<int> auxFaces;
    auxFaces.reserve(2);

    for (const auto& face : parentCell_to_face) {
        const auto& faceEntity =  Dune::cpgrid::EntityRep<1>(face.index(), true);
        const auto& faceTag = current_view_data_->face_tag_[faceEntity];
        // Add I_FACE false
        bool addIfalse = isNewBornOnEdge02 || isNewBornOnEdge04 || isNewBornOnEdge26 || isNewBornOnEdge46;
        if( addIfalse && (faceTag == 0)  && (!face.orientation()))  {
            auxFaces.push_back(face.index());
        }
        // Add J_FACE false
        bool addJfalse = isNewBornOnEdge01 || isNewBornOnEdge04 || isNewBornOnEdge15 || isNewBornOnEdge45;
        if( addJfalse && (faceTag == 1)  && (!face.orientation()))  {
            auxFaces.push_back(face.index());
        }
        // Add K_FACE false
        bool addKfalse = isNewBornOnEdge01 ||  isNewBornOnEdge13 || isNewBornOnEdge23 || isNewBornOnEdge02;
        if( addKfalse && (faceTag == 2) && (!face.orientation())) {
            auxFaces.push_back(face.index());
        }
        // Add I_FACE true
        bool addItrue = isNewBornOnEdge13 || isNewBornOnEdge15 || isNewBornOnEdge37 || isNewBornOnEdge57;
        if( addItrue && (faceTag == 0)  && (face.orientation()))  {
            auxFaces.push_back(face.index());
        }
        // Add J_FACE true
        bool addJtrue = isNewBornOnEdge23|| isNewBornOnEdge26 || isNewBornOnEdge37 || isNewBornOnEdge67;
        if( addJtrue && (faceTag == 1)  && (face.orientation()))  {
            auxFaces.push_back(face.index());
        }
        // Add K_FACE true
        bool addKtrue = isNewBornOnEdge45 || isNewBornOnEdge67 || isNewBornOnEdge46 || isNewBornOnEdge57;
        if(addKtrue && (faceTag == 2) && (face.orientation())) {
            auxFaces.push_back(face.index());
        }
    }
    return {auxFaces[0], auxFaces[1]};
}

int CpGrid::getParentFaceWhereNewRefinedCornerLaysOn(const std::array<int,3>& cells_per_dim,
                                                     int cornerIdxInLgr, int parentCellIdxOnLevel0) const
{
    assert(isRefinedNewBornCornerOnLgrBoundary(cells_per_dim, cornerIdxInLgr));

    const auto& parentCell_to_face = current_view_data_->cell_to_face_[cpgrid::EntityRep<0>(parentCellIdxOnLevel0, true)];
    if(parentCell_to_face.size()>6){
        OPM_THROW(std::logic_error, "The associted parent cell has more than six faces. Refinment/Adaptivity not supported yet.");
    }
    const auto& ijk = getRefinedCornerIJK(cells_per_dim, cornerIdxInLgr);

    bool isOnParentCell_I_FACEfalse_and_newBornCorn = ( (ijk[0] == 0) && ((ijk[1] % cells_per_dim[1] != 0) || (ijk[2] % cells_per_dim[2] !=0) ));
    bool isOnParentCell_I_FACEtrue_and_newBornCorn = ( (ijk[0] == cells_per_dim[0]) && ((ijk[1] % cells_per_dim[1] != 0) || (ijk[2] % cells_per_dim[2] !=0) ));
    bool isOnParentCell_J_FACEfalse_and_newBornCorn = ( (ijk[1] == 0) && ((ijk[0] % cells_per_dim[0] != 0) || (ijk[2] % cells_per_dim[2] !=0) ));
    bool isOnParentCell_J_FACEtrue_and_newBornCorn = ( (ijk[1] == cells_per_dim[1]) && ((ijk[0] % cells_per_dim[0] != 0) || (ijk[2] % cells_per_dim[2] !=0) ));
    bool isOnParentCell_K_FACEfalse_and_newBornCorn = ( (ijk[2] == 0) && ((ijk[1] % cells_per_dim[1] != 0) || (ijk[0] % cells_per_dim[0] !=0) ));
    bool isOnParentCell_K_FACEtrue_and_newBornCorn = ( (ijk[2] == cells_per_dim[2]) && ((ijk[1] % cells_per_dim[1] != 0) || (ijk[0] % cells_per_dim[0] !=0) ));
    
    for (const auto& face : parentCell_to_face) {
        const auto& faceEntity =  Dune::cpgrid::EntityRep<1>(face.index(), true);
        const auto& faceTag = current_view_data_->face_tag_[faceEntity];
        if (isOnParentCell_I_FACEfalse_and_newBornCorn && (faceTag == 0) && !face.orientation()) { // I_FACE false
            return face.index();
        }
        if (isOnParentCell_I_FACEtrue_and_newBornCorn && (faceTag == 0) && face.orientation()) { // I_FACE true
            return face.index();
        }
        if (isOnParentCell_J_FACEfalse_and_newBornCorn && (faceTag == 1) && !face.orientation()) { // J_FACE false
            return face.index();
        }
        if (isOnParentCell_J_FACEtrue_and_newBornCorn && (faceTag == 1) && face.orientation()) { // J_FACE true
            return face.index();
        }
        if (isOnParentCell_K_FACEfalse_and_newBornCorn && (faceTag == 2) && !face.orientation()) { // K_FACE false
            return face.index();
        }
        if (isOnParentCell_K_FACEtrue_and_newBornCorn && (faceTag == 2) && face.orientation()) { // K_FACE true
            return face.index();
        }
    }
    OPM_THROW(std::logic_error, "Cannot find parent face index where new refined corner lays on.");
}

int CpGrid::getParentFaceWhereNewRefinedFaceLaysOn(const std::array<int,3>& cells_per_dim,
                                                   int faceIdxInLgr,
                                                   const std::shared_ptr<cpgrid::CpGridData>& elemLgr_ptr,
                                                   int parentCellIdxOnLevel0) const
{
    assert(isRefinedFaceOnLgrBoundary(cells_per_dim, faceIdxInLgr, elemLgr_ptr));
    const auto& ijk = getRefinedFaceIJK(cells_per_dim, faceIdxInLgr, elemLgr_ptr);
    const auto& parentCell_to_face = current_view_data_->cell_to_face_[cpgrid::EntityRep<0>(parentCellIdxOnLevel0, true)];
    
    if(parentCell_to_face.size()>6){
        OPM_THROW(std::logic_error, "The associted parent cell has more than six faces. Refinment/Adaptivity not supported yet.");
    }

    // Order defined in Geometry::refine
    // K_FACES  (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) + i
    // I_FACES  (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1))
    //           + (i*cells_per_dim[1]*cells_per_dim[2]) + (k*cells_per_dim[1]) + j
    // J_FACES  (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2] +1))
    //                    + ((cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2])
    //                    + (j*cells_per_dim[0]*cells_per_dim[2]) + (i*cells_per_dim[2]) + k
    int refined_k_faces = cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1);
    int refined_i_faces = (cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2];
    int refined_j_faces = cells_per_dim[0]*(cells_per_dim[1]+1)*cells_per_dim[2];
    assert( faceIdxInLgr < refined_k_faces + refined_i_faces + refined_j_faces);

    for (const auto& face : parentCell_to_face) {
        const auto& faceEntity =  Dune::cpgrid::EntityRep<1>(face.index(), true);
        const auto& faceTag = current_view_data_->face_tag_[faceEntity];
        if (faceIdxInLgr <  refined_k_faces ) { // It's a K_FACE
            if ((ijk[2] == 0) && (faceTag == 2) && !face.orientation()) { // {K_FACE, false}
                return face.index();
            }
            if ((ijk[2] == cells_per_dim[2]) && (faceTag == 2) && face.orientation()) { // {K_FACE, true}
                return face.index();
            }
        }
        if ((faceIdxInLgr >= refined_k_faces) && (faceIdxInLgr < refined_k_faces + refined_i_faces)) { // It's I_FACE
            if ((ijk[0] == 0) && (faceTag == 0) && !face.orientation()) { // {I_FACE, false}
                return face.index();
            }
            if ((ijk[0] == cells_per_dim[0]) && (faceTag == 0) && face.orientation()) { // {I_FACE, true}
                return face.index();
            }
        }
        if (faceIdxInLgr >= refined_k_faces + refined_i_faces) {// It's J_FACE
            if ((ijk[1] == 0) && (faceTag == 1) && !face.orientation()) { // {J_FACE, false}
                return face.index();
            }
            if ((ijk[1] == cells_per_dim[1]) && (faceTag == 1) && face.orientation()) { // {J_FACE, true}
                return face.index();
            }
        }
    }
    OPM_THROW(std::logic_error, "Cannot find parent face index where the new refined face lays on.");
}

int CpGrid::replaceLgr1CornerIdxByLgr2CornerIdx(const std::array<int,3>& cells_per_dim, int cornerIdxLgr1) const
{   
    const auto& ijk = getRefinedCornerIJK(cells_per_dim, cornerIdxLgr1);
    // Order defined in Geometry::refine
    //  (j*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (i*(cells_per_dim[2]+1)) + k
    if (ijk[0] == cells_per_dim[0]) {
        return   (ijk[1]*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + ijk[2];
    }
    if (ijk[1] == cells_per_dim[1]) {
        return  (ijk[0]*(cells_per_dim[2]+1)) + ijk[2];
    }
    if (ijk[2] == cells_per_dim[2]) {
        return  (ijk[1]*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (ijk[0]*(cells_per_dim[2]+1));
    }
    else {
        OPM_THROW(std::logic_error, "Cannot convert corner index from one LGR to its neighboring LGR.");
    }
}

int CpGrid::replaceLgr1CornerIdxByLgr2CornerIdx(const std::array<int,3>& cells_per_dim, int cornerIdxLgr1, int elemLgr1, int parentFaceLastAppearanceIdx) const
{
    assert(newRefinedCornerLaysOnEdge(cells_per_dim, cornerIdxLgr1));
    const auto& faces = getParentFacesAssocWithNewRefinedCornLayingOnEdge(cells_per_dim, cornerIdxLgr1, elemLgr1);
    assert( (faces[0] == parentFaceLastAppearanceIdx) || (faces[1] == parentFaceLastAppearanceIdx));

    const auto& ijk = getRefinedCornerIJK(cells_per_dim, cornerIdxLgr1);
    const auto& parentCell_to_face = this->data_[0]->cell_to_face_[cpgrid::EntityRep<0>(elemLgr1, true)];

    if(parentCell_to_face.size()>6){
        OPM_THROW(std::logic_error, "The associted parent cell has more than six faces. Refinment/Adaptivity not supported yet.");
    }
    // Since lgr1 is associated with an element index smaller than "a last appearance lgr", then the only possibilities are
    // I_FACE true, J_FACE true, K_FACE true

    // Order defined in Geometry::refine
    //  (j*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (i*(cells_per_dim[2]+1)) + k

    for (const auto& face : parentCell_to_face) {
        const auto& faceEntity =  Dune::cpgrid::EntityRep<1>(face.index(), true);
        const auto& faceTag = current_view_data_->face_tag_[faceEntity];
        if (parentFaceLastAppearanceIdx == face.index() && face.orientation()) {
            if (faceTag == 0) { // I_FACE true
                // The same new born refined corner will have equal values of j and k, but i == 0 instead of cells_per_dim[0]
                return  (ijk[1]*(cells_per_dim[0]+1)*(cells_per_dim[2]+1))  + ijk[2];
            }
            if (faceTag == 1) { // J_FACE true
                // The same new born refined corner will have equal values of i and k, but j == 0 instead of cells_per_dim[1]
                return  (ijk[0]*(cells_per_dim[2]+1)) + ijk[2];
            }
            if (faceTag == 2) { // K_FACE true
                // The same new born refined corner will have equal values of i and j, but k == 0 instead of cells_per_dim[2]
                return  (ijk[1]*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (ijk[0]*(cells_per_dim[2]+1));
            }
        }
    }
    OPM_THROW(std::logic_error, "Cannot convert corner index from one LGR to its neighboring LGR.");
}

int  CpGrid::replaceLgr1FaceIdxByLgr2FaceIdx(const std::array<int,3>& cells_per_dim, int faceIdxInLgr1,
                                             const std::shared_ptr<Dune::cpgrid::CpGridData>& elemLgr1_ptr) const
{
    const auto& ijk = getRefinedFaceIJK(cells_per_dim, faceIdxInLgr1, elemLgr1_ptr);
    // lgr1 represents an element index < lgr2 (neighboring cells sharing a face with lgr1-element)
    // Order defined in Geometry::refine
    // K_FACES (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) + i
    // I_FACES  (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1))
    //           + (i*cells_per_dim[1]*cells_per_dim[2]) + (k*cells_per_dim[1]) + j
    // J_FACES  (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2] +1))
    //                    + ((cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2])
    //                    + (j*cells_per_dim[0]*cells_per_dim[2]) + (i*cells_per_dim[2]) + k
    if (ijk[0] == cells_per_dim[0]) {
        return  (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1)) + (ijk[2]*cells_per_dim[1]) + ijk[1];
    }
    if (ijk[1] == cells_per_dim[1]) {
        return  (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2] +1))
            + ((cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2]) + (ijk[0]*cells_per_dim[2]) + ijk[2];
    }
    if (ijk[2] == cells_per_dim[2]) {
        return  (ijk[1]*cells_per_dim[0]) + ijk[0];
    }
    else {
        OPM_THROW(std::logic_error, "Cannot convert face index from one LGR to its neighboring LGR.");
    }
}




} // namespace Dune


