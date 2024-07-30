//===========================================================================
//
// File: GridPartitioning.hpp
//
// Created: Mon Sep  7 10:09:13 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Bård Skaflestad     <bard.skaflestad@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2009, 2010, 2013 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.
  Copyright 2013       Dr. Markus Blatt - HPC-Simulation-Software & Services
  Copyright 2020, 2021 OPM-OP AS

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

#ifndef OPM_GRIDPARTITIONING_HEADER
#define OPM_GRIDPARTITIONING_HEADER

#include <vector>
#include <array>
#include <set>
#include <tuple>

#include <dune/common/parallel/mpihelper.hh>

#include <opm/grid/utility/OpmWellType.hpp>
#include <opm/grid/common/WellConnections.hpp>
#include <opm/grid/common/ZoltanGraphFunctions.hpp>

namespace Dune
{

    class CpGrid;

    struct OrderByFirst
    {
        bool operator()(const std::pair<int,int>& o, const std::pair<int,int>& v)
        {
            return o.first < v.first;
        }
    };

    /// Partition a CpGrid based on (ijk) coordinates, with splitting to ensure that each partition is connected.
    /// @param[in] grid the grid to partition
    /// @param[in] initial_split the number of parts in which to partition the grid, in each cardinal direction.
    ///                          Their product is the expected number of partitions produced.
    /// @param[out] num_part the resulting number of partitions. This may be lower than expected,
    ///                      because of inactive cells, or higher than expected,
    ///                      because of splits to ensure connectedness.
    /// @param[out] cell_part a vector containing, for each cell, its partition number
    void partition(const CpGrid& grid,
                   const std::array<int, 3>& initial_split,
                   int& num_part,
                   std::vector<int>& cell_part,
                   bool recursive = false,
                   bool ensureConnectivity = true);

    /// \brief Adds a layer of overlap cells to a partitioning.
    /// \param[in] grid The grid that is partitioned.
    /// \param[in] cell_part a vector containing each cells partition number.
    /// \param[out] cell_overlap a vector of sets that contains for each cell all
    ///             the partition numbers that it is an overlap cell of.
    /// \param[in] mypart The partition number of the processor.
    /// \param[in] all Whether to compute the overlap for all partions or just the
    ///            one associated by mypart.
    void addOverlapLayer(const CpGrid& grid,
                         const std::vector<int>& cell_part,
                         std::vector<std::set<int> >& cell_overlap,
                         int mypart, int overlapLayers, bool all=false);

    /// \brief Adds a layer of overlap cells to a partitioning.
    /// \param[in] grid The grid that is partitioned.
    /// \param[in] cell_part a vector containing each cells partition number.
    /// \param[inout] exportList List indices to export, each entry is a tuple
    /// of global index, process rank (to export to), attribute on remote.
    /// \param[inout] importList List indices to import, each entry is a tuple
    /// of global index, process rank (to import from), attribute here, local
    /// index here
    /// \param[in] cc The communication object
    /// \param[in] addCornerCells Switch for adding corner cells to overlap layer.
    /// \param[in] trans The transmissibilities on cell faces. When trans[i]==0, no overlap is added.
    /// \param[in] layer Number of overlap layers
    int addOverlapLayer(const CpGrid& grid, const std::vector<int>& cell_part,
                        std::vector<std::tuple<int,int,char>>& exportList,
                        std::vector<std::tuple<int,int,char,int>>& importList,
                        const Communication<Dune::MPIHelper::MPICommunicator>& cc,
                        bool addCornerCells, const double* trans, int layers = 1);

namespace cpgrid
{
#if HAVE_MPI
    /// \brief Creates lists as Zoltan would return for vanilla / user specified partitioning.
    ///
    /// \param grid The grid
    /// \param wells Pointer to vector with all possible wells (all perforations) of the problem.
    ///              nullptr is possible
    /// \param possibleFutureConnections Possible future connections of wells that might get added through an ACTIONX.
    ///                                  The grid will then be partitioned such that these connections are on the same
    ///                                  partition. If NULL, they will be neglected.
    /// \param transmissibilities C-array with transmissibilities or nullptr.
    /// \param parts The partitioning information. Contains for each cells the partition number (zero-based,
    ///              consecutive
    /// \return  A tuple consisting of a vector that contains for each local cell of the original grid the
    ///         the number of the process that owns it after repartitioning,
    ///         a vector containing a pair of name  and a boolean indicating whether this well has
    ///         perforated cells local to the process of all wells,
    ///         vector containing information for each exported cell (global id
    ///         of cell, process id to send to, attribute there), a vector containing
    ///         information for each imported cell (global index, process id that sends, attribute here, local index
    ///         here), and a WellConnections object containing information about the well connections
    ///         (if argument wells was not null and this is the root rank this will contain connections in
    ///          form of global indices)
    std::tuple<std::vector<int>, std::vector<std::pair<std::string,bool>>,
               std::vector<std::tuple<int,int,char> >,
               std::vector<std::tuple<int,int,char,int> >,
               WellConnections>
    createListsFromParts(const CpGrid& grid, const std::vector<cpgrid::OpmWellType> * wells,
                               const std::unordered_map<std::string, std::set<std::array<int,3>>>* possibleFutureConnections,
                               const double* transmissibilities, const std::vector<int>& parts,
                               bool allowDistributedWells, std::shared_ptr<cpgrid::CombinedGridWellGraph> gridAndWells = nullptr);

    /// \brief Creates a vanilla partitioning without a real loadbalancer
    ///
    /// The loadbalancing will take place on the cartesian grid.
    /// \param grid The grid
    /// \param wells Pointer to vector with all possible wells (all perforations) of the problem.
    ///              nullptr is possible
    /// \param possibleFutureConnections Possible future connections of wells that might get added through an ACTIONX.
    ///                                  The grid will then be partitioned such that these connections are on the same
    ///                                  partition. If NULL, they will be neglected.
    /// \param transmissibilities C-array with transmissibilities or nullptr.
    /// \return  A tuple consisting of a vector that contains for each local cell of the original grid the
    ///         the number of the process that owns it after repartitioning,
    ///         a vector containing a pair of name  and a boolean indicating whether this well has
    ///         perforated cells local to the process of all wells,
    ///         vector containing information for each exported cell (global id
    ///         of cell, process id to send to, attribute there), a vector containing
    ///         information for each imported cell (global index, process id that sends, attribute here, local index
    ///         here), and a WellConnections object containing information about the well connections
    ///         (if argument wells was not null and this is the root rank this will contain connections in
    ///          form of global indices)
    std::tuple<std::vector<int>, std::vector<std::pair<std::string,bool>>,
               std::vector<std::tuple<int,int,char> >,
               std::vector<std::tuple<int,int,char,int> >,
               WellConnections>
    vanillaPartitionGridOnRoot(const CpGrid& grid,
                               const std::vector<cpgrid::OpmWellType> * wells,
                               const std::unordered_map<std::string, std::set<std::array<int,3>>>* possibleFutureConnections,
                               const double* transmissibilities,
                               bool allowDistributedWells);
#endif
} // namespace cpgrid
} // namespace Dune


#endif // OPM_GRIDPARTITIONING_HEADER
