//===========================================================================
//
// File: CpGridDataTraits.hpp
//
// Created: Fri Mar 08 2023
//
// Author(s): Markus Blatt <markus@dr-blatt.de>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2014-2015 Dr. Blatt - HPC-Simulartion-Software & Services
  Copyright 2009-2023 Equinor ASA.

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

#ifndef OPM_CPGRIDDATATRAITS_HEADER
#define OPM_CPGRIDDATATRAITS_HEADER

#include <dune/common/parallel/mpihelper.hh>
#ifdef HAVE_DUNE_ISTL
#include <dune/istl/owneroverlapcopy.hh>
#endif

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 7)
#include <dune/common/parallel/communication.hh>
#include <dune/common/parallel/variablesizecommunicator.hh>
#else
#include <dune/common/parallel/collectivecommunication.hh>
#include <opm/grid/utility/VariableSizeCommunicator.hpp>
#endif

namespace Dune
{
namespace cpgrid
{
struct CpGridDataTraits
{
    /// \brief The type of the collective communication.
    using MPICommunicator = MPIHelper::MPICommunicator;

#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 7)
    using Communication = Dune::Communication<MPICommunicator>;
    using CollectiveCommunication = Dune::Communication<MPICommunicator>;
#else
    using CollectiveCommunication = Dune::CollectiveCommunication<MPICommunicator>;
    using Communication = Dune::CollectiveCommunication<MPICommunicator>;
#endif

#ifdef HAVE_DUNE_ISTL
    /// \brief The type of the set of the attributes
    using AttributeSet = Dune::OwnerOverlapCopyAttributeSet::AttributeSet;
#else
    /// \brief The type of the set of the attributes
    enum AttributeSet{owner, overlap, copy};
#endif

#if HAVE_MPI
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 7)
    /// \brief The type of the  Communicator.
    using Communicator = Dune::VariableSizeCommunicator<>;
#else
    /// \brief The type of the Communicator.
    using Communicator = Opm::VariableSizeCommunicator<>;
#endif

    /// \brief The type of the map describing communication interfaces.
    using InterfaceMap = Communicator::InterfaceMap;

    /// \brief type of OwnerOverlap communication for cells
    using CommunicationType = Dune::OwnerOverlapCopyCommunication<int,int>;

    /// \brief The type of the parallel index set
    using  ParallelIndexSet = typename CommunicationType::ParallelIndexSet;

    /// \brief The type of the remote indices information
    using RemoteIndices = Dune::RemoteIndices<ParallelIndexSet>;
#else
    using InterfaceMap = std::map<int, std::list<int> >;
#endif // HAVE_MPI
};

} // end namespace cpgrid
} // end space Opm

#endif // OPM_CPGRIDDATATRAITS_HEADER
