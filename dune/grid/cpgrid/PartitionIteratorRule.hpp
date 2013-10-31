//===========================================================================
//
// File: CpGrid.hpp
//
// Created: Fri Oct 31 2013
//
// Author(s): Markus Blatt <markus@dr-blatt.de>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2013 Dr. Blatt - HPC-Simulation & Services.
  Copyright 2013 Statoil ASA.

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
namespace Dune
{
namespace cpgrid
{
    /// A rule at what entities to stop
    ///
    template<PartitionIteratorType pitype>
    struct PartitionIteratorRule
    {
        enum {fullSet=false, emptySet=true};
        
        template<int codim>
        bool isInvalid(const Entity<codim>& e)
        {
            return true;
        }
    };
    
    template<>
    struct PartitionIteratorRule<Interior_Partition>
    {
        enum {fullSet=false, emptySet=false};
        template<int codim>
        bool isInvalid(const Entity<codim>& e)
        {
            if(e.partitionType()==InteriorEntity)
                return false;
            return true;
        }
    };
     
    template<>
    struct PartitionIteratorRule<InteriorBorder_Partition>
    {
        enum {fullSet=false, emptySet=false};
        template<int codim>
        bool isInvalid(const Entity<codim>& e)
        {
            if(e.partitionType()==InteriorEntity ||
               e.partitionType()==BorderEntity)
                return false;
            return true;
        }
    };
        
    template<>
    struct PartitionIteratorRule<Overlap_Partition>
    {
        enum {fullSet=false, emptySet=false};
        template<int codim>
        bool isInvalid(const Entity<codim>& e)
        {
            if(e.partitionType()==OverlapEntity)
                return false;
            return true;
        }
    };    
    
    template<>
    struct PartitionIteratorRule<OverlapFront_Partition>
        : public PartitionIteratorRule<Overlap_Partition>
    {
        // There are no front entities, therefore
        // we fall back to the behaviour of OverlapPartition
    };

    template<>
    struct PartitionIteratorRule<All_Partition>
    {
        enum {fullSet=true, emptySet=false};
        template<int codim>
        bool isInvalid(const Entity<codim>& e)
        {
            return false;
        }
    };

} // end namespace cpgrid
} // end namespace Dune
