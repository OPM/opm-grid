//===========================================================================
//
// File: CpGrid.hpp
//
// Created: Fri May 29 20:26:36 2009
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
Copyright 2009 SINTEF ICT, Applied Mathematics.
Copyright 2009 Statoil ASA.

This file is part of The Open Reservoir Simulator Project (OpenRS).

OpenRS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

OpenRS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with OpenRS.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPENRS_CPGRID_HEADER
#define OPENRS_CPGRID_HEADER


#include <string>
#include <map>

#include <dune/common/collectivecommunication.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/grid.hh>
#include <dune/common/timer.hh>

#include <dune/common/param/ParameterGroup.hpp>

#include "cpgrid/Entity.hpp"
#include "cpgrid/Geometry.hpp"
#include "cpgrid/Iterators.hpp"
#include "cpgrid/Indexsets.hpp"
#include "cpgrid/GridView.hpp"
#include "cpgrid/DefaultGeometryPolicy.hpp"

#include "common/SparseTable.hpp"


namespace Dune
{

    class CpGrid;

    ////////////////////////////////////////////////////////////////////////
    //
    //   CpGridTraits
    //
    ////////////////////////////////////////////////////////////////////////

    struct CpGridTraits
    {
	/// \brief The type that implementing the grid.
	typedef CpGrid Grid;

	/// \brief The type of the intersection at the leafs of the grid.
	typedef cpgrid::Intersection<CpGrid> LeafIntersection;
	/// \brief The type of the intersection at the levels of the grid.
	typedef cpgrid::Intersection<CpGrid> LevelIntersection;
	/// \brief The type of the intersection iterator at the leafs of the grid.
	typedef cpgrid::IntersectionIterator<CpGrid> LeafIntersectionIterator;
	/// \brief The type of the intersection iterator at the levels of the grid.
	typedef cpgrid::IntersectionIterator<CpGrid> LevelIntersectionIterator;

	/// \brief The type of the  hierarchic iterator.
	typedef cpgrid::HierarchicIterator HierarchicIterator;

	/// \brief Traits associated with a specific codim. 
	/// \tparam cd The codimension.
	template <int cd>
	struct Codim
	{
	    /// \brief The type of the geometry associated with the entity.
	    /// IMPORTANT: Codim<codim>::Geometry == Geometry<dim-codim,dimw>
	    typedef cpgrid::Geometry<3-cd, 3> Geometry;
	    /// \brief The type of the local geometry associated with the entity.
	    typedef cpgrid::Geometry<3-cd, 3> LocalGeometry;
	    /// \brief The type of the entity.
	    typedef cpgrid::Entity<cd, CpGrid> Entity;

	    /// \brief The type of the iterator over all level entities of this codim.
	    typedef cpgrid::Iterator<cd, All_Partition, CpGrid> LevelIterator;

	    /// \brief The type of the iterator over all leaf entities of this codim.
	    typedef cpgrid::Iterator<cd, All_Partition, CpGrid> LeafIterator;

	    /// \brief The type of the entity pointer for entities of this codim.
	    typedef cpgrid::EntityPointer<cd, CpGrid> EntityPointer;

	    /// \brief Traits associated with a specific grid partition type.
	    /// \tparam pitype The type of the grid partition.
	    template <PartitionIteratorType pitype>
	    struct Partition
	    {
		/// \brief The type of the iterator over the level entities of this codim on this partition.
		typedef cpgrid::Iterator<cd, pitype, CpGrid> LevelIterator;
		/// \brief The type of the iterator over the leaf entities of this codim on this partition.
		typedef cpgrid::Iterator<cd, pitype, CpGrid> LeafIterator;
	    };
	};

	/// \brief Traits associated with a specific grid partition type.
	/// \tparam pitype The type of the grid partition.
	template <PartitionIteratorType pitype>
	struct Partition
	{
	    /// \brief The type of the level grid view associated with this partition type.
// 	    typedef cpgrid::GridView<pitype> LevelGridView;
	    typedef Dune::GridView<DefaultLevelGridViewTraits<const CpGrid, pitype> > LevelGridView;

	    /// \brief The type of the leaf grid view associated with this partition type.
// 	    typedef cpgrid::GridView<pitype> LeafGridView;
	    typedef Dune::GridView<DefaultLeafGridViewTraits<const CpGrid, pitype> > LeafGridView;
	};

	/// \brief The type of the level index set.
	typedef cpgrid::IndexSet<CpGrid> LevelIndexSet;
	/// \brief The type of the leaf index set.
	typedef cpgrid::IndexSet<CpGrid> LeafIndexSet;
	/// \brief The type of the global id set.
	typedef cpgrid::IdSet<CpGrid> GlobalIdSet;
	/// \brief The type of the local id set.
	typedef cpgrid::IdSet<CpGrid> LocalIdSet;

	/// \brief The type of the collective communication.
	typedef Dune::CollectiveCommunication<CpGrid> CollectiveCommunication;
    };



    ////////////////////////////////////////////////////////////////////////
    //
    //   CpGridFamily
    //
    ////////////////////////////////////////////////////////////////////////

    struct CpGridFamily
    {
	typedef CpGridTraits Traits;
    };

    ////////////////////////////////////////////////////////////////////////
    //
    //   CpGrid
    //
    ////////////////////////////////////////////////////////////////////////

    /// \brief [<em> provides \ref Dune::Grid </em>]
    class CpGrid
	: public GridDefaultImplementation<3, 3, double, CpGridFamily >
    {
    public:
	/*
	// --- Enums ---
	enum { dimension = 3 };
	enum { dimensionworld = 3 };

	// --- Types ---

        /// The type used to store coordinates
        typedef double ctype;

	/// Traits for CpGrid
	/typedef CpGridTraits Traits;

	typedef cpgrid::GridView<All_Partition> LeafGridView;
	typedef cpgrid::GridView<All_Partition> LevelGridView;
	*/

	// --- Methods ---

	/// Default constructor
	CpGrid()
	    : index_set_(*this), id_set_(*this)
	{
	}


        /// Initialize the grid.
	void init(const parameter::ParameterGroup& param);

	/// Read the sintef legacy grid format ('topogeom').
	void readSintefLegacyFormat(const parameter::ParameterGroup& param);


        /// return grid name
        std::string name() const
        {
            return "CpGrid";
        }
    
        
        /// Return maximum level defined in this grid. Levels are numbered
        /// 0 ... maxlevel with 0 the coarsest level.
        int maxLevel() const
	{
            return 0;
        }
        


        /// Iterator to first entity of given codim on level
        template<int codim>
        typename Traits::template Codim<codim>::LevelIterator lbegin (int level) const
	{
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
            return cpgrid::Iterator<codim, All_Partition, CpGrid >(*this, 0);
        }
    
        
        /// one past the end on this level
        template<int codim>
        typename Traits::template Codim<codim>::LevelIterator lend (int level) const
	{
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
            return cpgrid::Iterator<codim,All_Partition, CpGrid >(*this, size(codim));
        }
        
        
        /// Iterator to first entity of given codim on level
        template<int codim, PartitionIteratorType PiType>
        typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lbegin (int level) const
	{
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
            return cpgrid::Iterator<codim,PiType, CpGrid >(*this, 0);
        }
        

        /// one past the end on this level
        template<int codim, PartitionIteratorType PiType>
        typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lend (int level) const
	{
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
            return cpgrid::Iterator<codim,PiType, CpGrid >(*this, size(codim));
        }


        /// Iterator to first leaf entity of given codim
        template<int codim>
        typename Traits::template Codim<codim>::LeafIterator leafbegin() const
	{
            return cpgrid::Iterator<codim,All_Partition, CpGrid >(*this, 0);
        }
        
    
        /// one past the end of the sequence of leaf entities
        template<int codim>
        typename Traits::template Codim<codim>::LeafIterator leafend() const
	{
            return cpgrid::Iterator<codim,All_Partition, CpGrid >(*this, size(codim));
        }
        
    
        /// Iterator to first leaf entity of given codim
        template<int codim, PartitionIteratorType PiType>
        typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafbegin() const
	{
            return cpgrid::Iterator<codim,PiType, CpGrid >(*this, 0);
        }
        
        
        /// one past the end of the sequence of leaf entities
        template<int codim, PartitionIteratorType PiType>
        typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafend() const
	{
            return cpgrid::Iterator<codim,PiType, CpGrid >(*this, size(codim));
        }
        

        /// \brief Number of grid entities per level and codim
        int size (int level, int codim) const
	{
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
            return size(codim);        
        }
        
        
        /// number of leaf entities per codim in this process
        int size (int codim) const
	{
	    switch (codim) {
	    case 0: return cell_to_face_.size();
	    case 1: return face_to_cell_.size();
            default: THROW("This grid has no entities of codim " << codim);
	    }
	    return -1;
        }
        
        
        /// number of entities per level, codim and geometry type in this process
        int size (int level, GeometryType type) const
	{
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
            return size(type);
        }
        
            
        /// number of leaf entities per codim and geometry type in this process
        int size (GeometryType) const
        {
	    THROW("size(level, type) not yet implemented");
            return -1;
        }
        
        
        /// \brief Access to the GlobalIdSet 
        const Traits::GlobalIdSet& globalIdSet() const
	{
            return id_set_;
        }
        
        
        /// \brief Access to the LocalIdSet 
        const Traits::LocalIdSet& localIdSet() const
	{
            return id_set_;
        }
        
        
        /// \brief Access to the LevelIndexSets 
        const Traits::LevelIndexSet& levelIndexSet(int level) const
        {
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
            return index_set_;
        }
        
        
        /// \brief Access to the LeafIndexSet 
        const Traits::LeafIndexSet& leafIndexSet() const
        {
            return index_set_;
        }


	/// \brief The leaf view
	LeafGridView leafView() const
	{
	    return LeafGridView(*this);
	}

	/*  No refinement implemented

        /// global refinement
        void globalRefine (int refCount)
        {
            hostgrid_->globalRefine(refCount);
        }
        
        /// \brief Mark entity for refinement
	///
	/// This only works for entities of codim 0.
	/// The parameter is currently ignored
	///
	/// \return <ul>
	/// <li> true, if marking was succesfull </li>
	/// <li> false, if marking was not possible </li>
	/// </ul>
	 
        bool mark(int refCount, const typename Traits::template Codim<0>::EntityPointer & e)
        {
            return hostgrid_->mark(refCount, getHostEntity<0>(*e));
        }
        
        /// \brief Return refinement mark for entity
	///
	/// \return refinement mark (1,0,-1)
	 
        int getMark(const typename Traits::template Codim<0>::EntityPointer & e) const
        {
            return hostgrid_->getMark(getHostEntity<0>(*e));
        }

        /// \todo Please doc me !
        bool preAdapt() {
            return hostgrid_->preAdapt();
        }
        
        
        /// Triggers the grid refinement process
        bool adapt()
        {
            return hostgrid_->adapt();
        }

        /// \brief Clean up refinement markers 
        void postAdapt() {
            return hostgrid_->postAdapt();
        }

	end of refinement section */


	/* No parallelism implemented

        /// \brief Size of the overlap on the leaf level 
        unsigned int overlapSize(int codim) const {
            return hostgrid_->overlapSize(codim);
        }
        
        
        /// \brief Size of the ghost cell layer on the leaf level 
        unsigned int ghostSize(int codim) const {
            return hostgrid_->ghostSize(codim);
        }
        
        
        /// \brief Size of the overlap on a given level 
        unsigned int overlapSize(int level, int codim) const {
            return hostgrid_->overlapSize(level,codim);
        }
        
        
        /// \brief Size of the ghost cell layer on a given level 
        unsigned int ghostSize(int level, int codim) const {
            return hostgrid_->ghostSize(level,codim);
        }
        
            
        /// \brief Distributes this grid over the available nodes in a distributed machine
	///
	/// \param minlevel The coarsest grid level that gets distributed
	/// \param maxlevel does currently get ignored
        void loadBalance(int strategy, int minlevel, int depth, int maxlevel, int minelement){
            DUNE_THROW(NotImplemented, "CpGrid::loadBalance()");
        }
        
        /// \brief The communication interface
	///  @param T: array class holding data associated with the entities
	///  @param P: type used to gather/scatter data in and out of the message buffer
	///  @param codim: communicate entites of given codim
	///  @param if: one of the predifined interface types, throws error if it is not implemented
	///  @param level: communicate for entities on the given level
	///
	///  Implements a generic communication function sending an object of type P for each entity
	///  in the intersection of two processors. P has two methods gather and scatter that implement
	///  the protocol. Therefore P is called the "protocol class".
        template<class T, template<class> class P, int codim>
        void communicate (T& t, InterfaceType iftype, CommunicationDirection dir, int level);
        
        /// The new communication interface.
	/// communicate objects for all codims on a given level        
        template<class DataHandle>
        void communicate (DataHandle& data, InterfaceType iftype, CommunicationDirection dir, int level) const
        {}
        
        template<class DataHandle>
        void communicate (DataHandle& data, InterfaceType iftype, CommunicationDirection dir) const
        {}        

	end of parallel section */
        
        /// dummy collective communication 
        const CollectiveCommunication& comm () const
        {
            return ccobj_;
        }

    private:

	template <int cd, class GridType>
	friend class cpgrid::Entity;

	template <class GridType>
	friend class cpgrid::Intersection;

	template <class GridType>
	friend class cpgrid::IndexSet;

        /// \todo Please doc me !
        CollectiveCommunication ccobj_;

	cpgrid::IndexSet<CpGrid> index_set_;

	cpgrid::IdSet<CpGrid> id_set_;

	cpgrid::GridView<All_Partition> view_;


	// Representing the topology
	cpgrid::OrientedEntityTable<0, 1> cell_to_face_;
	cpgrid::OrientedEntityTable<1, 0> face_to_cell_;

	typedef cpgrid::DefaultGeometryPolicy Geom;

	// Representing geometry
	Geom geometry_;
	template <int codim>
	const cpgrid::EntityVariable< cpgrid::Geometry<3 - codim, 3>, codim>& geomVector() const
	{
	    return geometry_.geomVector<codim>();
	}
	typedef FieldVector<ctype, 3> PointType;
	cpgrid::SignedEntityVariable<PointType, 1> face_normals_;

    }; // end Class CpGrid



    namespace Capabilities
    {
	/// \todo Please doc me !
	template <>
	struct hasEntity<CpGrid, 0>
	{
	    static const bool v = true;
	};
	/// \todo Please doc me !
	template <>
	struct hasEntity<CpGrid, 1>
	{
	    static const bool v = true;
	};
    
	/// \todo Please doc me !
	template <>
	struct isParallel<CpGrid>
	{
	    static const bool v = false;
	};
    
	/// \todo Please doc me !
	template <>
	struct isLevelwiseConforming<CpGrid>
	{
	    static const bool v = true;
	};

	/// \todo Please doc me !
// 	template <>
// 	struct isLeafwiseConforming<CpGrid>
// 	{
// 	    static const bool v = true;
// 	};

	/// \todo Please doc me !
	template <>
	struct hasHangingNodes<CpGrid>
	{
	    static const bool v = true;
	};

	/// \todo Please doc me !
	template <>
	struct hasBackupRestoreFacilities<CpGrid>
	{
	    static const bool v = false;
	};

	/// \todo Please doc me !
	template <>
	struct IsUnstructured<CpGrid>
	{
	    static const bool v = true;
	};

    }

} // namespace Dune


#endif // OPENRS_CPGRID_HEADER
