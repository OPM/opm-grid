//===========================================================================
//
// File: Indexsets.hpp
//
// Created: Fri May 29 23:30:01 2009
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
Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
Copyright 2009, 2010, 2022 Equinor ASA.

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

#ifndef OPM_INDEXSETS_HEADER
#define OPM_INDEXSETS_HEADER

#include <dune/geometry/type.hh>
#include <opm/common/ErrorMacros.hpp>
#include "GlobalIdMapping.hpp"
#include "Intersection.hpp"

#include <unordered_map>
namespace Dune
{
    namespace cpgrid
    {
    // forward declaration
    //class CpGridData;
        /// @brief
        /// @todo Doc me!
        /// @tparam
        class IndexSet
        {
            //friend class Dune::cpgrid::CpGridData;
        public:
            /// @brief
            /// @todo Doc me!
            typedef int IndexType;

            /** \brief Export the type of the entity used as parameter in the index(...) method */
            template <int cc>
            struct Codim
            {
              typedef cpgrid::Entity< cc > Entity;
            };

            /// @brief
            /// @todo Doc me!
            typedef std::vector<GeometryType> Types;

            /// @brief
            /// @todo Doc me!
            /// @param
            IndexSet(){}
            
            IndexSet(std::size_t numCells, std::size_t numPoints)
            {
                geom_types_[0].emplace_back(Dune::GeometryTypes::cube(3));
                geom_types_[3].emplace_back(Dune::GeometryTypes::cube(0));
                size_codim_map_[0] =  numCells; 
                size_codim_map_[3] = numPoints; 
            }

            /// \brief Destructor.
            ~IndexSet()
            {}

            /// @brief
            /// @todo Doc me!
            /// @param
            /// @return
            const Types& geomTypes(int codim) const
            {
                return geom_types_[codim];
            }

            /// @brief
            /// @todo Doc me!
            /// @param
            /// @return
            const Types& types(int codim) const
            {
                return geom_types_[codim];
            }

            /// @brief
            /// @todo Doc me!
            /// @param
            /// @return
            int size(GeometryType type) const
            {
                if (type.isCube()) {
                    return size(3 - type.dim());  // return grid_.size(type);
                } else {
                    return 0;
                }
            }


            /// @brief
            /// @todo Doc me!
            /// @param
            /// @return
            int size(int codim) const
            {
                return size_codim_map_[codim]; //grid_.size(codim)
            }


            /// @brief
            /// @todo Doc me!
            /// @tparam
            /// @return
            /// @param
            template<int cd>
            IndexType index(const cpgrid::Entity<cd>& e) const
            {
                return e.index();
            }

            /// @brief
            /// @todo Doc me!
            /// @tparam
            /// @return
            /// @param
            template<class EntityType>
            IndexType index(const EntityType& e) const
            {
                return e.index();
            }

            /// @brief
            /// @todo Doc me!
            /// @tparam
            /// @return
            /// @param
            template <int cc>
            IndexType subIndex(const cpgrid::Entity<0>& e, int i) const
            {
                return index(e.template subEntity<cc>(i));
            }

            /// @brief
            /// @todo Doc me!
            /// @tparam
            /// @return
            /// @param
            IndexType subIndex(const cpgrid::Entity<0>& e, int i, unsigned int cc) const;


            template<int codim>
	    IndexType subIndex(const cpgrid::Entity<codim>& /* e */, int /* i */, unsigned int /* cc */) const
	    {
	      DUNE_THROW(NotImplemented, "subIndex not implemented for codim"
			 << codim << "entities.");
	    }
            /// @brief
            /// @todo Doc me!
            /// @tparam
            /// @return
            /// @param
            template <class EntityType>
            bool contains(const EntityType& e) const
            {
                // return index(e) >= 0 && index(e) < grid_.size(EntityType::codimension); //EntityType::codimension == 0;
                return index(e) >= 0 && index(e) < this->size(EntityType::codimension);
            }

        private:
            // const CpGridData& grid_;
            Types geom_types_[4];
            std::array<int,4> size_codim_map_{0,0,0,0};
        };


        class IdSet
        {
            friend class ReversePointGlobalIdSet;
        public:
            typedef int IdType;

            IdSet(const CpGridData& grid)
                : grid_(grid)
            {
            }

            template<int cc>
            IdType id(const cpgrid::Entity<cc>& e) const
            {
                return computeId(e);
            }

            template<class EntityType>
            IdType id(const EntityType& e) const
            {
                return computeId(e);
            }

            /// return id of intersection (here face number)
            IdType id( const cpgrid::Intersection& intersection ) const
            {
                return intersection.id();
            }

            template<int cc>
            IdType subId(const cpgrid::Entity<0>& e, int i) const
            {
                return id(e.template subEntity<cc>(i));
            }

            IdType subId(const cpgrid::Entity<0>& e, int i, int cc) const;
        private:
            template<class EntityType>
            IdType computeId(const EntityType& e) const
            {
                IdType myId = 0;
                for( int c=0; c<EntityType::codimension; ++c )
                    myId += grid_.indexSet().size( c );
                return  myId + e.index();
            }
            const CpGridData& grid_;
        };


        class LevelGlobalIdSet : public GlobalIdMapping
        {
            friend class CpGridData;
            friend class ReversePointGlobalIdSet;
        public:
            typedef int IdType;

            void swap(std::vector<int>& cellMapping,
                      std::vector<int>& faceMapping,
                      std::vector<int>& pointMapping)
            {
                idSet_=nullptr;
                GlobalIdMapping::swap(cellMapping,
                                      faceMapping,
                                      pointMapping);
            }
            LevelGlobalIdSet(std::shared_ptr<const IdSet> ids, const CpGridData* view)  
                : idSet_(ids), view_(view)
            {}
            LevelGlobalIdSet()
                : idSet_(), view_()
            {}
            template<int codim>
            IdType id(const Entity<codim>& e) const
            {
                assert(view_ == e.pgrid_);
                return id(static_cast<const EntityRep<codim>&>(e));
            }
            template<int codim>
            IdType id(const EntityRep<codim>& e) const
            {
                if(idSet_)
                    return idSet_->id(e);
                else
                    return this->template getMapping<codim>()[e.index()];
            }

            template<int cc>
            IdType subId(const cpgrid::Entity<0>& e, int i) const
            {
                assert(view_ == e.pgrid_);
                return id(e.template subEntity<cc>(i));
            }

            IdType subId(const cpgrid::Entity<0>& e, int i, int cc) const;
        private:
            std::shared_ptr<const IdSet> idSet_;
            const CpGridData* view_;
        };

    /*!
     * \brief The global id set for Dune.
     *
     * You can pass it any entity of either the loadbalanced
     * or global grid that is stored on this process.
     */
    class GlobalIdSet
    {
    public:
        /// \brief The type of the id.
        using IdType = typename LevelGlobalIdSet::IdType;

        GlobalIdSet(const CpGridData& view);

        template<int codim>
        IdType id(const Entity<codim>& e) const
        {
            return levelIdSet(e.pgrid_).id(e);
        }

        template<int cc>
        IdType subId(const cpgrid::Entity<0>& e, int i) const
        {
            return levelIdSet(e.pgrid_).template subId<cc>(e, i);
        }

        IdType subId(const cpgrid::Entity<0>& e, int i, int cc) const;

        void insertIdSet(const CpGridData& view);
    private:
        /// \brief Get the correct id set of a level (global or distributed)
        const LevelGlobalIdSet& levelIdSet(const CpGridData* const data) const
        {
            auto candidate = idSets_.find(data);
            assert(candidate != idSets_.end());
            return *candidate->second;
        }
        /// \brief map of views onto idesets if the view.
        std::map<const CpGridData* const, std::shared_ptr<const LevelGlobalIdSet>> idSets_;
    };

    class ReversePointGlobalIdSet
    {
    public:
        ReversePointGlobalIdSet(const LevelGlobalIdSet& idSet)
        {
            if(idSet.idSet_)
            {
                grid_ = &(idSet.idSet_->grid_);
            }
            else
            {
                mapping_.reset(new std::unordered_map<int,int>);
                int localId = 0;
                for (const  auto& globalId: idSet.template getMapping<3>())
                    (*mapping_)[globalId] = localId++;
            }
        }
        int operator[](int i) const
        {
            if (mapping_)
            {
                return(*mapping_)[i];
            }
            else if (grid_)
            {
                return i - grid_->size(0) - grid_->size(1) - grid_->size(2);
            }

            OPM_THROW(std::runtime_error, "No grid or mapping. Should not be here!");
        }
        void release()
        {
            mapping_.reset(nullptr);
        }
    private:
        std::unique_ptr<std::unordered_map<int,int> > mapping_;
        const CpGridData* grid_ = nullptr;
    };

    } // namespace cpgrid
} // namespace Dune

#endif // OPM_INDEXSETS_HEADER
