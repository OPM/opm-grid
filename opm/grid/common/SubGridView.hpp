// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=2 sts=2:

// Note: this file is based on defaultgridview.hh from Dune, and is therefore
// licensed under the Dune license (GPLv2 + runtime exception),
// see https://dune-project.org/about/license/
// rather than the OPM license (GPLv3+)
// Copyright 2021 Dune contributors.
// Copyright 2021 SINTEF Digital, Mathematics and Cybernetics.

#ifndef OPM_SUBGRIDVIEW_HEADER
#define OPM_SUBGRIDVIEW_HEADER

#include <dune/common/exceptions.hh>
#include <dune/common/typetraits.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/gridview.hh>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <map>
#include <set>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace Dune
{

template <class GridImp>
class SubGridView;

template <class GridImp>
struct SubGridViewTraits {
    using GridViewImp = SubGridView<GridImp>;

    /** \brief type of the grid */
    using Grid = typename std::remove_const<GridImp>::type;

    /** \brief type of the index set */
    using IndexSet = typename Grid ::Traits ::LeafIndexSet;

    /** \brief type of the intersection */
    using Intersection = typename Grid ::Traits ::LeafIntersection;

    /** \brief type of the intersection iterator */
    using IntersectionIterator = typename Grid ::Traits ::LeafIntersectionIterator;

    /** \brief type of the collective communication */
    using CollectiveCommunication = typename Grid ::Traits ::CollectiveCommunication;

    template <int cd>
    struct Codim {
        using BaseIterator = typename Grid ::Traits ::template Codim<cd>::template Partition<All_Partition>::LeafIterator;

        using Entity = typename Grid ::Traits ::template Codim<cd>::Entity;
        using EntitySeed = typename Grid ::Traits ::template Codim<cd>::EntitySeed;

        using Geometry = typename Grid ::template Codim<cd>::Geometry;
        using LocalGeometry = typename Grid ::template Codim<cd>::LocalGeometry;

        /** \brief Define types needed to iterate over entities of a given partition type */
        template <PartitionIteratorType pit>
        struct Partition {
            /** \brief iterator over a given codim and partition type */
            using BaseIterator = typename Grid ::template Codim<cd>::template Partition<pit>::LeafIterator;
        };
    };

    enum { conforming = Capabilities ::isLeafwiseConforming<Grid>::v };
};


template <class GridImp>
class SubGridView
{
    using ThisType = SubGridView<GridImp>;

public:
    using Traits = SubGridViewTraits<GridImp>;

    /** \brief type of the grid */
    using Grid = typename Traits::Grid;

    /** \brief type of the index set */
    using IndexSet = typename Traits ::IndexSet;

    /** \brief type of the intersection */
    using Intersection = typename Traits ::Intersection;

    /** \brief type of the intersection iterator */
    using IntersectionIterator = typename Traits ::IntersectionIterator;

    /** \brief type of the collective communication */
    using CollectiveCommunication = typename Traits ::CollectiveCommunication;

    /** \brief Codim Structure */
    template <int cd>
    struct Codim : public Traits ::template Codim<cd> {
        using Entity = typename Traits::template Codim<cd>::Entity;
        class Iterator
        {
        public:
            Iterator(const SubGridView& view, std::size_t index)
                : view_(&view)
                , index_(index)
            {
            }
            const Entity& operator*() const
            {
                entity_ = this->view_->get(index_);
                return entity_;
            }
            const Entity* operator->() const
            {
                entity_ = this->view_->get(index_);
                return &entity_;
            }
            Iterator operator++()
            {
                ++index_;
                return *this;
            }
            Iterator operator++(int)
            {
                Iterator copy(*this);
                ++index_;
                return copy;
            }
            bool operator==(const Iterator& other) const
            {
                assert(view_ == other.view_);
                return index_ == other.index_;
            }
            bool operator!=(const Iterator& other) const
            {
                assert(view_ == other.view_);
                return index_ != other.index_;
            }
        private:
            const SubGridView* view_;
            std::size_t index_;
            mutable Entity entity_; // This may be low-performing for grids with large Entity objects.
        };
    };

    enum { conforming = Traits::conforming };
    enum { dimension = GridImp::dimension };

public:

    /// Construct a view of the codim 0 entities that can be constructed from the seeds input.
    ///
    /// The seeds input is moved from and will be in a valid but indeterminate state after the call.
    SubGridView(const Grid& grid,
                std::vector<typename Codim<0>::Entity::EntitySeed>&& seeds)
        : grid_(&grid)
        , subset_(std::move(seeds))
        , num_owned_(subset_.size())
    {
        // Add neighbouring not-owned entities to subset_
        using Seed = typename Codim<0>::Entity::EntitySeed;
        std::set<int> owned;
        std::map<int, Seed> neighbors;
        const auto& iset = grid_->leafIndexSet();
        const auto& leaf_view = grid_->leafGridView();
        for (const auto& seed : subset_) {
            // Add this entity to the set of owned indices.
            const auto& entity = grid_->entity(seed);
            owned.insert(iset.index(entity));
            // Iterating over all intersections, ...
            const auto end = leaf_view.iend(entity);
            for (auto it = leaf_view.ibegin(entity); it != end; ++it) {
                if (it->boundary()) {
                    continue;
                }
                if (it->neighbor()) {
                    const auto outside_entity = it->outside();
                    // ...for all neighbour entities, add to neighbors.
                    const int index = iset.index(outside_entity);
                    const Seed outside_seed = outside_entity.seed();
                    auto elem = std::make_pair(index, outside_seed);
                    neighbors.insert(elem/*{index, seed}*/);
                }
            }
        }
        // Now that owned is complete, we can eliminate any owned entries.
        std::map<int, Seed> unowned_neighbors;
        for (const auto& nb : neighbors) {
            if (owned.count(nb.first) == 0) {
                unowned_neighbors.insert(nb);
            }
        }
        subset_.resize(subset_.size() + unowned_neighbors.size());
        int count = num_owned_;
        for (const auto& [index, seed] : unowned_neighbors) {
            subset_[count] = seed;
            ++count;
        }
    }

    /** \brief obtain a const reference to the underlying hierarchic grid */
    const Grid& grid() const
    {
        assert(grid_);
        return *grid_;
    }

    /** \brief obtain the index set */
    const IndexSet& indexSet() const
    {
        throw std::logic_error("SubGridView::indexSet() not implemented.");
        return grid().leafIndexSet();
    }

    /** \brief obtain number of entities in a given codimension */
    int size(int codim) const
    {
        if (codim == 0) {
            return subset_.size();
        } else {
            throw std::logic_error("Not implemented.");
        }
    }

    /** \brief obtain number of entities with a given geometry type */
    int size(const GeometryType& type) const
    {
        throw std::logic_error("Not implemented.");
        return grid().size(type);
    }

    /** \brief obtain begin iterator for this view */
    template <int cd>
    typename Codim<cd>::Iterator begin() const
    {
        using Iterator = typename Codim<cd>::Iterator;
        return Iterator(*this, 0);
    }

    /** \brief obtain begin iterator for this view */
    // template <int cd, PartitionIteratorType pit>
    // typename Codim<cd>::template Partition<pit>::Iterator begin() const
    // {
    //     return grid().template leafbegin<cd, pit>();
    // }

    /** \brief obtain end iterator for this view */
    template <int cd>
    typename Codim<cd>::Iterator end() const
    {
        using Iterator = typename Codim<cd>::Iterator;
        return Iterator(*this, subset_.size());
    }

    /** \brief obtain end iterator for this view */
    // template <int cd, PartitionIteratorType pit>
    // typename Codim<cd>::template Partition<pit>::Iterator end() const
    // {
    //     return grid().template leafend<cd, pit>();
    // }

    /** \brief obtain begin intersection iterator with respect to this view */
    IntersectionIterator ibegin(const typename Codim<0>::Entity& entity) const
    {
        return entity.impl().ileafbegin();
    }

    /** \brief obtain end intersection iterator with respect to this view */
    IntersectionIterator iend(const typename Codim<0>::Entity& entity) const
    {
        return entity.impl().ileafend();
    }

    /** \brief obtain collective communication object */
    const CollectiveCommunication& comm() const
    {
        return grid().comm();
    }

    /** \brief Return size of the overlap region for a given codim on the grid view.  */
    int overlapSize(int codim) const
    {
        throw std::logic_error("Not implemented.");
        return grid().overlapSize(codim);
    }

    /** \brief Return size of the ghost region for a given codim on the grid view.  */
    int ghostSize(int codim) const
    {
        throw std::logic_error("Not implemented.");
        return grid().ghostSize(codim);
    }

    /** communicate data on this view */
    template <class DataHandleImp, class DataType>
    void
    communicate(CommDataHandleIF<DataHandleImp, DataType>& data, InterfaceType iftype, CommunicationDirection dir) const
    {
        return grid().communicate(data, iftype, dir);
    }

private:
    typename Codim<0>::Entity get(std::size_t ii) const
    {
        return grid_->entity(subset_[ii]);
    }
    const Grid* grid_;
    std::vector<typename Codim<0>::Entity::EntitySeed> subset_;
    const int num_owned_;
};

} // namespace Dune

#endif // OPM_SUBGRIDVIEW_HEADER
