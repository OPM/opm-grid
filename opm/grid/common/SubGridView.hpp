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
#include <vector>

namespace Dune
{

template <class GridImp>
class SubGridView;

template <class GridImp>
struct SubGridViewTraits {
    typedef SubGridView<GridImp> GridViewImp;

    /** \brief type of the grid */
    typedef typename std::remove_const<GridImp>::type Grid;

    /** \brief type of the index set */
    typedef typename Grid ::Traits ::LeafIndexSet IndexSet;

    /** \brief type of the intersection */
    typedef typename Grid ::Traits ::LeafIntersection Intersection;

    /** \brief type of the intersection iterator */
    typedef typename Grid ::Traits ::LeafIntersectionIterator IntersectionIterator;

    /** \brief type of the collective communication */
    typedef typename Grid ::Traits ::CollectiveCommunication CollectiveCommunication;

    template <int cd>
    struct Codim {
        typedef typename Grid ::Traits ::template Codim<cd>::template Partition<All_Partition>::LeafIterator BaseIterator;

        typedef typename Grid ::Traits ::template Codim<cd>::Entity Entity;
        typedef typename Grid ::Traits ::template Codim<cd>::EntitySeed EntitySeed;

        typedef typename Grid ::template Codim<cd>::Geometry Geometry;
        typedef typename Grid ::template Codim<cd>::LocalGeometry LocalGeometry;

        /** \brief Define types needed to iterate over entities of a given partition type */
        template <PartitionIteratorType pit>
        struct Partition {
            /** \brief iterator over a given codim and partition type */
            typedef typename Grid ::template Codim<cd>::template Partition<pit>::LeafIterator BaseIterator;
        };
    };

    enum { conforming = Capabilities ::isLeafwiseConforming<Grid>::v };
};


template <class GridImp>
class SubGridView
{
    typedef SubGridView<GridImp> ThisType;

public:
    typedef SubGridViewTraits<GridImp> Traits;

    /** \brief type of the grid */
    typedef typename Traits::Grid Grid;

    /** \brief type of the index set */
    typedef typename Traits ::IndexSet IndexSet;

    /** \brief type of the intersection */
    typedef typename Traits ::Intersection Intersection;

    /** \brief type of the intersection iterator */
    typedef typename Traits ::IntersectionIterator IntersectionIterator;

    /** \brief type of the collective communication */
    typedef typename Traits ::CollectiveCommunication CollectiveCommunication;

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
    {
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
        throw std::logic_error("Not implemented.");
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
};

} // namespace Dune

#endif // OPM_SUBGRIDVIEW_HEADER
