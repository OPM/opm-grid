/*
  Copyright 2014 Andreas Lauser

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

#ifndef OPM_INTERSECTION_MAPPER_HPP
#define OPM_INTERSECTION_MAPPER_HPP

#include <dune/grid/common/mcmgmapper.hh>

#include <cstdint>
#include <unordered_map>

namespace Dune {
/*!
 * \brief A class to map an arbitrary intersection to a unique index.
 *
 * This class only assumes that the element indices for the interior and exterior
 * elements of an intersection are unique for each intersection (i.e. that no "contact
 * area" between two elements is represented using more than a single intersection), and
 * that an intersection is seen from both cells. These two assumptions are fullfilled by
 * all known DUNE grid managers.
 *
 * On average, the mapping process is in O(1), but this cannot be guaranteed without
 * cooperation of the grid manager. Mean construction time of an object of this class
 * takes O(m*n), where 'm' is the number of elements in the grid view and 'n' represents
 * the maximum number of intersections per element.
 *
 * Note, that this class is _not_ specific to Dune::CpGrid!
 */
template <class GridView>
class IntersectionMapper
{
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, Dune::MCMGElementLayout> ElementMapper;

public:
    typedef typename GridView::Intersection Intersection;

    IntersectionMapper(const GridView& gridView)
        : gridView_(gridView)
        , elementMap_(gridView_)
    {
        update();
    }

    void update()
    {
        elementMap_.update();

        // reserve enough space in the map for a conforming hexahedron grid. It does not
        // matter too much if this number is an over- or underestimate: In the former
        // case, a bit of memory will be wasted, and in the latter, the code will become
        // slighly slower due to occasional resizes of the unordered_map...
        int numElems = gridView_.size(/*codim=*/0);
        isIdToIndex_.clear();
        isIdToIndex_.reserve(numElems*3 * 3/2);

        numIs_ = 0;

        auto elemIt = gridView_.template begin</*codim=*/0>();
        const auto& elemEndIt = gridView_.template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            auto isIt = gridView_.ibegin(*elemIt);
            const auto& isEndIt = gridView_.iend(*elemIt);
            for (; isIt != isEndIt; ++isIt) {
                std::uint64_t isId = isId_(*isIt);

                if (isIdToIndex_.count(isId) > 0)
                    continue;
                isIdToIndex_[isId] = numIs_;

                ++ numIs_;
            }
        }
    }

    int map(const Intersection& is) const
    { return isIdToIndex_.at(isId_(is)); }

    int size() const
    { return numIs_; }

private:
    std::uint64_t isId_(const Intersection& is) const
    {
        static const int elemIdxShift = 32; // bits

        int insideElemIdx = elementMap_.map(*is.inside());
        int outsideElemIdx =
            is.boundary()
            ? (elementMap_.size() + is.indexInInside())
            : elementMap_.map(*is.outside());

        int elem1Idx = std::min(insideElemIdx, outsideElemIdx);
        std::uint64_t elem2Idx = std::max(insideElemIdx, outsideElemIdx);

        return (elem2Idx<<elemIdxShift) + elem1Idx;
    }

    int numIs_;
    std::unordered_map<std::uint64_t, int> isIdToIndex_;
    const GridView gridView_;
    ElementMapper elementMap_;
};

} // namespace Dune

#endif // OPM_INTERSECTION_MAPPER_HPP
