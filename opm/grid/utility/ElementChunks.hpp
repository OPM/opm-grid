/*
  Copyright 2025 SINTEF Digital

  This file is part of the Open Porous Media project (OPM).

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
#ifndef OPM_ELEMENT_CHUNKS_HEADER
#define OPM_ELEMENT_CHUNKS_HEADER

#include <opm/grid/utility/createThreadIterators.hpp>

#include <cstddef>
#include <iterator>
#include <utility>
#include <vector>

namespace Opm
{

/// Class to simplify creating parallel yet performant threaded loops
/// over a grid.
///
///
/// The basic problem is that iteration over a Dune grid is not
/// random-access, and therefore not amenable to typical OpenMP
/// pragmas. The ThreadedEntityIterator class provided in
/// opm-simulators is a possibility, but it uses locking, and is
/// therefore not very efficient. The iteration over chunks
/// facilitated by this class is significantly faster in many cases.
///
/// Typical loop over grid without this facility:
/// for (const auto& elem : elements(gridview)) {
///     // Do something with elem
/// }
///
/// Typical OpenMP-threaded loop over grid using this class:
/// #pragma omp parallel for
/// for (const auto& chunk : ElementChunks(gridview, num_threads)) {
///     for (const auto& elem : chunk) {
///         // Do something with elem
///     }
/// }
///
/// The ElementChunks object stores a vector of iterators, so if you
/// have several such loops it can be a good idea to create the object
/// once instead of once per loop:
/// ElementChunks chunks(gridview, num_threads)
/// // First loop
/// #pragma omp parallel for
/// for (const auto& chunk : chunks) {
///     for (const auto& elem : chunk) {
///         // Do something with elem
///     }
/// }
/// // ...
/// // Second loop
/// #pragma omp parallel for
/// for (const auto& chunk : chunks) {
///     for (const auto& elem : chunk) {
///         // Do something else with elem
///     }
/// }
template <class GridView, class PartitionSet>
class ElementChunks
{
private:
    using Iter = decltype(std::begin(elements(std::declval<const GridView&>(), PartitionSet())));
    using Storage = std::vector<Iter>;
    using StorageIter = decltype(Storage().cbegin());

public:
    ElementChunks(const GridView& gv,
                  const PartitionSet included_partition,
                  const std::size_t num_chunks)
    {
        grid_chunk_iterators_
            = Opm::createChunkIterators(elements(gv, included_partition), gv.size(0), num_chunks);
    }

    struct Chunk
    {
        Chunk(const Iter& i1, const Iter& i2) : pi_(i1, i2) {}
        auto begin() const { return pi_.first; }
        auto end() const { return pi_.second; }
        std::pair<Iter, Iter> pi_;
    };

    struct ChunkIterator : public StorageIter
    {
        explicit ChunkIterator(const StorageIter& itit) : StorageIter(itit) {}
        Chunk operator*() const
        {
            const StorageIter it = *this;
            return Chunk{*it, *(it+1)};
        }
    };

    auto begin() const
    {
        return ChunkIterator(grid_chunk_iterators_.begin());
    }
    auto end() const
    {
        // Note the decrement, we want a StorageIterator pointing to
        // the last element in grid_chunk_iterators_, not the
        // one-beyond-last-element iterator.
        return ChunkIterator(--grid_chunk_iterators_.end());
    }
    auto size() const
    {
        return grid_chunk_iterators_.size() - 1;
    }
private:
    Storage grid_chunk_iterators_;
};


} // namespace Opm

#endif // OPM_ELEMENT_CHUNKS_HEADER
