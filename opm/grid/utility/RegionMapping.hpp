/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_REGIONMAPPING_HEADER_INCLUDED
#define OPM_REGIONMAPPING_HEADER_INCLUDED

#include <opm/grid/utility/IteratorRange.hpp>

#include <unordered_map>
#include <vector>

#if HAVE_MPI
#include <mpi.h>
#endif

namespace Opm
{

    /**
     * Forward and reverse mappings between cells and
     * regions/partitions (e.g., the ECLIPSE-style 'SATNUM',
     * 'PVTNUM', or 'EQUILNUM' arrays).
     *
     * \tparam Region Type of a forward region mapping.  Expected
     *                to provide indexed access through
     *                operator[]() as well as inner types
     *                'value_type', 'size_type', and
     *                'const_iterator'.
     */
    template < class Region = std::vector<int> >
    class RegionMapping {
    public:
        /**
         * Constructor.
         *
         * \param[in] reg Forward region mapping, restricted to
         *                active cells only.
         * \param[in] comm Global communicator to base region communicator on.
         */
        explicit
        RegionMapping(const Region& reg
#if HAVE_MPI
                      , MPI_Comm comm = MPI_COMM_WORLD
#endif
                     )
            : reg_(reg)
        {
#if HAVE_MPI
            rev_.init(reg_, comm);
#else
            rev_.init(reg_);
#endif
        }

        ~RegionMapping()
        {
#if HAVE_MPI
            if (rev_.comm_ != MPI_COMM_NULL)
                MPI_Comm_free(&rev_.comm_);
#endif
        }

        /**
         * Type of forward (cell-to-region) mapping result.
         * Expected to be an integer.
         */
        typedef typename Region::value_type RegionId;

        /**
         * Type of reverse (region-to-cell) mapping (element)
         * result.
         */
        typedef typename Region::size_type CellId;

        /**
         * Type of reverse region-to-cell range bounds and
         * iterators.
         */
        typedef typename std::vector<CellId>::const_iterator CellIter;

        using Range = iterator_range<CellIter>;

        /**
         * Compute region number of given active cell.
         *
         * \param[in] c Active cell
         * \return Region to which @c c belongs.
         */
        RegionId
        region(const CellId c) const { return reg_[c]; }

        const std::vector<RegionId>&
        activeRegions() const
        {
            return rev_.active;
        }

        /**
         * Extract active cells in particular region.
         *
         * \param[in] r Region number
         *
         * \return Range of active cells in region @c r.  Empty if @c r is
         * not an active region.
         */
        Range
        cells(const RegionId r) const {
            const auto id = rev_.binid.find(r);

            if (id == rev_.binid.end()) {
                // Region 'r' not an active region.  Return empty.
                return Range(rev_.c.end(), rev_.c.end());
            }

            const auto i = id->second;

            return Range(rev_.c.begin() + rev_.p[i + 0],
                         rev_.c.begin() + rev_.p[i + 1]);
        }

#if HAVE_MPI
        MPI_Comm comm() const { return rev_.comm_; }
#endif

    private:
        /**
         * Copy of forward region mapping (cell-to-region).
         */
        Region reg_;

        /**
         * Reverse mapping (region-to-cell).
         */
        struct {
            typedef typename std::vector<CellId>::size_type Pos;

            std::unordered_map<RegionId, Pos> binid;
            std::vector<RegionId>             active;

            std::vector<Pos>    p;   /**< Region start pointers */
            std::vector<CellId> c;   /**< Region cells */

#if HAVE_MPI
            MPI_Comm comm_;
#endif

            /**
             * Compute reverse mapping.  Standard linear insertion
             * sort algorithm.
             */
#if HAVE_MPI
            void
            init(const Region& reg, MPI_Comm comm)
#else
            void
            init(const Region& reg)
#endif
            {
                binid.clear();
                for (const auto& r : reg) {
                    ++binid[r];
                }

                p     .clear();  p.emplace_back(0);
                active.clear();
                {
                    Pos n = 0;
                    for (auto& id : binid) {
                        active.push_back(id.first);
                        p     .push_back(id.second);

                        id.second = n++;
                    }
                }
#if HAVE_MPI
                MPI_Comm_split(comm, reg.empty() ? MPI_UNDEFINED : 1,
                               0, &comm_);
#endif

                for (decltype(p.size()) i = 1, sz = p.size(); i < sz; ++i) {
                    p[0] += p[i];
                    p[i]  = p[0] - p[i];
                }

                assert (p[0] == static_cast<Pos>(reg.size()));

                c.resize(reg.size());
                {
                    CellId i = 0;
                    for (const auto& r : reg) {
                        auto& pos  = p[ binid[r] + 1 ];
                        c[ pos++ ] = i++;
                    }
                }

                p[0] = 0;
            }
        } rev_; /**< Reverse mapping instance */
    };

} // namespace Opm

#endif // OPM_REGIONMAPPING_HEADER_INCLUDED
