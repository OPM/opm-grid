/*
  Copyright 2017 SINTEF ICT, Applied Mathematics.
  Copyright 2017 Statoil ASA.

  This file is part of the Open Porous Media Project (OPM).

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

#ifndef OPM_REPAIRZCORN_HEADER_INCLUDED
#define OPM_REPAIRZCORN_HEADER_INCLUDED

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

/// \file
///
/// Facility for ensuring that the ZCORN vector (corner-point depth) of an
/// ECLIPSE corner-point description meets certain consistency requirements.
///
/// In particular, the ZCORN values should represent depth and therefore not
/// decrease along pillars.  Moreover, any cell's top corners should not be
/// below that cell's bottom corner on the same pillar and no cell's bottom
/// corner should be below that cell's lower neighbour's upper corner on the
/// same pillar.
///
/// The final two conditions may for instance occur as a result of round-off
/// when gridding software outputs ZCORN as \c float values while using \c
/// double values internally.
///
/// This code will fail if all cells not explicitly deactivated in the input
/// model are twisted, meaning that every such cell have at least one pillar
/// for which ZCORN increases over the cell and at least one pillar for
/// which ZCORN decreases.  This is hopefully a pathological case that does
/// not occur in real input decks.

namespace Opm { namespace UgGridHelpers {

    /// Facility for ensuring that the ZCORN vector of an ECLIPSE
    /// corner-point description meets certain consistency requirements.
    class RepairZCORN
    {
    public:
        /// Constructor.
        ///
        /// \tparam CartDims Representation of Cartesian model dimensions.
        ///     Must support \code operator[]() \endcode.  Typically \code
        ///     std::vector<int> \endcode or \code std::array<int, 3>
        ///     \endcode or similar (e.g., \code std::array<std::size_t, 3>
        ///     \endcode).
        ///
        /// \param[in] zcorn Original ZCORN values.  Consumed by \c
        ///    RepairZCORN object.  Pass a copy of the values if you need to
        ///    retain the original vector for comparison purposes.
        ///
        /// \param[in] actnum Explicit cell activation flag.  Empty input
        ///    treated as all cells active, otherwise standard ECL \c ACTNUM
        ///    array.
        ///
        /// \param[in] cartDims Model's Cartesian dimensions.  Must have at
        ///    least three elements with indices \c 0, \c 1, and \c 2.
        template <class CartDims>
        RepairZCORN(std::vector<double>&&   zcorn,
                    const std::vector<int>& actnum,
                    const CartDims&         cartDims)
            : active_   (actnum, cartDims)
            , zcorn_idx_(cartDims)
            , zcorn_    (std::move(zcorn))
        {
            this->ensureZCornIsDepth();
            this->ensureTopNotBelowBottom();
            this->ensureBottomNotBelowLowerTop(cartDims[2]);
        }

        /// Statistics about modified ZCORN values.
        struct ZCornChangeCount
        {
            /// Number of cells affected by operation.
            std::size_t cells{0};

            /// Total number of cell corners affected by operation.
            std::size_t corners{0};
        };

        /// Retrieve sanitized ZCORN values.
        ///
        /// Destroys internal object state (creates new ZCORN value object
        /// by moving internal state out of \code *this \endcode).
        std::vector<double> destructivelyGrabSanitizedValues()
        {
            return std::move(this->zcorn_);
        }

        /// Whether or not input ZCORN values were interpreted as elevations
        /// and switched to depth values.
        ///
        /// We effect the switching to depth values by reversing the signs
        /// of all elements--including those pertaining to deactivated cells
        /// (i.e., those for which \code ACTNUM == 0 \endcode).
        bool switchedToDepth() const
        {
            return this->switchedToDepth_;
        }

        /// Retrieve statistics about ZCORN values changed when ensuring
        /// that cells' top corners are not below their respective bottom
        /// corners (on the same pillar).
        const ZCornChangeCount& statTopBelowBottom() const
        {
            return this->topBelowBottom_;
        }

        /// Retrieve statistics about ZCORN values changed when ensuring
        /// that cells' bottom corners are not below their respective lower
        /// neighbour's top corners (on the same pillar).
        const ZCornChangeCount& statBottomBelowLowerTop() const
        {
            return this->bottomBelowLowerTop_;
        }

    private:
        /// Simplified mapping of model cells that are not explicitly
        /// deactivated.
        ///
        /// Linear cell IDs are relative to global (uncompressed) cell
        /// numbering.
        class ActiveCells
        {
        public:
            /// Constructor.
            ///
            /// Throws if size of activation flag array does not match model
            /// dimension.
            ///
            /// \tparam CartDims Representation of Cartesian model
            ///     dimensions.  Must support \code operator[]() \endcode.
            ///     Typically \code std::vector<int> \endcode or \code
            ///     std::array<int, 3> \endcode or similar (e.g., \code
            ///     std::array<std::size_t, 3> \endcode).
            ///
            /// \param[in] actnum Explicit cell activation flag.  Empty
            ///    input treated as all cells active, otherwise standard ECL
            ///    \c ACTNUM array.
            ///
            /// \param[in] cartDims Model's Cartesian dimensions.  Must have
            ///    at least three elements.
            template <class CartDims>
            ActiveCells(const std::vector<int>& actnum,
                        const CartDims&         cartDims)
                : nx_(cartDims[0])
                , ny_(cartDims[1])
                , nz_(cartDims[2])
            {
                const auto nglob = nx_ * ny_ * nz_;

                if (actnum.empty()) {
                    this->is_active_.resize(nglob, true);
                    this->acells_   .resize(nglob, 0);

                    std::iota(std::begin(this->acells_),
                              std::end  (this->acells_), 0);
                }
                else if (actnum.size() == nglob) {
                    this->is_active_.resize(nglob, false);
                    this->acells_   .reserve(nglob);

                    for (auto i = 0*nglob; i < nglob; ++i) {
                        if (actnum[i] != 0) {
                            this->acells_.push_back(i);
                            this->is_active_[i] = true;
                        }
                    }
                }
                else {
                    throw std::invalid_argument {
                        "ACTNUM vector does not match global size"
                    };
                }
            }

            /// Representation of Cartesian cell position.
            using IndexTuple = std::array<std::size_t, 3>;

            /// Retrieve cell position (\code (I, J, K) \endcode tuple) in
            /// global grid.
            ///
            /// \param[in] globCell Global cell ID of particular cell.
            ///
            /// \return Cartesian index position of \p globCell.
            IndexTuple getCellIJK(const std::size_t globCell) const
            {
                auto c = globCell;

                auto i = c % nx_;  c /= nx_;
                auto j = c % ny_;  c /= ny_;
                auto k = c % nz_;

                return {{ i, j, k }};
            }

            /// Retrieve global IDs of model's explicitly active cells.
            const std::vector<std::size_t>& activeGlobal() const
            {
                return this->acells_;
            }

            /// Retrieve total number of uncompressed cells in model.
            std::size_t numGlobalCells() const
            {
                return this->nx_ * this->ny_ * this->nz_;
            }

            /// Retrieve global ID of cell immediately below particular
            /// cell.
            ///
            /// Assumes that lateral cell positions are valid.
            ///
            /// \param[in] ijk Cartesian index position of single cell.
            ///
            /// \return Global ID of cell immediately below \p ijk.  If no
            ///    such cell exists, e.g., because the input cell is in the
            ///    bottom layer of the model or all cells below it are
            ///    explicitly deactivated, then function neighbourBelow()
            ///    returns \code static_cast<std::size_t>(-1) \endcode.
            std::size_t neighbourBelow(IndexTuple ijk) const
            {
                if (ijk[2] >= nz_ - 1) {
                    return -1;
                }

                ijk[2] += 1;

                auto below = this->globalCellIdx(ijk);
                while ((below < this->numGlobalCells()) &&
                       (! this->is_active_[below]))
                {
                    ijk[2] += 1;

                    below = this->globalCellIdx(ijk);
                }

                return below;
            }

        private:
            /// Number of cells in model's X direction.
            const std::size_t nx_;

            /// Number of cells in model's Y direction.
            const std::size_t ny_;

            /// Number of cells in model's Z direction.
            const std::size_t nz_;

            /// Predicate for whether or not a particular global cell is
            /// active (stores the result of \code ACTNUM != 0 \endcode for
            /// all global cells.
            std::vector<bool> is_active_;

            /// Global IDs of all of the model's active cells.
            std::vector<std::size_t> acells_;

            /// Translate Cartesian cell position to global cell ID.
            ///
            /// \param[in] Cartesian cell position.  Lateral positions (I,J)
            ///    assumed to be in valid range \code [0 .. nx_-1] \endcode
            ///    and \code [0 .. ny-1] \endcode respectively.
            std::size_t globalCellIdx(const IndexTuple& ijk) const
            {
                if (ijk[2] > nz_ - 1) { return -1; }

                return ijk[0] + nx_*(ijk[1] + ny_*ijk[2]);
            }
        };

        /// Layer of indirection to simplify extracting appropriate subsets
        /// of global ZCORN array.
        class ZCornIndex
        {
        public:
            /// Constructor
            template <class CartDims>
            ZCornIndex(const CartDims& cartDims)
                : nx_         (cartDims[0])
                , ny_         (cartDims[1])
                , nz_         (cartDims[2])
                , layerOffset_((2 * nx_) * (2 * ny_))
            {}

            /// Linear indices into global ZCORN array of cell's corners on
            /// particular pillar.
            struct PillarPointIDX
            {
                /// Linear ZCORN index of cell's top pillar point.
                std::size_t top;

                /// Linear ZCORN index of cell's bottom pillar point.
                std::size_t bottom;
            };

            /// Retrieve cell's top and bottom pillar points on all pillars.
            ///
            /// \param[in] ijk Cartesian position of particular cell.
            ///
            /// \return Linear ZCORN indices of cell \p ijk's top and bottom
            ///    pillar points on all of the cell's connecting pillars.
            template <class IndexTuple>
            std::array<PillarPointIDX, 4>
            pillarPoints(const IndexTuple& ijk) const
            {
                const auto start = this->getStartIDX(ijk);

                return {
                    this->p00(start),
                    this->p10(start),
                    this->p01(start),
                    this->p11(start)
                };
            }

        private:
            /// Number of cells in model's X direction.
            const std::size_t nx_;

            /// Number of cells in model's Y direction.
            const std::size_t ny_;

            /// Number of cells in model's Z direction.
            const std::size_t nz_;

            /// Number of values in a single ZCORN depth layer.
            const std::size_t layerOffset_;

            /// Compute linear ZCORN offset of current cell.
            ///
            /// \param[in] ijk Cartesian position of current cell.
            ///
            /// \return linear ZCORN offset of cell \p ijk.
            template <class IndexTuple>
            std::size_t getStartIDX(const IndexTuple& ijk) const
            {
                return 2*ijk[0] + 2*nx_*(2*ijk[1] + 2*ny_*(2 * ijk[2]));
            }

            /// Retrieve top and bottom pillar point indices of cell's lower
            /// left pillar.
            ///
            /// \param[in] start Linear ZCORN offset of current cell.
            ///
            /// \return Linear ZCORN indices of cell's top and pillar points
            ///    on cell's pillar (0,0).
            PillarPointIDX p00(const std::size_t start) const
            {
                return this->pillarpts(start, this->offset(0, 0));
            }

            /// Retrieve top and bottom pillar point indices of cell's lower
            /// right pillar.
            ///
            /// \param[in] start Linear ZCORN offset of current cell.
            ///
            /// \return Linear ZCORN indices of cell's top and pillar points
            ///    on cell's pillar (1,0).
            PillarPointIDX p10(const std::size_t start) const
            {
                return this->pillarpts(start, this->offset(1, 0));
            }

            /// Retrieve top and bottom pillar point indices of cell's upper
            /// left pillar.
            ///
            /// \param[in] start Linear ZCORN offset of current cell.
            ///
            /// \return Linear ZCORN indices of cell's top and pillar points
            ///    on cell's pillar (0,1).
            PillarPointIDX p01(const std::size_t start) const
            {
                return this->pillarpts(start, this->offset(0, 1));
            }

            /// Retrieve top and bottom pillar point indices of cell's upper
            /// right pillar.
            ///
            /// \param[in] start Linear ZCORN offset of current cell.
            ///
            /// \return Linear ZCORN indices of cell's top and pillar points
            ///    on cell's pillar (1,1).
            PillarPointIDX p11(const std::size_t start) const
            {
                return this->pillarpts(start, this->offset(1, 1));
            }

            /// Compute relative ZCORN offset of cell's local pillars.
            ///
            /// \param[in] i X-direction pillar offset of cell's pillar.
            ///    Must be 0 or 1.
            ///
            /// \param[in] j Y-direction pillar offset of cell's pillar.
            ///    Must be 0 or 1.
            ///
            /// \return Relative ZCORN offset of cell's local pillars.
            std::size_t offset(const std::size_t i, const std::size_t j) const
            {
                assert ((i == 0) || (i == 1));
                assert ((j == 0) || (j == 1));

                return i + j*2*nx_;
            }

            /// Retrieve cell's top and bottom pillar points on single
            /// pillar.
            ///
            /// \param[in] start Linear ZCORN offset of current cell.
            ///
            /// \param[in] off ZCORN offset of particular pillar relative to
            ///    cell's ZCORN indices.
            PillarPointIDX pillarpts(const std::size_t start,
                                     const std::size_t off) const
            {
                return {
                    start + off,
                    start + off + this->layerOffset_
                };
            }
        };

        /// Model's active cells
        const ActiveCells active_;

        /// Indirection map into model's ZCORN array.
        const ZCornIndex zcorn_idx_;

        /// Model's ZCORN array.  Subject to change.
        std::vector<double> zcorn_;

        /// Whether or not initial ZCORN values were interpreted as
        /// elevations (decreasing values for increasing layer index).
        bool switchedToDepth_{false};

        /// Statistics about TBB operation.
        ZCornChangeCount topBelowBottom_;

        /// Statistics about BBLT operation.
        ZCornChangeCount bottomBelowLowerTop_;

        /// Ensure ZCORN values are depths.
        ///
        /// Modifies \code this->zcorn_ \endcode.
        void ensureZCornIsDepth()
        {
            if (this->zcornIsElevation()) {
                for (auto& zc : this->zcorn_) {
                    zc = - zc;
                }

                this->switchedToDepth_ = true;
            }
        }

        /// Ensure ZCORN values honour restriction that cell's top corner is
        /// not below cell's bottom corner on same pillar for every cell.
        ///
        /// Modifies \code this->zcorn_ \endcode.
        void ensureTopNotBelowBottom()
        {
            for (const auto& globCell : this->active_.activeGlobal()) {
                this->ensureTopNotBelowBottom(globCell);
            }
        }

        /// Ensure ZCORN values honour restriction that cell's bottom corner
        /// is not below cell's lower neighbour's upper corner on same
        /// pillar for every cell.
        ///
        /// Modifies \code this->zcorn_ \endcode.
        ///
        /// \param[in] nz Number of layers in model.
        void ensureBottomNotBelowLowerTop(const std::size_t nz)
        {
            if (nz == 0) { return; }

            auto bottomLayer = [nz](const std::size_t layerID)
            {
                return layerID == (nz - 1);
            };

            for (const auto& globCell : this->active_.activeGlobal()) {
                const auto& ijk = this->active_.getCellIJK(globCell);

                if (bottomLayer(ijk[2])) { continue; }

                this->ensureCellBottomNotBelowLowerTop(ijk);
            }
        }

        /// Ensure ZCORN values honour restriction that cell's top corner is
        /// not below cell's bottom corner on same pillar.
        ///
        /// Modifies \code this->zcorn_ \endcode.
        ///
        /// \tparam CellIndex Integral type representing cell indices.
        ///
        /// \param[in] globCell Global ID of particular cell.
        template <class CellIndex>
        void ensureTopNotBelowBottom(const CellIndex globCell)
        {
            const auto cornerCnt0 = this->topBelowBottom_.corners;

            const auto ijk = this->active_.getCellIJK(globCell);

            for (const auto& pt : this->zcorn_idx_.pillarPoints(ijk)) {
                const auto zb = this->zcorn_[pt.bottom];
                auto&      zt = this->zcorn_[pt.top];

                if (zt > zb) {  // Top below bottom (ZCORN is depth)
                    zt = zb;

                    this->topBelowBottom_.corners += 1;
                }
            }

            this->topBelowBottom_.cells +=
                this->topBelowBottom_.corners > cornerCnt0;
        }

        /// Ensure ZCORN values honour restriction that cell's bottom corner
        /// is not below cell's lower neighbour's upper corner on same
        /// pillar.
        ///
        /// Modifies \code this->zcorn_ \endcode.
        ///
        /// \tparam IndexTuple Representation of Cartesian index tuples.
        ///    Typically \code ActiveCells::IndexTuple \endcode.
        ///
        /// \param[in] ijk Cartesian position of particular cell.
        template <class IndexTuple>
        void ensureCellBottomNotBelowLowerTop(const IndexTuple& ijk)
        {
            const auto below = this->active_.neighbourBelow(ijk);

            if (below >= this->active_.numGlobalCells()) {
                return;
            }

            const auto cornerCnt0 = this->bottomBelowLowerTop_.corners;

            const auto& up    = this->zcorn_idx_.pillarPoints(ijk);
            const auto  d_ijk = this->active_.getCellIJK(below);
            const auto& down  = this->zcorn_idx_.pillarPoints(d_ijk);

            for (auto n = up.size(), i = 0*n; i < n; ++i) {
                const auto zbu = this->zcorn_[up  [i].bottom];
                auto&      ztd = this->zcorn_[down[i].top];

                if (zbu > ztd) { // Bottom below lower top (ZCORN is depth)
                    ztd = zbu;

                    this->bottomBelowLowerTop_.corners += 1;
                }
            }

            this->bottomBelowLowerTop_.cells +=
                this->bottomBelowLowerTop_.corners > cornerCnt0;
        }

        /// Determine whether or not input ZCORN array represents elevations
        /// (values decreasing for increasing layer index) or depths (values
        /// increasing for increasing layer index).
        bool zcornIsElevation() const
        {
            auto all_signs = std::vector<int>{};
            all_signs.reserve(this->active_.numGlobalCells());

            for (const auto& globCell : this->active_.activeGlobal()) {
                all_signs.push_back(this->getZCornSign(globCell));
            }

            const int ignore = 0;

            // Elevation == ZCORN decreasing => all non-ignored cell signs
            // negative.
            return (first(all_signs, ignore) == -1)
                && allEqual(all_signs, ignore, -1);
        }

        /// Retrieve sign of single cell's ZCORN change.
        ///
        /// \tparam CellIndex Integral type representing cell indices.
        ///
        /// \param[in] globCell Global ID of particular cell.
        ///
        /// \return Sign of \p globCell's ZCORN change.  Positive (+1) if
        ///    ZCORN does not *DECREASE* along any of the cell's pillars,
        ///    zero (0) if ZCORN increases along some of the pillars and
        ///    decreases along others (or does not change at all), and
        ///    negative (-1) if ZCORN does not *INCREASE* along any of the
        ///    cell's pillars.
        template <typename CellIndex>
        int getZCornSign(const CellIndex globCell) const
        {
            auto sign = [](const double x) -> int
            {
                return (x > 0.0) - (x < 0.0);
            };

            const auto ijk = this->active_.getCellIJK(globCell);

            auto sgn = std::vector<int>{};  sgn.reserve(4);

            for (const auto& pt : this->zcorn_idx_.pillarPoints(ijk)) {
                const auto dz =
                    this->zcorn_[pt.bottom] - this->zcorn_[pt.top];

                sgn.push_back(sign(dz));
            }

            const int ignore = 0;

            if (! allEqual(sgn, ignore)) {
                return 0;
            }

            return sgn.front();
        }

        /// Check whether or not all elements of a vector are equal, while
        /// ignoring particular element value.
        ///
        /// \param[in] coll Sequence of elements.
        ///
        /// \param[in] ignore Special value that should be considered equal
        ///    (i.e., ignored) for comparison purposes.  Typically zero.
        ///
        /// \return Whether or not all non-ignored elements of \p coll have
        ///    the same value.
        bool allEqual(const std::vector<int>& coll,
                      const int               ignore) const
        {
            return this->allEqual(coll, ignore, ignore);
        }

        /// Check whether or not all elements of a vector are equal, while
        /// ignoring particular element value.
        ///
        /// \param[in] coll Sequence of elements.
        ///
        /// \param[in] ignore Special value that should be considered equal
        ///    (i.e., ignored) for comparison purposes.  Typically zero.
        ///
        /// \param[in] lookfor Search hint for which value to look for in
        ///    the element sequence.  If \code lookfor == ignore \endcode,
        ///    then allEqual() will first perform a linear scan over the
        ///    sequence to find the first non-ignored element and
        ///    subsequently compare all elements to this value.  Otherwise,
        ///    the elements will be compared directly to \p lookfor.  In
        ///    essence this parameter is an optimisation that enables
        ///    bypassing the initial value search if the desired value is
        ///    already known.
        ///
        /// \return Whether or not all non-ignored elements of \p coll have
        ///    the same value.
        bool allEqual(const std::vector<int>& coll,
                      const int               ignore,
                      const int               lookfor) const
        {
            const auto x0 = (lookfor != ignore)
                ? lookfor : first(coll, ignore);

            return std::all_of(std::begin(coll), std::end(coll),
                               [x0, ignore](const int xi)
                   {
                       return (xi == x0) || (xi == ignore);
                   });
        }

        /// Find first non-ignored element value in vector.
        ///
        /// \param[in] coll Sequence of elements.
        ///
        /// \param[in] ignore Special value that should be ignored for
        ///    comparison purposes.  Typically zero.
        ///
        /// \return First non-ignored element value in \p coll.  If there
        ///    are no non-ignored elements in \p coll, then first() returns
        ///    \p ignore.
        int first(const std::vector<int>& coll,
                  const int               ignore) const
        {
            auto e = std::end(coll);

            auto p = std::find_if(std::begin(coll), e,
                                  [ignore](const int xi)
                                  {
                                      return xi != ignore;
                                  });

            if (p == e) {
                return ignore;
            }

            return *p;
        }
    };

}} // namespace Opm::UgGridHelpers

#endif // OPM_REPAIRZCORN_HEADER_INCLUDED
