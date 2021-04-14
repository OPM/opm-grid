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

#ifndef OPM_MINPVPROCESSOR_HEADER_INCLUDED
#define OPM_MINPVPROCESSOR_HEADER_INCLUDED


#include <opm/grid/utility/ErrorMacros.hpp>
#include <opm/grid/utility/OpmParserIncludes.hpp>

#include <array>

namespace Opm
{

    /// \brief Transform a corner-point grid ZCORN field to account for MINPV processing.
    class MinpvProcessor
    {
    public:

        struct PinchPair {
            int cell1;
            int cell2;
            int removed_cell;

            PinchPair(int c1, int c2, int rm) :
                cell1(std::min(c1, c2)),
                cell2(std::max(c1, c2)),
                removed_cell(rm)
            {}
        };


        /// \brief Create a processor.
        /// \param[in]   nx   logical cartesian number of cells in I-direction
        /// \param[in]   ny   logical cartesian number of cells in J-direction
        /// \param[in]   nz   logical cartesian number of cells in K-direction
        MinpvProcessor(const int nx, const int ny, const int nz);
        /// Change zcorn so that it respects the minpv property.
        /// \param[in]       thickness thickness of the cell
        /// \param[in]       z_tolerance cells with thickness below z_tolerance will be bypassed in the minpv process.
        /// \param[in]       pv       pore volumes of all logical cartesian cells
        /// \param[in]       minpvv   minimum pore volume to accept a cell
        /// \param[in]       actnum   active cells, inactive cells are not considered
        /// \param[in]       mergeMinPVCells flag to determine whether cells below minpv
        /// should be included in the cell below
        /// \param[in, out]  zcorn    ZCORN array to be manipulated
        /// After processing, all cells that have lower pore volume than minpv
        /// will have the zcorn numbers changed so they are zero-thickness. Any
        /// cell below will be changed to include the deleted volume if mergeMinPCCells is true
        /// els the volume will be lost
        std::vector<PinchPair> process(const std::vector<double>& thickness,
                                       const double z_tolerance,
                                       const std::vector<double>& pv,
                                       const std::vector<double>& minpvv,
                                       const std::vector<int>& actnum,
                                       const bool mergeMinPVCells,
                                       double* zcorn,
                                       bool pinchNOGAP = false) const;
    private:
        std::array<int,8> cornerIndices(const int i, const int j, const int k) const;
        std::array<double, 8> getCellZcorn(const int i, const int j, const int k, const double* z) const;
        void setCellZcorn(const int i, const int j, const int k, const std::array<double, 8>& cellz, double* z) const;
        std::array<int, 3> dims_;
        std::array<int, 3> delta_;
    };

    inline MinpvProcessor::MinpvProcessor(const int nx, const int ny, const int nz) :
        dims_( {{nx,ny,nz}} ),
        delta_( {{1 , 2*nx , 4*nx*ny}} )
    { }



    inline std::vector<MinpvProcessor::PinchPair> MinpvProcessor::process(const std::vector<double>& thickness,
                                                                          const double z_tolerance,
                                                                          const std::vector<double>& pv,
                                                                          const std::vector<double>& minpvv,
                                                                          const std::vector<int>& actnum,
                                                                          const bool mergeMinPVCells,
                                                                          double* zcorn,
                                                                          bool pinchNOGAP) const
    {
        // Algorithm:
        // 1. Process each column of cells (with same i and j
        //    coordinates) from top (low k) to bottom (high k).
        // 2. For each cell 'c' visited, check if its pore volume
        //    pv[c] is less than minpvv[c] .
        // 3. If below the minpv threshold, move the lower four
        //    zcorn associated with the cell c to coincide with
        //    the upper four (so it becomes degenerate).
        // 4. Look for the next active cell by skipping
        //    inactive cells with thickness below the z_tolerance.
        // 5. If mergeMinPVcells:
        //    is true, the higher four zcorn associated with the cell below
        //    is moved to these values (so it gains the deleted volume).
        //    is false, a nnc is created between the cell above the removed
        //    cell and the cell below it. Note that the connection is only
        //    created if the cell below and above are active
        //    Inactive cells with thickness below z_tolerance and cells with porv<minpv
        //    are bypassed.
        // 6. If pinchNOGAP (only has an effect if mergeMinPVcells==false holds):
        //    is true active cells with porevolume less than minpvv will only be disregarded
        //    if their thickness is below z_tolerance and nncs will be created in this case.


        // return a list of the non-neighbor connection.
        std::vector<PinchPair> nnc;

        // Check for sane input sizes.
        const size_t log_size = dims_[0] * dims_[1] * dims_[2];
        if (pv.size() != log_size) {
            OPM_THROW(std::runtime_error, "Wrong size of PORV input, must have one element per logical cartesian cell.");
        }
        if (!actnum.empty() && actnum.size() != log_size) {
            OPM_THROW(std::runtime_error, "Wrong size of ACTNUM input, must have one element per logical cartesian cell.");
        }

        // Main loop.
        for (int kk = 0; kk < dims_[2]; ++kk) {
            for (int jj = 0; jj < dims_[1]; ++jj) {
                for (int ii = 0; ii < dims_[0]; ++ii) {
                    const int c = ii + dims_[0] * (jj + dims_[1] * kk);
                    if (pv[c] < minpvv[c] && (actnum.empty() || actnum[c])) {
                        // Move deeper (higher k) coordinates to lower k coordinates.
                        // i.e remove the cell
                        std::array<double, 8> cz = getCellZcorn(ii, jj, kk, zcorn);
                        for (int count = 0; count < 4; ++count) {
                            cz[count + 4] = cz[count];
                        }
                        setCellZcorn(ii, jj, kk, cz, zcorn);

                        // Find the next cell
                        int kk_iter = kk + 1;
                        if (kk_iter == dims_[2]) // we are at the end of the pillar.
                            continue;

                        int c_below = ii + dims_[0] * (jj + dims_[1] * (kk_iter));
                        // bypass inactive cells with thickness less then the tolerance
                        while ( ((actnum.empty() || !actnum[c_below]) && (thickness[c_below] <= z_tolerance))  ){
                            // move these cell to the posistion of the first cell to make the
                            // coordinates strictly sorted
                            setCellZcorn(ii, jj, kk_iter, cz, zcorn);
                            kk_iter ++;
                            if (kk_iter == dims_[2])
                                break;

                            c_below = ii + dims_[0] * (jj + dims_[1] * (kk_iter));
                        }

                        if (kk_iter == dims_[2]) // we have come to the end of the pillar.
                            continue;

                        // create nnc if false or merge the cells if true
                        if (!mergeMinPVCells) {

                            // We are at the top, so no nnc is created.
                            if (kk == 0)
                                continue;

                            int c_above = ii + dims_[0] * (jj + dims_[1] * (kk - 1));

                            // Bypass inactive cells with thickness below tolerance and active cells with volume below minpv
                            auto above_active = actnum.empty() || actnum[c_above];
                            auto above_inactive = actnum.empty() || !actnum[c_above]; // \todo Kept original, but should be !actnum.empty() && !actnum[c_above]
                            auto above_thin = thickness[c_above] < z_tolerance;
                            auto above_small_pv = pv[c_above] < minpvv[c_above];
                            if ((above_inactive && above_thin) || (above_active && above_small_pv
                                                                   && (!pinchNOGAP || above_thin) ) ) {
                                for (int topk = kk - 2; topk > 0; --topk) {
                                    c_above = ii + dims_[0] * (jj + dims_[1] * (topk));
                                    above_active = actnum.empty() || actnum[c_above];
                                    above_inactive = actnum.empty() || !actnum[c_above];
                                    auto above_significant_pv = pv[c_above] > minpvv[c_above];
                                    auto above_broad = thickness[c_above] > z_tolerance;
                                    // \todo if condition seems wrong and should be the negation of above?
                                    if ( (above_active && (above_significant_pv || (pinchNOGAP && above_broad) ) ) || (above_inactive && above_broad)) {
                                        break;
                                    }
                                }
                            }

                            // Bypass inactive cells with thickness below tolerance and active cells with volume below minpv
                            auto below_active = actnum.empty() || actnum[c_below];
                            auto below_inactive = actnum.empty() || !actnum[c_below]; // \todo Kept original, but should be !actnum.empty() && !actnum[c_below]
                            auto below_thin = thickness[c_below] < z_tolerance;
                            auto below_small_pv = pv[c_below] < minpvv[c];
                            if ((below_inactive && below_thin) || (below_active && below_small_pv
                                                                   && (!pinchNOGAP || below_thin ) ) ) {
                                for (int botk = kk_iter + 1; botk <  dims_[2]; ++botk) {
                                    c_below = ii + dims_[0] * (jj + dims_[1] * (botk));
                                    below_active = actnum.empty() || actnum[c_below];
                                    below_inactive = actnum.empty() || !actnum[c_below]; // \todo Kept original, but should be !actnum.empty() && !actnum[c_below]
                                    auto below_significant_pv = pv[c_below] > minpvv[c_below];
                                    auto below_broad = thickness[c_above] > z_tolerance;
                                    // \todo if condition seems wrong and should be the negation of above?
                                    if ( (below_active && (below_significant_pv || (pinchNOGAP && below_broad) ) ) || (below_inactive && below_broad)) {
                                        break;
                                    }
                                }
                            }

                            // Add a connection if the cell above and below is active and has porv > minpv
                            if ((actnum.empty() || (actnum[c_above] && actnum[c_below])) && pv[c_above] > minpvv[c_above] && pv[c_below] > minpvv[c_below]) {
                                    nnc.emplace_back(c_above, c_below, c);
                            }
                        } else {

                            // Set lower k coordinates of cell below to upper cells's coordinates.
                            // i.e fill the void using the cell below
                            std::array<double, 8> cz_below = getCellZcorn(ii, jj, kk_iter, zcorn);
                            for (int count = 0; count < 4; ++count) {
                                cz_below[count] = cz[count];
                            }
                            setCellZcorn(ii, jj, kk_iter, cz_below, zcorn);
                        }

                    }
                }

            }
        }
        return nnc;
    }



    inline std::array<int,8> MinpvProcessor::cornerIndices(const int i, const int j, const int k) const
    {
        const int ix = 2*(i*delta_[0] + j*delta_[1] + k*delta_[2]);
        std::array<int, 8> ixs = {{ ix,                         ix + delta_[0],
                                    ix + delta_[1],             ix + delta_[1] + delta_[0],
                                    ix + delta_[2],             ix + delta_[2] + delta_[0],
                                    ix + delta_[2] + delta_[1], ix + delta_[2] + delta_[1] + delta_[0] }};

        return ixs;
    }



    // Returns the eight z-values associated with a given cell.
    // The ordering is such that i runs fastest. That is, with
    // L = low and H = high:
    // {LLL, HLL, LHL, HHL, LLH, HLH, LHH, HHH }.
    inline std::array<double, 8> MinpvProcessor::getCellZcorn(const int i, const int j, const int k, const double* z) const
    {
        const std::array<int, 8> ixs = cornerIndices(i, j, k);
        std::array<double, 8> cellz;
        for (int count = 0; count < 8; ++count) {
            cellz[count] = z[ixs[count]];
        }
        return cellz;
    }



    inline void MinpvProcessor::setCellZcorn(const int i, const int j, const int k, const std::array<double, 8>& cellz, double* z) const
    {
        const std::array<int, 8> ixs = cornerIndices(i, j, k);
        for (int count = 0; count < 8; ++count) {
            z[ixs[count]] = cellz[count];
        }
    }



} // namespace Opm

#endif // OPM_MINPVPROCESSOR_HEADER_INCLUDED
