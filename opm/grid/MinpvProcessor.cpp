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
#include <config.h>
#include <opm/grid/MinpvProcessor.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <algorithm>

namespace Opm {


void MinpvProcessor::Result::add_nnc(int cell1, int cell2)
{
    auto key = std::min(cell1, cell2);
    auto value = std::max(cell1,cell2);

    this->nnc.insert({key, value});
}

MinpvProcessor::MinpvProcessor(const int nx, const int ny, const int nz) :
    dims_( {{nx,ny,nz}} ),
    delta_( {{1 , 2*nx , 4*nx*ny}} )
{ }

double MinpvProcessor::computeGap(const std::array<double,8>& coord_above,
                                  const std::array<double,8>& coord_below) const
{
    std::array<double, 4> vertical_gap;
    for (std::size_t i = 0; i < 4; ++i) {
        vertical_gap[i] = coord_below[i] - coord_above[4 + i];
        assert(vertical_gap[i] >= 0);
    }
    double min_val = *std::min_element(vertical_gap.begin(), vertical_gap.end());
    if (min_val < 1e-6){
        return 0;
    }

    return min_val;
}

MinpvProcessor::Result
MinpvProcessor::process(const std::vector<double>& thickness,
                        const double z_tolerance,
                        const double max_gap,
                        const std::vector<double>& pv,
                        const std::vector<double>& minpvv,
                        const std::vector<int>& actnum,
                        const bool mergeMinPVCells,
                        double* zcorn,
                        bool pinchNOGAP,
                        bool pinchOption4ALL,
                        const std::vector<double>& permz,
                        const std::function<double(int)>& multz,
                        const double tolerance_unique_points) const
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
    // 7. Default maximum gap allowed if option 3 in PINCH is omitted is 1e20 in any unit
    //    If pinch is not specified it is 1e20 in SI units, which should still behave similar
    //    to infinity.


    Result result;

    // Check for sane input sizes.
    const size_t log_size = dims_[0] * dims_[1] * dims_[2];
    if (pv.size() != log_size) {
        OPM_THROW(std::runtime_error, "Wrong size of PORV input, must have one element per logical cartesian cell.");
    }
    if (!actnum.empty() && actnum.size() != log_size) {
        OPM_THROW(std::runtime_error, "Wrong size of ACTNUM input, must have one element per logical cartesian cell.");
    }
    if (pinchOption4ALL && permz.empty()) {
        OPM_THROW(std::runtime_error,
                  "If option 4 of PINCH keyword is ALL, then the deck needs to specify PERMZ or PERMX");
        }

    // Main loop.
    for (int jj = 0; jj < dims_[1]; ++jj) {
        for (int ii = 0; ii < dims_[0]; ++ii) {
            for (int kk = 0; kk < dims_[2]; ++kk) {
                // For a corner case for option ALL
                // where one of the cells in-between has 0 transmissibility
                // we will omit the nnc
                bool option4ALLZero = false;
                const int c = ii + dims_[0] * (jj + dims_[1] * kk);
                bool c_active = actnum.empty() || actnum[c];
                bool c_thin = (thickness[c] <= z_tolerance);
                bool c_thin_inactive = !c_active && c_thin;
                bool c_low_pv_active = (pv[c] < minpvv[c] && c_active) || (c_thin && c_active);

                if (c_low_pv_active || c_thin_inactive) {
                    std::array<double, 8> cz = getCellZcorn(ii, jj, kk, zcorn);
                    // Cell is either inactive or made inactive due to MINPV

                    // Move deeper (higher k) coordinates to lower k coordinates.
                    // i.e remove the cell
                    for (int count = 0; count < 4; ++count) {
                        cz[count + 4] = cz[count];
                    }
                    setCellZcorn(ii, jj, kk, cz, zcorn);

                    if (c_low_pv_active) {
                        // Inactive due to MINPV mark it as removed
                        result.removed_cells.push_back(c);
                    }

                    if (kk == dims_[2] - 1) {
                        // this is cell at the bottom of the grid
                        // no neighbor below for an NNC.
                        continue;
                    }

                    // In the case of PinchNOGAP this cell must be thin to allow NNCs, if it was deactivated
                    // via PINCH, too.
                    // In addition skip NNC if PINCH option 4 is ALL and we know that Z transmissibilty will
                    // be zero because of multz or permz
                    bool nnc_allowed = (!c_low_pv_active || (!pinchNOGAP || thickness[c] <= z_tolerance))
                        && (!pinchOption4ALL || (permz[c] != 0.0 && multz(c) != 0.0) );

                    if (pinchOption4ALL)
                    {
                        option4ALLZero = option4ALLZero || (!permz.empty() && permz[c] == 0) || multz(c) == 0;
                    }

                    // Find the next cell below
                    int kk_iter = kk + 1;

                    int c_below = ii + dims_[0] * (jj + dims_[1] * (kk_iter));
                    bool active = actnum.empty() || actnum[c_below];
                    bool thin = (thickness[c_below] <= z_tolerance);
                    bool thin_inactive = !active && thin;
                    bool low_pv_active = ((pv[c_below] < minpvv[c_below]) && active) || (thin && active);


                    while ( (thin_inactive || low_pv_active) && kk_iter < dims_[2] )
                    {
                        // bypass inactive cells with thickness less then the tolerance
                        if (thin_inactive)
                        {
                            // move these cell to the position of the first cell to make the
                            // coordinates strictly sorted
                            setCellZcorn(ii, jj, kk_iter, cz, zcorn);
                        }
                        if (low_pv_active)
                        {
                            // In the case of PichNOGAP this cell must be thin to allow NNCs, too.
                            nnc_allowed = nnc_allowed && (!pinchNOGAP || thin);
                            // Cell is made inactive due to MINPV
                            // It might make sense to always proceed as in the else branch,
                            // but we try to keep changes due to refactoring smalle here and
                            // mimic the old approach
                            if (mergeMinPVCells)
                            {
                                // original algorithm would have extended this cells before
                                // the collapsing. Doing the same.
                                setCellZcorn(ii, jj, kk_iter, cz, zcorn);
                            }
                            else
                            {
                                // original algorithm collapses the unextended cell
                                cz = getCellZcorn(ii, jj, kk_iter, zcorn);
                                for (int count = 0; count < 4; ++count) {
                                    cz[count + 4] = cz[count];
                                }
                                setCellZcorn(ii, jj, kk_iter, cz, zcorn);
                            }
                            result.removed_cells.push_back(c_below);
                        }
                        // Skip NNC if PINCH option 4 is ALL and we know that Z transmissibilty will
                        // be zero because of multz or permz
                        nnc_allowed = nnc_allowed &&
                            (!pinchOption4ALL || (permz[c_below] != 0.0 && multz(c_below) != 0.0));

                        if (pinchOption4ALL) {
                            option4ALLZero = option4ALLZero || (!permz.empty() && permz[c_below] == 0) || multz(c_below) == 0;
                        }

                        // move to next lower cell
                        kk_iter = kk_iter + 1;
                        if (kk_iter == dims_[2])
                        {
                            break;
                        }

                        c_below = ii + dims_[0] * (jj + dims_[1] * (kk_iter));
                        active = actnum.empty() || actnum[c_below];
                        thin = (thickness[c_below] <= z_tolerance);
                        thin_inactive = (!actnum.empty() && !actnum[c_below]) && thin;
                        low_pv_active = ((pv[c_below] < minpvv[c_below]) && active) || (thin && active);
                    }

                    // create nnc if false or merge the cells if true
                    if (mergeMinPVCells){// && c_low_pv_active) {
                        // try to make a topological connected grid
                        // in Flow this is currently called only for edge_conformal grids
                        // however zcorn inbetween is not modified to make zcorn sorted
                        // Set lower k coordinates of cell below to upper cells's coordinates.
                        // i.e fill the void using the cell below
                        if (kk==0  || kk_iter == dims_[2]) {
                            kk = kk_iter;
                            continue;
                        }
                        //bottom cell not active, hence no nnc is created
                        if (!actnum.empty() && !actnum[c_below]) {
                            kk = kk_iter;
                            continue;
                        }

                        std::array<double, 8> cz_below = getCellZcorn(ii, jj, kk_iter, zcorn);
                        for (int count = 0; count < 4; ++count) {
                            cz_below[count] = cz[count];
                        }

                        setCellZcorn(ii, jj, kk_iter, cz_below, zcorn);
                    }
                    else
                    {

                        // No top or bottom cell, so no nnc is created.
                        if (kk == 0 || kk_iter == dims_[2]) {
                            kk = kk_iter;
                            continue;
                        }
                        // bottom cell not active, hence no nnc is created
                        if (!actnum.empty() && !actnum[c_below]) {
                            kk = kk_iter;
                            continue;
                        }

                        // Bypass inactive cells with thickness below tolerance and
                        // active cells with volume below minpv
                        int k_above = kk-1;
                        int c_above = ii + dims_[0] * (jj + dims_[1] * (kk-1));
                        auto above_active = actnum.empty() || actnum[c_above];
                        auto above_inactive = !actnum.empty() && !actnum[c_above];
                        auto above_thin = thickness[c_above] < z_tolerance;
                        auto above_small_pv = (pv[c_above] < minpvv[c_above]) ||  above_thin;

                        if ((above_inactive && above_thin) || (above_active && above_small_pv
                                                               && (!pinchNOGAP || above_thin) ) ) {
                            for (k_above = kk - 2; k_above > 0; --k_above) {
                                c_above = ii + dims_[0] * (jj + dims_[1] * (k_above));
                                above_active = actnum.empty() || actnum[c_above];
                                above_inactive = !actnum.empty() && !actnum[c_above];
                                auto above_significant_pv = !((pv[c_above] < minpvv[c_above]) || (thickness[c_above] < z_tolerance));
                                auto above_broad = thickness[c_above] > z_tolerance;

                                // \todo if condition seems wrong and should be the negation of above?
                                if ( (above_active && (above_significant_pv || (pinchNOGAP && above_broad) ) ) || (above_inactive && above_broad)) {
                                    break;
                                }

                                nnc_allowed = nnc_allowed &&
                                    (!pinchOption4ALL || (permz[c_above] != 0.0 && multz(c_above) != 0.0) );

                                if (pinchOption4ALL) {
                                    option4ALLZero =  option4ALLZero || (!permz.empty() && permz[c_above] == 0.0) || multz(c_above) == 0.0;
                                }
                            }
                        }

                        // Allow nnc only of total thickness of pinched out cells is below threshold.
                        // and sum of gaps is below threshold
                        const std::array<double, 8> cz_below = getCellZcorn(ii, jj, kk_iter, zcorn);
                        const std::array<double, 8> cz_above = getCellZcorn(ii, jj, k_above, zcorn);
                        // top cell might not have been inspected for option 4 ALL before
                        option4ALLZero = option4ALLZero || (!permz.empty() && permz[c_above] == 0.0) || multz(c_above) == 0.0;
                        nnc_allowed = nnc_allowed && (computeGap(cz_above, cz_below) < max_gap) && (!pinchOption4ALL || !option4ALLZero) ;

                        //bool
                        above_small_pv = (pv[c_above] < minpvv[c_above]) || (thickness[c_above] < z_tolerance);
                        bool below_small_pv = (pv[c_below] < minpvv[c_below]) || (thickness[c_below] < z_tolerance);
                        if ( nnc_allowed &&
                             (actnum.empty() || (actnum[c_above] && actnum[c_below])) &&
                             !(above_small_pv) && !(below_small_pv) ){
                            result.add_nnc(c_above, c_below);
                        }
                        kk = kk_iter;
                    }
                }
                else
                {
                    if (kk < dims_[2] - 1 && (actnum.empty() || actnum[c]) && !((pv[c] < minpvv[c]) || (thickness[c] < z_tolerance)) &&
                        multz(c) != 0.0)
                    {
                        // Check whether there is a gap to the neighbor below whose thickness is less
                        // than MAX_GAP. In that case we need to create an NNC if there is a gap between the two cells.
                        int kk_below = kk + 1;
                        int c_below = ii + dims_[0] * (jj + dims_[1] * kk_below);

                        if ( (actnum.empty() || actnum[c_below])
                             && 
                             !((pv[c_below] < minpvv[c_below])  || (thickness[c_below] < z_tolerance) ) ) 
                        {
                            // Check MAX_GAP threshold
                            std::array<double, 8> cz = getCellZcorn(ii, jj, kk, zcorn);
                            std::array<double, 8> cz_below = getCellZcorn(ii, jj, kk_below, zcorn);
                            bool vertically_connected = true; // If true a connection will be there anyway -> Skip NNC

                            for(int i = 0; i < 4; ++i) {
                                vertically_connected = vertically_connected && std::abs(cz_below[i] - cz[4+i])
                                    <= tolerance_unique_points;
                            }

                            if (!vertically_connected && computeGap(cz, cz_below) < max_gap) {
                                result.add_nnc(c, c_below);
                            }
                        }
                    }
                }
            }
        }
    }

    return result;
}

std::array<int,8>
MinpvProcessor::cornerIndices(const int i, const int j, const int k) const
{
    const int ix = 2*(i*delta_[0] + j*delta_[1] + k*delta_[2]);
    std::array<int, 8> ixs = {{ ix,                         ix + delta_[0],
                                ix + delta_[1],             ix + delta_[1] + delta_[0],
                                ix + delta_[2],             ix + delta_[2] + delta_[0],
                                ix + delta_[2] + delta_[1], ix + delta_[2] + delta_[1] + delta_[0] }};

    return ixs;
}

std::array<double, 8>
MinpvProcessor::getCellZcorn(const int i, const int j,
                             const int k, const double* z) const
{
    const std::array<int, 8> ixs = cornerIndices(i, j, k);
    std::array<double, 8> cellz;
    for (int count = 0; count < 8; ++count) {
        cellz[count] = z[ixs[count]];
    }
    return cellz;
}

void MinpvProcessor::setCellZcorn(const int i, const int j, const int k,
                             const std::array<double, 8>& cellz, double* z) const
{
    const std::array<int, 8> ixs = cornerIndices(i, j, k);
    for (int count = 0; count < 8; ++count) {
        z[ixs[count]] = cellz[count];
    }
}

} // namespace Opm
