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

#include <array>
#include <cstddef>
#include <functional>
#include <map>
#include <numeric>
#include <vector>

namespace Opm
{

    /// \brief Transform a corner-point grid ZCORN field to account for MINPV processing.
    class MinpvProcessor
    {
    public:

        struct Result {
            std::vector<std::size_t> removed_cells;
            std::map<int,int> nnc;

            void add_nnc(int cell1, int cell2);
        };


        /// \brief Create a processor.
        /// \param[in]   nx   logical cartesian number of cells in I-direction
        /// \param[in]   ny   logical cartesian number of cells in J-direction
        /// \param[in]   nz   logical cartesian number of cells in K-direction
        MinpvProcessor(const int nx, const int ny, const int nz);
        /// Change zcorn so that it respects the minpv property.
        /// \param[in]       thickness thickness of the cell
        /// \param[in]       z_tolerance cells with thickness below z_tolerance will be bypassed in the minpv process.
        /// \param[in]       max_gap  Maximum gap of pinched out layers allowed when creating NNCs
        /// \param[in]       pv       pore volumes of all logical cartesian cells
        /// \param[in]       minpvv   minimum pore volume to accept a cell
        /// \param[in]       actnum   active cells, inactive cells are not considered
        /// \param[in]       mergeMinPVCells flag to determine whether cells below minpv
        /// should be included in the cell below
        /// \param[in, out]  zcorn    ZCORN array to be manipulated
        /// \param[in]       pinchNOGAP whether PINCH has option NOGAP (default false)
        /// \param[in]       pinchOption4ALL Whether option 4 of PINCH is set to all. In that case we do not create NNCs if
        ///                             zero transmissibility in Z direction betwenn the cells is expected.
        ///                             (default: false)
        /// \param[in]       permz     Cell permeability in Z direction.
        /// \oaram[in]       multz     transmissiblity multipliers in Z direction
        /// After processing, all cells that have lower pore volume than minpv
        /// will have the zcorn numbers changed so they are zero-thickness. Any
        /// cell below will be changed to include the deleted volume if mergeMinPCCells is true
        /// els the volume will be lost
        /// \param[in]       tolerance_unique_points Tolerance used to identify points based on their cooridinates.
        Result process(const std::vector<double>& thickness,
                       const double z_tolerance,
                       const double max_gap,
                       const std::vector<double>& pv,
                       const std::vector<double>& minpvv,
                       const std::vector<int>& actnum,
                       const bool mergeMinPVCells,
                       double* zcorn,
                       const bool pinchNOGAP = false,
                       const bool pinchOption4ALL = false,
                       const std::vector<double>& permz = {},
                       const std::function<double(int)>& multZ = [](int){ return 0;},
                       const double tolerance_unique_points = 0) const;
    private:
        double computeGap(const std::array<double,8>& coord_above, const std::array<double,8>& coord_below) const;
        std::array<int,8> cornerIndices(const int i, const int j, const int k) const;
        // Returns the eight z-values associated with a given cell.
        // The ordering is such that i runs fastest. That is, with
        // L = low and H = high:
        // {LLL, HLL, LHL, HHL, LLH, HLH, LHH, HHH }.
        std::array<double, 8> getCellZcorn(const int i, const int j, const int k, const double* z) const;
        void setCellZcorn(const int i, const int j, const int k, const std::array<double, 8>& cellz, double* z) const;
        std::array<int, 3> dims_;
        std::array<int, 3> delta_;
    };

} // namespace Opm

#endif // OPM_MINPVPROCESSOR_HEADER_INCLUDED
