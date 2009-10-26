//===========================================================================
//
// File: ReservoirPropertyCommon.hpp
//
// Created: Mon Oct 26 08:23:31 2009
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
  Copyright 2009 SINTEF ICT, Applied Mathematics.
  Copyright 2009 Statoil ASA.

  This file is part of The Open Reservoir Simulator Project (OpenRS).

  OpenRS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OpenRS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OpenRS.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPENRS_RESERVOIRPROPERTYCOMMON_HEADER
#define OPENRS_RESERVOIRPROPERTYCOMMON_HEADER

#include <dune/common/Units.hpp>
#include <dune/grid/common/EclipseGridParser.hpp>
#include <dune/solvers/common/Matrix.hpp>

namespace Dune
{

    /// @brief A property class for incompressible two-phase flow.
    /// @tparam dim the dimension of the space, used for giving permeability tensors the right size.
    template <int dim, class RPImpl, class RockType>
    class ReservoirPropertyCommon
    {
    public:
        /// @brief Tensor type for read-only access to permeability.
        typedef ImmutableCMatrix PermTensor;
        /// @brief Tensor type to be used for holding copies of permeability tensors.
        typedef OwnCMatrix       MutablePermTensor;
        /// @brief Tensor type for read and write access to permeability.
        typedef SharedCMatrix    SharedPermTensor;

        /// @brief The number of phases 
        enum { NumberOfPhases = 2 };

        /// @brief Default constructor.
        ReservoirPropertyCommon();

        /// @brief Initialize from a grdecl file.
        /// @param parser the parser holding the grdecl data.
	/// @param global_cell the mapping from cell indices to the logical
	///                    cartesian indices of the grdecl file.
        void init(const EclipseGridParser& parser,
                  const std::vector<int>& global_cell,
                  const std::string* rock_list_filename = 0);

        /// @brief Initialize a uniform reservoir.
        /// @param num_cells number of cells in the grid.
        /// @param uniform_poro the uniform porosity.
        /// @param uniform_perm the uniform (scalar) permeability.
        void init(const int num_cells,
                  const double uniform_poro = 0.2,
                  const double uniform_perm = 100.0*prefix::milli*unit::darcy);

	/// @brief Viscosity of first (water) phase.
	/// @return the viscosity value.
        double viscosityFirstPhase() const;

	/// @brief Viscosity of second (oil) phase.
	/// @return the viscosity value.
        double viscositySecondPhase() const;

	/// @brief Density of first (water) phase.
	/// @return the density value.
        double densityFirstPhase() const;

	/// @brief Density of second (oil) phase.
	/// @return the density value.
        double densitySecondPhase() const;

        /// @brief Read-access to porosity.
        /// @param cell_index index of a grid cell.
        /// @return porosity value of the cell.
        double porosity(int cell_index) const;

        /// @brief Read-access to permeability.
        /// @param cell_index index of a grid cell.
        /// @return permeability value of the cell.
        PermTensor permeability(int cell_index) const;

        /// @brief Read- and write-access to permeability. Use with caution.
        /// @param cell_index index of a grid cell.
        /// @return permeability value of the cell.
	SharedPermTensor permeabilityModifiable(int cell_index);

        /// @brief Densities for both phases.
	/// @tparam Vector a class with size() and operator[].
        /// @param cell_index index of a grid cell (not used).
        /// @param[out] density the phase densities.
	///                     Expected to be of size 2 before (and after) the call.
        template<class Vector>
        void phaseDensity(int /*cell_index*/, Vector& density) const;

        /// @brief Difference of densities.
        /// @return densityFirstPhase() - densitySecondPhase()
        double densityDifference() const;

        /// @brief A factor useful in cfl computations.
        /// @return the cfl factor.
        double cflFactor() const;

        /// @brief A factor useful in gravity cfl computations.
        /// @return the gravity cfl factor.
        double cflFactorGravity() const;

        /// @brief Capillary pressure.
        /// @param cell_index index of a grid cell.
	/// @param saturation a saturation value.
        /// @return capillary pressure at the given cell and saturation.
        double capillaryPressure(int cell_index, double saturation) const;

    protected:
	// Methods
        void assignPorosity(const EclipseGridParser& parser,
                            const std::vector<int>& global_cell);
        void assignPermeability(const EclipseGridParser& parser,
                                const std::vector<int>& global_cell);
        void assignRockTable(const EclipseGridParser& parser,
                             const std::vector<int>& global_cell);
        void readRocks(const std::string& rock_list_file);

	// Supporting Barton/Nackman trick (also known as the curiously recurring template pattern).
	RPImpl& asImpl();

	// Data members.
        std::vector<double>        porosity_;
        std::vector<double>        permeability_;
        std::vector<unsigned char> permfield_valid_;
        double density1_;
        double density2_;
        double viscosity1_;
        double viscosity2_;
        double cfl_factor_;
        double cfl_factor_gravity_;
        std::vector<RockType> rock_;
        std::vector<int> cell_to_rock_;
    };


} // namespace Dune

#include "ReservoirPropertyCommon_impl.hpp"

#endif // OPENRS_RESERVOIRPROPERTYCOMMON_HEADER
