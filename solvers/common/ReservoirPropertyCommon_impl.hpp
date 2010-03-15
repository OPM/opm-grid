//===========================================================================
//
// File: ReservoirPropertyCommon_impl.hpp
//
// Created: Mon Oct 26 08:29:09 2009
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
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.

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

#ifndef OPENRS_RESERVOIRPROPERTYCOMMON_IMPL_HEADER
#define OPENRS_RESERVOIRPROPERTYCOMMON_IMPL_HEADER


#include <fstream>
#include <boost/static_assert.hpp>
#include <boost/array.hpp>
#include <dune/grid/common/EclipseGridInspector.hpp>

namespace Dune
{

    namespace {

        /// @brief
        ///    Classify and verify a given permeability specification
        ///    from a structural point of view.  In particular, we
        ///    verify that there are no off-diagonal permeability
        ///    components such as @f$k_{xy}@f$ unless the
        ///    corresponding diagonal components are known as well.
        ///
        /// @param parser [in]
        ///    An Eclipse data parser capable of answering which
        ///    permeability components are present in a given input
        ///    deck.
        ///
        /// @return
        ///    An enum value with the following possible values:
        ///        ScalarPerm     only one component was given.
        ///        DiagonalPerm   more than one component given.
        ///        TensorPerm     at least one cross-component given.
        ///        None           no components given.
        ///        Invalid        invalid set of components given.
        PermeabilityKind classifyPermeability(const EclipseGridParser& parser)
        {
            const bool xx = parser.hasField("PERMX" );
            const bool xy = parser.hasField("PERMXY");
            const bool xz = parser.hasField("PERMXZ");

            const bool yx = parser.hasField("PERMYX");
            const bool yy = parser.hasField("PERMY" );
            const bool yz = parser.hasField("PERMYZ");

            const bool zx = parser.hasField("PERMZX");
            const bool zy = parser.hasField("PERMZY");
            const bool zz = parser.hasField("PERMZ" );

            int num_comp = xx + xy + xz + yx + yy + yz + zx + zy + zz;
            int num_cross_comp = xy + xz + yx + yz + zx + zy;
            PermeabilityKind retval = None;
            if (num_cross_comp) {
                retval = TensorPerm;
            } else {
                if (num_comp == 1) {
                    retval = ScalarPerm;
                } else if (num_comp >= 2) {
                    retval = DiagonalPerm;
                }
            }

            bool ok = true;
            if (num_comp) {
                // At least one tensor component specified on input.
                // Verify that any remaining components are OK from a
                // structural point of view.  In particular, there
                // must not be any cross-components (e.g., k_{xy})
                // unless the corresponding diagonal component (e.g.,
                // k_{xx}) is present as well...
                //
                ok =        xx || !(xy || xz || yx || zx) ;
                ok = ok && (yy || !(yx || yz || xy || zy));
                ok = ok && (zz || !(zx || zy || xz || yz));
            }
            if (!ok) {
                retval = Invalid;
            }

            return retval;
        }


        /// @brief
        ///    Copy isotropic (scalar) permeability to other diagonal
        ///    components if the latter have not (yet) been assigned a
        ///    separate value.  Specifically, this function assigns
        ///    copies of the @f$i@f$ permeability component (e.g.,
        ///    'PERMX') to the @f$j@f$ and @f$k@f$ permeability (e.g.,
        ///    'PERMY' and 'PERMZ') components if these have not
        ///    previously been assigned.
        ///
        /// @param kmap
        ///    Permeability indirection map.  In particular @code
        ///    kmap[i] @endcode is the index (an integral number in
        ///    the set [1..9]) into the permeability tensor
        ///    representation of function @code fillTensor @endcode
        ///    which represents permeability component @code i
        ///    @endcode.
        ///
        /// @param [in] i
        /// @param [in] j
        /// @param [in] k
        void setScalarPermIfNeeded(boost::array<int,9>& kmap,
                                   int i, int j, int k)
        {
            if (kmap[j] == 0) { kmap[j] = kmap[i]; }
            if (kmap[k] == 0) { kmap[k] = kmap[i]; }
        }


        /// @brief
        ///   Extract pointers to appropriate tensor components from
        ///   input deck.  The permeability tensor is, generally,
        ///   @code
        ///        [ kxx  kxy  kxz ]
        ///    K = [ kyx  kyy  kyz ]
        ///        [ kzx  kzy  kzz ]
        ///   @endcode
        ///   We store these values in a linear array using natural
        ///   ordering with the column index cycling the most rapidly.
        ///   In particular we use the representation
        ///   @code
        ///        [  0    1    2    3    4    5    6    7    8  ]
        ///    K = [ kxx, kxy, kxz, kyx, kyy, kyz, kzx, kzy, kzz ]
        ///   @endcode
        ///   Moreover, we explicitly enforce symmetric tensors by
        ///   assigning
        ///   @code
        ///     3     1       6     2       7     5
        ///    kyx = kxy,    kzx = kxz,    kzy = kyz
        ///   @endcode
        ///   However, we make no attempt at enforcing positive
        ///   definite tensors.
        ///
        /// @param [in]  parser
        ///    An Eclipse data parser capable of answering which
        ///    permeability components are present in a given input
        ///    deck as well as retrieving the numerical value of each
        ///    permeability component in each grid cell.
        ///
        /// @param [out] tensor
        /// @param [out] kmap
        PermeabilityKind fillTensor(const EclipseGridParser&                 parser,
                                    std::vector<const std::vector<double>*>& tensor,
                                    boost::array<int,9>&                     kmap)
        {
            PermeabilityKind kind = classifyPermeability(parser);
            if (kind == Invalid) {
                THROW("Invalid set of permeability fields given.");
            }
            ASSERT (tensor.size() == 1);
            for (int i = 0; i < 9; ++i) { kmap[i] = 0; }

            enum { xx, xy, xz,    // 0, 1, 2
                   yx, yy, yz,    // 3, 4, 5
                   zx, zy, zz };  // 6, 7, 8

            // -----------------------------------------------------------
            // 1st row: [kxx, kxy, kxz]
            if (parser.hasField("PERMX" )) {
                kmap[xx] = tensor.size();
                tensor.push_back(&parser.getFloatingPointValue("PERMX" ));

                setScalarPermIfNeeded(kmap, xx, yy, zz);
            }
            if (parser.hasField("PERMXY")) {
                kmap[xy] = kmap[yx] = tensor.size();  // Enforce symmetry.
                tensor.push_back(&parser.getFloatingPointValue("PERMXY"));
            }
            if (parser.hasField("PERMXZ")) {
                kmap[xz] = kmap[zx] = tensor.size();  // Enforce symmetry.
                tensor.push_back(&parser.getFloatingPointValue("PERMXZ"));
            }

            // -----------------------------------------------------------
            // 2nd row: [kyx, kyy, kyz]
            if (parser.hasField("PERMYX")) {
                kmap[yx] = kmap[xy] = tensor.size();  // Enforce symmetry.
                tensor.push_back(&parser.getFloatingPointValue("PERMYX"));
            }
            if (parser.hasField("PERMY" )) {
                kmap[yy] = tensor.size();
                tensor.push_back(&parser.getFloatingPointValue("PERMY" ));

                setScalarPermIfNeeded(kmap, yy, zz, xx);
            }
            if (parser.hasField("PERMYZ")) {
                kmap[yz] = kmap[zy] = tensor.size();  // Enforce symmetry.
                tensor.push_back(&parser.getFloatingPointValue("PERMYZ"));
            }

            // -----------------------------------------------------------
            // 3rd row: [kzx, kzy, kzz]
            if (parser.hasField("PERMZX")) {
                kmap[zx] = kmap[xz] = tensor.size();  // Enforce symmetry.
                tensor.push_back(&parser.getFloatingPointValue("PERMZX"));
            }
            if (parser.hasField("PERMZY")) {
                kmap[zy] = kmap[yz] = tensor.size();  // Enforce symmetry.
                tensor.push_back(&parser.getFloatingPointValue("PERMZY"));
            }
            if (parser.hasField("PERMZ" )) {
                kmap[zz] = tensor.size();
                tensor.push_back(&parser.getFloatingPointValue("PERMZ" ));

                setScalarPermIfNeeded(kmap, zz, xx, yy);
            }
            return kind;
        }

    } // anonymous namespace




    // ----- Methods of ReservoirPropertyCommon -----




    template <int dim, class RPImpl, class RockType>
    ReservoirPropertyCommon<dim, RPImpl, RockType>::ReservoirPropertyCommon()
#if 1
        : density1_  (1013.9*unit::kilogram/unit::cubic(unit::meter)),
          density2_  ( 834.7*unit::kilogram/unit::cubic(unit::meter)),
          viscosity1_(   1.0*prefix::centi*unit::Poise),
          viscosity2_(   3.0*prefix::centi*unit::Poise),
#else
        : density1_  (1000.0*unit::kilogram/unit::cubic(unit::meter)),
          density2_  (1000.0*unit::kilogram/unit::cubic(unit::meter)),
          viscosity1_(   1.0*prefix::centi*unit::Poise),
          viscosity2_(   1.0*prefix::centi*unit::Poise),
#endif
          permeability_kind_(Invalid)
    {
    }


    template <int dim, class RPImpl, class RockType>
    void ReservoirPropertyCommon<dim, RPImpl, RockType>::init(const EclipseGridParser& parser,
                                                              const std::vector<int>& global_cell,
                                                              const double perm_threshold,
                                                              const std::string* rock_list_filename,
                                                              const bool use_jfunction_scaling,
                                                              const double sigma,
                                                              const double theta)
    {
        // This code is mostly copied from ReservoirPropertyCommon::init(...).
        BOOST_STATIC_ASSERT(dim == 3);

        permfield_valid_.assign(global_cell.size(),
                                std::vector<unsigned char>::value_type(0));

        assignPorosity    (parser, global_cell);
        assignPermeability(parser, global_cell, perm_threshold);
        assignRockTable   (parser, global_cell);

        if (rock_list_filename) {
            readRocks(*rock_list_filename);
        }

        // Added section. This is a hack, because not all rock classes
        // may care about J-scaling. They still have to implement
        // setUseJfunctionScaling() and setSigmaAndTheta(), though the
        // latter may throw if called for a rock where it does not make sense.
        int num_rocks = rock_.size();
        for (int i = 0; i < num_rocks; ++i) {
            rock_[i].setUseJfunctionScaling(use_jfunction_scaling);
            if (use_jfunction_scaling) {
                rock_[i].setSigmaAndTheta(sigma, theta);
            }
        }
        // End of added section.

        asImpl().computeCflFactors();
    }



    template <int dim, class RPImpl, class RockType>
    void ReservoirPropertyCommon<dim, RPImpl, RockType>::init(const int num_cells,
                                                              const double uniform_poro,
                                                              const double uniform_perm)
    {
        permfield_valid_.assign(num_cells, std::vector<unsigned char>::value_type(1));
        porosity_.assign(num_cells, uniform_poro);
        permeability_.assign(dim*dim*num_cells, 0.0);
        for (int i = 0; i < num_cells; ++i) {
            SharedPermTensor K = permeabilityModifiable(i);
            for (int dd = 0; dd < dim; ++dd) {
                K(dd, dd) = uniform_perm;
            }
        }
        cell_to_rock_.assign(num_cells, 0);
        asImpl().computeCflFactors();
    }

    template <int dim, class RPImpl, class RockType>
    void ReservoirPropertyCommon<dim, RPImpl, RockType>::setViscosities(double v1, double v2)
    {
        viscosity1_ = v1;
        viscosity2_ = v2;
    }

    template <int dim, class RPImpl, class RockType>
    void ReservoirPropertyCommon<dim, RPImpl, RockType>::setDensities(double d1, double d2)
    {
        density1_ = d1;
        density2_ = d2;
    }

    template <int dim, class RPImpl, class RockType>
    double ReservoirPropertyCommon<dim, RPImpl, RockType>::viscosityFirstPhase() const
    {
        return viscosity1_;
    }


    template <int dim, class RPImpl, class RockType>
    double ReservoirPropertyCommon<dim, RPImpl, RockType>::viscositySecondPhase() const
    {
        return viscosity2_;
    }


    template <int dim, class RPImpl, class RockType>
    double ReservoirPropertyCommon<dim, RPImpl, RockType>::densityFirstPhase() const
    {
        return density1_;
    }


    template <int dim, class RPImpl, class RockType>
    double ReservoirPropertyCommon<dim, RPImpl, RockType>::densitySecondPhase() const
    {
        return density2_;
    }


    template <int dim, class RPImpl, class RockType>
    double ReservoirPropertyCommon<dim, RPImpl, RockType>::porosity(int cell_index) const
    {
        return porosity_[cell_index];
    }


    template <int dim, class RPImpl, class RockType>
    typename ReservoirPropertyCommon<dim, RPImpl, RockType>::PermTensor
    ReservoirPropertyCommon<dim, RPImpl, RockType>::permeability(int cell_index) const
    {
        ASSERT (permfield_valid_[cell_index]);

        const PermTensor K(dim, dim, &permeability_[dim*dim*cell_index]);
        return K;
    }


    template <int dim, class RPImpl, class RockType>
    typename ReservoirPropertyCommon<dim, RPImpl, RockType>::SharedPermTensor
    ReservoirPropertyCommon<dim, RPImpl, RockType>::permeabilityModifiable(int cell_index)
    {
        // Typically only used for assigning synthetic perm values.
        SharedPermTensor K(dim, dim, &permeability_[dim*dim*cell_index]);

        // Trust caller!
        permfield_valid_[cell_index] = std::vector<unsigned char>::value_type(1);

        return K;
    }


    template <int dim, class RPImpl, class RockType>
    template<class Vector>
    void ReservoirPropertyCommon<dim, RPImpl, RockType>::phaseDensities(int /*cell_index*/, Vector& density) const
    {
        ASSERT (density.size() >= NumberOfPhases);
        density[0] = densityFirstPhase();
        density[1] = densitySecondPhase();
    }


    template <int dim, class RPImpl, class RockType>
    double ReservoirPropertyCommon<dim, RPImpl, RockType>::densityDifference() const
    {
        return density1_ - density2_;
    }


    template <int dim, class RPImpl, class RockType>
    double ReservoirPropertyCommon<dim, RPImpl, RockType>::cflFactor() const
    {
        return cfl_factor_;
    }


    template <int dim, class RPImpl, class RockType>
    double ReservoirPropertyCommon<dim, RPImpl, RockType>::cflFactorGravity() const
    {
        return cfl_factor_gravity_;
    }


    template <int dim, class RPImpl, class RockType>
    double ReservoirPropertyCommon<dim, RPImpl, RockType>::capillaryPressure(int cell_index, double saturation) const
    {
        if (rock_.size() > 0) {
            int r = cell_to_rock_[cell_index];
            return rock_[r].capPress(permeability(cell_index), porosity(cell_index), saturation);
        } else {
            // HACK ALERT!
            // Use zero capillary pressure if no known rock table exists.
            return 0.0;
        }
    }


    template <int dim, class RPImpl, class RockType>
    void ReservoirPropertyCommon<dim, RPImpl, RockType>::writeSintefLegacyFormat(const std::string& grid_prefix) const
    {
        int num_cells = porosity_.size();
        // Write porosity.
	{
            std::string filename = grid_prefix + "-poro.dat";
	    std::ofstream file(filename.c_str());
	    if (!file) {
		THROW("Could not open file " << filename);
	    }
            file << num_cells << '\n';
            std::copy(porosity_.begin(), porosity_.end(), std::ostream_iterator<double>(file, "\n"));
	}
        // Write permeability.
	{
            std::string filename = grid_prefix + "-perm.dat";
	    std::ofstream file(filename.c_str());
	    if (!file) {
		THROW("Could not open file " << filename);
	    }
            file << num_cells << '\n';
            switch (permeability_kind_) {
            case TensorPerm:
                std::copy(permeability_.begin(), permeability_.end(), std::ostream_iterator<double>(file, "\n"));
                break;
            case DiagonalPerm:
                for (int c = 0; c < num_cells; ++c) {
                    int index = c*dim*dim;
                    for (int dd = 0; dd < dim; ++dd) {
                        file << permeability_[index + (dim + 1)*dd] << ' ';
                    }
                    file << '\n';
                }
                break;
            case ScalarPerm:
            case None: // Treated like a scalar permeability.
                for (int c = 0; c < num_cells; ++c) {
                    file << permeability_[c*dim*dim] << '\n';
                }
                break;
            default:
                THROW("Cannot write invalid permeability.");
            }
	}
    }


    // ------ Private methods ------


    template <int dim, class RPImpl, class RockType>
    RPImpl& ReservoirPropertyCommon<dim, RPImpl, RockType>::asImpl()
    {
	return static_cast<RPImpl&>(*this);
    }




    template <int dim, class RPImpl, class RockType>
    void ReservoirPropertyCommon<dim, RPImpl, RockType>::assignPorosity(const EclipseGridParser& parser,
                                                         const std::vector<int>& global_cell)
    {
        porosity_.assign(global_cell.size(), 1.0);

        if (parser.hasField("PORO")) {
            const std::vector<double>& poro = parser.getFloatingPointValue("PORO");

            for (int c = 0; c < int(porosity_.size()); ++c) {
                porosity_[c] = poro[global_cell[c]];
            }
        }
    }




    template <int dim, class RPImpl, class RockType>
    void ReservoirPropertyCommon<dim, RPImpl, RockType>::assignPermeability(const EclipseGridParser& parser,
                                                                            const std::vector<int>& global_cell,
                                                                            double perm_threshold)
    {
        EclipseGridInspector insp(parser);
        boost::array<int, 3> dims = insp.gridSize();
        int num_global_cells = dims[0]*dims[1]*dims[2];
        ASSERT (num_global_cells > 0);

        permeability_.assign(dim * dim * global_cell.size(), 0.0);

        std::vector<const std::vector<double>*> tensor;
        tensor.reserve(10);

        const std::vector<double> zero(num_global_cells, 0.0);
        tensor.push_back(&zero);

        BOOST_STATIC_ASSERT(dim == 3);
        boost::array<int,9> kmap;
        permeability_kind_ = fillTensor(parser, tensor, kmap);

        // Assign permeability values only if such values are
        // given in the input deck represented by 'parser'.  In
        // other words: Don't set any (arbitrary) default values.
        // It is infinitely better to experience a reproducible
        // crash than subtle errors resulting from a (poorly
        // chosen) default value...
        //
        if (tensor.size() > 1) {
            const int nc  = global_cell.size();
            int       off = 0;

            for (int c = 0; c < nc; ++c, off += dim*dim) {
                SharedPermTensor K(dim, dim, &permeability_[off]);
                int       kix  = 0;
                const int glob = global_cell[c];

                for (int i = 0; i < dim; ++i) {
                    for (int j = 0; j < dim; ++j, ++kix) {
                        K(i,j) = unit::convert::from((*tensor[kmap[kix]])[glob],
                                                     prefix::milli*unit::darcy);
                    }
                    K(i,i) = std::max(K(i,i), perm_threshold);
                }

                permfield_valid_[c] = std::vector<unsigned char>::value_type(1);
            }
        }
    }




    template <int dim, class RPImpl, class RockType>
    void ReservoirPropertyCommon<dim, RPImpl, RockType>::assignRockTable(const EclipseGridParser& parser,
                                                          const std::vector<int>& global_cell)
    {
        const int nc = global_cell.size();

        cell_to_rock_.assign(nc, 0);

        if (parser.hasField("SATNUM")) {
            const std::vector<int>& satnum = parser.getIntegerValue("SATNUM");

            for (int c = 0; c < nc; ++c) {
                // Note: SATNUM is FORTRANish, ranging from 1 to n, therefore we subtract one.
                cell_to_rock_[c] = satnum[global_cell[c]] - 1;
            }
        }
    }




    template <int dim, class RPImpl, class RockType>
    void ReservoirPropertyCommon<dim, RPImpl, RockType>::readRocks(const std::string& rock_list_file)
    {
        std::ifstream rl(rock_list_file.c_str());
        if (!rl) {
            THROW("Could not open file " << rock_list_file);
        }
        int num_rocks = -1;
        rl >> num_rocks;
        ASSERT(num_rocks >= 1);
	rock_.resize(num_rocks);
	std::string dir(rock_list_file.begin(), rock_list_file.begin() + rock_list_file.find_last_of('/') + 1);
        for (int i = 0; i < num_rocks; ++i) {
            std::string spec;
	    while (spec.empty()) {
		std::getline(rl, spec);
	    }
            rock_[i].read(dir, spec);
        }
    }




} // namespace Dune


#endif // OPENRS_RESERVOIRPROPERTYCOMMON_IMPL_HEADER
