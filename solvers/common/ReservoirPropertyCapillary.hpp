//===========================================================================
//
// File: ReservoirPropertyCapillary.hpp
//
// Created: Fri Jul  3 12:28:48 2009
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

#ifndef OPENRS_RESERVOIRPROPERTYCAPILLARY_HEADER
#define OPENRS_RESERVOIRPROPERTYCAPILLARY_HEADER

#include <fstream>
#include <boost/static_assert.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/common/Units.hpp>
#include <dune/grid/common/EclipseGridParser.hpp>
#include <dune/grid/common/EclipseGridInspector.hpp>

#include <dune/solvers/common/NonuniformTableLinear.hpp>
#include <dune/solvers/common/Matrix.hpp>


namespace Dune
{

    namespace {

        /// @brief
        ///    Verify that a given permeability specification is sound
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
        ///    Whether or not the input permeability is sound from a
        ///    structural point of view.  If no permeability is
        ///    specified on input, we deem that the tensor is ok.
        bool structurallyReasonableTensor(const EclipseGridParser& parser)
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

            bool ret;
            if (xx || xy || xz ||
                yx || yy || yz ||
                zx || zy || zz) {
                // At least one tensor component specified on input.
                // Verify that any remaining components are OK from a
                // structural point of view.  In particular, there
                // must not be any cross-components (e.g., k_{xy})
                // unless the corresponding diagonal component (e.g.,
                // k_{xx}) is present as well...
                //
                ret =         xx || !(xy || xz || yx || zx) ;
                ret = ret && (yy || !(yx || yz || xy || zy));
                ret = ret && (zz || !(zx || zy || xz || yz));
            } else {
                // No permeability components.  We deem this officially OK!
                ret = true;
            }

            return ret;
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
        ///    permeabilty component in each grid cell.
        ///
        /// @param [out] tensor
        /// @param [out] kmap
        void fillTensor(const EclipseGridParser&                 parser,
                        std::vector<const std::vector<double>*>& tensor,
                        boost::array<int,9>&                     kmap)
        {
            ASSERT (structurallyReasonableTensor(parser));
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
        }
    }

    /// A property class for incompressible two-phase flow.
    /// @brief
    /// @todo Doc me!
    /// @tparam
    template <int dim>
    class ReservoirPropertyCapillary
    {
    public:
        /// @brief
        /// @todo Doc me!
        typedef ImmutableCMatrix PermTensor;
        /// @brief
        /// @todo Doc me!
        typedef OwnCMatrix       MutablePermTensor;
        /// @brief
        /// @todo Doc me!
        typedef SharedCMatrix    SharedPermTensor;

        /// @brief
        /// @todo Doc me!
        enum { NumberOfPhases = 2 };

        /// @brief
        /// @todo Doc me!
        ReservoirPropertyCapillary()
            : density1_  (1013.9*unit::kilogram/unit::cubic(unit::meter)),
              density2_  ( 834.7*unit::kilogram/unit::cubic(unit::meter)),
              viscosity1_(   1.0*prefix::centi*unit::Poise),
              viscosity2_(   3.0*prefix::centi*unit::Poise)
        {
        }

        /// @brief
        /// @todo Doc me!
        /// @param
        void init(const EclipseGridParser& parser,
                  const std::vector<int>& global_cell,
                  const std::string* rock_list_filename = 0)
        {
            BOOST_STATIC_ASSERT(dim == 3);

            permfield_valid_.assign(global_cell.size(),
                                    std::vector<unsigned char>::value_type(0));

            assignPorosity    (parser, global_cell);
            assignPermeability(parser, global_cell);
            assignRockTable   (parser, global_cell);

            if (rock_list_filename) {
                readRocks(*rock_list_filename);
            }

            computeCflFactors();
        }

        /// @brief
        /// @todo Doc me!
        /// @param
        void init(const int num_cells,
                  const double uniform_poro = 0.2,
                  const double uniform_perm = 100.0*prefix::milli*unit::darcy)
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
            computeCflFactors();
        }

        double viscosityFirstPhase() const
        {
            return viscosity1_;
        }

        double viscositySecondPhase() const
        {
            return viscosity2_;
        }

        double densityFirstPhase() const
        {
            return density1_;
        }

        double densitySecondPhase() const
        {
            return density2_;
        }

        /// @brief
        /// @todo Doc me!
        /// @param
        /// @return
        double porosity(int cell_index) const
        {
            return porosity_[cell_index];
        }
        /// @brief
        /// @todo Doc me!
        /// @param
        /// @return
        PermTensor permeability(int cell_index) const
        {
            ASSERT (permfield_valid_[cell_index]);

            const PermTensor K(dim, dim, &permeability_[dim*dim*cell_index]);
            return K;
        }
        /// @brief
        /// @todo Doc me!
        /// @param
        /// @return
        SharedPermTensor permeabilityModifiable(int cell_index)
        {
            // Typically only used for assigning synthetic perm values.
            SharedPermTensor K(dim, dim, &permeability_[dim*dim*cell_index]);

            // Trust caller!
            permfield_valid_[cell_index] = std::vector<unsigned char>::value_type(1);

            return K;
        }
                /// @brief
        /// @todo Doc me!
        /// @param
        /// @return
        double mobilityFirstPhase(int cell_index, double saturation) const
        {
            return relPermFirstPhase(cell_index, saturation) / viscosity1_;
        }
        /// @brief
        /// @todo Doc me!
        /// @param
        /// @return
        double mobilitySecondPhase(int cell_index, double saturation) const
        {
            return relPermSecondPhase(cell_index, saturation) / viscosity2_;
        }
        /// @brief
        /// @todo Doc me!
        /// @param
        /// @return
        double totalMobility(int cell_index, double saturation) const
        {
            double l1 = mobilityFirstPhase(cell_index, saturation);
            double l2 = mobilitySecondPhase(cell_index, saturation);
            return l1 + l2;
        }
        /// @brief
        /// @todo Doc me!
        /// @param
        /// @return
	template <class MatrixType>
        void anisoTotalMobility(int cell_index, double saturation, MatrixType& total_mobility) const
        {
	    double scalar_totmob = totalMobility(cell_index, saturation);
	    // For now, just make an identity tensor.
	    for (int row = 0; row < dim; ++row) {
		for (int col = 0; col < dim; ++col) {
		    total_mobility(row, col) = 0.0;
		}
		total_mobility(row, row) = scalar_totmob;
	    }
        }
        /// @brief
        /// @todo Doc me!
        /// @param
        /// @return
        double fractionalFlow(int cell_index, double saturation) const
        {
            double l1 = mobilityFirstPhase(cell_index, saturation);
            double l2 = mobilitySecondPhase(cell_index, saturation);
            return l1/(l1 + l2);
        }
        /// @brief
        /// @todo Doc me!
        /// @param
        /// @return
        template<class Vector>
        void phaseDensity(int /*cell_index*/, Vector& density) const
        {
            ASSERT (density.size() >= NumberOfPhases);
            density[0] = densityFirstPhase();
            density[1] = densitySecondPhase();
        }
        /// @brief
        /// @todo Doc me!
        /// @param
        template<class Vector>
        void phaseMobility(int cell_index, double sat, Vector& mob) const
        {
            ASSERT (mob.size() >= NumberOfPhases);
            mob[0] = mobilityFirstPhase(cell_index, sat);
            mob[1] = mobilitySecondPhase(cell_index, sat);
        }
        /// @brief
        /// @todo Doc me!
        /// @return
        double densityDifference() const
        {
            return density1_ - density2_;
        }
        /// @brief
        /// @todo Doc me!
        /// @return
        double cflFactor() const
        {
            return cfl_factor_;
        }
        /// @brief
        /// @todo Doc me!
        /// @param
        /// @return
        double cflFactorGravity() const
        {
            return cfl_factor_gravity_;
        }
        /// @brief
        /// @todo Doc me!
        /// @param
        /// @return
        double capillaryPressure(int cell_index, double saturation) const
        {
            if (rock_.size() > 0) {
                // p_c = J\frac{\sigma \cos \theta}{\sqrt{k/\phi}}
                double sigma_cos_theta = 1.0; // An approximation.
                double perm = trace(permeability(cell_index))/double(dim);
                double poro = porosity(cell_index);
                double sqrt_k_phi = std::sqrt(perm/poro);
                int r = cell_to_rock_[cell_index];
                return rock_[r].Jfunc_(saturation)
                    *sigma_cos_theta/sqrt_k_phi;
            } else {
                // HACK ALERT!
                // Use zero capillary pressure if no known rock table exists.
                return 0.0;
            }
        }

    private:
        double relPermFirstPhase(int cell_index, double saturation) const
        {
            if (rock_.size() > 0) {
                const int region = cell_to_rock_[cell_index];
                ASSERT (region < int(rock_.size()));
                return rock_[region].krw_(saturation);
            } else {
                // HACK ALERT!
                // Use quadratic rel-perm if no known rock table exists.
                return saturation * saturation;
            }
        }
        double relPermSecondPhase(int cell_index, double saturation) const
        {
            if (rock_.size() > 0) {
                const int region = cell_to_rock_[cell_index];
                ASSERT (region < int(rock_.size()));
                return rock_[region].kro_(saturation);
            } else {
                // HACK ALERT!
                // Use quadratic rel-perm if no known rock table exists.
                return (1 - saturation) * (1 - saturation);
            }
        }
        void cflMobs(int rock, double s, double& mob_first, double& mob_gravity) const
        {
            if (rock == -1) {
                // No rock dependency, we might just as well use the first cell.
                const int cell_index = 0;
                double l1 = mobilityFirstPhase(cell_index, s);
                double l2 = mobilitySecondPhase(cell_index, s);
                mob_first = l1/(l1 + l2);
                mob_gravity = l1*l2/(l1 + l2);
            } else {
                double krw = rock_[rock].krw_(s);
                double kro = rock_[rock].kro_(s);
                double l1 = krw/viscosity1_;
                double l2 = kro/viscosity2_;
                mob_first = l1/(l1 + l2);
                mob_gravity = l1*l2/(l1 + l2);
            }
            
            // This is a hack for now, we should make this rock-dependant,
            // for the multi-rock case.
            const int cell_index = 0;
            double l1 = mobilityFirstPhase(cell_index, s);
            double l2 = mobilitySecondPhase(cell_index, s);
            mob_first = l1/(l1 + l2);
            mob_gravity = l1*l2/(l1 + l2);
        }
        std::pair<double, double> computeSingleRockCflFactors(int rock) const
        {
            const int N = 257;
            double delta = 1.0/double(N - 1);
            double last_m1, last_mg;
            double max_der1 = -1e100;
            double max_derg = -1e100;
            cflMobs(rock, 0.0, last_m1, last_mg);
            for (int i = 1; i < N; ++i) {
                double s = double(i)*delta;
                double m1, mg;
                cflMobs(rock, s, m1, mg);
                double est_deriv_m1 = std::fabs(m1 - last_m1)/delta;
                double est_deriv_mg = std::fabs(mg - last_mg)/delta;
                max_der1 = std::max(max_der1, est_deriv_m1);
                max_derg = std::max(max_derg, est_deriv_mg);
                last_m1 = m1;
                last_mg = mg;
            }
            return std::make_pair(1.0/max_der1, 1.0/max_derg);
        }
        void computeCflFactors()
        {
            if (rock_.empty()) {
                std::pair<double, double> fac = computeSingleRockCflFactors(-1);
                cfl_factor_ = fac.first;
                cfl_factor_gravity_ = fac.second;
            } else {
                cfl_factor_ = 1e100;
                cfl_factor_gravity_ = 1e100;
                for (int r = 0; r < int(rock_.size()); ++r) {
                    std::pair<double, double> fac = computeSingleRockCflFactors(r);
                    cfl_factor_ = std::min(cfl_factor_, fac.first);
                    cfl_factor_gravity_ = std::min(cfl_factor_gravity_, fac.second);
                }
            }
        }

        void assignPorosity(const EclipseGridParser& parser,
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

        void assignPermeability(const EclipseGridParser& parser,
                                const std::vector<int>& global_cell)
        {
            int num_global_cells = -1;
            if (parser.hasField("SPECGRID")) {
		const SPECGRID& sgr = parser.getSpecGrid();
		num_global_cells = 
		    sgr.dimensions[0] * sgr.dimensions[1] * sgr.dimensions[2];
            }
            ASSERT (num_global_cells > 0);

            permeability_.assign(dim * dim * global_cell.size(), 0.0);

            std::vector<const std::vector<double>*> tensor;
            tensor.reserve(10);

            const std::vector<double> zero(num_global_cells, 0.0);
            tensor.push_back(&zero);

            BOOST_STATIC_ASSERT(dim == 3);
            boost::array<int,9> kmap;
            fillTensor(parser, tensor, kmap);

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
                    }

                    permfield_valid_[c] = std::vector<unsigned char>::value_type(1);
                }
            }
        }

        void assignRockTable(const EclipseGridParser& parser,
                             const std::vector<int>& global_cell)
        {
            const int nc = global_cell.size();

            cell_to_rock_.assign(nc, 0);

            if (parser.hasField("SATNUM")) {
                const std::vector<int>& satnum = parser.getIntegerValue("SATNUM");

                for (int c = 0; c < nc; ++c) {
                    cell_to_rock_[c] = satnum[global_cell[c]] - 1;
                }
            }
        }

        struct Rock {
            typedef utils::NonuniformTableLinear<double> TabFunc;
            TabFunc krw_;
            TabFunc kro_;
            TabFunc Jfunc_;
            TabFunc invJfunc_;
        };

        void readRocks(const std::string& rock_list_file)
        {
            std::ifstream rl(rock_list_file.c_str());
            if (!rl) {
                THROW("Could not open file " << rock_list_file);
            }
            int num_rocks = -1;
            rl >> num_rocks;
            ASSERT(num_rocks >= 1);
            for (int i = 0; i < num_rocks; ++i) {
                std::string rockname;
                rl >> rockname;
                std::string rockfilename = rock_list_file;
                rockfilename.replace(rockfilename.begin() + rockfilename.find_last_of('/') + 1,
                                     rockfilename.end(), rockname);
                std::ifstream rock_stream(rockfilename.c_str());
                if (!rock_stream) {
                    THROW("Could not open file " << rockfilename);
                }
                rock_.push_back(readStatoilFormat(rock_stream));
            }
        }

        Rock readStatoilFormat(std::istream& is)
        {
            std::string firstline;
            std::getline(is, firstline);
            typedef FieldVector<double, 4> Data;
            std::istream_iterator<Data> start(is);
            std::istream_iterator<Data> end;
            std::vector<Data> data(start, end);
            if (!is.eof()) {
                THROW("Reading stopped but we're not at eof: something went wrong reading data.");
            }
            std::vector<double> svals, krw, kro, Jfunc;
            for (int i = 0; i < int(data.size()); ++i) {
                svals.push_back(data[i][0]);
                krw.push_back(data[i][1]);
                kro.push_back(data[i][2]);
                Jfunc.push_back(data[i][3]);
            }
            typedef typename Rock::TabFunc TFun;
            Rock r;
            r.krw_ = TFun(svals, krw);
            r.kro_ = TFun(svals, kro);
            r.Jfunc_ = TFun(svals, Jfunc);
            std::vector<double> invJfunc(Jfunc);
            std::reverse(invJfunc.begin(), invJfunc.end());
            std::vector<double> invsvals(svals);
            std::reverse(invsvals.begin(), invsvals.end());
            r.invJfunc_ = TFun(invJfunc, invsvals);
            return r;
        }

        std::vector<double>        porosity_;
        std::vector<double>        permeability_;
        std::vector<unsigned char> permfield_valid_;

        double density1_;
        double density2_;
        double viscosity1_;
        double viscosity2_;
        double cfl_factor_;
        double cfl_factor_gravity_;
        std::vector<Rock> rock_;
        std::vector<int> cell_to_rock_;
    };



} // namespace Dune


#endif // OPENRS_RESERVOIRPROPERTYCAPILLARY_HEADER
