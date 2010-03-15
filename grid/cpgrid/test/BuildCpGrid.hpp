//===========================================================================
//
// File: BuildCpGrid.hpp
//
// Created: Mon Aug 24 14:12:25 2009
//
// Author(s): Bård Skaflestad     <bard.skaflestad@sintef.no>
//            Atgeirr F Rasmussen <atgeirr@sintef.no>
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

#ifndef OPENRS_BUILDCPGRID_HEADER
#define OPENRS_BUILDCPGRID_HEADER

#include <algorithm>
#include <cmath>
#include <fstream>
#include <ostream>
#include <string>
#include <vector>

#include <dune/common/ErrorMacros.hpp>
#include <dune/grid/preprocess/preprocess.h>
#include <dune/grid/CpGrid.hpp>

namespace Dune {
    class CoordinateArray {
    public:
        CoordinateArray(int nx, int ny, int nz)
            : nx_(nx), ny_(ny), nz_(nz),
              data_((nx+1) * (ny+1) * (nz+1))
        {}
        double& operator()(int i, int j, int k)
        {
            ASSERT ((0 <= i) && (i <= nx_));
            ASSERT ((0 <= j) && (j <= ny_));
            ASSERT ((0 <= k) && (k <= nz_));
            return data_[i + (nx_+1)*(j + (ny_+1)*k)];
        }
    private:
        int nx_, ny_, nz_;
        std::vector<double> data_;
    };


    class BuildCpGrid {
    public:
        BuildCpGrid(int nx, int ny, int nz)
            : nx_(nx), ny_(ny), nz_(nz),
              X_ (nx,      ny,      nz),
              Y_ (nx,      ny,      nz),
              Z_ (nx,      ny,      nz)
        {}

        void build(Dune::CpGrid& g,
                   double        z_tol = 0.0,
                   bool          write_grdecl = true)
        {
            std::vector<double> coord;
            this->coord(coord);

            std::vector<double> zcorn;
            this->zcorn(zcorn);
            setZCORN(zcorn);

            std::vector<int> actnum(Nx() * Ny() * Nz(), 1);
            setACTNUM(actnum);

            if (write_grdecl) {
                std::ofstream GRDECL(fileName().c_str());
                writeGRDECL(GRDECL, coord, zcorn, actnum);
            }

            struct grdecl grdecl;
            grdecl.dims[0] = Nx();
            grdecl.dims[1] = Ny();
            grdecl.dims[2] = Nz();
            grdecl.coord   = &coord [0];
            grdecl.zcorn   = &zcorn [0];
            grdecl.actnum  = &actnum[0];

            g.processEclipseFormat(grdecl, z_tol);
        }

    protected:
        int Nx() const { return nx_; }
        int Ny() const { return ny_; }
        int Nz() const { return nz_; }

        double& X(int i, int j, int k) { return X_(i, j, k); }
        double& Y(int i, int j, int k) { return Y_(i, j, k); }
        double& Z(int i, int j, int k) { return Z_(i, j, k); }

        virtual void        setACTNUM(std::vector<int>& /* actnum */)
        {
            // All active by default
        }

        virtual void        setZCORN (std::vector<double>& zcorn) = 0;
        virtual std::string fileName ()                           = 0;

    private:
        const int nx_, ny_, nz_;
        CoordinateArray X_, Y_, Z_;

        void coord(std::vector<double>& coord)
        {
            coord.clear();
            coord.resize(6 * (nx_+1) * (ny_+1));

            double* c = &coord[0];
            for (int j = 0; j < ny_+1; ++j)
            for (int i = 0; i < nx_+1; ++i) {
                *c++ = X_(i, j,  0 ); // X-top
                *c++ = Y_(i, j,  0 ); // Y-top
                *c++ = Z_(i, j,  0 ); // Z-top
                *c++ = X_(i, j, nz_); // X-bot
                *c++ = Y_(i, j, nz_); // Y-bot
                *c++ = Z_(i, j, nz_); // Z-bot
            }
        }

        void zcorn(std::vector<double>& zcorn)
        {
            zcorn.clear();
            zcorn.resize((2*nx_) * (2*ny_) * (2*nz_));

            double* z = &zcorn[0];
            for (int k = 1; k < 2*nz_ + 1; ++k)
            for (int j = 1; j < 2*ny_ + 1; ++j)
            for (int i = 1; i < 2*nx_ + 1; ++i) {
                *z++ = Z_(i/2, j/2, k/2);
            }
        }

        template<typename charT, class traits>
        void writeGRDECL(std::basic_ostream<charT,traits>& grdecl,
                         const std::vector<double>&        coord ,
                         const std::vector<double>&        zcorn ,
                         const std::vector<int>&           actnum)
        {
            ASSERT (int(coord .size()) == 6 * (nx_+1) * (ny_+1));
            ASSERT (int(zcorn .size()) == (2*nx_) * (2*ny_) * (2*nz_));
            ASSERT (int(actnum.size()) == nx_ * ny_ * nz_);

            grdecl.precision(15);
            grdecl << "SPECGRID\n"
                   << nx_ << ' '
                   << ny_ << ' '
                   << nz_ << " 1 F\n/\n\n"
                   << "COORD\n";
            grdecl.setf(std::ios::scientific | std::ios::showpos);
            for (int i = 0; i < 6 * (nx_ + 1) * (ny_ + 1); ++i) {
                grdecl << coord[i] << (((i + 1) % 6 == 0) ? "\n" : " ");
            }
            grdecl << "/\n\nZCORN\n";
            for (int i = 0; i < (2*nx_) * (2*ny_) * (2*nz_); ++i) {
                grdecl << zcorn[i] << (((i + 1) % 8 == 0) ? "\n" : " ");
            }
            grdecl << "/\n\nACTNUM\n";
            grdecl.unsetf(std::ios::scientific | std::ios::showpos);
            for (int i = 0; i < nx_ * ny_ * nz_; ++i) {
                grdecl << actnum[i] << (((i + 1) % 8 == 0) ? "\n" : " ");
            }
            grdecl << "/\n";
        }
    };


    class SimpleFault : public BuildCpGrid {
    public:
        SimpleFault(int    nx, int    ny, int    nz,
                    double hx, double hy, double hz, double drop)
            : BuildCpGrid(nx, ny, nz),
              hx_(hx), hy_(hy), hz_(hz), drop_(drop)
        {
            const double pi = 3.14159265358979323846264338327950288;
            for (int k = 0; k < nz+1; ++k) {
                const double zeta = double(k) / (nz + 1);
                for (int j = 0; j < ny+1; ++j) {
                    const double eta = double(j) / (ny + 1);
                    for (int i = 0; i < nx+1; ++i) {
                        const double xi = double(i) / (nx + 1);

                        // Make box.
                        X(i,j,k) = xi   * hx_;
                        Y(i,j,k) = eta  * hy_;
                        Z(i,j,k) = zeta * hz_;

                        // Add perturbation.
                        X(i,j,k) += 0.2 * (0.5-std::abs(xi-0.5)) * (zeta - 0.5) * hx_;

                        const double x = X(i,j,k) / hx_;
                        Z(i,j,k) -= (0.050*std::sin(pi * ( x  - 0.5)) +
                                     0.075*std::sin(pi * (eta + 2*x))) * hz_;
                    }
                }
            }

            std::vector<double> z_srt(nz + 1);
            for (int j = 0; j < ny+1; ++j) {
                for (int i = 0; i < nx+1; ++i) {
                    for (int k = 0; k < nz+1; ++k) {
                        z_srt[k] = Z(i,j,k);
                    }

                    std::sort(z_srt.begin(), z_srt.end());

                    for (int k = 0; k < nz+1; ++k) {
                        Z(i,j,k) = z_srt[k];
                    }
                }
            }
        }


    private:
        const double hx_, hy_, hz_;
        const double drop_;

        void setZCORN(std::vector<double>& zcorn)
        {
            const int nx = Nx(), ny = Ny(), nz = Nz();
            const int imin = 2 * (nx / 2);

            for (int k = 0   ; k < 2*nz; ++k)
            for (int j = 0   ; j < 2*ny; ++j)
            for (int i = imin; i < 2*nx; ++i) {
                zcorn[i + 2*nx*(j + 2*ny*k)] += drop_;
            }
        }

        std::string fileName() { return "SimpleFault.grdecl"; }
    };


    class SlopingFault : public BuildCpGrid {
    public:
        SlopingFault(int    nx, int    ny, int    nz,
                     double hx, double hy, double hz, double drop)
            : BuildCpGrid(nx, ny, nz),
              hx_(hx), hy_(hy), hz_(hz), drop_(drop)
        {
            for (int k = 0; k < nz+1; ++k) {
                const double zeta = double(k) / (nz + 1);
                for (int j = 0; j < ny+1; ++j) {
                    const double eta = double(j) / (ny + 1);
                    for (int i = 0; i < nx+1; ++i) {
                        const double xi = double(i) / (nx + 1);

                        // Make box.
                        X(i,j,k) = xi   * hx_;
                        Y(i,j,k) = eta  * hy_;
                        Z(i,j,k) = zeta * hz_;

                        // Add perturbation.
                        X(i,j,k) += 0.2 * (0.5-std::abs(xi-0.5)) * (zeta-0.5) * hx_;
                    }
                }
            }
        }

    private:
        const double hx_, hy_, hz_;
        const double drop_;

        void setZCORN(std::vector<double>& zcorn)
        {
            const int nx = Nx(), ny = Ny(), nz = Nz();
            const int imin = 2 * (nx / 2);

            for (int k = 0   ; k < 2*nz; ++k)
            for (int j = 0   ; j < 2*ny; ++j)
            for (int i = imin; i < 2*nx; ++i) {
                zcorn[i + 2*nx*(j + 2*ny*k)] += drop_;
            }
        }

        std::string fileName() { return "SlopingFault.grdecl"; }
    };

} // namespace Dune
#endif // OPENRS_BUILDCPGRID_HEADER
