//===========================================================================
//
// File: DefaultGeometryPolicy.hpp
//
// Created: Tue Jun  2 16:23:01 2009
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
Copyright 2009, 2010, 2022 Equinor ASA.

This file is part of The Open Porous Media project  (OPM).

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

#ifndef OPM_DEFAULTGEOMETRYPOLICY_HEADER
#define OPM_DEFAULTGEOMETRYPOLICY_HEADER

#include "EntityRep.hpp"

namespace Dune
{

    class CpGrid;

    namespace cpgrid
    {
        template<int mydim, int dim>
        class Geometry;
        /// @brief
        /// @todo Doc me!
        class DefaultGeometryPolicy
        {
            friend class CpGridData;
            template<int mydim, int dim>
            friend class Geometry;
            friend class ::Dune::CpGrid;
        public:
            /// @brief
            /// @todo Doc me
            DefaultGeometryPolicy()
            {
            }

            /// @brief
            /// @todo Doc me
            /// @param
            DefaultGeometryPolicy(const EntityVariable<cpgrid::Geometry<3, 3>, 0>& cell_geom,
                                  const EntityVariable<cpgrid::Geometry<2, 3>, 1>& face_geom,
                                  const EntityVariable<cpgrid::Geometry<0, 3>, 3>& point_geom)
                : cell_geom_(cell_geom), face_geom_(face_geom), point_geom_(point_geom)
            {
            }

            /// @brief
            /// @todo Doc me!
            /// @tparam
            /// @param
            /// @return
            template <int codim>
            const EntityVariable<cpgrid::Geometry<3 - codim, 3>, codim>& geomVector() const
            {
                static_assert(codim != 2, "");
                return geomVector(std::integral_constant<int,codim>());
            }
            
            /// \brief Get cell geometry
            const EntityVariable<cpgrid::Geometry<3, 3>, 0>& geomVector(const std::integral_constant<int, 0>&) const
            {
                return cell_geom_;
            }
            /// \brief Get cell geometry
            EntityVariable<cpgrid::Geometry<3, 3>, 0>& geomVector(const std::integral_constant<int, 0>&)
            {
                return cell_geom_;
            }
            /// \brief Get face geometry
            const EntityVariable<cpgrid::Geometry<2, 3>, 1>& geomVector(const std::integral_constant<int, 1>&) const
            {
                return face_geom_;
            }
            /// \brief Get face geometry
            EntityVariable<cpgrid::Geometry<2, 3>, 1>& geomVector(const std::integral_constant<int, 1>&)
            {
                return face_geom_;
            }

            /// \brief Get point geometry
            template<int codim>
            const EntityVariable<cpgrid::Geometry<0, 3>, 3>& geomVector(const std::integral_constant<int, codim>&) const
            {
                static_assert(codim==3, "Codim has to be 3");
                return point_geom_;
            }/// \brief Get point geometry
            template<int codim>
            EntityVariable<cpgrid::Geometry<0, 3>, 3>& geomVector(const std::integral_constant<int, codim>&)
            {
                static_assert(codim==3, "Codim has to be 3");
                return point_geom_;
            }
            
        private:
            EntityVariable<cpgrid::Geometry<3, 3>, 0> cell_geom_;
            EntityVariable<cpgrid::Geometry<2, 3>, 1> face_geom_;
            EntityVariable<cpgrid::Geometry<0, 3>, 3> point_geom_;
        };



    } // namespace cpgrid
} // namespace Dune


#endif // OPM_DEFAULTGEOMETRYPOLICY_HEADER
