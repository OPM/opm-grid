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

#include <config.h>

#include <opm/grid/cpgrid/DefaultGeometryPolicy.hpp>
#include <opm/grid/cpgrid/Geometry.hpp>

#include <memory>

namespace Dune::cpgrid
{

DefaultGeometryPolicy::DefaultGeometryPolicy()
    : cell_geom_ptr_(std::make_shared<EntityVariable<Geometry<3, 3>, 0>>()),
      face_geom_ptr_(std::make_shared<EntityVariable<Geometry<2, 3>, 1>>()),
      point_geom_ptr_(std::make_shared<EntityVariable<Geometry<0, 3>, 3>>())
{
}

DefaultGeometryPolicy::DefaultGeometryPolicy(const EntityVariable<Geometry<3, 3>, 0>& cell_geom,
                                             const EntityVariable<Geometry<2, 3>, 1>& face_geom,
                                             const EntityVariable<Geometry<0, 3>, 3>& point_geom)
    : cell_geom_ptr_(std::make_shared<EntityVariable<Geometry<3, 3>, 0>>(cell_geom)),
      face_geom_ptr_(std::make_shared<EntityVariable<Geometry<2, 3>, 1>>(face_geom)),
      point_geom_ptr_(std::make_shared<EntityVariable<Geometry<0, 3>, 3>>(point_geom))
{
}

} // namespace Dune::cpgrid
