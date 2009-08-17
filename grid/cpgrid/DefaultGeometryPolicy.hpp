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

#ifndef OPENRS_DEFAULTGEOMETRYPOLICY_HEADER
#define OPENRS_DEFAULTGEOMETRYPOLICY_HEADER

#include <boost/mpl/if.hpp>
#include "Geometry.hpp"
#include "Entity.hpp"

namespace Dune
{
    namespace cpgrid
    {

	struct GetCellGeom;
	struct GetFaceGeom;
	struct GetPointGeom;

	/// @brief
	/// @todo Doc me!
	class DefaultGeometryPolicy
	{
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
		BOOST_STATIC_ASSERT(codim != 2);
		typedef typename boost::mpl::if_c<codim == 0, GetCellGeom, 
		    typename boost::mpl::if_c<codim == 1, GetFaceGeom, GetPointGeom>::type >::type selector;
		return selector::value(*this);
	    }
	private:
	    friend class GetCellGeom;
	    friend class GetFaceGeom;
	    friend class GetPointGeom;
	    EntityVariable<cpgrid::Geometry<3, 3>, 0> cell_geom_;
	    EntityVariable<cpgrid::Geometry<2, 3>, 1> face_geom_;
	    EntityVariable<cpgrid::Geometry<0, 3>, 3> point_geom_;
	};

	/// @brief
	/// @todo Doc me!
	struct GetCellGeom
	{
	    
	    /// @brief
	    /// @todo Doc me!
	    /// @tparam
	    /// @param
	    /// @return
	    static const EntityVariable<cpgrid::Geometry<3, 3>, 0>&
	    value(const DefaultGeometryPolicy& geom)
	    {
		return geom.cell_geom_;
	    }
	};

	/// @brief
	/// @todo Doc me!
	struct GetFaceGeom
	{
	    /// @brief
	    /// @todo Doc me!
	    /// @tparam
	    /// @param
	    /// @return
	    static const EntityVariable<cpgrid::Geometry<2, 3>, 1>&
	    value(const DefaultGeometryPolicy& geom)
	    {
		return geom.face_geom_;
	    }
	};

	/// @brief
	/// @todo Doc me!
	struct GetPointGeom
	{
	    /// @brief
	    /// @todo Doc me!
	    /// @tparam
	    /// @param
	    /// @return 
	    static const EntityVariable<cpgrid::Geometry<0, 3>, 3>&
	    value(const DefaultGeometryPolicy& geom)
	    {
		return geom.point_geom_;
	    }
	};



    } // namespace cpgrid
} // namespace Dune


#endif // OPENRS_DEFAULTGEOMETRYPOLICY_HEADER
