/*===========================================================================
//
// File: facetopology.h
//
// Created: Fri Jun 19 08:47:10 2009
//
// Author: Jostein R. Natvig <Jostein.R.Natvig@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================*/

/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.

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

#ifndef OPM_FACETOPOLOGY_HEADER
#define OPM_FACETOPOLOGY_HEADER

#include <stdbool.h>

struct processed_grid;

void findconnections(bool edge_conformal,
                     int n,
                     int *pts[4],
                     int *intersectionlist,
                     int *work,
                     struct processed_grid *out);

#endif /* OPM_FACETOPOLOGY_HEADER */

/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
