//===========================================================================
//
// File: fortran.hpp
//
// Created: Sun Jun 21 18:50:37 2009
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

#ifndef OPENRS_FORTRAN_HEADER
#define OPENRS_FORTRAN_HEADER

// Fortran name mangling
#ifdef F77_MANGLE
#undef F77_MANGLE
#endif
#define F77_MANGLE(NAME)   NAME ## _

#ifdef F77_NAME
#undef F77_NAME
#endif
#define F77_NAME(lcase,UCASE) F77_MANGLE(lcase)

#ifdef F77_CHARACTER_TYPE
#undef F77_CHARACTER_TYPE
#define F77_CHARACTER_TYPE const char*, int

#ifdef F77_CHARACTER
#undef F77_CHARACTER
#define F77_CHARACTER(c) &c, 1

#endif // OPENRS_FORTRAN_HEADER
