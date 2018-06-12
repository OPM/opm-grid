/*
  Copyright 2017 IRIS AS

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
#ifndef OPM_GRID_OPMPARSERINCLUDES_HEADER_INCLUDED
#define OPM_GRID_OPMPARSERINCLUDES_HEADER_INCLUDED

#if HAVE_ECL_INPUT
#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Deck/DeckItem.hpp>
#include <opm/parser/eclipse/Deck/DeckKeyword.hpp>
#include <opm/parser/eclipse/Deck/DeckRecord.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/ConnectionSet.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well.hpp>

namespace Dune {
    namespace cpgrid {
        typedef Opm::Well OpmWellType;
        typedef Opm::EclipseState OpmEclipseStateType;
    }
}
#else // #if HAVE_ECL_INPUT

namespace Dune {
    namespace cpgrid {
        typedef int OpmWellType;
        typedef int OpmEclipseStateType;
    }
}

#endif // #if HAVE_ECL_INPUT

#endif // #ifndef OPM_GRID_OPMPARSERINCLUDES_HEADER_INCLUDED
