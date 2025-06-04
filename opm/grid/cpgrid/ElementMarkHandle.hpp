//===========================================================================
//
// File: ElementMarkHandle.hpp
//
// Created: June 2 2025
//
// Author(s): Antonella Ritorto <antonella.ritorto@opm-op.com>
//            Markus Blatt      <markus.blatt@opm-op.com>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2025 Equinor ASA
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

#ifndef OPM_ELEMENTMARKHANDLE_HEADER
#define OPM_ELEMENTMARKHANDLE_HEADER

#include <opm/grid/cpgrid/Entity.hpp>

namespace
{
/// \brief Handle for assignment of element marks (0: do nothing, 1: refine, -1: coarse - not supported yet).
struct ElementMarkHandle {

    using DataType = int;

    ElementMarkHandle(std::vector<DataType>& winningMark)
        : winningMark_(winningMark)
    {}

    bool fixedSize(std::size_t, std::size_t)
    {
        // For each element, gather/scatter its mark (0, 1, or -1).
        return true;
    }

    bool contains(std::size_t, std::size_t codim)
    {
        // Only communicate values attached to cells.
        return codim == 0;
    }

    template <class T> // T = Entity<0>
    std::size_t size([[maybe_unused]] const T& element)
    {
        return 1;
    }

    // Gather element mark (0: do nothing, 1: refine, -1 coarse - not supported yet)
    template <class B, class T> // T = Entity<0>
    void gather(B& buffer, const T& element)
    {
        buffer.write( winningMark_[element.index()] );
    }

    // Scatter element mark. Rewrite mark to the maximum.
    template <class B, class T> // T = Entity<0>
    void scatter(B& buffer, const T& element, [[maybe_unused]] std::size_t size)
    {
        DataType tmp_mark;
        buffer.read(tmp_mark);
        winningMark_[element.index()] = std::max(winningMark_[element.index()], tmp_mark);
    }

private:
    std::vector<DataType>& winningMark_;
};
} // namespace
#endif
