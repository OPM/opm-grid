/*
Copyright 2014 Statoil ASA.
Copyright 2014 Dr. Markus Blatt - HPC-Simulation-Software & Services

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
#ifndef OPM_GLOBALIDMAPPING_HEADER
#define OPM_GLOBALIDMAPPING_HEADER
#include <vector>
namespace Dune
{
namespace cpgrid
{
/// \brief Class managing the mappings of local indices to global ids.
///
///
class GlobalIdMapping
{
public:
    /// \brief Swap data for initialization
    /// \param cellMapping A vector with global id of index i at position i.
    /// \param faceMapping A vector with global id of index i at position i.
    /// \param pointMapping A vector with global id of index i at position i.

    void swap(std::vector<int>& cellMapping,
              std::vector<int>& faceMapping,
              std::vector<int>& pointMapping)
    {
        cellMapping_.swap(cellMapping);
        faceMapping_.swap(faceMapping);
        pointMapping_.swap(pointMapping);
    }
    /// \brief Get the vector with the mappings for a codimension
    /// \tparam codim The codimension.
    template<int codim>
    std::vector<int>& getMapping()
    {
        static_assert(codim == 0 || codim == 1 || codim==3,
                      "Mappings only available for codimension 0, 1, and 3");
        if(codim==0)
            return cellMapping_;
        if(codim==1)
            return faceMapping_;
        return pointMapping_;
    }

    /// \brief Get the vector with the mappings for a codimension
    /// \tparam codim The codimension.
    template<int codim>
    const std::vector<int>& getMapping() const
    {
        static_assert(codim == 0 || codim == 1 || codim==3,
                      "Mappings only available for codimension 0, 1, and 3");
        if(codim==0)
            return cellMapping_;
        if(codim==1)
            return faceMapping_;
        return pointMapping_;
    }
protected:
    /// \brief A vector containing the global id of cell with index i at position i.
    std::vector<int> cellMapping_;
    /// \brief A vector containing the global id of face with index i at position i.
    std::vector<int> faceMapping_;
    /// \brief A vector containing the global id of point with index i at position i.
    std::vector<int> pointMapping_;
};
}
}

#endif
