//===========================================================================
//
// File: DataHandleWrappers.hpp
//
// Created: Mon Nov 4 2013
//
// Author(s): Markus Blatt <markus@dr-blatt.de>
//
// $Date$
//
// $Revision$
//
//===========================================================================
/**
  Copyright 2019 Equinor AS.

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
#ifndef OPM_DATAHANDLEWRAPPERS_HEADER
#define OPM_DATAHANDLEWRAPPERS_HEADER

#include <array>
#include <vector>

#include "OrientedEntityTable.hpp"
#include "EntityRep.hpp"

namespace Dune
{
namespace cpgrid
{

/// \brief A data handle to send data attached to faces via cell communication
///
/// With it we can use the cell communication to also communicate data attached
/// to faces.
/// \warning As we send all faces of cell most of the faces will be send twice
/// which will temporarily waste some space and bandwidth.
///
/// \tparam Handle The type of the data handle to wrap. It must gather and scatter
///         only for codim-1 entities.
/// \warning The wrapped data handle must not used the last argument of
/// its' scatter method since this numer will be incorrect in most cases!
template<class Handle>
struct FaceViaCellHandleWrapper
{
    using DataType = typename Handle::DataType;
    using C2FTable = OrientedEntityTable<0, 1>;

    /// \brief Constructs the data handle
    ///
    /// \param handle Handle object wrapped and used for actual gather/scatter
    /// \param c2pGather Table to determine points when gathering
    /// \param c2p Table to determine points when scattering
    FaceViaCellHandleWrapper(Handle& handle,
                             const C2FTable& c2fGather,
                             const C2FTable& c2f)
        : handle_(handle), c2fGather_(c2fGather), c2f_(c2f)
    {}

    bool fixedSize(int, int)
    {
        return false; // as the faces per cell differ
    }
    template<class T>
    typename std::enable_if<T::codimension != 0, std::size_t>::type
    size(const T&)
    {
        OPM_THROW(std::logic_error, "This should never throw! We only know sizes for cells");
        return 1;
    }
    std::size_t size(const EntityRep<0>& t)
    {
        const auto& faces = c2fGather_[t];
        std::size_t size{};
        for (const auto& face : faces)
        {
            size += handle_.size(face);
        }
        return size;
    }
    bool contains(std::size_t dim, std::size_t codim)
    {
        return dim==3 && codim == 0;
    }
    template<class B, class T>
    typename std::enable_if<T::codimension != 0, void>::type
    gather(B&, const T&)
    {
        OPM_THROW(std::logic_error, "This should never throw! We can only gather cell values and indicate that");
    }
    template<class B>
    void gather(B& buffer, const EntityRep<0>& t)
    {
        const auto& faces = c2fGather_[t];
        for (const auto& face : faces)
        {
            handle_.gather(buffer, face);
        }
    }
    template<class B, class T>
    typename std::enable_if<T::codimension != 0, void>::type
    scatter(B&, const T&, std::size_t)
    {
        OPM_THROW(std::logic_error, "This should never throw! We can only gather cell values and indicate that");
    }
    template<class B>
    void scatter(B& buffer, const EntityRep<0>& t, std::size_t)
    {
        const auto& faces = c2f_[t];
        for (const auto& face : faces)
        {
            // Note that the size (last parameter) is not correct here.
            // Therefore this handle needs to know how many data items
            // to expect. Not usable outside of CpGrid.
            handle_.scatter(buffer, face, 1);
        }
    }
private:
    Handle& handle_;
    const C2FTable& c2fGather_, c2f_;
};

struct PointViaCellWarner
{
    static bool printWarn;
    void warn();
};

/// \brief A data handle to send data attached to points via cell communication
///
/// With it we can use the cell communication to also communicate data attached
/// to faces.
/// \warning As we send all points of cell most of the points will be send six times
/// which will temporarily waste some space and bandwidth.
///
/// \tparam Handle The type of the data handle to wrap. It must gather and scatter
///         only for codim-3 entities.
template<class Handle>
struct PointViaCellHandleWrapper : public PointViaCellWarner
{
    using DataType = typename Handle::DataType;
    using C2PTable = std::vector< std::array<int,8> >;

    /// \brief Constructs the data handle
    ///
    /// \param handle Handle object wrapped and used for actual gather/scatter
    /// \param c2pGather Table to determine points when gathering
    /// \param c2p Table to determine points when scattering
    PointViaCellHandleWrapper(Handle& handle,
                             const C2PTable& c2pGather,
                             const C2PTable& c2p)
        : handle_(handle), c2pGather_(c2pGather), c2p_(c2p)
    {}
    bool fixedSize(int i, int j)
    {
        if( ! handle_.fixedSize(i, j))
        {
            this->warn();
        }
        return handle_.fixedSize(i, j);
    }
    template<class T>
    typename std::enable_if<T::codimension != 0, std::size_t>::type
    size(const T&)
    {
        OPM_THROW(std::logic_error, "This should never throw! We only know sizes for cells");
        return 1;
    }
    std::size_t size(const EntityRep<0>& t)
    {
        const auto& points = c2pGather_[t.index()];
        std::size_t size{};
        for (const auto& point : points)
        {
            size += handle_.size(EntityRep<3>(point, true));
        }
        return size;
    }
    bool contains(std::size_t dim, std::size_t codim)
    {
        return dim==3 && codim == 0;
    }
    template<class B, class T>
    typename std::enable_if<T::codimension != 0, void>::type
    gather(B&, const T&)
    {
        OPM_THROW(std::logic_error, "This should never throw! We can only gather cell values and indicate that");
    }
    template<class B>
    void gather(B& buffer, const EntityRep<0>& t)
    {
        const auto& points = c2pGather_[t.index()];
        for (const auto& point : points)
        {
            handle_.gather(buffer, EntityRep<3>(point, true));
        }
    }
    template<class B, class T>
    typename std::enable_if<T::codimension != 0, void>::type
    scatter(B&, const T&, std::size_t)
    {
        OPM_THROW(std::logic_error, "This should never throw! We can only gather cell values and indicate that");
    }
    template<class B>
    void scatter(B& buffer, const EntityRep<0>& t, std::size_t s)
    {
        const auto& points = c2p_[t.index()];
        for (const auto& point : points)
        {
            handle_.scatter(buffer, EntityRep<3>(point, true), s/8);
        }
    }
private:
    Handle& handle_;
    const C2PTable& c2pGather_, c2p_;
};

} // end namespace cpgrid
} // end namespace Dune
#endif
