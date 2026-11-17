/*
  Copyright 2026 Equinor ASA

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
#ifndef DUNE_POLYHEDRALGRID_ENTITY2INDEXDATAHANDLE_HH
#define DUNE_POLYHEDRALGRID_ENTITY2INDEXDATAHANDLE_HH

#include <cstddef>
#include <iterator>
#include <vector>

namespace Dune
{

  /// \brief Buffer for the data of one entity, as used by the data handles
  ///        when data is moved between the global and the distributed view.
  template< class T >
  class PolyhedralGridMoveBuffer
  {
  public:
    typedef typename std::vector< T >::const_iterator const_iterator;

    /// An empty buffer that can hold the given number of values to write.
    explicit PolyhedralGridMoveBuffer ( const std::size_t size )
      : buffer_( size ), position_( 0 )
    {}

    /// A buffer to read the values in the given range from.
    template< class Iterator >
    PolyhedralGridMoveBuffer ( Iterator begin, Iterator end )
      : buffer_( begin, end ), position_( 0 )
    {}

    void read ( T& value ) { value = buffer_[ position_++ ]; }
    void write ( const T& value ) { buffer_[ position_++ ] = value; }

    const_iterator begin () const { return buffer_.begin(); }
    const_iterator end () const { return buffer_.end(); }

  private:
    std::vector< T > buffer_;
    std::size_t position_;
  };

  /// \brief Adapts a Dune data handle, which works on entities, to the
  ///        index based interface of VariableSizeCommunicator.
  ///
  /// \tparam Grid        The grid.
  /// \tparam DataHandle  A Dune::CommDataHandleIF compatible data handle.
  /// \tparam codim       The codimension of the entities that are communicated.
  template< class Grid, class DataHandle, int codim >
  class PolyhedralGridEntity2IndexDataHandle
  {
  public:
    typedef typename DataHandle::DataType DataType;

    PolyhedralGridEntity2IndexDataHandle ( const Grid& grid, DataHandle& data )
      : grid_( grid ), data_( data )
    {}

    bool fixedSize ()
    {
      return data_.fixedSize( Grid::dimension, codim );
    }

    std::size_t size ( const std::size_t i )
    {
      return data_.size( entity( i ) );
    }

    template< class B >
    void gather ( B& buffer, const std::size_t i )
    {
      data_.gather( buffer, entity( i ) );
    }

    template< class B >
    void scatter ( B& buffer, const std::size_t i, const std::size_t s )
    {
      data_.scatter( buffer, entity( i ), s );
    }

  private:
    typename Grid::template Codim< codim >::Entity entity ( const std::size_t i ) const
    {
      typedef typename Grid::template Codim< codim >::EntitySeed EntitySeed;
      return grid_.entity( EntitySeed( static_cast< int >( i ) ) );
    }

    const Grid& grid_;
    DataHandle& data_;
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_ENTITY2INDEXDATAHANDLE_HH
