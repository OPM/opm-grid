#ifndef DUNE_POLYHEDRALGRID_ITERATOR_HH
#define DUNE_POLYHEDRALGRID_ITERATOR_HH

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/entityiterator.hh>

#include <dune/grid/polyhedralgrid/entitypointer.hh>

namespace Dune
{

  // PolyhedralGridIterator
  // --------------

  template< class Grid, class HostIterator >
  class PolyhedralGridIterator
  : public PolyhedralGridEntityPointer< Grid, HostIterator >
  {
    typedef PolyhedralGridIterator< Grid, HostIterator > This;
    typedef PolyhedralGridEntityPointer< Grid, HostIterator > Base;

  protected:
    using Base::hostIterator_;
    using Base::releaseEntity;

    typedef typename Base::ExtraData ExtraData;

  public:
    PolyhedralGridIterator ( ExtraData data, const HostIterator &hostIterator )
    : Base( data, hostIterator )
    {}

    /** \brief increment */
    void increment ()
    {
      ++hostIterator_;
      releaseEntity();
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_ITERATOR_HH
