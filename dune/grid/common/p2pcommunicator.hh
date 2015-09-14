#ifndef DUNE_COMMUNICATOR_HEADER_INCLUDED
#define DUNE_COMMUNICATOR_HEADER_INCLUDED

#include <vector>
#include <set>
#include <map>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parallel/collectivecommunication.hh>

// the following implementation is only available in case MPI is available
#if HAVE_MPI
#include <dune/common/parallel/mpicollectivecommunication.hh>
#endif

namespace Dune
{
  class SimpleMessageBuffer
  {
    typedef std::vector< char >  BufferType;

    mutable BufferType buffer_;
    const double factor_;
    mutable size_t pos_;
public:
    /** \brief constructor taking memory reserve estimation factor (default is 1.1, i.e. 10% over estimation )
     */
    SimpleMessageBuffer( const double factor = 1.1 )
      : buffer_(), factor_( factor )
    {
      resetReadPosition();
    }

    /** \brief clear the buffer */
    void clear() { BufferType().swap( buffer_ ); resetReadPosition(); }
    /** \brief reset read position of buffer to beginning */
    void resetReadPosition() { pos_ = 0 ; }
    /** \brief return size of buffer */
    size_t size() const { return buffer_.size(); }

    /** \brief reserve memory for 'size' entries  */
    void reserve( const size_t size )
    {
        buffer_.reserve( size );
    }

    /** \brief resize buffer to 'size' entries  */
    void resize( const size_t size )
    {
        buffer_.resize( size );
    }

    /** \brief write value to buffer, value must implement the operator= correctly (i.e. no internal pointers etc.) */
    template <class T>
    void write( const T& value )
    {
      // union to access bytes in value
      const size_t tsize = sizeof( T );
      union { T value; char bytes[ tsize ]; } convert;
      convert.value = value;
      size_t pos  = buffer_.size();
      const size_t sizeNeeded = pos + tsize ;
      // reserve with some 10% overestimation
      if( buffer_.capacity() < sizeNeeded )
      {
        reserve( size_t(factor_ * sizeNeeded) ) ;
      }
      // resize to size need to store value
      buffer_.resize( sizeNeeded );
      for(unsigned int i=0; i<tsize; ++i )
      {
        buffer_[ pos++ ] = convert.bytes[ i ] ;
      }
    }

    /** \brief read value from buffer, value must implement the operator= correctly (i.e. no internal pointers etc.) */
    template <class T>
    void read( T& value ) const
    {
      // read bytes from stream and store in value
      const size_t tsize = sizeof( T );
      union { T value; char bytes[ tsize ]; } convert;
      assert( pos_ + tsize <= buffer_.size() );
      for( unsigned int i=0; i<tsize; ++i )
      {
        convert.bytes[ i ] = buffer_[ pos_++ ];
      }
      value = convert.value;
    }

    /** \brief return pointer to buffer and size for use with MPI functions */
    std::pair< char* , int > buffer() const
    {
      return std::make_pair( buffer_.data(), int(buffer_.size()) );
    }
  };

  /** \brief Point-2-Point communicator for exchange messages between processes */
  template <class MsgBuffer, class Comm = MPIHelper::MPICommunicator >
  class Point2PointCommunicator : public CollectiveCommunication< Comm >
  {
  protected:
    typedef CollectiveCommunication< Comm > BaseType;
    typedef Point2PointCommunicator< MsgBuffer, Comm > ThisType;

    // starting message tag
    static const int messagetag = 234;

    typedef std::map< int, int > linkage_t;
    typedef std::vector< int >   vector_t;

    linkage_t  _sendLinkage ;
    linkage_t  _recvLinkage ;

    vector_t   _sendDest ;
    vector_t   _recvSource ;

    mutable vector_t   _recvBufferSizes;
    mutable bool       _recvBufferSizesComputed;

  public :
    using BaseType :: rank;
    using BaseType :: size;

    // export type of message buffer
    typedef MsgBuffer MessageBufferType ;

    /* \brief data handle interface that needs to be implemented for use with some of
     * the exchange methods */
    class DataHandleInterface
    {
    protected:
      DataHandleInterface () {}
    public:
      virtual ~DataHandleInterface () {}
      virtual void   pack( const int link, MessageBufferType& os ) = 0 ;
      virtual void unpack( const int link, MessageBufferType& os ) = 0 ;
      // should contain work that could be done between send and receive
      virtual void localComputation () {}
    };

  public :
    // default constructor
    Point2PointCommunicator() : BaseType() { removeLinkage(); }
    // constructor taking CollectiveCommunication
    Point2PointCommunicator( const BaseType& comm ) : BaseType( comm ) { removeLinkage(); }

    // return new tag number for the exchange messages
    static int getMessageTag( const unsigned int icrement )
    {
      static int tag = messagetag + 2 ;
      // increase tag counter
      const int retTag = tag;
      tag += icrement ;
      // the MPI standard guaratees only up to 2^15-1
      // this needs to be revised for the all-to-all communication
      if( tag < 0 ) // >= 32767 )
      {
        // reset tag to initial value
        tag = messagetag + 2 ;
      }
      return retTag;
    }

    // return new tag number for the exchange messages
    static int getMessageTag()
    {
      return getMessageTag( 1 );
    }

    inline void computeDestinations( const linkage_t& linkage, vector_t& dest );
    inline void insertRequest( const std::set< int >& sendLinks, const std::set< int >& recvLinks );

  public:
    /** \brief return number of processes we will send data to */
    inline int sendLinks () const { return _sendLinkage.size(); }
    /** \brief return number of processes we will receive data from */
    inline int recvLinks () const { return _recvLinkage.size(); }

    /** \brief return vector containing possible recv buffer sizes */
    const vector_t& recvBufferSizes() const { return _recvBufferSizes; }

    /** \brief return send link number for a given send rank number */
    inline int sendLink (const int rank) const
    {
      assert (_sendLinkage.end () != _sendLinkage.find (rank)) ;
      return (* _sendLinkage.find (rank)).second ;
    }

    /** \brief return recv link number for a given recv rank number */
    inline int recvLink (const int rank) const
    {
      assert (_recvLinkage.end () != _recvLinkage.find (rank)) ;
      return (* _recvLinkage.find (rank)).second ;
    }

    /** \brief return vector containing all process numbers we will send to */
    const std::vector< int > &sendDest   () const { return _sendDest; }
    /** \brief return vector containing all process numbers we will receive from */
    const std::vector< int > &recvSource () const { return _recvSource; }

    /** \brief remove stored linkage */
    inline void removeLinkage () ;

    /** \brief exchange message buffers with peers defined by inserted linkage */
    virtual std::vector< MessageBufferType > exchange (const std::vector< MessageBufferType > &) const;

    /** \brief exchange data with peers, handle defines pack and unpack of data */
    virtual void exchange ( DataHandleInterface& ) const;

    /** \brief exchange data with peers, handle defines pack and unpack of data,
     *  if receive buffers are known from previous run and have not changed
     *  communication could be faster */
    virtual void exchangeCached ( DataHandleInterface& ) const;
  };

} // namespace Dune

// include inline implementation
#include "p2pcommunicator_impl.hh"

#endif // #ifndef DUNE_COMMUNICATOR_HEADER_INCLUDED
