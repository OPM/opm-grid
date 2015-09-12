#ifndef DUNE_POINT2POINTCOMMUNICATOR_IMPL_HEADER_INCLUDED
#define DUNE_POINT2POINTCOMMUNICATOR_IMPL_HEADER_INCLUDED

#include <iostream>

namespace Dune
{

  template <class MsgBuffer, class Comm >
  inline void
  Point2PointCommunicator< MsgBuffer, Comm >::
  removeLinkage()
  {
    _sendLinkage.clear();
    _recvLinkage.clear();
    _sendDest.clear();
    _recvSource.clear();

    // clear previously stored buffer sizes
    _recvBufferSizes.clear();
    _recvBufferSizesComputed = false ;
  }

  template <class MsgBuffer, class Comm >
  inline void
  Point2PointCommunicator< MsgBuffer, Comm >::
  computeDestinations( const linkage_t& linkage, vector_t& dest )
  {
    typedef linkage_t::const_iterator const_iterator ;
    dest.resize( linkage.size() );
    const const_iterator linkageEnd = linkage.end ();
    for( const_iterator i = linkage.begin (); i != linkageEnd; ++i )
    {
      dest[ (*i).second ] = (*i).first;
    }
  }

  template <class MsgBuffer, class Comm >
  inline void
  Point2PointCommunicator< MsgBuffer, Comm >::
  insertRequest( const std::set< int >& sendLinks, const std::set< int >& recvLinks )
  {
    // remove old linkage
    removeLinkage();

    const int me = rank ();

    {
      typedef std::map< int, int >::iterator iterator ;
      typedef std::set< int >::const_iterator const_iterator;

      const iterator sendEnd = _sendLinkage.end ();
      const iterator recvEnd = _recvLinkage.end ();
      const const_iterator sendLinksEnd = sendLinks.end ();
      int sendLink = 0 ;
      int recvLink = 0 ;
      for (const_iterator i = sendLinks.begin (); i != sendLinksEnd; ++i )
      {
        const int rank = (*i);
        // if rank was not inserted, insert with current link number
        if( rank != me && (_sendLinkage.find ( rank ) == sendEnd ) )
        {
          _sendLinkage.insert( std::make_pair( rank, sendLink++) );
        }
      }

      const const_iterator recvLinksEnd = recvLinks.end ();
      for (const_iterator i = recvLinks.begin (); i != recvLinksEnd; ++i )
      {
        const int rank = (*i);
        // if rank was not inserted, insert with current link number
        if( rank != me && (_recvLinkage.find ( rank ) == recvEnd ) )
        {
          _recvLinkage.insert( std::make_pair( rank, recvLink++) );
        }
      }
    }

    // compute send destinations
    computeDestinations( _sendLinkage, _sendDest );

    // compute send destinations
    computeDestinations( _recvLinkage, _recvSource );
  }


  //////////////////////////////////////////////////////////////////////////
  // non-blocking communication object
  // this class is defined here since it contains MPI information
  //////////////////////////////////////////////////////////////////////////
#if HAVE_MPI

#ifndef NDEBUG
// this is simply to avoid warning of unused variables
#define MY_INT_TEST int test =
#else
#define MY_INT_TEST
#endif

  template < class P2PCommunicator >
  class NonBlockingExchangeImplementation
  {
    typedef P2PCommunicator  P2PCommunicatorType ;
    const P2PCommunicatorType& _p2pCommunicator;

    const int _sendLinks;
    const int _recvLinks;
    const int _tag;

    MPI_Request* _sendRequest;
    MPI_Request* _recvRequest;

    const bool _recvBufferSizesKnown;
    bool _needToSend ;

    // no copying
    NonBlockingExchangeImplementation( const NonBlockingExchangeImplementation& );

    // return vector of send requests for number of send links is positive
    MPI_Request* createSendRequest() const
    {
      return ( _sendLinks > 0 ) ? new MPI_Request [ _sendLinks ] : 0;
    }

    // return vector of recv requests when
    // number of recv links is positive and symmetric is true
    MPI_Request* createRecvRequest( const bool recvBufferSizesKnown ) const
    {
      return ( _recvLinks > 0 && recvBufferSizesKnown ) ? new MPI_Request [ _recvLinks ] : 0;
    }

    // call cast operator on CollectiveCommunication to retreive MPI_Comm
    MPI_Comm mpiCommunicator() const { return static_cast< MPI_Comm > (_p2pCommunicator); }

  public:
    typedef typename P2PCommunicatorType :: DataHandleInterface DataHandleInterface;
    typedef typename P2PCommunicatorType :: MessageBufferType   MessageBufferType;

    NonBlockingExchangeImplementation( const P2PCommunicatorType& p2pComm,
                                       const int tag,
                                       const bool recvBufferSizesKnown = false )
      : _p2pCommunicator( p2pComm ),
        _sendLinks( _p2pCommunicator.sendLinks() ),
        _recvLinks( _p2pCommunicator.recvLinks() ),
        _tag( tag ),
        _sendRequest( createSendRequest() ),
        _recvRequest( createRecvRequest( recvBufferSizesKnown ) ),
        _recvBufferSizesKnown( recvBufferSizesKnown ),
        _needToSend( true )
    {
      // make sure every process has the same tag
      int mytag = tag ;
      assert ( mytag == _p2pCommunicator.max( mytag ) );
    }

    NonBlockingExchangeImplementation( const P2PCommunicatorType& p2pComm,
                                       const int tag,
                                       const std::vector< MessageBufferType > & in )
      : _p2pCommunicator( p2pComm ),
        _sendLinks( _p2pCommunicator.sendLinks() ),
        _recvLinks( _p2pCommunicator.recvLinks() ),
        _tag( tag ),
        _sendRequest( createSendRequest() ),
        _recvRequest( createRecvRequest( false ) ),
        _recvBufferSizesKnown( false ),
        _needToSend( false )
    {
      // make sure every process has the same tag
      int mytag = tag ;
      assert ( mytag == _p2pCommunicator.max( mytag ) );

      assert ( _sendLinks == int( in.size() ) );
      sendImpl( in );
    }

    /////////////////////////////////////////
    //  interface methods
    /////////////////////////////////////////
    ~NonBlockingExchangeImplementation()
    {
      if( _sendRequest )
      {
        delete [] _sendRequest;
        _sendRequest = 0;
      }

      if( _recvRequest )
      {
        delete [] _recvRequest;
        _recvRequest = 0;
      }
    }

    // virtual methods
    void send( const std::vector< MessageBufferType > & in ) { sendImpl( in ); }
    std::vector< MessageBufferType > receive() { return receiveImpl(); }

    //////////////////////////////////////////
    // implementation
    //////////////////////////////////////////

    // send data implementation
    void sendImpl( const std::vector< MessageBufferType > & osv )
    {
      // get mpi communicator
      MPI_Comm comm = mpiCommunicator();

      // get vector with destinations
      const std::vector< int >& sendDest = _p2pCommunicator.sendDest();

      // send data
      for (int link = 0; link < _sendLinks; ++link)
      {
        sendLink( sendDest[ link ], _tag, osv[ link ], _sendRequest[ link ], comm );
      }

      // set send info
      _needToSend = false ;
    }

    // receive data without buffer given
    std::vector< MessageBufferType > receiveImpl ()
    {
      // create vector of empty streams
      std::vector< MessageBufferType > out ( _recvLinks );
      receiveImpl( out );
      return out;
    }

    // receive data implementation with given buffers
    void receiveImpl ( std::vector< MessageBufferType >& out, DataHandleInterface* dataHandle = 0)
    {
      // do nothing if number of links is zero
      if( (_recvLinks + _sendLinks) == 0 ) return;

      // get mpi communicator
      MPI_Comm comm = mpiCommunicator();

      // get vector with destinations
      const std::vector< int >& recvSource = _p2pCommunicator.recvSource();

      // check whether out vector has more than one stream
      const bool useFirstStreamOnly = (out.size() == 1) ;
#ifdef NDEBUG
      {
        const int outSize = out.size();
        // check for all links messages
        for (int link = 0; link < outSize; ++link )
        {
          // contains no written data
          assert ( out[ link ].notReceived() );
        }
      }
#endif

      // flag vector holding information about received links
      std::vector< bool > linkNotReceived( _recvLinks, true );

      // count noumber of received messages
      int numReceived = 0;
      while( numReceived < _recvLinks )
      {
        // check for all links messages
        for (int link = 0; link < _recvLinks; ++link )
        {
          // if message was not received yet, check again
          if( linkNotReceived[ link ] )
          {
            // get appropriate object stream
            MessageBufferType& osRecv = useFirstStreamOnly ? out[ 0 ] : out[ link ];

            // check whether a message was completely received
            // if message was received the unpack data
            if( probeAndReceive( comm, recvSource[ link ], _tag, osRecv ) )
            {
              // if data handle was given do unpack
              if( dataHandle ) dataHandle->unpack( link, osRecv );

              // mark link as received
              linkNotReceived[ link ] = false ;

              // increase number of received messages
              ++ numReceived;
            }
          }
        }
      }

      // if send request exists, i.e. some messages have been sent
      if( _sendRequest )
      {
        // wait until all processes are done with receiving
        MY_INT_TEST MPI_Waitall ( _sendLinks, _sendRequest, MPI_STATUSES_IGNORE);
        assert (test == MPI_SUCCESS);
      }
    }

    // receive data implementation with given buffers
    void unpackRecvBufferSizeKnown( std::vector< MessageBufferType >& out, DataHandleInterface& dataHandle )
    {
      // do nothing if number of links is zero
      if( _recvLinks == 0 ) return;

      // flag vector holding information about received links
      std::vector< bool > linkNotReceived( _recvLinks, true );

      // count noumber of received messages
      int numReceived = 0;
      while( numReceived < _recvLinks )
      {
        // check for all links messages
        for (int link = 0; link < _recvLinks; ++link )
        {
          // if message was not received yet, check again
          if( linkNotReceived[ link ] )
          {
            assert( _recvRequest );
            // check whether message was received, and if unpack data
            if( receivedMessage( _recvRequest[ link ], out[ link ] ) )
            {
              // if data handle was given do unpack
              dataHandle.unpack( link, out[ link ] );

              // mark link as received
              linkNotReceived[ link ] = false ;
              // increase number of received messages
              ++ numReceived;
            }
          }
        }
      }

      // if send request exists, i.e. some messages have been sent
      if( _sendRequest )
      {
        // wait until all processes are done with receiving
        MY_INT_TEST MPI_Waitall ( _sendLinks, _sendRequest, MPI_STATUSES_IGNORE);
        assert (test == MPI_SUCCESS);
      }
    }

    // receive data implementation with given buffers
    void send( std::vector< MessageBufferType >& osSend,
               DataHandleInterface& dataHandle )
    {
      std::vector< MessageBufferType > osRecv;
      send( osSend, osRecv, dataHandle );
    }

    // receive data implementation with given buffers
    void send( std::vector< MessageBufferType >& osSend,
               std::vector< MessageBufferType >& osRecv,
               DataHandleInterface& dataHandle )
    {
      if( _needToSend )
      {
        // get mpi communicator
        MPI_Comm comm = mpiCommunicator();

        // get vector with destinations
        const std::vector< int >& sendDest = _p2pCommunicator.sendDest();

        // send data
        for (int link = 0; link < _sendLinks; ++link)
        {
          // pack data
          dataHandle.pack( link, osSend[ link ] );

          // send data
          sendLink( sendDest[ link ], _tag, osSend[ link ], _sendRequest[ link ], comm );
        }

        // set send info
        _needToSend = false ;
      }

      // resize receive buffer if in symmetric mode
      if( _recvBufferSizesKnown )
      {
        // get mpi communicator
        MPI_Comm comm = mpiCommunicator();

        osRecv.resize( _recvLinks );

        // get vector with destinations
        const std::vector< int >& recvSource = _p2pCommunicator.recvSource();
        const std::vector< int >& recvBufferSizes = _p2pCommunicator.recvBufferSizes();

        // send data
        for (int link = 0; link < _recvLinks; ++link)
        {
          // send data
          const int bufferSize = recvBufferSizes[ link ];

          // post receive if in symmetric mode
          assert( _recvRequest );
          assert( &_recvRequest[ link ] );
          postReceive( recvSource[ link ], _tag, bufferSize, osRecv[ link ], _recvRequest[ link ], comm );
        }
      }
    }

    // receive data implementation with given buffers
    void receive( DataHandleInterface& dataHandle )
    {
      // do work that can be done between send and receive
      dataHandle.localComputation() ;

      // create receive message buffers
      std::vector< MessageBufferType > out( 1 );
      // receive data
      receiveImpl( out, &dataHandle );
    }

    // receive data implementation with given buffers
    void exchange( DataHandleInterface& dataHandle )
    {
      const int recvLinks = _p2pCommunicator.recvLinks();
      // do nothing if number of links is zero
      if( (recvLinks + _sendLinks) == 0 ) return;

      // send message buffers, we need several because of the
      // non-blocking send routines, send might not be finished
      // when we start recieving
      std::vector< MessageBufferType > osSend ;
      std::vector< MessageBufferType > osRecv ;

      // if data was noy send yet, do it now
      if( _needToSend )
      {
        // resize message buffer vector
        osSend.resize( _sendLinks );

        // send data
        send( osSend, osRecv, dataHandle );
      }

      // now receive data
      if( _recvBufferSizesKnown )
        unpackRecvBufferSizeKnown( osRecv, dataHandle );
      else
        receive( dataHandle );
    }

  protected:
    int sendLink( const int dest, const int tag,
                  const MessageBufferType& msgBuffer, MPI_Request& request, MPI_Comm& comm )
    {
      // buffer = point to mem and size
      std::pair< char*, int > buffer = msgBuffer.buffer();

      MY_INT_TEST MPI_Isend ( buffer.first, buffer.second, MPI_BYTE, dest, tag, comm, &request );
      assert (test == MPI_SUCCESS);

      return buffer.second;
    }

    void postReceive( const int source, const int tag, const int bufferSize,
                      MessageBufferType& msgBuffer, MPI_Request& request, MPI_Comm& comm )
    {
      // reserve memory for receive buffer
      msgBuffer.resize( bufferSize );

      // get buffer and size
      std::pair< char*, int > buffer = msgBuffer.buffer();

      // MPI receive (non-blocking)
      {
        MY_INT_TEST MPI_Irecv ( buffer.first, buffer.second, MPI_BYTE, source, tag, comm, & request);
        assert (test == MPI_SUCCESS);
      }
    }

    // does receive operation for one link
    bool receivedMessage( MPI_Request& request, MessageBufferType& buffer )
    {
#ifndef NDEBUG
      // for checking whether the buffer size is correct
      MPI_Status status ;
#endif
      // msg received, 0 or 1
      int received = 0;

      // if receive of message is finished, unpack
      MPI_Test( &request, &received,
#ifndef NDEBUG
                &status
#else
                MPI_STATUS_IGNORE // ignore status in non-debug mode for performance reasons
#endif
              );

#ifndef NDEBUG
      if( received )
      {
        int checkBufferSize = -1;
        MPI_Get_count ( & status, MPI_BYTE, &checkBufferSize );
        if( checkBufferSize != buffer.size() )
          std::cout << "Buffer sizes don't match: "  << checkBufferSize << " " << buffer.size() << std::endl;
        assert( checkBufferSize == buffer.size() );
      }
#endif
      return bool(received);
    }

    // does receive operation for one link
    bool probeAndReceive( MPI_Comm& comm,
                          const int source,
                          const int tag,
                          MessageBufferType& osRecv )
    {
      // corresponding MPI status
      MPI_Status status;

      // msg available, 0 or 1
      // available does not mean already received
      int available = 0;

      // check for any message with tag (nonblocking)
      MPI_Iprobe( source, tag, comm, &available, &status );

      // receive message if available flag is true
      if( available )
      {
        // this should be the same, otherwise we got an error
        assert ( source == status.MPI_SOURCE );

        // length of message
        int bufferSize = -1;

        // get length of message
        {
          MY_INT_TEST MPI_Get_count ( &status, MPI_BYTE, &bufferSize );
          assert (test == MPI_SUCCESS);
        }

        // reserve memory
        osRecv.resize( bufferSize );

        // get buffer
        std::pair< char*, int > buffer = osRecv.buffer();

        // MPI receive (blocking)
        {
          MY_INT_TEST MPI_Recv ( buffer.first, buffer.second, MPI_BYTE, status.MPI_SOURCE, tag, comm, & status);
          assert (test == MPI_SUCCESS);
        }

        return true ; // received
      }
      return false ;  // not yet received
    }
  }; // end NonBlockingExchangeImplementation

#undef MY_INT_TEST
#endif // #if HAVE_MPI

  // --exchange

  template <class MsgBuffer, class Comm >
  inline void
  Point2PointCommunicator< MsgBuffer, Comm >::
  exchange( DataHandleInterface& handle ) const
  {
    assert( _recvBufferSizes.empty () );
#if HAVE_MPI
    NonBlockingExchangeImplementation< ThisType > nonBlockingExchange( *this, getMessageTag() );
    nonBlockingExchange.exchange( handle );
#endif
  }

  // --exchange
  template <class MsgBuffer, class Comm >
  inline std::vector< MsgBuffer >
  Point2PointCommunicator< MsgBuffer, Comm >::
  exchange( const std::vector< MessageBufferType > & in ) const
  {
#if HAVE_MPI
    // note: for the non-blocking exchange the message tag
    // should be different each time to avoid MPI problems
    NonBlockingExchangeImplementation< ThisType > nonBlockingExchange( *this, getMessageTag(), in );
    return nonBlockingExchange.receiveImpl();
#else
    // don nothing when MPI is not found
    return in;
#endif
  }

  // --exchange
  template <class MsgBuffer, class Comm >
  inline void
  Point2PointCommunicator< MsgBuffer, Comm >::
  exchangeCached( DataHandleInterface& handle ) const
  {
#if HAVE_MPI
    if( ! _recvBufferSizesComputed )
    {
      const int nSendLinks = sendLinks();
      std::vector< MsgBuffer > buffers( nSendLinks );
      // pack all data
      for( int link=0; link<nSendLinks; ++link )
      {
        handle.pack( link, buffers[ link ] );
      }
      // exchange data
      buffers = exchange( buffers );
      const int nRecvLinks = recvLinks();
      // unpack all data
      for( int link=0; link<nRecvLinks; ++link )
      {
        handle.unpack( link, buffers[ link ] );
      }
      // store receive buffer sizes
      _recvBufferSizes.resize( nRecvLinks );
      for( int link=0; link<nRecvLinks; ++link )
      {
        _recvBufferSizes[ link ] = buffers[ link ].size();
      }
      _recvBufferSizesComputed = true ;
    }
    else
    {
      NonBlockingExchangeImplementation< ThisType > nonBlockingExchange( *this, getMessageTag(), _recvBufferSizesComputed );
      nonBlockingExchange.exchange( handle );
    }
#endif
  }

} // namespace Dune
#endif // #ifndef DUNE_POINT2POINTCOMMUNICATOR_IMPL_HEADER_INCLUDED
