#undef NDEBUG
#include <config.h>

// Warning suppression for Dune includes.
#include <opm/common/utility/platform_dependent/disable_warnings.h>

#include <dune/grid/common/p2pcommunicator.hh>

// Re-enable warnings.
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <iostream>

void testBuffer()
{
  int iVal = 4;
  double dVal = M_LN2;

  long int lVal = std::numeric_limits< long int >::max() - 4;

  std::vector< double > values( 5 );
  for( int i=0; i<5; ++i )
    values[ i ] = 1.0/double(i+1);

  Dune::SimpleMessageBuffer buffer;

  buffer.write( iVal );
  buffer.write( dVal );
  buffer.write( lVal );
  for( int i=0; i<5; ++i )
    buffer.write( values[ i ] );

  assert( buffer.size() == (sizeof(int) + sizeof(double) + sizeof(long int)) + sizeof(double) * 5  );

  for( int i=0; i<3; ++i )
  {
    buffer.resetReadPosition();

    int iCheck = -1;
    buffer.read( iCheck );
    assert( iVal == iCheck );

    double dCheck = -1.;
    buffer.read( dCheck );
    assert( dVal == dCheck );

    long int lCheck = -1;
    buffer.read( lCheck );
    assert( lVal == lCheck );
    std::vector< double > checkValues( 5 );
    for( int j=0; j<5; ++j )
      buffer.read( checkValues[ j ] );
#ifndef NDEBUG
    for( int j=0; j<5; ++j )
      assert( std::abs(values[ j ] - checkValues[ j ] ) < 1e-12 );
#endif
  }
}

typedef Dune :: Point2PointCommunicator< Dune :: SimpleMessageBuffer > P2PCommunicatorType;

class DataHandle : public P2PCommunicatorType :: DataHandleInterface
{
  const P2PCommunicatorType& comm_;
  const bool output_ ;
public:
  typedef typename P2PCommunicatorType :: MessageBufferType MessageBufferType ;
  DataHandle( const P2PCommunicatorType& comm, const bool output )
    : comm_( comm ), output_( output ) {}

  void pack( const int /* link */, MessageBufferType& buffer )
  {
    int bsize = comm_.size() - comm_.rank();
    buffer.write( bsize );
    for( int r=comm_.rank(); r<comm_.size(); ++r )
      buffer.write( r );
  }

  void unpack( const int /* link */, MessageBufferType& buffer )
  {
    int bsize = -1;
    buffer.read( bsize );
    if( output_ )
    {
      std::cout << "Handle: Received bsize = " << bsize << std::endl;
    }
    for( int r=0; r<bsize; ++r )
    {
      int rr = -1;
      buffer.read( rr );
      if( output_ )
        std::cout << rr << std::endl;
    }
  }
};

void testCommunicator( const bool output )
{
  typedef typename P2PCommunicatorType :: MessageBufferType MessageBufferType ;

  P2PCommunicatorType comm;

  const int size = comm.size();
  const int rank = comm.rank();

  std::set<int> send;
  send.insert( rank < size-1 ? rank+1 : 0 );
  std::set<int> recv;
  recv.insert( rank > 0 ? rank-1 : size-1 );
  if( rank > 0 )
    send.insert( 0 );
  if( rank == 0 )
  {
    for( int i=1; i<size; ++i )
      recv.insert( i );
  }

  comm.insertRequest( send, recv );

  const int sendLinks = comm.sendLinks();
  std::vector< MessageBufferType > sendBuffers( sendLinks );

  for( int i=0; i<sendLinks; ++i )
  {
    int bsize = size-rank;
    sendBuffers[ i ].write( bsize );
    for( int r=rank; r<size; ++r )
      sendBuffers[ i ].write( r );
  }

  // exchange buffers
  std::vector< MessageBufferType > recvBuffers = comm.exchange( sendBuffers );

  const int recvLinks = comm.recvLinks();
  assert( int(recvBuffers.size()) == recvLinks );
  for( int i=0; i<recvLinks; ++i )
  {
    int bsize = -1;
    recvBuffers[ i ].read( bsize );
    if( output )
    {
      std::cout << "Received bsize = " << bsize << std::endl;
    }
    for( int r=0; r<bsize; ++r )
    {
      int rr = -1;
      recvBuffers[ i ].read( rr );
      if( output )
        std::cout << rr << std::endl;
    }
  }

  {
    // use handle to perform the same operations as above
    DataHandle handle( comm, output );
    comm.exchange( handle );
  }

  for( int i=0; i<5; ++i )
  {
    // use handle to perform the same operations as above
    DataHandle handle( comm, output );
    comm.exchangeCached( handle );
  }
}

int main(int argc, char** argv)
{
  // initialize MPI
  Dune::MPIHelper::instance( argc, argv );
  // test buffer
  testBuffer();
  // test communication, needs to be run with more than 1 core to be effective
  testCommunicator( false );
  return 0;
}
