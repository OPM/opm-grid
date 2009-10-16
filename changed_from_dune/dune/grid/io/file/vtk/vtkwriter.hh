// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// $Id: vtkwriter.hh 10931 2009-08-10 16:42:23Z bska $

#ifndef DUNE_VTKWRITER_HH
#define DUNE_VTKWRITER_HH

#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <vector>
#include <list>

#include <dune/common/deprecated.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/iteratorfacades.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/genericreferenceelements.hh>
#include <dune/grid/common/gridenums.hh>

#include "b64enc.hh"

/** @file
    @author Peter Bastian, Christian Engwer
    @brief Provides file i/o for the visualization toolkit
*/

/**
   @todo put vtk io intro here ...

   details and examples regarding the VTK file format can be found here:

   http://www.geophysik.uni-muenchen.de/intranet/it-service/applications/paraview/vtk-file-formats/
   (not available any more)

   http://www.geophysik.uni-muenchen.de/~moder/Paraview/VTK_File_Formats.php
   (alternative)
*/

namespace Dune
{
    /** \brief options for VTK output
        \ingroup VTK */
  struct VTKOptions
  {
    enum OutputType {
      /** @brief Output to the file is in ascii. */
      ascii,
      /** @brief Output to the file is inline binary. */
      binary, 
      /** @brief Ouput is appended binary to the file. */
      binaryappended
      // /** @brief Output to the file is compressed inline binary. */
      // binarycompressed, 
      // /** @brief Ouput is compressed and appended to the file. */
      // compressedappended
    };
    enum DataMode {
      /** @brief Output conforming data. */
      conforming,
      /** @brief Output non conforming data. */
      nonconforming
    };
  };


  // map type to name in data array
  template<class T>
  struct VTKTypeNameTraits {
    std::string operator () (){
      return "";
    }
  };

  template<>
  struct VTKTypeNameTraits<char> {
    std::string operator () () {
      return "Int8";
    }
    typedef int PrintType;
  };

  template<>
  struct VTKTypeNameTraits<unsigned char> {
    std::string operator () () {
      return "UInt8";
    }   
    typedef int PrintType;
  };

  template<>
  struct VTKTypeNameTraits<short> {
    std::string operator () () {
      return "Int16";
    }   
    typedef short PrintType;
  };

  template<>
  struct VTKTypeNameTraits<unsigned short> {
    std::string operator () () {
      return "UInt16";
    }   
    typedef unsigned short PrintType;
  };

  template<>
  struct VTKTypeNameTraits<int> {
    std::string operator () () {
      return "Int32";
    }   
    typedef int PrintType;
  };

  template<>
  struct VTKTypeNameTraits<unsigned int> {
    std::string operator () () {
      return "UInt32";
    }   
    typedef unsigned int PrintType;
  };

  template<>
  struct VTKTypeNameTraits<float> {
    std::string operator () () {
      return "Float32";
    }   
    typedef float PrintType;
  };

  template<>
  struct VTKTypeNameTraits<double> {
    std::string operator () () {
      return "Float64";
    }   
    typedef double PrintType;
  };


  /** \brief A base class for grid functions with any return type and dimension
      \ingroup VTK
      
      Trick : use double as return type
  */
  template <class Grid>
  class VTKFunction
  {
  public:
    typedef typename Grid::ctype ctype;
    enum { dim = Grid::dimension };
    typedef typename Grid::template Codim< 0 >::Entity Entity;
    
    //! return number of components (1 for scalar valued functions, 3 for
    //! vector valued function in 3D etc.)
    virtual int ncomps () const = 0;
    
    //! evaluate single component comp in the entity e at local coordinates xi
    /*! Evaluate the function in an entity at local coordinates.
      @param[in]  comp   number of component to be evaluated
      @param[in]  e      reference to grid entity of codimension 0
      @param[in]  xi     point in local coordinates of the reference element of e
      \return            value of the component
    */
    virtual double evaluate (int comp, const Entity& e, const Dune::FieldVector<ctype,dim>& xi) const = 0;
    
    //! get name
    virtual std::string name () const = 0;
    
    //! virtual destructor
    virtual ~VTKFunction () {}
  };

  /** 
   * @brief Writer for the ouput of grid functions in the vtk format.
   * @ingroup VTK
   *
   * Writes arbitrary grid functions (living on cells or vertices of a grid)
   * to a file suitable for easy visualization with 
   * <a href="http://public.kitware.com/VTK/">The Visualization Toolkit (VTK)</a>.
   */
  template< class GridView >
  class VTKWriter {
    template<int dim>
    struct P0Layout
    {
      bool contains (Dune::GeometryType gt)
        {
            return gt.dim()==dim;
        }
    };
    
    template<int dim>
    struct P1Layout
    {
      bool contains (Dune::GeometryType gt)
        {
            return gt.dim()==0;
        }
    };    

    // extract types
    typedef typename GridView::Grid Grid;
    typedef typename Grid::ctype DT;
    enum { n = GridView::dimension };
    enum { w = GridView::dimensionworld };

    typedef typename GridView::template Codim< 0 >::Entity Cell;
    typedef typename GridView::template Codim< n >::Entity Vertex;
    typedef Cell Entity;
    
    typedef typename GridView::IndexSet IndexSet;
    
    static const PartitionIteratorType VTK_Partition = InteriorBorder_Partition;
    
    typedef typename GridView::template Codim< 0 >
      ::template Partition< VTK_Partition >::Iterator
      GridCellIterator;
    typedef typename GridView::template Codim< n >
      ::template Partition< VTK_Partition >::Iterator
      GridVertexIterator;
    
    typedef MultipleCodimMultipleGeomTypeMapper< GridView, P1Layout > VertexMapper;

  public:
    typedef Dune::VTKFunction<Grid> VTKFunction;
    typedef shared_ptr< Dune::VTKFunction<Grid> > VTKFunctionPtr;

  protected:
    typedef typename std::list<VTKFunctionPtr>::const_iterator FunctionIterator;
    
    class CellIterator : public GridCellIterator
    {
    public:
      CellIterator(const GridCellIterator & x) : GridCellIterator(x) {};
      const FieldVector<DT,n> position() const
        {
            return GenericReferenceElements<DT,n>::general((*this)->type()).position(0,0);
        }
    };

    CellIterator cellBegin() const
    {
      return gridView_.template begin< 0, VTK_Partition >();
    }

    CellIterator cellEnd() const
    {
      return gridView_.template end< 0, VTK_Partition >();
    }
    
    class VertexIterator :
      public ForwardIteratorFacade<VertexIterator, const Entity>//Entity, Entity&, int>
    {
      GridCellIterator git;
      GridCellIterator gend;
      VTKOptions::DataMode datamode;
      int index;
      const VertexMapper & vertexmapper;
      std::vector<bool> visited;
      const std::vector<int> & number;
      int offset;
    protected:
      void basicIncrement ()
      {
        if( git == gend )
          return;
        ++index;
        const int numCorners = git->template count< n >();
        if( index == numCorners )
        {
          offset += numCorners;
          index = 0;

          ++git;
          while( (git != gend) && (git->partitionType() != InteriorEntity) )
            ++git;
        }
      }
    public:
      VertexIterator(const GridCellIterator & x,
                     const GridCellIterator & end,
                     const VTKOptions::DataMode & dm,
                     const VertexMapper & vm,
                     const std::vector<int> & num) :
        git(x), gend(end), datamode(dm), index(0),
        vertexmapper(vm), visited(vm.size(), false),
        number(num), offset(0)
        {
          if (datamode == VTKOptions::conforming && git != gend)
            visited[vertexmapper.template map(*git,index,n)] = true;
        };
      void increment ()
        {
          switch (datamode)
          {
          case VTKOptions::conforming:
            while(visited[vertexmapper.template map(*git,index,n)])
            {
              basicIncrement();
              if (git == gend) return;
            }
            visited[vertexmapper.template map(*git,index,n)] = true;
            break;
          case VTKOptions::nonconforming:
            basicIncrement();
            break;
          }
       }
      bool equals (const VertexIterator & cit) const
        {
          return git == cit.git
            && index == cit.index && datamode == cit.datamode;
        }
      const Entity& dereference() const
        {
          return *git;
        }
      int id () const
        {
          switch (datamode)
          {
          case VTKOptions::conforming:
            return
              number[vertexmapper.template map(*git,renumber(*git,index),n)];
          case VTKOptions::nonconforming:
            return offset + renumber(*git,index);
          default:
            DUNE_THROW(IOError,"VTKWriter: unsupported DataMode" << datamode);
          }
        }
      int localindex () const
        {
          return index;
        }
      const FieldVector<DT,n> & position () const
        {
          return GenericReferenceElements<DT,n>::general(git->type()).position(index,n);
        }
    };    
    
    VertexIterator vertexBegin () const
    {
      return VertexIterator( gridView_.template begin< 0, VTK_Partition >(),
                             gridView_.template end< 0, VTK_Partition >(),
                             datamode, *vertexmapper, number );
    }

    VertexIterator vertexEnd () const
    {
      return VertexIterator( gridView_.template end< 0, VTK_Partition >(),
                             gridView_.template end< 0, VTK_Partition >(),
                             datamode, *vertexmapper, number );
    }
    
    class CornerIterator :
      public ForwardIteratorFacade<CornerIterator, Entity, Entity&, int>
    {
      GridCellIterator git;
      GridCellIterator gend;
      VTKOptions::DataMode datamode;
      int index;
      const VertexMapper & vertexmapper;
      std::vector<bool> visited;
      const std::vector<int> & number;
      int offset;
    protected:
      void basicIncrement ()
      {
        if( git == gend )
          return;
        ++index;
        const int numCorners = git->template count< n >();
        if( index == numCorners )
        {
          offset += numCorners;
          index = 0;

          ++git;
          while( (git != gend) && (git->partitionType() != InteriorEntity) )
            ++git;
        }
      }
    public:
      CornerIterator(const GridCellIterator & x,
                     const GridCellIterator & end,
                     const VTKOptions::DataMode & dm,
                     const VertexMapper & vm,
                     const std::vector<int> & num) :
        git(x), gend(end), datamode(dm), index(0),
        vertexmapper(vm),
        number(num), offset(0) {};
      void increment ()
        {
          basicIncrement();
        }
      bool equals (const CornerIterator & cit) const
        {
          return git == cit.git
            && index == cit.index && datamode == cit.datamode;
        }
      Entity& dereference() const
        {
          return *git;
        }
      int id () const
        {
          switch (datamode)
          {
          case VTKOptions::conforming:
            return
              number[vertexmapper.template map(*git,renumber(*git,index),n)];
          case VTKOptions::nonconforming:
            return offset + renumber(*git,index);
          default:
            DUNE_THROW(IOError,"VTKWriter: unsupported DataMode" << datamode);
          }
        }
      int localindex () const
        {
          return index;
        }
    };
    
    CornerIterator cornerBegin () const
    {
      return CornerIterator( gridView_.template begin< 0, VTK_Partition >(),
                             gridView_.template end< 0, VTK_Partition >(),
                             datamode, *vertexmapper, number );
    }
    
    CornerIterator cornerEnd () const
    {
      return CornerIterator( gridView_.template end< 0, VTK_Partition >(),
                             gridView_.template end< 0, VTK_Partition >(),
                             datamode, *vertexmapper, number );
    }    
        
  private:
      /** \brief take a vector and interpret it as cell data
          \ingroup VTK
      */
    template<class V>
    class P0VectorWrapper : public VTKFunction  
    {
      typedef MultipleCodimMultipleGeomTypeMapper< GridView, P0Layout > VM0;
    public:
      //! return number of components
      virtual int ncomps () const
        {
          return 1;
        }

      //! evaluate
      virtual double evaluate (int comp, const Entity& e, const Dune::FieldVector<DT,n>& xi) const
        {
          return v[mapper.map(e)];
        }
      
      //! get name
      virtual std::string name () const
        {
          return s;
        }

      //! construct from a vector and a name
      P0VectorWrapper ( const GridView &gridView, const V &v_, std::string s_ )
      : g( gridView.grid() ),
        is( gridView.indexSet() ),
        v( v_ ),
        s( s_ ),
        mapper( gridView )
      {
        if (v.size()!=(unsigned int)mapper.size())
          DUNE_THROW(IOError,"VTKWriter::P0VectorWrapper: size mismatch");
      }

      virtual ~P0VectorWrapper() {}
      
    private:
      const Grid& g;
      const IndexSet &is;
      const V& v;
      std::string s;
      VM0 mapper;
    };

      /** \brief take a vector and interpret it as vertex data
          \ingroup VTK
          \tparam V Data container type
      */
    template<class V>
    class P1VectorWrapper : public VTKFunction  
    {
      typedef MultipleCodimMultipleGeomTypeMapper< GridView, P1Layout > VM1;

    public:
      //! return number of components
      virtual int ncomps () const
        {
          return 1;
        }

      //! evaluate
      virtual double evaluate (int comp, const Entity& e, const Dune::FieldVector<DT,n>& xi) const
        {
          double min=1E100;
          int imin=-1;
          Dune::GeometryType gt = e.type();
          for (int i=0; i<e.template count<n>(); ++i)
          {
            Dune::FieldVector<DT,n> 
              local = Dune::GenericReferenceElements<DT,n>::general(gt).position(i,n);
            local -= xi;
            if (local.infinity_norm()<min)
            {
              min = local.infinity_norm();
              imin = i;
            }
          }
          return v[mapper.template map(e,imin,n)];
        }
      
      //! get name
      virtual std::string name () const
        {
          return s;
        }

      //! construct from a vector and a name
      P1VectorWrapper ( const GridView &gridView, const V &v_, std::string s_ )
      : g( gridView.grid() ),
        is( gridView.indexSet() ),
        v( v_ ),
        s( s_ ),
        mapper( gridView )
      {
        if (v.size()!=(unsigned int)mapper.size())
          DUNE_THROW(IOError,"VTKWriter::P1VectorWrapper: size mismatch");
      }

      virtual ~P1VectorWrapper() {}
      
    private:
      const Grid& g;
      const IndexSet &is;
      const V& v;
      std::string s;
      VM1 mapper;
    };

  public:
    /**
     * @brief Construct a VTKWriter working on a specific GridView.
     * 
     * 
     * @param gridView The gridView the grid functions live on. (E. g. a LevelGridView.)
     * @param dm The data mode.
     */
    explicit VTKWriter ( const GridView &gridView,
                         VTKOptions::DataMode dm = VTKOptions::conforming )
    : gridView_( gridView ),
      datamode( dm )
    {
      indentCount = 0;
      numPerLine = 4*3; //should be a multiple of 3 !
    }

    /**
     * @brief Add a grid function that lives on the cells of the grid to the visualization.
     * @param p Dune:shared_ptr to the function to visualize
     */
    void addCellData (const VTKFunctionPtr & p)
      {
        celldata.push_back(p);
      }

    /**
     * @brief Add a grid function that lives on the cells of the grid to the visualization.
     * @param p The function to visualize.  The VTKWriter object will take
     *          ownership of the VTKFunction *p and delete it when it's done.
     */
      void addCellData (VTKFunction* p) // DUNE_DEPRECATED
      {
        celldata.push_back(VTKFunctionPtr(p));
      }

    /**
     * @brief Add a grid function (represented by container) that lives on the cells of 
     * the grid to the visualization.
     *
     * The container has to have random access via operator[] (e. g. std::vector). The 
     * value of the grid function for an arbitrary element
     * will be accessed by calling operator[] with the id of the element.
     *
     * @param v The container with the values of the grid function for each cell.
     * @param name A name to identify the grid function.
     */
    template<class V>
    void addCellData (const V& v, std::string name)
      {
        VTKFunction* p = new P0VectorWrapper< V >( gridView_, v, name );
        celldata.push_back(p);
      }

    /**
     * @brief Add a grid function that lives on the vertices of the grid to the visualization.
     * @param p The function to visualize.  The VTKWriter object will take
     *          ownership of the VTKFunction *p and delete it when it's done.
     */
      void addVertexData (VTKFunction* p) // DUNE_DEPRECATED
      {
        vertexdata.push_back(VTKFunctionPtr(p));
      }

    /**
     * @brief Add a grid function that lives on the vertices of the grid to the visualization.
     * @param p Dune:shared_ptr to the function to visualize
     */
    void addVertexData (const VTKFunctionPtr & p)
      {
        vertexdata.push_back(p);
      }

    /**
     * @brief Add a grid function (represented by container) that lives on the vertices of the 
     * grid to the visualization output.
     *
     * The container has to have random access via operator[] (e. g. std::vector). The value 
     * of the grid function for an arbitrary element
     * will be accessed by calling operator[] with the id of the element.
     *
     * @param v The container with the values of the grid function for each cell.
     * @param name A name to identify the grid function.
     */
    template<class V>
    void addVertexData (const V& v, std::string name)
      {
        VTKFunction* p = new P1VectorWrapper< V >( gridView_, v, name );
        vertexdata.push_back(p);
      }

    //! clear list of registered functions
    void clear ()
      {
        celldata.clear();
        vertexdata.clear();
      }

    //! destructor
    virtual ~VTKWriter ()
      {
        this->clear();
      }

    /** \brief write output (interface might change later)
     *
     *  \param[in]  name  basic name to write (may not contain a path)
     *  \param[in]  type  type of output (e.g,, ASCII) (optional)
     */
    std::string write ( const std::string &name,
                        VTKOptions::OutputType type = VTKOptions::ascii )
    {
      return write( name, type, gridView_.comm().rank(), gridView_.comm().size() );
    }

    /** \brief write output (interface might change later)
     *
     *  \param[in]  name        basic name to write (may not contain a path)
     *  \param[in]  path        path to data output
     *  \param[in]  extendpath  path keyword for each process
     *  \param[in]  type        type of output (e.g,, ASCII) (optional)
     */
    std::string pwrite ( const char* name,  const char* path, const char* extendpath,
                         VTKOptions::OutputType type = VTKOptions::ascii )
    {
      return pwrite( name, path, extendpath, type, gridView_.comm().rank(), gridView_.comm().size() );
    }

  protected:
    std::string write ( const std::string &name,
                        VTKOptions::OutputType type,
                        const int commRank,
                        const int commSize )
    {
      // make data mode visible to private functions
      outputtype = type;

      // reset byte counter for binary appended output
      bytecount = 0;

      // generate filename for process data
      std::ostringstream pieceName;
      if( commSize > 1 )
      {
        pieceName << "s" << std::setfill( '0' ) << std::setw( 4 ) << commSize << ":";
        pieceName << "p" << std::setfill( '0' ) << std::setw( 4 ) << commRank << ":";
      }
      pieceName << name << (GridView::dimension > 1 ? ".vtu" : ".vtp");

      // write process data
      std::ofstream file;
      file.open( pieceName.str().c_str(), std::ios::binary );
      writeDataFile( file );
      file.close();

      // for serial jobs we're done here
      if( commSize == 1 )
        return pieceName.str();

      // synchronize processes
      gridView_.comm().barrier();

      // generate name of parallel header
      std::ostringstream parallelName;
      parallelName << "s" << std::setfill( '0' ) << std::setw( 4 ) << commSize << ":";
      parallelName << name << (GridView::dimension > 1 ? ".pvtu" : ".pvtp");

      // on process 0: write out parallel header
      if( commRank == 0 )
      {
        file.open( parallelName.str().c_str() );
        writeParallelHeader( file, name.c_str(), ".", commSize );
        file.close();
      }

      // synchronize processes
      gridView_.comm().barrier();
      return parallelName.str();
    }

    //! write output; interface might change later
    std::string pwrite ( const char* name,  const char* path, const char* extendpath,
                         VTKOptions::OutputType ot,
                         const int commRank,
                         const int commSize )
      {
        // make data mode visible to private functions
        outputtype=ot;

        // reset byte counter for binary appended output
        bytecount = 0;

        // do some magic because paraview can only cope with relative pathes to piece files
        std::ofstream file;
        char piecepath[ 4096 ];
        char relpiecepath[ 4096 ];
        int n=strlen(path);
        int m=strlen(extendpath);
        if (n>0 && path[0]=='/' && path[n-1]=='/')
        {
          // 1) path is an absolute path to the directory where the pvtu file will be placed
          // 2) extendpath is an absolute path from "/" where the pieces are placed
          // 3) pieces are addressed relative in the pvtu files
          if (m==0)
          {
            // write pieces to root :-)
            piecepath[0] = '/';
            piecepath[1] = '\0';
          }
          else
          {
            // make piecepath absolute with trailing "/"
            char *p=piecepath;
            if (extendpath[0]!='/')
            {
              *p = '/';
              p++;
            }
            for (int i=0; i<m; i++)
            {
              *p = extendpath[i];
              p++;
            }
            if (*(p-1)!='/')
            {
              *p = '/';
              p++;
            }
            *p = '\0';
          }
          // path and piecepath are either "/" or have leading and trailing /
          // count slashes in path 
          int k=0;
          const char *p=path;
          while (*p!='\0')
          {
            if (*p=='/') k++;
            p++;
          }
          char *pp = relpiecepath;
          if (k>1)
          {
            for (int i=0; i<k; i++)
            {
              *pp='.'; pp++; *pp='.'; pp++; *pp='/'; pp++; 
            }
          }
          // now copy the extendpath
          for (int i=0; i<m; i++)
          {
            if (i==0 && extendpath[i]=='/') continue;
            *pp = extendpath[i];
            pp++;
          }
          if ( pp!=relpiecepath && (*(pp-1)!='/') )
          {
            *pp = '/';
            pp++;
          }
          *pp = '\0';
        }
        else
        {
          // 1) path is a relative path to the directory where pvtu files are placed
          // 2) extendpath is relative to where the pvtu files are and there the pieces are placed
          if (n==0 || m==0)
            sprintf(piecepath,"%s%s",path,extendpath);
          else
          {
            // both are non-zero
            if (path[n-1]!='/' && extendpath[0]!='/')
              sprintf(piecepath,"%s/%s",path,extendpath);
            else
              sprintf(piecepath,"%s%s",path,extendpath);
          }
          // the pieces are relative to the pvtu files
          sprintf(relpiecepath,"%s",extendpath);
        }
        char fullname[ 8192 ];
        if (GridView::dimension>1)
          sprintf(fullname,"%s/s%04d:p%04d:%s.vtu",piecepath, commSize, commRank, name);
        else
          sprintf(fullname,"%s/s%04d:p%04d:%s.vtp",piecepath, commSize, commRank, name);
        file.open(fullname,std::ios::binary);
        writeDataFile(file);
        file.close();
        gridView_.comm().barrier();
        if( commRank  ==0 )
        {
          if (GridView::dimension>1)
            sprintf(fullname,"%s/s%04d:%s.pvtu",path, commSize, name);
          else
            sprintf(fullname,"%s/s%04d:%s.pvtp",path, commSize, name);
          file.open(fullname);
          writeParallelHeader(file,name,relpiecepath, commSize );
          file.close();
        }
        gridView_.comm().barrier();
        return fullname;
      }

  protected:
    enum VTKGeometryType
    {
      vtkLine = 3,
      vtkTriangle = 5,
      vtkQuadrilateral = 9,
      vtkTetrahedron = 10,
      vtkHexahedron = 12,
      vtkPrism = 13,
      vtkPyramid = 14
    };
    
    //! mapping from GeometryType to VTKGeometryType
    static VTKGeometryType vtkType(const Dune::GeometryType & t)
      {
        if (t.isLine())
          return vtkLine;
        if (t.isTriangle())
          return vtkTriangle;
        if (t.isQuadrilateral())
          return vtkQuadrilateral;
        if (t.isTetrahedron())
          return vtkTetrahedron;
        if (t.isPyramid())
          return vtkPyramid;
        if (t.isPrism())
          return vtkPrism;
        if (t.isHexahedron())
          return vtkHexahedron;
        DUNE_THROW(IOError,"VTKWriter: unsupported GeometryType " << t);
      }

  private:
    //! write header file in parallel case to stream
    void writeParallelHeader ( std::ostream& s, const char* piecename, const char* piecepath,
                               const int commSize )
      {
        // xml header
        s << "<?xml version=\"1.0\"?>\n";

        // VTKFile
        s << "<VTKFile type=\"P" << getTypeString()
          << "\" version=\"0.1\" byte_order=\"" << getEndiannessString() << "\">\n";
        indentUp();

        // PUnstructuredGrid
        indent(s); 
        s << "<P" << getTypeString() << " GhostLevel=\"0\">\n";
        indentUp();

        // PPointData
        indent(s); s << "<PPointData";
        for (FunctionIterator it=vertexdata.begin(); it!=vertexdata.end(); ++it)
          if ((*it)->ncomps()==1)
          {
            s << " Scalars=\"" << (*it)->name() << "\"" ;
            break;
          }
        for (FunctionIterator it=vertexdata.begin(); it!=vertexdata.end(); ++it)
          if ((*it)->ncomps()>1)
          {
            s << " Vectors=\"" << (*it)->name() << "\"" ;
            break;
          }
        s << ">\n";
        indentUp();
        for (FunctionIterator it=vertexdata.begin(); it!=vertexdata.end(); ++it)
        {
          indent(s); s << "<PDataArray type=\"Float32\" Name=\"" << (*it)->name() << "\" ";
          s << "NumberOfComponents=\"" << ((*it)->ncomps()>1?3:1) << "\" ";
          s << "format=\"" << getFormatString() << "\"/>\n";
        }
        indentDown();
        indent(s); s << "</PPointData>\n";

        // PCellData
        indent(s); s << "<PCellData";
        for (FunctionIterator it=celldata.begin(); it!=celldata.end(); ++it)
          if ((*it)->ncomps()==1)
          {
            s << " Scalars=\"" << (*it)->name() << "\"" ;
            break;
          }
        for (FunctionIterator it=celldata.begin(); it!=celldata.end(); ++it)
          if ((*it)->ncomps()>1)
          {
            s << " Vectors=\"" << (*it)->name() << "\"" ;
            break;
          }
        s << ">\n";
        indentUp();
        for (FunctionIterator it=celldata.begin(); it!=celldata.end(); ++it)
        {
          indent(s); s << "<PDataArray type=\"Float32\" Name=\"" << (*it)->name() << "\" ";
          s << "NumberOfComponents=\"" << ((*it)->ncomps()>1?3:1) << "\" ";
          s << "format=\"" << getFormatString() << "\"/>\n";
        }
        indentDown();
        indent(s); s << "</PCellData>\n";

        // PPoints
        indent(s); s << "<PPoints>\n";
        indentUp();
        indent(s); s << "<PDataArray type=\"Float32\" Name=\"Coordinates\" NumberOfComponents=\"" << "3" << "\" ";
        s << "format=\"" << getFormatString() << "\"/>\n";
        indentDown();
        indent(s); s << "</PPoints>\n";

        // Pieces
        for( int i = 0; i < commSize; ++i )
        {
          char fullname[ 4096 ];
          if (GridView::dimension>1)
            sprintf(fullname,"%s/s%04d:p%04d:%s.vtu",piecepath, commSize, i,piecename);
          else
            sprintf(fullname,"%s/s%04d:p%04d:%s.vtp",piecepath, commSize, i,piecename);
          indent(s); s << "<Piece Source=\"" << fullname << "\"/>\n";
        }

        // /PUnstructuredGrid
        indentDown();
        indent(s); 
        s << "</P" << getTypeString() << ">\n";

        // /VTKFile
        indentDown();
        s << "</VTKFile>\n";

        s.flush();
      }

    //! write data file to stream
    void writeDataFile (std::ostream& s)
      {
        // xml header
        s << "<?xml version=\"1.0\"?>\n";

        // VTKFile
        s << "<VTKFile type=\"" << getTypeString()
          << "\" version=\"0.1\" byte_order=\"" << getEndiannessString() << "\">\n";
        indentUp();

        // Grid characteristics
        vertexmapper = new VertexMapper( gridView_ );
        if (datamode == VTKOptions::conforming)
        {
          number.resize(vertexmapper->size());
          for (std::vector<int>::size_type i=0; i<number.size(); i++) number[i] = -1;
        }
        countEntities(nvertices, ncells, ncorners);

        // UnstructuredGrid
        indent(s);
        s << "<" << getTypeString() << ">\n";
        indentUp();

        // Piece
        indent(s);
        if (n>1)
          s << "<Piece NumberOfPoints=\"" << nvertices << "\" NumberOfCells=\"" << ncells << "\">\n";
        else
          s << "<Piece NumberOfPoints=\"" << nvertices << "\""
            << " NumberOfVerts=\"0\""
            << " NumberOfLines=\"" << ncells << "\">" 
            << " NumberOfPolys=\"0\"\n";
        indentUp();

        // PointData
        writeVertexData(s);

        // CellData
        writeCellData(s);

        // Points
        writeGridPoints(s);

        // Cells
        writeGridCells(s);

        // /Piece
        indentDown();
        indent(s); s << "</Piece>\n";

        // /UnstructuredGrid
        indentDown();
        indent(s); 
        s << "</" << getTypeString() << ">\n";

        // write appended binary dat section
        if (outputtype==VTKOptions::binaryappended)
          writeAppendedData(s);

        // /VTKFile
        indentDown();
        s << "</VTKFile>\n";

        s.flush();
        
        delete vertexmapper; number.clear();
      }

  protected:
    std::string getEndiannessString() const
      {
          short i = 1;
          if (reinterpret_cast<char*>(&i)[1] == 1)
              return "BigEndian";
          else
              return "LittleEndian";
      }
      
    std::string getFormatString() const
      {
          if (outputtype==VTKOptions::ascii)
              return "ascii";
          if (outputtype==VTKOptions::binary)
              return "binary";
          if (outputtype==VTKOptions::binaryappended)
              return "appended";
          DUNE_THROW(IOError, "VTKWriter: unsupported OutputType" << outputtype);
      }
      
    std::string getTypeString() const
      {
          if (n==1)
              return "PolyData";
          else
              return "UnstructuredGrid";
      }
      
    //! count the vertices, cells and corners
    virtual void countEntities(int &nvertices, int &ncells, int &ncorners)
      {
        nvertices = 0;
        ncells = 0;
        ncorners = 0;
        for (CellIterator it=cellBegin(); it!=cellEnd(); ++it)
        {
          ncells++;
          for (int i=0; i<it->template count<n>(); ++i)
          {
            ncorners++;
            if (datamode == VTKOptions::conforming)
            {
              int alpha = vertexmapper->template map(*it,i,n);
              if (number[alpha]<0)
                number[alpha] = nvertices++;
            }
            else
            {
              nvertices++;
            }
          }
        }
      }

    //! write cell data
    virtual void writeCellData (std::ostream& s)
      {
        indent(s); s << "<CellData";
        for (FunctionIterator it=celldata.begin(); it!=celldata.end(); ++it)
          if ((*it)->ncomps()==1)
          {
            s << " Scalars=\"" << (*it)->name() << "\"" ;
            break;
          }
        for (FunctionIterator it=celldata.begin(); it!=celldata.end(); ++it)
          if ((*it)->ncomps()>1)
          {
            s << " Vectors=\"" << (*it)->name() << "\"" ;
            break;
          }
        s << ">\n";
        indentUp();
        for (FunctionIterator it=celldata.begin(); it!=celldata.end(); ++it)
        {
          VTKDataArrayWriter<float> *p = makeVTKDataArrayWriter<float>(s, (*it)->name().c_str(), (*it)->ncomps(), (*it)->ncomps()*ncells);
          for (CellIterator i=cellBegin(); i!=cellEnd(); ++i)
          {
            for (int j=0; j<(*it)->ncomps(); j++)
              p->write((*it)->evaluate(j,*i,i.position()));
			//vtk file format: a vector data always should have 3 comps(with 3rd comp = 0 in 2D case)
			if((*it)->ncomps()==2)
			  p->write(0.0);
          }
          delete p;
        }
        indentDown();
        indent(s); s << "</CellData>\n";
        s.flush();
      }

    //! write vertex data
    virtual void writeVertexData (std::ostream& s)
      {
        indent(s); s << "<PointData";
        for (FunctionIterator it=vertexdata.begin(); it!=vertexdata.end(); ++it)
          if ((*it)->ncomps()==1)
          {
            s << " Scalars=\"" << (*it)->name() << "\"" ;
            break;
          }
        for (FunctionIterator it=vertexdata.begin(); it!=vertexdata.end(); ++it)
          if ((*it)->ncomps()>1)
          {
            s << " Vectors=\"" << (*it)->name() << "\"" ;
            break;
          }
        s << ">\n";
        indentUp();
        for (FunctionIterator it=vertexdata.begin(); it!=vertexdata.end(); ++it)
        {
          VTKDataArrayWriter<float> *p = makeVTKDataArrayWriter<float>(s, (*it)->name().c_str(), (*it)->ncomps(), (*it)->ncomps()*nvertices);
          for (VertexIterator vit=vertexBegin(); vit!=vertexEnd(); ++vit)
          {
            for (int j=0; j<(*it)->ncomps(); j++)
              p->write((*it)->evaluate(j,*vit,vit.position()));
			//vtk file format: a vector data always should have 3 comps(with 3rd comp = 0 in 2D case)
			if((*it)->ncomps()==2)
			  p->write(0.0);
          }
          delete p;
        }
        indentDown();
        indent(s); s << "</PointData>\n";
        s.flush();
      }

    //! write the positions of vertices
    virtual void writeGridPoints (std::ostream& s)
      {
        indent(s); s << "<Points>\n";
        indentUp();

        VTKDataArrayWriter<float> *p = makeVTKDataArrayWriter<float>(s, "Coordinates", 3, 3*nvertices);
        VertexIterator vEnd = vertexEnd();
        for (VertexIterator vit=vertexBegin(); vit!=vEnd; ++vit)
        {
          int dimw=w;
          for (int j=0; j<std::min(dimw,3); j++)
              p->write(vit->geometry().corner(vit.localindex())[j]);
          for (int j=std::min(dimw,3); j<3; j++)
            p->write(0.0);
        }
        delete p;
      
        indentDown();
        indent(s); s << "</Points>\n";
        s.flush();
      }

    //! write the connectivity array
    virtual void writeGridCells (std::ostream& s)
      {
        indent(s); 
        if (n>1)
          s << "<Cells>\n";
        else
          s << "<Lines>\n";
        indentUp();

        // connectivity
        VTKDataArrayWriter<int> *p1 = makeVTKDataArrayWriter<int>(s, "connectivity", 1, ncorners);
        for (CornerIterator it=cornerBegin(); it!=cornerEnd(); ++it)
          p1->write(it.id());
        delete p1;

        // offsets
        VTKDataArrayWriter<int> *p2 = makeVTKDataArrayWriter<int>(s, "offsets", 1, ncells);
        {
          int offset = 0;
          for (CellIterator it=cellBegin(); it!=cellEnd(); ++it)
          {
            offset += it->template count<n>();
            p2->write(offset);
          }
        }
        delete p2;

        // types
        if (n>1)
        {
          VTKDataArrayWriter<unsigned char> *p3 = makeVTKDataArrayWriter<unsigned char>(s, "types", 1, ncells);
          for (CellIterator it=cellBegin(); it!=cellEnd(); ++it)
          {
            int vtktype = vtkType(it->type());
            p3->write(vtktype);
          }
          delete p3;
        }

        indentDown();
        indent(s); 
        if (n>1)
          s << "</Cells>\n";
        else
          s << "</Lines>\n";
        s.flush();
      }

    //! write the appended data sections
    virtual void writeAppendedData (std::ostream& s)
      {
        indent(s); s << "<AppendedData encoding=\"raw\">\n";
        indentUp();
        indent(s); s << "_"; // indicates start of binary data

        SimpleStream stream(s);

        // write length before each data block
        unsigned int blocklength;

        // point data     
        for (FunctionIterator it=vertexdata.begin(); it!=vertexdata.end(); ++it)
        {
		  
		  blocklength = nvertices * (*it)->ncomps() * sizeof(float);
		  //vtk file format: a vector data always should have 3 comps(with 3rd comp = 0 in 2D case)
		  if((*it)->ncomps()==2)
			blocklength = nvertices * (3) * sizeof(float);
          stream.write(blocklength);
          std::vector<bool> visited(vertexmapper->size(), false);
          for (VertexIterator vit=vertexBegin(); vit!=vertexEnd(); ++vit)
          {
            for (int j=0; j<(*it)->ncomps(); j++)
            {
              float data = (*it)->evaluate(j,*vit,vit.position());
              stream.write(data);
            }
			//vtk file format: a vector data always should have 3 comps(with 3rd comp = 0 in 2D case)
			if((*it)->ncomps()==2)
            {
			  float data=0.0;
			  stream.write(data);
            }
          }
        }

        // cell data
        for (FunctionIterator it=celldata.begin(); it!=celldata.end(); ++it)
        {
          blocklength = ncells * (*it)->ncomps() * sizeof(float);
          stream.write(blocklength);
          for (CellIterator i=cellBegin(); i!=cellEnd(); ++i)
          {
            for (int j=0; j<(*it)->ncomps(); j++)
            {
              float data = (*it)->evaluate(j,*i,i.position());
              stream.write(data);
            }
			//vtk file format: a vector data always should have 3 comps(with 3rd comp = 0 in 2D case)
			if((*it)->ncomps()==2)
            {
			  float data=0.0;
			  stream.write(data);
            }
          }
        }
      
        // point coordinates
        blocklength = nvertices * 3 * sizeof(float);
        stream.write(blocklength);
        std::vector<bool> visited(vertexmapper->size(), false);
        for (VertexIterator vit=vertexBegin(); vit!=vertexEnd(); ++vit)
        {
          int dimw=w;
          float data;
          for (int j=0; j<std::min(dimw,3); j++)
          {
              data = vit->geometry().corner(vit.localindex())[j];
            stream.write(data);
          }
          data = 0;
          for (int j=std::min(dimw,3); j<3; j++)
            stream.write(data);
        }
      
        // connectivity
        blocklength = ncorners * sizeof(unsigned int);
        stream.write(blocklength);
        for (CornerIterator it=cornerBegin(); it!=cornerEnd(); ++it)
        {
          stream.write(it.id());
        }

        // offsets
        blocklength = ncells * sizeof(unsigned int);
        stream.write(blocklength);
        {
          int offset = 0;
          for (CellIterator it=cellBegin(); it!=cellEnd(); ++it)
          {
            offset += it->template count<n>();
            stream.write(offset);
          }
        }

        // cell types
        if (n>1)
        {
          blocklength = ncells * sizeof(unsigned char);
          stream.write(blocklength);
          for (CellIterator it=cellBegin(); it!=cellEnd(); ++it)
          {
            unsigned char vtktype = vtkType(it->type());
            stream.write(vtktype);
          }
        }

        s << std::endl;
        indentDown();
        indent(s); s << "</AppendedData>\n";
        s.flush();
      }
    
    //! base class for data array writers
    template<class T>
    class VTKDataArrayWriter
    {
    public:
      //! write one data element
      virtual void write (T data) = 0;
      //! virtual destructor
      virtual ~VTKDataArrayWriter () {}
    };

  private:
    //! a streaming writer for data array tags, uses ASCII inline format
    template<class T>
    class VTKAsciiDataArrayWriter : public VTKDataArrayWriter<T>
    {
    public:
      //! make a new data array writer
      VTKAsciiDataArrayWriter (std::ostream& theStream, std::string name, int ncomps) 
        : s(theStream), counter(0), numPerLine(12)
        {
          VTKTypeNameTraits<T> tn;
          s << "<DataArray type=\"" << tn() << "\" Name=\"" << name << "\" ";
		  //vtk file format: a vector data always should have 3 comps(with 3rd comp = 0 in 2D case)
          if (ncomps>3)
            DUNE_THROW(IOError, "VTKWriter does not support more than 3 components");
          s << "NumberOfComponents=\"" << (ncomps>1?3:1) << "\" ";
          s << "format=\"ascii\">\n";
        }

      //! write one data element to output stream
      void write (T data)
        {
          typedef typename VTKTypeNameTraits<T>::PrintType PT;
          s << (PT) data << " ";
          counter++;
          if (counter%numPerLine==0) s << "\n";
        }

      //! finish output; writes end tag
      ~VTKAsciiDataArrayWriter ()
        {
          if (counter%numPerLine!=0) s << std::endl;
          s << "</DataArray>\n";
        }

    private:
      std::ostream& s;
      int counter;
      int numPerLine;
    };

    // a streaming writer for data array tags, uses binary inline format
    template<class T>
    class VTKBinaryDataArrayWriter : public VTKDataArrayWriter<T>
    {
    public:
      //! make a new data array writer
      VTKBinaryDataArrayWriter (std::ostream& theStream, std::string name, int ncomps, int nitems) 
        : s(theStream)
        {
          VTKTypeNameTraits<T> tn;
          ncomps = (ncomps>1?3:1);
          s << "<DataArray type=\"" << tn() << "\" Name=\"" << name << "\" ";
		  //vtk file format: a vector data always should have 3 comps(with 3rd comp = 0 in 2D case)
          if (ncomps>3)
            DUNE_THROW(IOError, "VTKWriter does not support more than 3 components");
          s << "NumberOfComponents=\"" << ncomps << "\" ";
          s << "format=\"binary\">\n";
          // reset chunk
          chunk.txt.read(0,0);
          // store size
          unsigned long int size = ncomps*nitems*sizeof(T);
          b64enc(size);
          flush();
        }

      //! write one data element to output stream
      void write (T data)
        {
          b64enc(data);
        }

      //! finish output; writes end tag
      ~VTKBinaryDataArrayWriter ()
        {
            flush();
            s << "\n</DataArray>\n";
            s.flush();
        }

    private:
      template <class X>
      void b64enc(X & data)
        {
            char* p = reinterpret_cast<char*>(&data);
            for (size_t len = sizeof(X); len > 0; len--,p++)
            {
                chunk.txt.put(*p);
                if (chunk.txt.size == 3)
                {
                    chunk.data.write(obuf);
                    s.write(obuf,4);
                }
            }
        }

      void flush()
        {
            if (chunk.txt.size > 0)
            {
                chunk.data.write(obuf);
                s.write(obuf,4);
            }
        }
        
      std::ostream& s;
      b64chunk chunk;
      char obuf[4];
    };

    //! a streaming writer for data array tags, uses binary appended format
    template<class T>
    class VTKBinaryAppendedDataArrayWriter : public VTKDataArrayWriter<T>
    {
    public:
      //! make a new data array writer
      VTKBinaryAppendedDataArrayWriter (std::ostream& theStream, std::string name, int ncomps, unsigned int& bc) 
        : s(theStream),bytecount(bc)
        {
          VTKTypeNameTraits<T> tn;
          s << "<DataArray type=\"" << tn() << "\" Name=\"" << name << "\" ";
		  //vtk file format: a vector data always should have 3 comps(with 3rd comp = 0 in 2D case)
          if (ncomps>3)
            DUNE_THROW(IOError, "VTKWriter does not support more than 3 components");
          s << "NumberOfComponents=\"" << (ncomps>1?3:1) << "\" ";
          s << "format=\"appended\" offset=\""<< bytecount << "\" />\n";
          bytecount += 4; // header
        }

      //! write one data element to output stream
      void write (T data)
        {
          bytecount += sizeof(T);
        }

    private:
      std::ostream& s;
      unsigned int& bytecount;
    };

  protected:
    /** @brief Make a VTKDataArrayWriter with new
     *
     * @param s           The stream to write to
     * @param name        The name of the vtk array
     * @param components  The number of components of the vector
     * @param totallength the total number of entries, i.e. components*vectors
     */
    template<class T>
    VTKDataArrayWriter<T> *makeVTKDataArrayWriter(std::ostream &s,
                                                  const char *name,
                                                  unsigned int components,
                                                  unsigned int totallength)
    {
      switch(outputtype) {
      case VTKOptions::ascii:
        return new VTKAsciiDataArrayWriter<T>(s, name, components);
      case VTKOptions::binary:
        return new VTKBinaryDataArrayWriter<T>(s, name, components, totallength);
      case VTKOptions::binaryappended:
        return new VTKBinaryAppendedDataArrayWriter<T>(s, name, components, bytecount);
      }
      DUNE_THROW(IOError, "VTKWriter: unsupported OutputType" << outputtype);
    }

    //! write out data in binary
    class SimpleStream
    {
    public:
      //! make a new stream
      SimpleStream (std::ostream& theStream)
        : s(theStream)
        {}
      //! write data to stream
      template<class T>
      void write (T data)
        {
          char* p = reinterpret_cast<char*>(&data);
          s.write(p,sizeof(T));
        }
    private:
      std::ostream& s;
    };

    //! increase indentation level
    void indentUp ()
      {
        indentCount++;
      }

    //! decrease indentation level
    void indentDown ()
      {
        indentCount--;
      }

    //! write indentation to stream
    void indent (std::ostream& s)
      {
        for (int i=0; i<indentCount; i++) 
          s << "  ";
      }

    //! renumber VTK -> Dune
    static int renumber (const GeometryType &t, int i)
      {
        static const int quadRenumbering[4] = {0,1,3,2};
        static const int cubeRenumbering[8] = {0,1,3,2,4,5,7,6};
        static const int prismRenumbering[6] = {0,2,1,3,5,4};
        static const int pyramidRenumbering[6] = {0,2,1,3,5,4};
        if (t.isQuadrilateral())
          return quadRenumbering[i];
        if (t.isPyramid())
          return pyramidRenumbering[i];
        if (t.isPrism())
          return prismRenumbering[i];
        if (t.isHexahedron())
          return cubeRenumbering[i];
        return i;
      }
    static int renumber (const Entity& e, int i)
      {
        return renumber(e.type(), i);
      }

    // the list of registered functions
    std::list<VTKFunctionPtr> celldata;
    std::list<VTKFunctionPtr> vertexdata;

  private: 
    // the grid
    GridView gridView_;

    // indend counter
    int indentCount;
    int numPerLine;

    // temporary grid information
  protected:
    int ncells;
    int nvertices;
    int ncorners;
  private:
    VertexMapper* vertexmapper;
    std::vector<int> number;
    VTKOptions::DataMode datamode;
  protected:
    VTKOptions::OutputType outputtype;
  private:
    unsigned int bytecount;
  };

  /** \brief VTKWriter on the leaf grid
      \ingroup VTK
   */
  template< class Grid >
  class LeafVTKWriter
  : public VTKWriter< typename Grid::LeafGridView >
  {
    typedef VTKWriter< typename Grid::LeafGridView > Base;

  public:
      /** \brief Construct a VTK writer for the leaf level of a given grid */
    explicit LeafVTKWriter ( const Grid &grid,
                             VTKOptions::DataMode dm = VTKOptions::conforming ) DUNE_DEPRECATED
    : Base( grid.leafView(), dm )
    {}
  };

  /** \brief VTKWriter on a given level grid
      \ingroup VTK
   */
  template< class Grid >
  class LevelVTKWriter
  : public VTKWriter< typename Grid::LevelGridView >
  {
    typedef VTKWriter< typename Grid::LevelGridView > Base;

  public:
      /** \brief Construct a VTK writer for a certain level of a given grid */
    LevelVTKWriter ( const Grid &grid, int level,
                     VTKOptions::DataMode dm = VTKOptions::conforming ) DUNE_DEPRECATED
    : Base( grid.levelView( level ), dm )
    {}
  };
  
}

#endif
