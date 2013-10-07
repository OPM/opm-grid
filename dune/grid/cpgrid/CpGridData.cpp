#include"config.h"
#include"CpGridData.hpp"
#include"Intersection.hpp"
#include"Entity.hpp"
#include"OrientedEntityTable.hpp"
#include"Indexsets.hpp"
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/grid/common/GridPartitioning.hpp>
#include <algorithm>

namespace Dune
{
namespace cpgrid
{
CpGridData::CpGridData()
    : index_set_(), local_id_set_(),
    ccobj_(Dune::MPIHelper::getCommunicator())
{}
CpGridData::CpGridData(CpGrid& grid)
  : index_set_(new IndexSet(*this)),   local_id_set_(new IdSet(*this)),
    ccobj_(Dune::MPIHelper::getCommunicator()), use_unique_boundary_ids_(false)
{}
CpGridData::~CpGridData()
{
    if(index_set_) delete index_set_;
    if(local_id_set_) delete local_id_set_;
}

void CpGridData::computeUniqueBoundaryIds()
{
    // Perhaps we should make available a more comprehensive interface
    // for EntityVariable, so that we don't have to build a separate
    // vector and assign() to unique_boundary_ids_ at the end.
    int num_faces = face_to_cell_.size();
    std::vector<int> ids(num_faces, 0);
    int count = 0;
    for (int i = 0; i < num_faces; ++i) {
        cpgrid::EntityRep<1> face(i, true);
        if (face_to_cell_[face].size() == 1) {
            // It's on the boundary.
            // Important! Since boundary ids run from 1 to n,
            // we use preincrement instead of postincrement below.
            ids[i] = ++count;
        }
    }
    unique_boundary_ids_.assign(ids.begin(), ids.end());
#ifdef VERBOSE
    std::cout << "computeUniqueBoundaryIds() gave all boundary intersections\n"
              << "unique boundaryId()s ranging from 1 to " << count << std::endl;
#endif
}

int CpGridData::size(int codim) const
{
    switch (codim) {
    case 0: return cell_to_face_.size();
    case 1: return 0;
    case 2: return 0;
    case 3: return geomVector<3>().size();
    default: return 0;
    }
}


 // A functor that counts existent entries and renumbers them.
struct CountExistent
{
    void operator()(int& i)
    {
        if(i < std::numeric_limits<int>::max())
            count++;
    }
    int count;
};


std::unique_ptr<CpGridData>  CpGridData::distributeGlobalGrid(CpGrid& grid,
                                                              CpGridData& view_data,
                          std::vector<int>& cell_part,
                          Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator>& ccobj) const
{
    int my_rank=ccobj.rank();
    //int size=ccobj.size();
#ifdef DEBUG
    typedef std::vector<int>::const_iterator Iterator;
    for(Iterator i=cell_part.begin(); i!= cell_part.end(); ++i)
        if(*i>=size)
            OPM_THROW(std::runtime_error, "rank for cell is too big");
#endif
    // vector with the set of ranks that
    std::vector<std::set<int> > overlap; 
    
    //Dune::OwnerOverlapCopyCommunication<int,Dune::ParallelLocalIndex> ownerOverlap;
    typedef OwnerOverlapCopyAttributeSet::AttributeSet AttributeSet;
    typedef Dune::ParallelIndexSet<int,Dune::ParallelLocalIndex<AttributeSet> > IndexSet;
    overlap.resize(cell_part.size());
    addOverlapLayer(grid, cell_part, overlap, my_rank, false);
    // count number of cells
    struct CellCounter
    {
        int myrank;
        int count;
        typedef std::set<int>::const_iterator NeighborIterator;
        std::set<int> neighbors;
        typedef Dune::ParallelLocalIndex<AttributeSet> Index;
        IndexSet* indexset;
        std::vector<int> global2local;
        /**
         * @brief Adds an index with flag owner to the index set.
         *
         * An index is added if the rank is our rank.
         * @param i The global index
         * @param rank The rank that owns the index.
         */
        void operator() (int i, int rank)
        {
            if(rank==myrank)
            {
                global2local.push_back(count);
                indexset->add(i, Index(++count, OwnerOverlapCopyAttributeSet::owner, true));
            }
            else
                global2local.push_back(std::numeric_limits<int>::max());
        }
        /**
         * @brief Adds an index with flag overlap to the index set.
         * @param i The global index.
         * @param ov The set of ranks where the index is in the overlap
         * region.
         * @param rank The rank that owns this index.
         */
        void operator() (int i, const std::set<int>& ov, int rank)
        {
            std::set<int>::const_iterator iter=
                ov.find(myrank);
            if(iter!=ov.end())
            {
                global2local.push_back(count);
                indexset->add(i, Index(++count, OwnerOverlapCopyAttributeSet::overlap, true));
                neighbors.insert(rank);
            }
            else
                global2local.push_back(std::numeric_limits<int>::max());
            
        }
    } cell_counter;
    cell_counter.myrank=my_rank;
    cell_counter.count=0;
    cell_counter.global2local.reserve(cell_part.size());
    IndexSet indexset;
    cell_counter.indexset=&indexset;
    // set up the index set.
    cell_counter.indexset->beginResize();
    typedef std::vector<std::set<int> >::const_iterator OIterator;
    std::vector<int>::const_iterator ci=cell_part.begin();
    for(OIterator end=overlap.end(), begin=overlap.begin(), i=begin; i!=end; ++i, ++ci)
    {
        if(i->size())
            cell_counter(i-begin, *i, *ci);
        else
            cell_counter(i-begin, *ci);
    }
    cell_counter.indexset->endResize();
    // setup the remote indices.
    //typedef Dune::OwnerOverlapCopyCommunication<int,Dune::ParallelLocalIndex>::RemoteIndices RemoteIndices;
    typedef Dune::RemoteIndices<IndexSet> RemoteIndices;
    typedef RemoteIndexListModifier<RemoteIndices::ParallelIndexSet, RemoteIndices::Allocator,
                                    false> Modifier;
    typedef RemoteIndices::RemoteIndex RemoteIndex;
    RemoteIndices remoteIndices;

    // Create a map of ListModifiers
    { //extra scope to call destructor of the Modifiers
        std::map<int,Modifier> modifiers;
        for(CellCounter::NeighborIterator n=cell_counter.neighbors.begin(), end=cell_counter.neighbors.end();
            n != end; ++n)
            modifiers.insert(std::make_pair(*n, remoteIndices.getModifier<false,false>(*n)));
        // Insert remote indices. For each entry in the index set, see wether there are overlap occurences and add them.
        for(IndexSet::const_iterator i=indexset.begin(), end=indexset.end();
            i!=end; ++i)
        {
            std::set<int>::const_iterator iter=overlap[i->local()].find(my_rank);
            if(iter!=overlap[i->local()].end())
                modifiers[cell_part[i->local()]]
                    .insert(RemoteIndex(OwnerOverlapCopyAttributeSet::owner,&(*i)));
        }
    }

    // We can identify existing cells with the help of the index set.
    // Now we need to compute the existing faces and points. Either exist
    // if they are reachable from an existing cell.
    // We use std::numeric_limits<int>::max() to indicate non-existent entities.
    std::vector<int> face_indicator(view_data.geometry_.geomVector<1>().size(),
                                    std::numeric_limits<int>::max());
    std::vector<int> point_indicator(view_data.geometry_.geomVector<3>().size(), 
                                     std::numeric_limits<int>::max());
    for(IndexSet::iterator i=indexset.begin(), end=indexset.end();
            i!=end; ++i)
    {
        typedef boost::iterator_range<EntityRep<1>*>::iterator RowIter;
        int row_index=i->global();
        // Somehow g++-4.4 does not find functions of father even if we
        // change inheritance of OrientedEntityTable to public.
        // Therfore we use an ugly task to base class here.
        Opm::SparseTable<EntityRep<1> >& c2f=view_data.cell_to_face_;
        for(RowIter f=c2f[row_index].begin(), fend=c2f[i->global()].end();
            f!=fend; ++f)
        {
            int findex=f->index();
            --face_indicator[findex];
            // points reachable from a cell exist, too.
            for(auto p=view_data.face_to_point_[findex].begin(), 
                    pend=view_data.face_to_point_[findex].end(); p!=pend; ++p)
            {
                assert(p>=0);
                --point_indicator[*p];
            }
        }
    }
    // Set up the new topology arrays
    EntityVariable<cpgrid::Geometry<3, 3>, 0> cell_geom;
    std::vector<cpgrid::Geometry<3, 3> > tmp_cell_geom(indexset.size());
    auto global_cell_geom=view_data.geomVector<0>();
    std::vector<int> global_cell;
    global_cell.resize(indexset.size());

    // Copy the existing cells.
    for(auto i=indexset.begin(), end=indexset.end(); i!=end; ++i)
    {
        tmp_cell_geom[i->local()]=static_cast<std::vector<cpgrid::Geometry<3, 3> >&>(global_cell_geom)[i->global()];
        global_cell[i->local()]=view_data.global_cell_[i->global()];
    }
    static_cast<std::vector<cpgrid::Geometry<3, 3> >&>(cell_geom).swap(tmp_cell_geom);

    CountExistent exCount;
    // count the existing faces, renumber, and allocate space.
    exCount.count = 0;
    std::for_each(face_indicator.begin(), face_indicator.end(), CountExistent());
    EntityVariable<cpgrid::Geometry<2, 3>, 1> face_geom;
    std::vector<cpgrid::Geometry<2, 3> > tmp_face_geom(exCount.count);
    EntityVariable<enum face_tag, 1> face_tag;
    std::vector<enum face_tag> tmp_face_tag(exCount.count);
    SignedEntityVariable<PointType, 1> face_normals;
    std::vector<PointType> tmp_face_normals(exCount.count);
    const std::vector<cpgrid::Geometry<2, 3> >& global_face_geom=view_data.geomVector<1>();
    const std::vector<enum face_tag>& global_face_tag=view_data.face_tag_;
    const std::vector<PointType>& global_face_normals=view_data.face_normals_;
    
    // Now copy the face geometries that do exist.
    auto f = tmp_face_geom.begin();
    auto ft = tmp_face_tag.begin();
    auto fn = tmp_face_normals.begin();
    for(auto begin=face_indicator.begin(), fi=begin,  end=face_indicator.end(); fi!=end; 
        ++fi)
    {
        if(*fi<std::numeric_limits<int>::max())
        {
            *f=global_face_geom[fi-begin];
            *ft=global_face_tag[fi-begin];
            *fn=global_face_normals[fi-begin];
            ++f; ++ft; ++fn;
        }
    }
    static_cast<std::vector<PointType>&>(face_normals).swap(tmp_face_normals);
    static_cast<std::vector<cpgrid::Geometry<2, 3> >&>(face_geom).swap(tmp_face_geom);
    static_cast<std::vector<enum face_tag>&>(face_tag).swap(tmp_face_tag);

    // Count the existing points and allocate space
    exCount.count=0;
    exCount = std::for_each(point_indicator.begin(), point_indicator.end(),exCount);
    std::vector<cpgrid::Geometry<0, 3> > tmp_point_geom(exCount.count);
    EntityVariable<cpgrid::Geometry<0, 3>, 3> point_geom;
    const std::vector<cpgrid::Geometry<0, 3> >& global_point_geom=view_data.geomVector<3>();
    
    // Now copy the point geometries that do exist.
    auto p= tmp_point_geom.begin();
    for(auto begin=point_indicator.begin(), pi=begin,  end=point_indicator.end(); pi!=end;
        ++pi)
    {
        if(*pi<std::numeric_limits<int>::max())
        {
            *p=global_point_geom[pi-begin];
            ++p;
        }
    }
    // swap the underlying vectors to get data into point_geom
    static_cast<std::vector<cpgrid::Geometry<0, 3> >&>(point_geom).swap(tmp_point_geom);
    
    // Copy the vectors to geometry. There is no other way currently.
    DefaultGeometryPolicy geometry=cpgrid::DefaultGeometryPolicy(cell_geom, face_geom, 
                                            point_geom);

    // Create the topology information. This is stored in sparse matrix like data structures.
    // First conunt the size of the nonzeros of the cell_to_face data.
    int data_size=0;
    for(auto i=indexset.begin(), end=indexset.end();
        i!=end; ++i)
    {
        Opm::SparseTable<EntityRep<1> >& c2f=view_data.cell_to_face_;
        data_size+=c2f.rowSize(i->global());
    }
    
    //- cell_to_face_ : extract owner/overlap rows from cell_to_face_
    // Construct the sparse matrix like data structure.
    OrientedEntityTable<0, 1> cell_to_face;
    cell_to_face.reserve(indexset.size(), data_size);

    for(auto i=indexset.begin(), end=indexset.end();
        i!=end; ++i)
    {
        typedef boost::iterator_range<EntityRep<1>*>::iterator RowIter;
        Opm::SparseTable<EntityRep<1> >& c2f=view_data.cell_to_face_;
        auto row=c2f[i->global()];
        // create the new row, i.e. copy orientation and use new face indicator.
        std::vector<EntityRep<1> > new_row(row.size());
        std::vector<EntityRep<1> >::iterator  ncell=new_row.begin();
        for(RowIter cell=row.begin(), cend=row.end(); cell!=cend; ++cell, ++ncell)
            ncell->setValue(face_indicator[cell->index()], cell->orientation());
        // Append the new row to the matrix
        cell_to_face.appendRow(new_row.begin(), new_row.end());
    }

    // Calculate the number of nonzeros needed for the face_to_cell sparse matrix
    // To speed things up, we only use an upper limit here.
    data_size=0;
    Opm::SparseTable<EntityRep<0> >& f2c=view_data.face_to_cell_;
    for(auto f=face_indicator.begin(), fend=face_indicator.end(); f!=fend; ++f)
        if(*f<std::numeric_limits<int>::max())
            data_size += f2c.rowSize(*f);

    OrientedEntityTable<1, 0> face_to_cell;
    face_to_cell.reserve(f2c.size(), data_size);
    
    //- face_to cell_ : extract rows that connect to an existent cell
    std::vector<int> cell_indicator(view_data.cell_to_face_.size(),
                                    std::numeric_limits<int>::max());
    for(auto i=indexset.begin(), end=indexset.end(); i!=end; ++i)
        cell_indicator[i->global()]=i->local();
    
    for(auto begin=face_indicator.begin(), f=begin, fend=face_indicator.end(); f!=fend; ++f)
    {
        if(*f<std::numeric_limits<int>::max())
        {
            // face does exist
            std::vector<EntityRep<0> > new_row;
            auto old_row = f2c[f-begin];
            new_row.reserve(old_row.size());
            // push back connected existent cells.
            // for those cells we use the new cell_indicator and copy the 
            // orientation of the old cell.
            for(auto cell = old_row.begin(), cend=old_row.end();cell!=cend; ++cell)
            {
                if(cell_indicator[cell->index()]<std::numeric_limits<int>::max())
                    new_row.push_back(EntityRep<0>(cell_indicator[cell->index()], 
                                                   cell->orientation()));
            }
            face_to_cell.appendRow(new_row.begin(), new_row.end());
        }
    }
    
    // Compute the number of non zeros of the face_to_point matrix.
    data_size=0;
    for(auto f=face_indicator.begin(), fend=face_indicator.end(); f!=fend; ++f)
        if(*f<std::numeric_limits<int>::max())
            data_size += view_data.face_to_point_.rowSize(*f);

    
    Opm::SparseTable<int> face_to_point;
    face_to_point.reserve(view_data.face_to_point_.size(), data_size);

    //- face_to_point_ : extract row associated with existing faces_
    for(auto begin=face_indicator.begin(), f=begin, fend=face_indicator.end(); f!=fend; ++f)
    {
        if(*f<std::numeric_limits<int>::max())
        {
            // face does exist
            std::vector<int> new_row;
            auto old_row = view_data.face_to_point_[f-begin];
            new_row.reserve(old_row.size());
            for(auto point = old_row.begin(), pend=old_row.end(); point!=pend; ++point)
            {
                assert(point_indicator[*point]<std::numeric_limits<int>::max());
                new_row.push_back(point_indicator[*point]);
            }
            face_to_point.appendRow(new_row.begin(), new_row.end());
        }
    }
    
    // \TODO - logical_cartesian_size ????
    // - unique_boundary_ids_ : extract the ones that correspond existent faces
    EntityVariable<int, 1> unique_boundary_ids;
    
    if(view_data.unique_boundary_ids_.size())
    {
        // Unique boundary ids are inherited from the global grid.
        auto id=view_data.unique_boundary_ids_.begin();
        unique_boundary_ids.reserve(view_data.face_to_cell_.size());
        for(auto f=face_indicator.begin(), fend=face_indicator.end(); f!=fend; ++f, ++id)
        {
            if(*f<std::numeric_limits<int>::max())
            {
                unique_boundary_ids.push_back(*id);
            }
        }
    }

    // TODO setup index set and id set!!!
    
    std::unique_ptr<CpGridData> grid_data;
    return grid_data;
}

} // end namespace cpgrid
} // end namespace Dune
