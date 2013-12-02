#include"config.h"
#include <algorithm>
#include <map>
#include <vector>
#include"CpGridData.hpp"
#include"Intersection.hpp"
#include"Entity.hpp"
#include"OrientedEntityTable.hpp"
#include"Indexsets.hpp"
#include"PartitionTypeIndicator.hpp"
#include <dune/grid/common/GridPartitioning.hpp>
#include <dune/common/parallel/remoteindices.hh>
#include <dune/common/enumset.hh>
#include <opm/core/utility/SparseTable.hpp>

namespace Dune
{
namespace cpgrid
{


CpGridData::CpGridData(const CpGridData& g)
  : partition_type_indicator_(new PartitionTypeIndicator(*this)), ccobj_(g.ccobj_)
{}

CpGridData::CpGridData()
    : index_set_(new IndexSet(*this)), local_id_set_(new IdSet(*this)),
      global_id_set_(new GlobalIdSet(local_id_set_)), partition_type_indicator_(new PartitionTypeIndicator(*this)),
      ccobj_(Dune::MPIHelper::getCommunicator()), use_unique_boundary_ids_(false)
{}

#if HAVE_MPI
CpGridData::CpGridData(MPI_Comm comm)
    : index_set_(new IndexSet(*this)), local_id_set_(new IdSet(*this)),
      global_id_set_(new GlobalIdSet(local_id_set_)), partition_type_indicator_(new PartitionTypeIndicator(*this)),
      ccobj_(comm), use_unique_boundary_ids_(false)
{}
#endif

CpGridData::CpGridData(CpGrid& grid)
  : index_set_(new IndexSet(*this)),   local_id_set_(new IdSet(*this)),
    global_id_set_(new GlobalIdSet(local_id_set_)),  partition_type_indicator_(new PartitionTypeIndicator(*this)),
    ccobj_(Dune::MPIHelper::getCommunicator()), use_unique_boundary_ids_(false)
{}


template<class InterfaceMap>
void freeInterfaces(InterfaceMap& map)
{
    typedef typename InterfaceMap::iterator Iter;
    for(Iter i=map.begin(), end=map.end(); i!=end; ++i)
    {
        i->second.first.free();
        i->second.second.free();
    }
}

template<class InterfaceMap>
void freeInterfaces(tuple<InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap>&
                    interfaces)
{
    freeInterfaces(get<0>(interfaces));
    freeInterfaces(get<1>(interfaces));
    freeInterfaces(get<2>(interfaces));
    freeInterfaces(get<3>(interfaces));
    freeInterfaces(get<4>(interfaces));
}

CpGridData::~CpGridData()
{
    freeInterfaces(face_interfaces);
    freeInterfaces(point_interfaces);
    delete index_set_;
    delete local_id_set_;
    delete global_id_set_;
    delete partition_type_indicator_;
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
    CountExistent() : count() {}

    void operator()(int& i)
    {
        if(i < std::numeric_limits<int>::max())
            count++;
    }
    int count;
};

struct AssignAndIncrement
{
    AssignAndIncrement() : i_(){}
    void operator()(int& val){ if(val<std::numeric_limits<int>::max()) val=i_++; }
    int i_;
} assigner;

/** 
 * @brief Counts the number of global ids and sets them up.
 * @param indicator A vector indicating whether an entity exists.
 * @param ids A vector to the the global ids in.
 * @param idSet The idSet of the global grid.
 * @return the number of entities that exist.
 */
template<int codim>
int setupAndCountGlobalIds(const std::vector<int>& indicator, std::vector<int>& ids,
                           const IdSet& idSet)
{
    int count = std::count_if(indicator.begin(),
                              indicator.end(), 
                              std::bind2nd(std::less<int>(),
                                           std::numeric_limits<int>::max()));
    ids.resize(count);
    typedef typename std::vector<int>::const_iterator VIter;
    for(VIter ibegin=indicator.begin(), i=ibegin, iend= indicator.end(); 
        i!=iend; ++i)
    {
        ids[*i]=idSet.id(EntityRep<codim>(i-ibegin,true));
    }
    return count;
}

template<class T>
struct GetRowType
{};

template<class T>
struct GetRowType<Opm::SparseTable<T> >
{
    typedef typename Opm::SparseTable<T>::row_type type;
};
template<class E, class A>
struct GetRowType<std::vector<E,A> >
{
    typedef typename std::vector<E,A>::value_type type;
};

PartitionType getPartitionType(const PartitionTypeIndicator& p, const EntityRep<1>& f,
                               const CpGridData&)
{
    return p.getPartitionType(f);
}

PartitionType getPartitionType(const PartitionTypeIndicator& p, int i,
                               const CpGridData& grid)
{
    return p.getPartitionType(Entity<3>(grid, i, true));
}

int getIndex(const int* i)
{
    return *i;
}

template<class T>
int getIndex(T i)
{
    return i->index();
}

template<class T>
struct AttributeDataHandle
{
    typedef std::pair<int,char> DataType;
    
    AttributeDataHandle(int rank, const PartitionTypeIndicator& indicator,
                        std::vector<std::map<int, char> >& vals,
                        const T& cell_to_entity,
                        const CpGridData& grid)
        : rank_(rank), indicator_(indicator), vals_(vals),
        c2e_(cell_to_entity), grid_(grid)
    {}
    bool fixedsize()
    {
        return true;
    }
    std::size_t size(std::size_t i)
    {
        return c2e_[i].size();
    }
    template<class B>
    void gather(B buffer, std::size_t i)
    {
        typedef typename GetRowType<T>::type::const_iterator RowIter;
        for(RowIter f=c2e_[i].begin(), fend=c2e_[i].end();
            f!=fend; ++f)
        {
            char t=getPartitionType(indicator_, *f, grid_);
            buffer.write(std::make_pair(rank_, t));
        }
    }

    template<class B>
    void scatter(B buffer, std::size_t i, std::size_t s)
    {
        typedef typename GetRowType<T>::type::const_iterator RowIter;
        for(RowIter f=c2e_[i].begin(), fend=c2e_[i].end();
            f!=fend; ++f, --s)
        {
            std::pair<int,char> rank_attr;
            buffer.read(rank_attr);
            vals_[getIndex(f)].insert(rank_attr);
        }
    }
    int rank_;
    const PartitionTypeIndicator& indicator_;
    std::vector<std::map<int, char> >& vals_;
    const T& c2e_;
    const CpGridData& grid_;
};
 

template<class T, class Functor, class FromSet, class ToSet>
struct InterfaceFunctor
{
    InterfaceFunctor(std::map<int,std::pair<T,T> >& m)
        : map_(m)
    {}
    void operator()(int rank, std::size_t index, PartitionType mine, PartitionType other)
    {
        if(from.contains(mine) && to.contains(other))
            func(map_[rank].first, index);
        if(from.contains(other) && to.contains(mine))
            func(map_[rank].second, index);
    }
    std::map<int,std::pair<T,T> >& map_;
    FromSet from;
    ToSet   to;
    Functor func;
};

struct InterfaceIncrementor
{
    template<class T>
    void operator()(T& t, std::size_t)
    {
        ++t;
    }    
};

struct InterfaceAdder
{
    void operator()(InterfaceInformation& info, std::size_t index)
    {
        info.add(index);
    }
};

template<class Tuple>
struct InterfaceTupleFunctor
{
    InterfaceTupleFunctor(Tuple& t)
    : t_(t)
    {}
    
    void operator()(int rank, std::size_t index, PartitionType mine, PartitionType other)
    {
        get<0>(t_)(rank, index, mine, other);
        get<1>(t_)(rank, index, mine, other);
        get<2>(t_)(rank, index, mine, other);
        get<3>(t_)(rank, index, mine, other);
        get<4>(t_)(rank, index, mine, other);
    }
    Tuple& t_;
};

struct Converter
{
    typedef EnumItem<PartitionType, InteriorEntity> Interior;
    typedef EnumItem<PartitionType, BorderEntity> Border;
    typedef EnumItem<PartitionType, OverlapEntity> Overlap;
    typedef EnumItem<PartitionType, FrontEntity> Front;
    typedef Combine<Interior,Border> InteriorBorder;
    typedef Combine<Overlap,Front> OverlapFront;

    typedef tuple<InteriorBorder, InteriorBorder,        Overlap     , Overlap, 
                  AllSet<PartitionType> > SourceTuple;
    typedef tuple<InteriorBorder, AllSet<PartitionType>, OverlapFront, AllSet<PartitionType>,
                  AllSet<PartitionType> > DestinationTuple;
};

/**
 * \brief Reserves space for an interface.
 * \tparam i The index of the interface.
 * \param sizes A vector with the sizes of the interface for each neighboring rank.
 * \param interface The communication interfaces.
 */
template<std::size_t i, class InterfaceMap>
void reserve(const std::vector<std::map<int,std::pair<std::size_t,std::size_t> > >& sizes,
             tuple<InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap>&
             interfaces)
{
    typename InterfaceMap::iterator iiter=get<i>(interfaces).begin();
    typedef typename std::map<int,std::pair<std::size_t,std::size_t> >::const_iterator Iter;
    for(Iter iter=sizes[i].begin(), end =sizes[i].end(); iter!=end; ++iter, ++iiter)
    {
        iiter->second.first.reserve(iter->second.first);
        iiter->second.second.reserve(iter->second.second);
    }
}

/**
 * \brief Reserves space for the interfaces.
 * \param sizes A vector with the sizes of the interface for each neighboring rank.
 * \param interface The communication interfaces.
 */
template<class InterfaceMap>
void reserve(const std::vector<std::map<int,std::pair<std::size_t,std::size_t> > >& sizes,
             tuple<InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap>&
             interfaces)
{
    reserve<0>(sizes, interfaces);
    reserve<1>(sizes, interfaces);
    reserve<2>(sizes, interfaces);
    reserve<3>(sizes, interfaces);
    reserve<4>(sizes, interfaces);
}

/**
 * \brief A functor that calculates the size of the interface.
 * \tparam i The indentifier of the interface.
 */
template<std::size_t i>
struct SizeFunctor : 
    public InterfaceFunctor<std::size_t, InterfaceIncrementor,
                            typename tuple_element<i,typename Converter::SourceTuple>::type,
                            typename tuple_element<i,typename Converter::DestinationTuple>::type>
{
    typedef InterfaceFunctor<std::size_t, InterfaceIncrementor,
                             typename tuple_element<i,typename Converter::SourceTuple>::type,
                             typename tuple_element<i,typename Converter::DestinationTuple>::type>
    Base;
    SizeFunctor(std::map<int,std::pair<std::size_t,std::size_t> >& m)
        :Base(m)
    {}
};

/**
 * \brief A functor that adds indices to the interface.
 * \tparam i The indentifier of the interface.
 */
template<std::size_t i>
struct AddFunctor :
    public InterfaceFunctor<InterfaceInformation, InterfaceAdder,
                             typename tuple_element<i,typename Converter::SourceTuple>::type,
                             typename tuple_element<i,typename Converter::DestinationTuple>::type>
{
    typedef InterfaceFunctor<InterfaceInformation, InterfaceAdder,
                             typename tuple_element<i,typename Converter::SourceTuple>::type,
                             typename tuple_element<i,typename Converter::DestinationTuple>::type>
    Base;
    AddFunctor(std::map<int,std::pair<InterfaceInformation,InterfaceInformation> >& m)
        : Base(m)
    {}
};

class FacePartitionTypeIterator
{
public:
    FacePartitionTypeIterator(const PartitionTypeIndicator* part)
    : indicator_(part), index_()
    {}
    void operator++()
    {
        ++index_;
    }
    PartitionType operator*()
    {
        return indicator_->getFacePartitionType(index_);
    }
private:
    const PartitionTypeIndicator* indicator_;
    int index_;
};


/**
 * \brief Applies a functor the each pair of the interface.
 * \tparam Functor The type of the functor to apply.
 * \param attributes[in] A vector that contains for each index a map from other
 * process ranks to the attribute there.
 * \param my_attributes[in] A vector with the attributes of each index on this process.
 * \param func The functor.
 */
template<class Functor, class T>
void iterate_over_attributes(std::vector<std::map<int,char> >& attributes,
                             T my_attribute_iter, Functor& func)
{
    typedef typename std::vector<std::map<int,char> >::const_iterator Iter;
    for(Iter begin=attributes.begin(), i=begin, end=attributes.end(); i!=end; 
        ++i, ++my_attribute_iter)
    {
        typedef typename std::map<int,char>::const_iterator MIter;
        for(MIter m=i->begin(), mend=i->end(); m!=mend; ++m)
        {
            func(m->first,i-begin,PartitionType(*my_attribute_iter), PartitionType(m->second));
        }
    }
}

        
/**
 * \brief Creates the communication interface for either faces or points.
 * \param attributes[in] A vector that contains for each index a map from other
 * process ranks to the attribute there.
 * \param my_attributes[in] A vector with the attributes of each index on this process.
 * \param[out] interfaces The tuple with the interface maps for communication.
 */
template<class InterfaceMap,class T>
void createInterfaces(std::vector<std::map<int,char> >& attributes,
                      T partition_type_iterator,
                      tuple<InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap>&
                      interfaces)
{
    // calculate sizes
    std::vector<std::map<int,std::pair<std::size_t,std::size_t> > > sizes(5);
    typedef tuple<SizeFunctor<0>,SizeFunctor<1>,SizeFunctor<2>,SizeFunctor<3>,
                  SizeFunctor<4> > SizeTuple;
    SizeTuple
        size_functor_tuple = make_tuple(SizeFunctor<0>(sizes[0]),
                                        SizeFunctor<1>(sizes[1]),
                                        SizeFunctor<2>(sizes[2]),
                                        SizeFunctor<3>(sizes[3]),
                                        SizeFunctor<4>(sizes[4]));
    InterfaceTupleFunctor<SizeTuple> size_functor(size_functor_tuple);
    iterate_over_attributes(attributes, partition_type_iterator, size_functor);
    // reserve space
    reserve(sizes, interfaces);
    // add indices to the interface
    typedef tuple<AddFunctor<0>,AddFunctor<1>,AddFunctor<2>,AddFunctor<3>,
                  AddFunctor<4> > AddTuple;
    AddTuple
        add_functor_tuple(AddFunctor<0>(get<0>(interfaces)),
                          AddFunctor<1>(get<1>(interfaces)),
                          AddFunctor<2>(get<2>(interfaces)),
                          AddFunctor<3>(get<3>(interfaces)),
                          AddFunctor<4>(get<4>(interfaces)));
    InterfaceTupleFunctor<AddTuple> add_functor(add_functor_tuple);
    iterate_over_attributes(attributes, partition_type_iterator, add_functor);
                                
}


void CpGridData::distributeGlobalGrid(const CpGrid& grid,
                                      const CpGridData& view_data,
                                      const std::vector<int>& cell_part)
{
    Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator>& ccobj=ccobj_;
    std::unique_ptr<CpGridData> grid_data(new CpGridData);
    int my_rank=ccobj.rank();
#ifdef DEBUG
    int size=ccobj.size();
    typedef std::vector<int>::const_iterator Iterator;
    for(Iterator i=cell_part.begin(); i!= cell_part.end(); ++i)
        if(*i>=size)
            OPM_THROW(std::runtime_error, "rank for cell is too big");
#endif
    // vector with the set of ranks that
    std::vector<std::set<int> > overlap; 
    
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
        ParallelIndexSet* indexset;
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
                indexset->add(i, Index(count++, AttributeSet::owner, true));
            }
            else
            {
                neighbors.insert(rank);
                global2local.push_back(std::numeric_limits<int>::max());
            }
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
                indexset->add(i, Index(count++, AttributeSet::overlap, true));
                neighbors.insert(rank);
            }
            else
                global2local.push_back(std::numeric_limits<int>::max());
            
        }
    } cell_counter;
    cell_counter.myrank=my_rank;
    cell_counter.count=0;
    cell_counter.global2local.reserve(cell_part.size());
    cell_counter.indexset=&cell_indexset_;
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
    typedef RemoteIndexListModifier<RemoteIndices::ParallelIndexSet, RemoteIndices::Allocator,
                                    false> Modifier;
    typedef RemoteIndices::RemoteIndex RemoteIndex;
    RemoteIndices cell_remote_indices(cell_indexset_, cell_indexset_, ccobj_);


    // Create a map of ListModifiers
    if(cell_counter.neighbors.size()){ //extra scope to call destructor of the Modifiers
        std::map<int,Modifier> modifiers;
        for(CellCounter::NeighborIterator n=cell_counter.neighbors.begin(), end=cell_counter.neighbors.end();
            n != end; ++n)
            modifiers.insert(std::make_pair(*n, cell_remote_indices.getModifier<false,false>(*n)));
        // Insert remote indices. For each entry in the index set, see wether there are overlap occurences and add them.
        for(ParallelIndexSet::const_iterator i=cell_indexset_.begin(), end=cell_indexset_.end();
            i!=end; ++i)
        {
            std::set<int>::const_iterator iter=overlap[i->local()].find(my_rank);
            if(iter!=overlap[i->local()].end())
                modifiers[cell_part[i->local()]]
                    .insert(RemoteIndex(AttributeSet::owner,&(*i)));
        }
    }
    else
    {
        // Force update of the sync counter in the remote indices.
        Modifier* mod=new Modifier(cell_remote_indices.getModifier<false,false>(0));
        delete mod;
    }

    // We can identify existing cells with the help of the index set.
    // Now we need to compute the existing faces and points. Either exist
    // if they are reachable from an existing cell.
    // We use std::numeric_limits<int>::max() to indicate non-existent entities.
    std::vector<int> face_indicator(view_data.geometry_.geomVector<1>().size(),
                                    std::numeric_limits<int>::max());
    std::vector<int> point_indicator(view_data.geometry_.geomVector<3>().size(), 
                                     std::numeric_limits<int>::max());
    for(ParallelIndexSet::iterator i=cell_indexset_.begin(), end=cell_indexset_.end();
            i!=end; ++i)
    {
        typedef boost::iterator_range<const EntityRep<1>*>::iterator RowIter;
        int row_index=i->global();
        // Somehow g++-4.4 does not find functions of father even if we
        // change inheritance of OrientedEntityTable to public.
        // Therfore we use an ugly cast to base class here.
        const Opm::SparseTable<EntityRep<1> >& c2f=view_data.cell_to_face_;
        for(RowIter f=c2f[row_index].begin(), fend=c2f[row_index].end();
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
    
       
    // renumber  face and point indicators
    std::for_each(face_indicator.begin(), face_indicator.end(), AssignAndIncrement());
    std::for_each(point_indicator.begin(), point_indicator.end(), AssignAndIncrement());

    std::vector<int> map2GlobalFaceId;
    int noExistingFaces = setupAndCountGlobalIds<1>(face_indicator, map2GlobalFaceId,
                                                    *view_data.local_id_set_);
    std::vector<int> map2GlobalPointId;
    int noExistingPoints = setupAndCountGlobalIds<3>(point_indicator, map2GlobalPointId,
                                                     *view_data.local_id_set_);
    std::vector<int> map2GlobalCellId(cell_indexset_.size());
    for(ParallelIndexSet::const_iterator i=cell_indexset_.begin(), end=cell_indexset_.end();
        i!=end; ++i)
    {
        map2GlobalCellId[i->local()]=view_data.local_id_set_->id(EntityRep<0>(i->global(), true));
    }

    global_id_set_->swap(map2GlobalCellId, map2GlobalFaceId, map2GlobalPointId);
    // Set up the new topology arrays
    EntityVariable<cpgrid::Geometry<3, 3>, 0> cell_geom;
    std::vector<cpgrid::Geometry<3, 3> > tmp_cell_geom(cell_indexset_.size());
    auto global_cell_geom=view_data.geomVector<0>();
    std::vector<int> global_cell;
    global_cell.resize(cell_indexset_.size());

    // Copy the existing cells.
    for(auto i=cell_indexset_.begin(), end=cell_indexset_.end(); i!=end; ++i)
    {
        tmp_cell_geom[i->local()]=static_cast<std::vector<cpgrid::Geometry<3, 3> >&>(global_cell_geom)[i->global()];
        global_cell[i->local()]=view_data.global_cell_[i->global()];
    }
    static_cast<std::vector<cpgrid::Geometry<3, 3> >&>(cell_geom).swap(tmp_cell_geom);

    // count the existing faces, renumber, and allocate space.
    EntityVariable<cpgrid::Geometry<2, 3>, 1> face_geom;
    std::vector<cpgrid::Geometry<2, 3> > tmp_face_geom(noExistingFaces);
    std::vector<enum face_tag> tmp_face_tag(noExistingFaces);
    std::vector<PointType> tmp_face_normals(noExistingFaces);
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
    static_cast<std::vector<PointType>&>(face_normals_).swap(tmp_face_normals);
    static_cast<std::vector<cpgrid::Geometry<2, 3> >&>(face_geom).swap(tmp_face_geom);
    static_cast<std::vector<enum face_tag>&>(face_tag_).swap(tmp_face_tag);

    // Count the existing points and allocate space
    std::vector<cpgrid::Geometry<0, 3> > tmp_point_geom(noExistingPoints);
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
    geometry_=cpgrid::DefaultGeometryPolicy(cell_geom, face_geom, 
                                            point_geom);

    // Create the topology information. This is stored in sparse matrix like data structures.
    // First conunt the size of the nonzeros of the cell_to_face data.
    int data_size=0;
    for(auto i=cell_indexset_.begin(), end=cell_indexset_.end();
        i!=end; ++i)
    {
        const Opm::SparseTable<EntityRep<1> >& c2f=view_data.cell_to_face_;
        data_size+=c2f.rowSize(i->global());
    }
    
    //- cell_to_face_ : extract owner/overlap rows from cell_to_face_
    // Construct the sparse matrix like data structure.
    //OrientedEntityTable<0, 1> cell_to_face;
    cell_to_face_.reserve(cell_indexset_.size(), data_size);
    cell_to_point_.reserve(cell_indexset_.size());
    
    for(auto i=cell_indexset_.begin(), end=cell_indexset_.end();
        i!=end; ++i)
    {
        typedef boost::iterator_range<const EntityRep<1>*>::const_iterator RowIter;
        const Opm::SparseTable<EntityRep<1> >& c2f=view_data.cell_to_face_;
        auto row=c2f[i->global()];
        // create the new row, i.e. copy orientation and use new face indicator.
        std::vector<EntityRep<1> > new_row(row.size());
        std::vector<EntityRep<1> >::iterator  nface=new_row.begin();
        for(RowIter face=row.begin(), fend=row.end(); face!=fend; ++face, ++nface)
            nface->setValue(face_indicator[face->index()], face->orientation());
        // Append the new row to the matrix
        cell_to_face_.appendRow(new_row.begin(), new_row.end());
        for(int j=0; j<8; ++j)
            cell_to_point_[i->local()][j]=point_indicator[view_data.cell_to_point_[i->global()][j]];
    }

    // Calculate the number of nonzeros needed for the face_to_cell sparse matrix
    // To speed things up, we only use an upper limit here.
    data_size=0;
    const Opm::SparseTable<EntityRep<0> >& f2c=view_data.face_to_cell_;
    for(auto f=face_indicator.begin(), fend=face_indicator.end(); f!=fend; ++f)
        if(*f<std::numeric_limits<int>::max())
            data_size += f2c.rowSize(*f);

    face_to_cell_.reserve(f2c.size(), data_size);
    
    //- face_to cell_ : extract rows that connect to an existent cell
    std::vector<int> cell_indicator(view_data.cell_to_face_.size(),
                                    std::numeric_limits<int>::max());
    for(auto i=cell_indexset_.begin(), end=cell_indexset_.end(); i!=end; ++i)
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
                // The statement below results in all faces having two neighbours
                // except for those at the domain boundary.
                // Note that along the front partition there are invalid neighbours
                // marked with index std::numeric_limits<int>::max()
                new_row.push_back(EntityRep<0>(cell_indicator[cell->index()], 
                                                   cell->orientation()));
            }
            face_to_cell_.appendRow(new_row.begin(), new_row.end());
        }
    }
    
    // Compute the number of non zeros of the face_to_point matrix.
    data_size=0;
    for(auto f=face_indicator.begin(), fend=face_indicator.end(); f!=fend; ++f)
        if(*f<std::numeric_limits<int>::max())
            data_size += view_data.face_to_point_.rowSize(*f);

    face_to_point_.reserve(view_data.face_to_point_.size(), data_size);

    //- face_to_point__ : extract row associated with existing faces_
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
            face_to_point_.appendRow(new_row.begin(), new_row.end());
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
    // Compute the partition type for cell
    partition_type_indicator_->cell_indicator_.resize(cell_indexset_.size());
    for(ParallelIndexSet::const_iterator i=cell_indexset_.begin(), end=cell_indexset_.end();
            i!=end; ++i)
    {
        partition_type_indicator_->cell_indicator_[i->local()]=
            i->local().attribute()==AttributeSet::owner?
            InteriorEntity:OverlapEntity;
    }
    
    // Compute partition type for points
    // We initialize all point with interior. Then we loop over the faces. If a face is of 
    // type border, then the type of the point is overwritten with border. In the other cases
    // we set type of the point to the one of the face as long as the type of the point is 
    // not boder.
    partition_type_indicator_->point_indicator_.resize(geometry_.geomVector<3>().size(),
                                                      InteriorEntity);
    for(int i=0; i<face_to_point_.size(); ++i)
    {
        for(auto p=view_data.face_to_point_[i].begin(), 
                pend=view_data.face_to_point_[i].end(); p!=pend; ++p)
        {
            if(partition_type_indicator_->point_indicator_[*p]!=BorderEntity)
                partition_type_indicator_->point_indicator_[*p]=
                    partition_type_indicator_->getFacePartitionType(i);
        }
    }

    // Compute the interface information for cells
    get<InteriorBorder_All_Interface>(cell_interfaces)
        .build(cell_remote_indices, EnumItem<AttributeSet, AttributeSet::owner>(),
               AllSet<AttributeSet>());
    get<Overlap_OverlapFront_Interface>(cell_interfaces)
        .build(cell_remote_indices, EnumItem<AttributeSet, AttributeSet::overlap>(),
               EnumItem<AttributeSet, AttributeSet::overlap>());
     get<Overlap_All_Interface>(cell_interfaces)
        .build(cell_remote_indices, EnumItem<AttributeSet, AttributeSet::overlap>(),
                                 AllSet<AttributeSet>());
    get<All_All_Interface>(cell_interfaces)
        .build(cell_remote_indices, AllSet<AttributeSet>(), AllSet<AttributeSet>());

    // Now we use the all_all communication of the cells to compute which faces and points
    // are also present on other processes and with what attribute.
    std::vector<std::map<int,char> > face_attributes(noExistingFaces);
    AttributeDataHandle<Opm::SparseTable<EntityRep<1> > > 
        face_handle(ccobj_.rank(), *partition_type_indicator_,
                    face_attributes, static_cast<Opm::SparseTable<EntityRep<1> >&>(cell_to_face_),
                    *this);
    Dune::VariableSizeCommunicator<> comm(get<All_All_Interface>(cell_interfaces));
    comm.forward(face_handle);
    createInterfaces(face_attributes, FacePartitionTypeIterator(partition_type_indicator_),
                     face_interfaces);
    std::vector<std::map<int,char> >().swap(face_attributes);
    std::vector<std::map<int,char> > point_attributes(noExistingPoints);
    AttributeDataHandle<std::vector<std::array<int,8> > >
        point_handle(ccobj_.rank(), *partition_type_indicator_,
                     point_attributes, cell_to_point_, *this);
    comm.forward(point_handle);
    createInterfaces(point_attributes, partition_type_indicator_->point_indicator_.begin(),
                     point_interfaces);
    
    
}

} // end namespace cpgrid
} // end namespace Dune
