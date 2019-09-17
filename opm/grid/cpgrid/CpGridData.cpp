#include"config.h"
#include <algorithm>
#include <map>
#include <vector>
#include"CpGridData.hpp"
#include"DataHandleWrappers.hpp"
#include"Intersection.hpp"
#include"Entity.hpp"
#include"OrientedEntityTable.hpp"
#include"Indexsets.hpp"
#include"PartitionTypeIndicator.hpp"

// Warning suppression for Dune includes.
#include <opm/grid/utility/platform_dependent/disable_warnings.h>

#include <opm/grid/common/GridPartitioning.hpp>
#include <dune/common/parallel/remoteindices.hh>
#include <dune/common/enumset.hh>
#include <opm/grid/utility/SparseTable.hpp>

#include <opm/grid/utility/platform_dependent/reenable_warnings.h>
#include <opm/grid/CpGrid.hpp>

namespace Dune
{
namespace cpgrid
{


CpGridData::CpGridData(const CpGridData& g)
    : index_set_(new IndexSet(*this)), local_id_set_(new IdSet(*this)),
      global_id_set_(new GlobalIdSet(local_id_set_)), partition_type_indicator_(new PartitionTypeIndicator(*this)), ccobj_(g.ccobj_)
{
#if HAVE_MPI
    ccobj_=CollectiveCommunication(MPI_COMM_SELF);
    cell_interfaces_=std::make_tuple(Interface(ccobj_),Interface(ccobj_),Interface(ccobj_),Interface(ccobj_),Interface(ccobj_));
#endif
}

CpGridData::CpGridData()
    : index_set_(new IndexSet(*this)), local_id_set_(new IdSet(*this)),
      global_id_set_(new GlobalIdSet(local_id_set_)), partition_type_indicator_(new PartitionTypeIndicator(*this)),
      ccobj_(Dune::MPIHelper::getCommunicator()), use_unique_boundary_ids_(false)
{
#if HAVE_MPI
    ccobj_=CollectiveCommunication(MPI_COMM_SELF);
    cell_interfaces_=std::make_tuple(Interface(ccobj_),Interface(ccobj_),Interface(ccobj_),Interface(ccobj_),Interface(ccobj_));
#endif
}

#if HAVE_MPI
CpGridData::CpGridData(MPI_Comm comm)
    : index_set_(new IndexSet(*this)), local_id_set_(new IdSet(*this)),
      global_id_set_(new GlobalIdSet(local_id_set_)), partition_type_indicator_(new PartitionTypeIndicator(*this)),
      ccobj_(comm), use_unique_boundary_ids_(false)
{
    cell_interfaces_=std::make_tuple(Interface(ccobj_),Interface(ccobj_),Interface(ccobj_),Interface(ccobj_),Interface(ccobj_));
}
#endif

CpGridData::CpGridData(CpGrid&)
  : index_set_(new IndexSet(*this)),   local_id_set_(new IdSet(*this)),
    global_id_set_(new GlobalIdSet(local_id_set_)),  partition_type_indicator_(new PartitionTypeIndicator(*this)),
    ccobj_(Dune::MPIHelper::getCommunicator()), use_unique_boundary_ids_(false)
{
#if HAVE_MPI
    //ccobj_=CollectiveCommunication(MPI_COMM_SELF);
    cell_interfaces_=std::make_tuple(Interface(ccobj_),Interface(ccobj_),Interface(ccobj_),Interface(ccobj_),Interface(ccobj_));
#endif
}

#if HAVE_MPI
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
void freeInterfaces(std::tuple<InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap>&
                    interfaces)
{
    freeInterfaces(std::get<0>(interfaces));
    freeInterfaces(std::get<1>(interfaces));
    freeInterfaces(std::get<2>(interfaces));
    freeInterfaces(std::get<3>(interfaces));
    freeInterfaces(std::get<4>(interfaces));
}
#endif

CpGridData::~CpGridData()
{
#if HAVE_MPI
    // code deactivated, because users cannot access face indices and therefore
    // communication on faces makes no sense!
    //freeInterfaces(face_interfaces_);
    freeInterfaces(point_interfaces_);
#endif
    delete index_set_;
    delete local_id_set_;
    delete global_id_set_;
    delete partition_type_indicator_;
}

void CpGridData::populateGlobalCellIndexSet()
{
    cell_indexset_.beginResize();
    for (int index = 0, end = size(0); index != end ; ++index){
        cell_indexset_.add(global_id_set_->id(EntityRep<0>(index, true)),
                           ParallelIndexSet::LocalIndex(index, AttributeSet::owner, true));
    }
    cell_indexset_.endResize();
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

#if HAVE_MPI

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
        if(*i<std::numeric_limits<int>::max())
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

/// \brief Handle for face tag, normal and boundary id
struct FaceTagNormalBIdHandle
{
    using DataType = std::tuple<face_tag,FieldVector<double, 3>,int>;
    using TagContainer = EntityVariable<enum face_tag, 1>;
    using BIdContainer = EntityVariable<int, 1>;
    using NormalContainer = SignedEntityVariable<FieldVector<double, 3>, 1>;

    FaceTagNormalBIdHandle(const TagContainer& gatherTags, const NormalContainer& gatherNormals, const BIdContainer& gatherBIds,
                           TagContainer& scatterTags,  NormalContainer& scatterNormals, BIdContainer& scatterBIds)
        : gatherTags_(gatherTags), gatherNormals_(gatherNormals), gatherBIds_(gatherBIds),
          scatterTags_(scatterTags), scatterNormals_(scatterNormals), scatterBIds_(scatterBIds)
    {}
    bool fixedsize(int, int)
    {
        return true;
    }
    bool contains(std::size_t dim, std::size_t codim)
    {
        return dim==3 && codim == 1;
    }
    template<class T>
    std::size_t size(const T&)
    {
        return 1;
    }
    template<class B, class T>
    void gather(B& buffer, const T& t)
    {
        if ( t.orientation() )
        {
            buffer.write(DataType(gatherTags_[t], gatherNormals_[t],
                                  gatherBIds_[t]));
        }
        else
        {
            auto normal = -gatherNormals_[t];
            buffer.write(DataType(gatherTags_[t], normal,
                                  gatherBIds_[t]));
        }
    }
    template<class B, class T>
    void scatter(B& buffer, T& t, std::size_t )
    {
        DataType tmp;
        buffer.read(tmp);
        scatterTags_[t]=std::get<0>(tmp);
        scatterNormals_.get(t.index()) = std::get<1>(tmp);
        scatterBIds_[t] = std::get<2>(tmp);
    }
private:
    const TagContainer& gatherTags_;
    const NormalContainer& gatherNormals_;
    const BIdContainer& gatherBIds_;
    TagContainer& scatterTags_;
    NormalContainer& scatterNormals_;
    BIdContainer& scatterBIds_;
};
/// \brief Handle for face tag, normal and boundary id
struct FaceTagNormalHandle
{
    using DataType = std::tuple<face_tag,FieldVector<double, 3> >;
    using TagContainer = EntityVariable<enum face_tag, 1>;
    using NormalContainer = SignedEntityVariable<FieldVector<double, 3>, 1>;

    FaceTagNormalHandle(const TagContainer& gatherTags, const NormalContainer& gatherNormals,
                        TagContainer& scatterTags,  NormalContainer& scatterNormals)
        : gatherTags_(gatherTags), gatherNormals_(gatherNormals),
          scatterTags_(scatterTags), scatterNormals_(scatterNormals)
    {}
    bool fixedsize(int, int)
    {
        return true;
    }
    bool contains(std::size_t dim, std::size_t codim)
    {
        return dim==3 && codim == 1;
    }
    template<class T>
    std::size_t size(const T&)
    {
        return 4; // 3 coordinates + 1 volume
    }
    template<class B, class T>
    void gather(B& buffer, const T& t)
    {
        if ( t.orientation() )
        {
            buffer.write(std::make_tuple(gatherTags_[t], gatherNormals_[t]));
        }
        else
        {
            auto normal = - gatherNormals_[t];
            buffer.write(std::make_tuple(gatherTags_[t], normal));
        }
    }
    template<class B, class T>
    void scatter(B& buffer, T& t, std::size_t )
    {
        DataType tmp;
        buffer.read(tmp);
        scatterTags_[t]=std::get<0>(tmp);
        scatterNormals_.get(t.index()) = std::get<1>(tmp);
    }
private:
    const TagContainer& gatherTags_;
    const NormalContainer& gatherNormals_;
    TagContainer& scatterTags_;
    NormalContainer& scatterNormals_;
};
struct PointGeometryHandle
{
    using DataType = double;
    using Container = EntityVariable<cpgrid::Geometry<0, 3>, 3>;

    PointGeometryHandle(const Container& gatherCont, Container& scatterCont)
        : gatherPoints_(gatherCont), scatterPoints_(scatterCont)
    {}
    bool fixedsize(int, int)
    {
        return true;
    }
    bool contains(std::size_t dim, std::size_t codim)
    {
        return dim==3 && codim == 3;
    }
    template<class T>
    std::size_t size(const T&)
    {
        return 3;
    }
    template<class B, class T>
    void gather(B& buffer, const T& t)
    {
        auto& geom = gatherPoints_[t];
        for (int i = 0; i < 3; i++)
            buffer.write(geom.center()[i]);
    }
    template<class B, class T>
    void scatter(B& buffer, const T& t, std::size_t )
    {
        using Vector = typename Geometry<0, 3>::GlobalCoordinate;
        Vector pos;

        for (int i = 0; i < 3; i++)
            buffer.read(pos[i]);
        scatterPoints_[t] = Geometry<0, 3>(pos);
    }
private:
    const Container& gatherPoints_;
    Container& scatterPoints_;
};

struct FaceGeometryHandle
{
    using DataType = double;
    using Geom = Geometry<2, 3>;
    using Container = EntityVariable<Geom, 1>;

    FaceGeometryHandle(const Container& gatherCont, Container& scatterCont)
        : gatherPoints_(gatherCont), scatterPoints_(scatterCont)
    {}
    bool fixedsize(int, int)
    {
        return true;
    }
    bool contains(std::size_t dim, std::size_t codim)
    {
        return dim==3 && codim == 1;
    }
    template<class T>
    std::size_t size(const T&)
    {
        return 4; // 3 coordinates + 1 volume
    }
    template<class B, class T>
    void gather(B& buffer, const T& t)
    {
        auto& geom = gatherPoints_[t];
        for (int i = 0; i < 3; i++)
            buffer.write(geom.center()[i]);
        buffer.write(geom.volume());
    }
    template<class B, class T>
    void scatter(B& buffer, const T& t, std::size_t )
    {
        using Vector = typename Geom::GlobalCoordinate;
        Vector pos;
        double vol;

        for (int i = 0; i < 3; i++)
            buffer.read(pos[i]);

        buffer.read(vol);
        scatterPoints_[t] = Geom(pos, vol);
    }
private:
    const Container& gatherPoints_;
    Container& scatterPoints_;
};


struct CellGeometryHandle
{
    using DataType = double;
    using Geom = Geometry<3, 3>;
    using Container = EntityVariable<Geom, 0>;

    CellGeometryHandle(const Container& gatherCont, Container& scatterCont,
                       const EntityVariable<cpgrid::Geometry<0, 3>, 3>& pointGeom,
                       const std::vector< std::array<int,8> >& cell2Points)
        : gatherCont_(gatherCont), scatterCont_(scatterCont), pointGeom_(pointGeom),
          cell2Points_(cell2Points)
    {}
    bool fixedsize(int, int)
    {
        return true;
    }
    bool contains(std::size_t dim, std::size_t codim)
    {
        return dim==3 && codim == 0;
    }
    template<class T>
    std::size_t size(const T&)
    {
        return 4; // 3 coordinates + 1 volume
    }
    template<class B, class T>
    typename std::enable_if<T::codimension != 0, void>::type
    gather(B&, const T&)
    {
        OPM_THROW(std::logic_error, "This should never throw!");
    }
    template<class B>
    void gather(B& buffer, const EntityRep<0>& t)
    {
        auto& geom = gatherCont_[t];
        for (int i = 0; i < 3; i++)
            buffer.write(geom.center()[i]);
        buffer.write(geom.volume());
    }
    template<class B, class T>
    typename std::enable_if<T::codimension != 0, void>::type
    scatter(B&, const T&, std::size_t)
    {
        OPM_THROW(std::logic_error, "This should never throw!");
    }
    template<class B>
    void scatter(B& buffer, const EntityRep<0>& t, std::size_t )
    {
        using Vector = typename Geom::GlobalCoordinate;
        Vector pos;
        double vol;

        for (int i = 0; i < 3; i++)
            buffer.read(pos[i]);

        buffer.read(vol);
        scatterCont_[t] = Geom(pos, vol, pointGeom_, cell2Points_[t.index()].data());
    }
private:
    const Container& gatherCont_;
    Container& scatterCont_;
    const EntityVariable<cpgrid::Geometry<0, 3>, 3>& pointGeom_;
    const std::vector< std::array<int,8> >& cell2Points_;
};

struct Cell2PointsDataHandle
{
    using DataType = int;
    using Vector = std::vector<std::array<int,8> >;
    Cell2PointsDataHandle(const Vector& globalCell2Points,
                          const GlobalIdSet& globalIds,
                          Vector& localCell2Points,
                          std::vector<int>& flatGlobalPoints)
        : globalCell2Points_(globalCell2Points), globalIds_(globalIds),
          localCell2Points_(localCell2Points), flatGlobalPoints_(flatGlobalPoints)
    {}
    bool fixedsize(int, int)
    {
        return true;
    }
    bool contains(std::size_t dim, std::size_t codim)
    {
        return dim==3 && codim == 0;
    }
    template<class T>
    std::size_t size(const T&)
    {
        return 8;
    }
    template<class B, class T>
    void gather(B& buffer, const T& t)
    {
        std::size_t i = t.index();
        assert(i < globalCell2Points_.size());
        const auto& points = globalCell2Points_[i];
        std::for_each(points.begin(), points.end(),
                      [&buffer, this](const int& point){
                          buffer.write(globalIds_.id(EntityRep<3>(point, true)));});
    }
    template<class B, class T>
    void scatter(B& buffer, const T& t, std::size_t )
    {
        auto i = t.index();
        auto& points = localCell2Points_[i];
        std::for_each(points.begin(), points.end(),
                      [&buffer, this](int& point){
                          buffer.read(point);
                          this->flatGlobalPoints_.push_back(point);
                      });
    }
private:
    const Vector& globalCell2Points_;
    const GlobalIdSet& globalIds_;
    Vector& localCell2Points_;
    std::vector<int>& flatGlobalPoints_;
};

template<class Table, int from>
struct RowSizeDataHandle
{
    using DataType = int;
    RowSizeDataHandle(const Table& global,
                      std::vector<int>& noEntries)
        : global_(global), noEntries_(noEntries)
    {}
    bool fixedsize(int, int)
    {
        return true;
    }
    template<class T>
    std::size_t size(const T&)
    {
        return 1;
    }
    bool contains(std::size_t dim, std::size_t codim)
    {
        return dim==3 && codim == from;
    }
    template<class B, class T>
    typename std::enable_if<T::codimension != from, void>::type
    gather(B&, const T&)
    {
        OPM_THROW(std::logic_error, "This should never throw. Only available for other codims.");
    }
    template<class B>
    void gather(B& buffer, const EntityRep<from>& t)
    {
        buffer.write(global_.rowSize(t));
    }
    template<class B, class T>
    void scatter(B& buffer, const T& t, std::size_t)
    {
        buffer.read(noEntries_[t.index()]);
    }
private:
    const Table& global_;
    std::vector<int>& noEntries_;
};

template<int from>
struct SparseTableEntity
{
    SparseTableEntity(const Opm::SparseTable<int>& table)
        : table_(table)
    {}
    int rowSize(const EntityRep<from>& index) const
    {
        return table_.rowSize(index.index());
    }
private:
    const Opm::SparseTable<int>& table_;
};

struct SparseTableDataHandle
{
    using Table = Opm::SparseTable<int>;
    using DataType = int;
    static constexpr int from = 1;
    SparseTableDataHandle(const Table& global,
                          const GlobalIdSet& globalIds,
                          Table& local,
                          const std::map<int,int>& global2Local)
        : global_(global), globalIds_(globalIds), local_(local), global2Local_(global2Local)
    {}
    bool fixedsize(int, int)
    {
        return false;
    }
    bool contains(std::size_t dim, std::size_t codim)
    {
        return dim==3 && codim == from;
    }
    template<class T>
    std::size_t size(const T& t)
    {
        return global_.rowSize(t.index());
    }
    template<class B, class T>
    void gather(B& buffer, const T& t)
    {
        const auto& entries = global_[t.index()];
        std::for_each(entries.begin(), entries.end(), [&buffer, this](const DataType& i){buffer.write(globalIds_.id(EntityRep<3>(i, true)));});
    }
    template<class B, class T>
    void scatter(B& buffer, const T& t, std::size_t )
    {
        const auto& entries = local_[t.index()];
        for (auto&& point : entries)
        {
            int i{};
            buffer.read(i);
            if ( point != std::numeric_limits<int>::max() )
            {
                // face already processed
                continue;
            }
            auto candidate = global2Local_.find(i);
            assert(candidate != global2Local_.end());
            point = candidate->second;
        }
        OPM_THROW(std::logic_error, "This should never throw!");
    }
    private:
    const Table& global_;
    const GlobalIdSet& globalIds_;
    Table& local_;
    const std::map<int,int>& global2Local_;
};

template<int from, int to>
struct OrientedEntityTableDataHandle
{
    using DataType = int;
    using Table = OrientedEntityTable<from, to>;
    using ToEntity = typename Table::ToType;
    using FromEntity = typename Table::FromType;
    OrientedEntityTableDataHandle(const Table& global, Table& local,
                                  const GlobalIdSet* globalIds = nullptr)
        : global_(global), local_(local), globalIds_(globalIds)
    {}
    bool fixedsize(int, int)
    {
        return false;
    }
    bool contains(std::size_t dim, std::size_t codim)
    {
        return dim==3 && codim == from;
    }
    template<class T>
    std::size_t size(const T&)
    {
        OPM_THROW(std::logic_error, "This should never throw!");
        return 0;
    }
    std::size_t size(const FromEntity& t)
    {
        return global_.rowSize(t);
    }
    template<class B, class T>
    void gather(B&, const T&)
    {
        OPM_THROW(std::logic_error, "This should never throw!");
    }
    template<class B>
    void gather(B& buffer, const FromEntity& t)
    {
        const auto& entries = global_.row(t);
        if (globalIds_)
        {
            std::for_each(entries.begin(), entries.end(),
                          [&buffer, this](const ToEntity& i){
                              int id = globalIds_->id(i);
                              if (!i.orientation())
                                  id = ~id;
                              buffer.write(id);});
        }
        else
        {
            std::for_each(entries.begin(), entries.end(), [&buffer](const ToEntity& i){buffer.write(i.signedIndex());});
        }
    }
    template<class B, class T>
    void scatter(B&, const T&, std::size_t)
    {
        OPM_THROW(std::logic_error, "This should never throw!");
    }
private:
    const Table& global_;
    Table local_;
    const GlobalIdSet* globalIds_;
};

class C2FDataHandle
    : public OrientedEntityTableDataHandle<0,1>
{
public:
    C2FDataHandle(const Table& global, const GlobalIdSet& globalIds, Table& local,
                  std::vector<int>& unsignedGlobalFaceIds)
        : OrientedEntityTableDataHandle<0,1>(global, local, &globalIds),
          unsignedGlobalFaceIds_(unsignedGlobalFaceIds)
    {}
    template<class B>
    void scatter(B& buffer, const FromEntity& t, std::size_t s)
    {
        auto& entries = local_.row(t);
        int i{};
        for (auto&& entry : entries)
        {
            buffer.read(i);
            if ( i < 0 )
            {
                entry = ToEntity(~i, false);
                unsignedGlobalFaceIds_.push_back(~i);
            }
            else
            {
                entry = ToEntity(i, true);
                unsignedGlobalFaceIds_.push_back(i);
            }
        }
    }
    template<class B, class T>
    void scatter(B&, const T&, std::size_t)
    {
        OPM_THROW(std::logic_error, "This should never throw!");
    }
private:
    std::vector<int>& unsignedGlobalFaceIds_;
};

template<class IndexSet>
class F2CDataHandle
    : public OrientedEntityTableDataHandle<1,0>
{
public:
    F2CDataHandle(const Table& global, Table& local,
                  const IndexSet& global2Local)
        : OrientedEntityTableDataHandle<1,0>(global, local),
          global2Local_(global2Local)
    {}
    template<class B>
    void scatter(B& buffer, const FromEntity& t, std::size_t s)
    {
        auto& entries = local_.row(t);
        int i{};
        for (auto&& entry : entries)
        {
            buffer.read(i);
            if (entry.index() != std::numeric_limits<int>::max())
            {
                // face already processed, continue to save map lookup
                continue;
            }
            bool orientation = true;
            if ( i < 0 )
            {
                i = ~i;
                orientation = false;
            }
            using IndexPair = typename IndexSet::IndexPair;
            auto candidate = std::lower_bound(global2Local_.begin(), global2Local_.end(),
                                          IndexPair(i),
                                          [](const IndexPair& p1,
                                             const IndexPair& p2){
                                              return p1.global() < p2.global();});
            if (candidate == global2Local_.end() || i != candidate->global())
            {
                // mark cell as being stored elsewhere
                i = std::numeric_limits<int>::max();
            }
            else
            {
                i = candidate->local();
            }
            entry = ToEntity(i, orientation);
        }
    }
private:
    const IndexSet& global2Local_;
};

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
    void gather(B& buffer, std::size_t i)
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
    void scatter(B& buffer, std::size_t i, std::size_t s)
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
        std::get<0>(t_)(rank, index, mine, other);
        std::get<1>(t_)(rank, index, mine, other);
        std::get<2>(t_)(rank, index, mine, other);
        std::get<3>(t_)(rank, index, mine, other);
        std::get<4>(t_)(rank, index, mine, other);
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

    typedef std::tuple<InteriorBorder, InteriorBorder,        Overlap     , Overlap,
                       AllSet<PartitionType> > SourceTuple;
    typedef std::tuple<InteriorBorder, AllSet<PartitionType>, OverlapFront, AllSet<PartitionType>,
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
             std::tuple<InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap>&
             interfaces)
{
    typedef typename std::map<int,std::pair<std::size_t,std::size_t> >::const_iterator Iter;
    const std::map<int,std::pair<std::size_t,std::size_t> >& sizeMap=sizes[i];
    InterfaceMap& interfaceMap=std::get<i>(interfaces);

    for(Iter iter=sizeMap.begin(), end =sizeMap.end(); iter!=end; ++iter)
    {
        std::pair<InterfaceInformation,InterfaceInformation>& interface=interfaceMap[iter->first];
        interface.first.reserve(iter->second.first);
        interface.second.reserve(iter->second.second);
    }
}

/**
 * \brief Reserves space for the interfaces.
 * \param sizes A vector with the sizes of the interface for each neighboring rank.
 * \param interface The communication interfaces.
 */
template<class InterfaceMap>
void reserve(const std::vector<std::map<int,std::pair<std::size_t,std::size_t> > >& sizes,
             std::tuple<InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap>&
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
                            typename std::tuple_element<i,typename Converter::SourceTuple>::type,
                            typename std::tuple_element<i,typename Converter::DestinationTuple>::type>
{
    typedef InterfaceFunctor<std::size_t, InterfaceIncrementor,
                             typename std::tuple_element<i,typename Converter::SourceTuple>::type,
                             typename std::tuple_element<i,typename Converter::DestinationTuple>::type>
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
                             typename std::tuple_element<i,typename Converter::SourceTuple>::type,
                             typename std::tuple_element<i,typename Converter::DestinationTuple>::type>
{
    typedef InterfaceFunctor<InterfaceInformation, InterfaceAdder,
                             typename std::tuple_element<i,typename Converter::SourceTuple>::type,
                             typename std::tuple_element<i,typename Converter::DestinationTuple>::type>
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
                      std::tuple<InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap>&
                      interfaces)
{
    // calculate sizes
    std::vector<std::map<int,std::pair<std::size_t,std::size_t> > > sizes(5);
    typedef std::tuple<SizeFunctor<0>,SizeFunctor<1>,SizeFunctor<2>,SizeFunctor<3>,
                       SizeFunctor<4> > SizeTuple;

    SizeTuple
    size_functor_tuple = std::make_tuple(SizeFunctor<0>(sizes[0]),
                                         SizeFunctor<1>(sizes[1]),
                                         SizeFunctor<2>(sizes[2]),
                                         SizeFunctor<3>(sizes[3]),
                                         SizeFunctor<4>(sizes[4]));
    InterfaceTupleFunctor<SizeTuple> size_functor(size_functor_tuple);
    iterate_over_attributes(attributes, partition_type_iterator, size_functor);
    // reserve space
    reserve(sizes, interfaces);
    // add indices to the interface
    typedef std::tuple<AddFunctor<0>,AddFunctor<1>,AddFunctor<2>,AddFunctor<3>,
                       AddFunctor<4> > AddTuple;
    AddTuple
        add_functor_tuple(AddFunctor<0>(std::get<0>(interfaces)),
                          AddFunctor<1>(std::get<1>(interfaces)),
                          AddFunctor<2>(std::get<2>(interfaces)),
                          AddFunctor<3>(std::get<3>(interfaces)),
                          AddFunctor<4>(std::get<4>(interfaces)));
    InterfaceTupleFunctor<AddTuple> add_functor(add_functor_tuple);
    iterate_over_attributes(attributes, partition_type_iterator, add_functor);

}

#endif // #if HAVE_MPI

void CpGridData::computeGeometry(CpGrid& grid,
                                 const DefaultGeometryPolicy&  globalGeometry,
                                 const OrientedEntityTable<0, 1>& globalCell2Faces,
                                 const std::vector< std::array<int,8> >& globalCell2Points,
                                 DefaultGeometryPolicy& geometry,
                                 const OrientedEntityTable<0, 1>& cell2Faces,
                                 const std::vector< std::array<int,8> >& cell2Points)
{
    FaceGeometryHandle faceGeomHandle(globalGeometry.geomVector(std::integral_constant<int,1>()),
                                      geometry.geomVector(std::integral_constant<int,1>()));
    FaceViaCellHandleWrapper<FaceGeometryHandle>
        wrappedFaceGeomHandle(faceGeomHandle, globalCell2Faces, cell2Faces);
    grid.scatterData(wrappedFaceGeomHandle);

    PointGeometryHandle pointGeomHandle(globalGeometry.geomVector(std::integral_constant<int,3>()),
                                             geometry.geomVector(std::integral_constant<int,3>()));
    PointViaCellHandleWrapper<PointGeometryHandle>
        wrappedPointGeomHandle(pointGeomHandle, globalCell2Points, cell2Points);
    grid.scatterData(wrappedPointGeomHandle);

    CellGeometryHandle cellGeomHandle(globalGeometry.geomVector(std::integral_constant<int,0>()),
                                      geometry.geomVector(std::integral_constant<int,0>()),
                                      geometry.geomVector(std::integral_constant<int,3>()),
                                      cell2Points);
    grid.scatterData(cellGeomHandle);
}

void computeFace2Point(CpGrid& grid,
                       const OrientedEntityTable<0, 1>& globalCell2Faces,
                       const GlobalIdSet& globalIds,
                       const OrientedEntityTable<0, 1>& cell2Faces,
                       const Opm::SparseTable<int>& globalFace2Points,
                       Opm::SparseTable<int>& face2Points,
                       const std::map<int,int>& global2local,
                       std::size_t noFaces)
{
    std::vector<int> rowSizes(noFaces);
    using EntityTable = SparseTableEntity<1>;
    using RowSizeDataHandle = RowSizeDataHandle<EntityTable, 1>;
    EntityTable wrappedGlobal(globalFace2Points);
    RowSizeDataHandle rowSizeHandle(wrappedGlobal, rowSizes);
    FaceViaCellHandleWrapper<RowSizeDataHandle>
        wrappedSizeHandle(rowSizeHandle, globalCell2Faces, cell2Faces);
    grid.scatterData(wrappedSizeHandle);
    face2Points.allocate(rowSizes.begin(), rowSizes.end());
    // Use entity with index INT_MAX to mark unprocessed row entries
    for (int row = 0, size = face2Points.size(); row < size; ++row)
    {
        for (auto&& point : face2Points[row])
        {
            point = std::numeric_limits<int>::max();
        }
    }
    SparseTableDataHandle handle(globalFace2Points, globalIds, face2Points, global2local);
    FaceViaCellHandleWrapper<SparseTableDataHandle>
        wrappedHandle(handle, globalCell2Faces, cell2Faces);
    grid.scatterData(wrappedHandle);
}

template<class IndexSet>
void computeFace2Cell(CpGrid& grid,
                      const OrientedEntityTable<0, 1>& globalCell2Faces,
                      const OrientedEntityTable<0, 1>& cell2Faces,
                      const OrientedEntityTable<1, 0>& globalFace2Cells,
                      OrientedEntityTable<1, 0>& face2Cells,
                      IndexSet& global2local)
{
    std::vector<int> rowSizes(global2local.size());
    using Table = OrientedEntityTable<1,0>;
    RowSizeDataHandle<Table,1> rowSizeHandle(globalFace2Cells, rowSizes);
    FaceViaCellHandleWrapper<RowSizeDataHandle<Table,1> > wrappedSizeHandle(rowSizeHandle,
                                                                        globalCell2Faces, cell2Faces);
    grid.scatterData(wrappedSizeHandle);
    face2Cells.allocate(rowSizes.begin(), rowSizes.end());
    // Use entity with index INT_MAX to mark unprocessed row entries
    for (int row = 0, size = face2Cells.size(); row < size; ++row)
    {
        for (auto&& face : face2Cells.row(EntityRep<1>(row, true)))
        {
            face = EntityRep<0>(std::numeric_limits<int>::max(), true);
        }
    }
    F2CDataHandle<IndexSet> entryHandle(globalFace2Cells, face2Cells, global2local);
    FaceViaCellHandleWrapper<F2CDataHandle<IndexSet> > wrappedEntryHandle(entryHandle,
                                                                          globalCell2Faces, cell2Faces);
    grid.scatterData(wrappedEntryHandle);
#ifndef NDEBUG
    for (int row = 0, size = face2Cells.size(); row < size; ++row)
    {
        bool oneValid = false;
        for (auto&& face : face2Cells.row(EntityRep<1>(row, true)))
        {
            oneValid = oneValid || face.index() != std::numeric_limits<int>::max();
        }
        assert(oneValid);
    }
#endif
}


std::map<int,int> computeCell2Face(CpGrid& grid,
                                    const OrientedEntityTable<0, 1>& globalCell2Faces,
                                    const GlobalIdSet& globalIds,
                                    OrientedEntityTable<0, 1>& cell2Faces,
                                    std::vector<int>& map2Global,
                                    std::size_t noCells)
{
    std::vector<int> rowSizes(noCells);
    using Table = OrientedEntityTable<0,1>;
    RowSizeDataHandle<Table,0> rowSizeHandle(globalCell2Faces, rowSizes);
    grid.scatterData(rowSizeHandle);
    cell2Faces.allocate(rowSizes.begin(), rowSizes.end());
    map2Global.reserve((noCells*6)*1.1);
    C2FDataHandle handle(globalCell2Faces, globalIds, cell2Faces,
                         map2Global);
    grid.scatterData(handle);
    // make map2Global a map from local index to global id
    std::sort(map2Global.begin(),map2Global.end());
    auto newEnd = std::unique(map2Global.begin(),map2Global.end());
    map2Global.resize(newEnd - map2Global.begin());
    // Convert face ids to local ones
    std::map<int, int> map2Local;
    auto current_global = map2Global.begin();
    int localId = 0;
    // \todo improve since we are inserting values by sorted keys.
    std::generate_n(std::inserter(map2Local, map2Local.begin()),
                    newEnd - current_global,
                    [&localId, &current_global](){ return std::make_pair(*(current_global++), localId++); });
    // translate global to local ids
    for (int row = 0, size = cell2Faces.size(); row < size; ++row)
    {
        for (auto&& face : cell2Faces.row(EntityRep<0>(row, true)))
        {
            auto index = face.signedIndex();
            if (index < 0 )
            {
                face = EntityRep<1>(map2Local[~index], false);
            }
            else
            {
                face = EntityRep<1>(map2Local[index], true);
            }
        }
    }
    return map2Local;
}


std::map<int,int> computeCell2Point(CpGrid& grid,
                                    const std::vector<std::array<int,8> >& globalCell2Points,
                                    const GlobalIdSet& globalIds,
                                    std::vector<std::array<int,8> >& cell2Points,
                                    std::vector<int>& map2Global,
                                    std::size_t noCells)
{
    cell2Points.resize(noCells);
    map2Global.reserve(noCells*8);
    Cell2PointsDataHandle handle(globalCell2Points, globalIds, cell2Points,
                                 map2Global);
    grid.scatterData(handle);
    // make map2Global a map from local index to global id
    std::sort(map2Global.begin(),map2Global.end());
    auto newEnd = std::unique(map2Global.begin(),map2Global.end());
    map2Global.resize(newEnd - map2Global.begin());
    // Convert point ids to local ones
    std::map<int, int> map2Local;
    auto current_global = map2Global.begin();
    int localId = 0;
    // \todo improve since we are inserting values by sorted keys.
    std::generate_n(std::inserter(map2Local, map2Local.begin()),
                    newEnd - current_global,
                    [&localId, &current_global](){ return std::make_pair(*(current_global++), localId++); });
    for (auto&& points : cell2Points)
    {
        for (auto&& point : points)
        {
            point = map2Local[point];
        }
    }
    return map2Local;
}


void CpGridData::distributeGlobalGrid(CpGrid& grid,
                                      const CpGridData& view_data,
                                      const std::vector<int>& /* cell_part */)
{
#if HAVE_MPI
    // setup the remote indices.
    cell_remote_indices_.setIndexSets(cell_indexset_, cell_indexset_, ccobj_);
    cell_remote_indices_.template rebuild<false>(); // We could probably also compute this on our own, like before?

    // We can identify existing cells with the help of the index set.
    // Now we need to compute the existing faces and points. Either exist
    // if they are reachable from an existing cell.
    // We use std::numeric_limits<int>::max() to indicate non-existent entities.
    std::vector<int> map2GlobalFaceId;
    std::vector<int> map2GlobalPointId;
    std::map<int,int> point_indicator =
        computeCell2Point(grid, view_data.cell_to_point_, *view_data.global_id_set_, cell_to_point_,
                          map2GlobalPointId, cell_indexset_.size());

    // create global ids array for cells. The parallel index set uses the global id
    // as the global index.
    std::vector<int> map2GlobalCellId(cell_indexset_.size());
    for(ParallelIndexSet::const_iterator i=cell_indexset_.begin(), end=cell_indexset_.end();
        i!=end; ++i)
    {
        map2GlobalCellId[i->local()]=i->global();
    }

    std::map<int,int> face_indicator =
        computeCell2Face(grid, view_data.cell_to_face_, *view_data.global_id_set_, cell_to_face_,
                         map2GlobalFaceId, cell_indexset_.size());

    auto noExistingPoints = map2GlobalPointId.size();
    auto noExistingFaces = map2GlobalFaceId.size();

    global_id_set_->swap(map2GlobalCellId, map2GlobalFaceId, map2GlobalPointId);

    computeFace2Cell(grid, view_data.cell_to_face_, cell_to_face_,
                     view_data.face_to_cell_, face_to_cell_, cell_indexset_);
    computeFace2Point(grid,  view_data.cell_to_face_, *view_data.global_id_set_, cell_to_face_,
                      view_data.face_to_point_, face_to_point_, point_indicator,
                      noExistingFaces);


    logical_cartesian_size_=view_data.logical_cartesian_size_;

    // Set up the new topology arrays
    geometry_.geomVector(std::integral_constant<int,1>()).resize(noExistingFaces);
    geometry_.geomVector(std::integral_constant<int,0>()).resize(cell_to_face_.size());
    geometry_.geomVector(std::integral_constant<int,3>()).resize(noExistingPoints);

    computeGeometry(grid, view_data.geometry_, view_data.cell_to_face_, view_data.cell_to_point_,
                    geometry_, cell_to_face_, cell_to_point_);

    global_cell_.resize(cell_indexset_.size());

    // Copy the existing cells.
    for (auto i = cell_indexset_.begin(), end = cell_indexset_.end(); i != end; ++i)
    {
        global_cell_[i->local()]=view_data.global_cell_[i->global()];
    }

    // Scatter face tags, normals, and boundary ids.
    auto noBids = view_data.unique_boundary_ids_.size();
    bool hasBids = ccobj_.max(noBids);
    face_tag_.resize(noExistingFaces);
    face_normals_.resize(noExistingFaces);

    if (hasBids)
    {
        unique_boundary_ids_.resize(noExistingFaces);
        FaceTagNormalBIdHandle faceHandle(view_data.face_tag_, view_data.face_normals_, view_data.unique_boundary_ids_,
                                          face_tag_, face_normals_, unique_boundary_ids_);
        FaceViaCellHandleWrapper<FaceTagNormalBIdHandle>
        wrappedFaceHandle(faceHandle, view_data.cell_to_face_, cell_to_face_);
        grid.scatterData(wrappedFaceHandle);
    }
    else
    {
        FaceTagNormalHandle faceHandle(view_data.face_tag_, view_data.face_normals_,
                                       face_tag_, face_normals_);
        FaceViaCellHandleWrapper<FaceTagNormalHandle>
        wrappedFaceHandle(faceHandle, view_data.cell_to_face_, cell_to_face_);
        grid.scatterData(wrappedFaceHandle);
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
    // We initialize all points with interior. Then we loop over the faces. If a face is of
    // type border, then the type of the point is overwritten with border. In the other cases
    // we set the type of the point to the one of the face as long as the type of the point is
    // not border.
    partition_type_indicator_->point_indicator_.resize(geometry_.geomVector<3>().size(),
                                                       OverlapEntity);
    for(int i=0; i<face_to_point_.size(); ++i)
    {
        for(auto p=face_to_point_[i].begin(),
                pend=face_to_point_[i].end(); p!=pend; ++p)
        {
            PartitionType new_type=partition_type_indicator_->getFacePartitionType(i);
            PartitionType old_type=PartitionType(partition_type_indicator_->point_indicator_[*p]);
            if(old_type==InteriorEntity)
            {
                if(new_type!=OverlapEntity)
                    partition_type_indicator_->point_indicator_[*p]=new_type;
            }
            if(old_type==OverlapEntity)
                partition_type_indicator_->point_indicator_[*p]=new_type;
            if(old_type==FrontEntity && new_type==BorderEntity)
                partition_type_indicator_->point_indicator_[*p]=new_type;
        }
    }

    // Compute the interface information for cells
    std::get<InteriorBorder_All_Interface>(cell_interfaces_)
        .build(cell_remote_indices_, EnumItem<AttributeSet, AttributeSet::owner>(),
               AllSet<AttributeSet>());
    std::get<Overlap_OverlapFront_Interface>(cell_interfaces_)
        .build(cell_remote_indices_, EnumItem<AttributeSet, AttributeSet::copy>(),
               EnumItem<AttributeSet, AttributeSet::copy>());
    std::get<Overlap_All_Interface>(cell_interfaces_)
        .build(cell_remote_indices_, EnumItem<AttributeSet, AttributeSet::copy>(),
                                 AllSet<AttributeSet>());
    std::get<All_All_Interface>(cell_interfaces_)
        .build(cell_remote_indices_, AllSet<AttributeSet>(), AllSet<AttributeSet>());

    // Now we use the all_all communication of the cells to compute which faces and points
    // are also present on other processes and with what attribute.
    const auto& all_all_cell_interface = std::get<All_All_Interface>(cell_interfaces_);

    // Work around a bug/deadlock in DUNE <=2.5.1 which happens if the
    // buffer cannot hold all data that needs to be send.
    // https://gitlab.dune-project.org/core/dune-common/merge_requests/416
    // For this we calculate an upper barrier of the number of
    // data items to be send manually and use it to construct a
    // VariableSizeCommunicator with sufficient buffer.
    std::size_t max_entries = 0;
    for (const auto& pair: all_all_cell_interface.interfaces() )
    {
        using std::max;
        max_entries = max(max_entries, pair.second.first.size());
        max_entries = max(max_entries, pair.second.second.size());
    }
    Dune::VariableSizeCommunicator<> comm(all_all_cell_interface.communicator(),
                                          all_all_cell_interface.interfaces(),
                                          max_entries*8*sizeof(int));
    /*
      // code deactivated, because users cannot access face indices and therefore
      // communication on faces makes no sense!
    std::vector<std::map<int,char> > face_attributes(noExistingFaces);
    AttributeDataHandle<Opm::SparseTable<EntityRep<1> > >
        face_handle(ccobj_.rank(), *partition_type_indicator_,
                    face_attributes, static_cast<Opm::SparseTable<EntityRep<1> >&>(cell_to_face_),
                    *this);
    if( std::get<All_All_Interface>(cell_interfaces_).interfaces().size() )
    {
        comm.forward(face_handle);
    }
    createInterfaces(face_attributes, FacePartitionTypeIterator(partition_type_indicator_),
                     face_interfaces_);
    std::vector<std::map<int,char> >().swap(face_attributes);
    */
    std::vector<std::map<int,char> > point_attributes(noExistingPoints);
    AttributeDataHandle<std::vector<std::array<int,8> > >
        point_handle(ccobj_.rank(), *partition_type_indicator_,
                     point_attributes, cell_to_point_, *this);
    if( static_cast<const Dune::Interface&>(std::get<All_All_Interface>(cell_interfaces_))
        .interfaces().size() )
    {
        comm.forward(point_handle);
    }
    createInterfaces(point_attributes, partition_type_indicator_->point_indicator_.begin(),
                     point_interfaces_);
#else // #if HAVE_MPI
    static_cast<void>(grid);
    static_cast<void>(view_data);
    static_cast<void>(cell_part);
    static_cast<void>(overlap_layers);
#endif
}

} // end namespace cpgrid
} // end namespace Dune
