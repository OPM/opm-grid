#include"CpGridData.hpp"
#include"Intersection.hpp"
#include"Entity.hpp"
#include"Indexsets.hpp"

namespace Dune
{
namespace cpgrid
{
CpGridData::CpGridData()
    : index_set_(), local_id_set_()
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
} // end namespace cpgrid
} // end namespace Dune
