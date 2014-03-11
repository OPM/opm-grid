#ifdef HAVE_CONFIG_H
#include"config.h"
#endif
#include"EntityRep.hpp"
#include"Intersection.hpp"
#include"Entity.hpp"

namespace Dune
{
namespace cpgrid
{

Intersection::Intersection(const CpGridData& grid, const EntityRep<0>& cell, int subindex, bool update_now)

		: pgrid_(&grid),
		  index_(cell.index()),
		  subindex_(subindex),
		  faces_of_cell_(grid.cell_to_face_[cell]),
		  global_geom_(cpgrid::Entity<1>(grid, faces_of_cell_[subindex_]).geometry()),
// 		  in_inside_geom_(global_geom_.center()
// 				  - cpgrid::Entity<0>(grid, index_).geometry().center(),
// 				  global_geom_.volume()),
		  nbcell_(cell.index()), // Init to self, which is invalid.
		  is_on_boundary_(false)
            {
                assert(index_ >= 0);
		if (update_now) {
		    update();
		}
            }
int Intersection::boundaryId() const
            {
                int ret = 0;
                if (boundary()) {
		    if (pgrid_->uniqueBoundaryIds()) {
			// Use the unique boundary ids.
			OrientedEntityTable<0,1>::ToType face = faces_of_cell_[subindex_];
			ret = pgrid_->unique_boundary_ids_[face];
		    } else {
			// Use the face tag based ids, i.e. 1-6 for i-, i+, j-, j+, k-, k+.
			typedef OrientedEntityTable<0,1>::ToType Face;
			const Face& f = faces_of_cell_[subindex_];
			const bool normal_is_in = !f.orientation();
			enum face_tag tag = pgrid_->face_tag_[f];

			switch (tag) {
			case LEFT:
			    //                   LEFT : RIGHT
			    ret = normal_is_in ? 1    : 2; // min(I) : max(I)
			    break;
			case BACK:
			    //                   BACK : FRONT
			    ret = normal_is_in ? 3    : 4; // min(J) : max(J)
			    break;
			case TOP:
			    // Note: TOP at min(K) as 'z' measures *depth*.
			    //                   TOP  : BOTTOM
			    ret = normal_is_in ? 5    : 6; // min(K) : max(K)
			    break;
			}
		    }
                }
                return ret;
            }
int Intersection::boundarySegmentIndex() const
            {
                // Since this is almost the same that we did for
                // 'unique boundary ids' we use those numbers, although since
                // they are 1-based and not 0-based we must be careful.
                if (!boundary()) {
                    OPM_THROW(std::runtime_error, "Cannot call boundarySegmentIndex() on non-boundaries.");
                }
                assert(!pgrid_->unique_boundary_ids_.empty());
                // Use the unique boundary ids (subtract 1).
                const EntityRep<1>& face = faces_of_cell_[subindex_];
                return pgrid_->unique_boundary_ids_[face] - 1;
            }
void Intersection::update()
	    {
                const EntityRep<1>& face = faces_of_cell_[subindex_];
		//global_geom_ = cpgrid::Entity<1>(*pgrid_, face).geometry();
                global_geom_ = pgrid_->geometry_.geomVector<1>()[face];
                OrientedEntityTable<1,0>::row_type cells_of_face = pgrid_->face_to_cell_[face];
		is_on_boundary_ = cells_of_face.size() == 1;
                // Wether there is no nother nbcell for this intersection
                // i.e. either this on the boundary or a front intersection
                bool has_no_nbcell = is_on_boundary_ ||
                    cells_of_face[0].index()==std::numeric_limits<int>::max() ||
                    cells_of_face[1].index()==std::numeric_limits<int>::max();
                if (has_no_nbcell) {
		    nbcell_ = index_; // self is invalid value
		} else {
                    assert(cells_of_face.size() == 2);
                    if (cells_of_face[0].index() == index_) {
                        nbcell_ = cells_of_face[1].index();
                    } else {
                        nbcell_ = cells_of_face[0].index();
                    }
// 		    in_outside_geom_ = LocalGeometry(global_geom_.center()
// 						     - outside().geometry().center(),
// 						     global_geom_.volume());
		}
            }

void Intersection::increment()
	    {
		++subindex_;
		if (subindex_ < faces_of_cell_.size()) {
		    update();
		}
	    }


Intersection::EntityPointer Intersection::inside() const
{
    return EntityPointer(*pgrid_, index_, true);
}

Intersection::EntityPointer Intersection::outside() const
{
    return EntityPointer(*pgrid_, nbcell(), true);
}
} // end namespace cpgrid
} // end namespace Dune 
