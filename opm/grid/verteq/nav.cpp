// Copyright (C) 2013 Uni Research AS
// This file is licensed under the GNU General Public License v3.0
#include <opm/grid/verteq/nav.hpp>
#include <cstdlib>
#include <iostream>
using namespace std;

const Dir Dir::DEC (0);
const Dir Dir::INC (1);

const Dim2D Dim2D::X (0);
const Dim2D Dim2D::Y (1);
const Dim3D Dim3D::Z (2);

template <typename Dim> Side <Dim>
Side <Dim>::from_tag (int tag) {
	// direction is the minor bits in the enumeration
	const div_t bits = div (tag, Dir::COUNT);
	const int dir_val = bits.rem;
	const int dim_val = bits.quot;
	return Side <Dim> (Dim (dim_val), Dir (dir_val));
}

// template instantiation to satisfy linker
template Side <Dim2D> Side <Dim2D>::from_tag (int);
template Side <Dim3D> Side <Dim3D>::from_tag (int);

// these needs to be initialized here instead of in the header
// because vector::resize takes a reference to the data and not
// a value as a parameter (to avoid copying)
const int Cart2D::NO_ELEM = -1;
const int Cart2D::NO_FACE = -1;
const int Cart2D::NO_NODE = -1;

// enumeration of all possible sides
template <>
const Side <Dim2D> Side <Dim2D>::ALL[] = {
	Side (Dim2D::X, Dir::DEC), // I-
	Side (Dim2D::X, Dir::INC), // I+
	Side (Dim2D::Y, Dir::DEC), // J-
	Side (Dim2D::Y, Dir::INC)  // J+
};

template <>
const Side <Dim3D> Side <Dim3D>::ALL[] = {
	Side (Dim2D::X, Dir::DEC), // I-
	Side (Dim2D::X, Dir::INC), // I+
	Side (Dim2D::Y, Dir::DEC), // J-
	Side (Dim2D::Y, Dir::INC), // J+
	Side (Dim3D::Z, Dir::DEC), // K-
	Side (Dim3D::Z, Dir::INC)  // K+
};

// sides that exists as standalone constants
const Side3D UP   (Dim3D::Z, Dir::DEC);
const Side3D DOWN (Dim3D::Z, Dir::INC);

// print carriers (for debugging)
ostream& operator << (ostream& os, const Coord2D& c) {
	return (os << '(' << c.i () << ',' << c.j () << ')'); // e.g. "(3,2)"
}

ostream& operator << (ostream& os, const Coord3D& c) {
	// e.g. "(3,2,5)"
	return (os << '(' << c.i () << ',' << c.j () << ',' << c.k () << ')');
}

static const char DIR_NAMES[] = {'-', '+'};
static const char DIM_NAMES[] = {'I', 'J', 'K'};

ostream& operator << (ostream& os, const Dir& d) {
	return (os << DIR_NAMES[d.val]); // e.g. '-'
}

ostream& operator << (ostream& os, const Dim2D& d) {
	return (os << DIM_NAMES[d.val]); // e.g. 'I'
}

ostream& operator << (ostream& os, const Side2D& s) {
	return (os << s.dim () << s.dir ()); // e.g. "I-"
}

ostream& operator << (ostream& os, const Side3D& s) {
	return (os << s.dim () << s.dir ()); // e.g. "I-"
}

ostream& operator << (ostream& os, const Corn3D& c) {
	// e.g. "(I-,J+)"
	return (os << '(' << DIM_NAMES[Dim3D::X.val] << c.i () << ','
										<< DIM_NAMES[Dim3D::Y.val] << c.j () << ','
										<< DIM_NAMES[Dim3D::Z.val] << c.k () << ')');
}
