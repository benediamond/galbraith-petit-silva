#include <array>

#pragma once

#include "extension.h"

using namespace std;
using namespace LiDIA;

class torsion_basis {
  public:
    int l;

    field_extension ext;              // contains K and L
    elliptic_curve<gf_element> E_ext; // weirdly necessary when reading in points... revisit
    array<point<gf_element>, 2> P;    // use the same name again...
    // should any/all of these be references?

    torsion_basis() {}
    torsion_basis(const int, const field_extension &, const elliptic_curve<gf_element> &,
                  const array<point<gf_element>, 2> &);
};

istream &operator>>(istream &in, torsion_basis &P);
ostream &operator<<(ostream &out, const torsion_basis &P);
