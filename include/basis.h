#include <array>

#pragma once

#include "extension.h"

class torsion_basis {
  public:
    int l;

    field_extension ext; // contains K and L
    LiDIA::elliptic_curve<LiDIA::gf_element> E_ext;
    // weirdly necessary when reading in points... revisit
    std::array<LiDIA::point<LiDIA::gf_element>, 2> P; // use the same name again...
    // should any/all of these be references?

    torsion_basis() {}
    torsion_basis(const int, const field_extension &,
                  const LiDIA::elliptic_curve<LiDIA::gf_element> &,
                  const std::array<LiDIA::point<LiDIA::gf_element>, 2> &);
};

std::istream &operator>>(std::istream &in, torsion_basis &P);
std::ostream &operator<<(std::ostream &out, const torsion_basis &P);
