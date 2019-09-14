#include <array>

#pragma once

#include "basis.h"
#include "step.h"

class maximal_order {    // keeps track of info required for the endomorphism theta.
    step step;           // includes psi, phi and omega
    LiDIA::gf_element u; // Silverman, III, Prop. 1.4

  public:
    LiDIA::bigint c;

    LiDIA::bigint p;
    LiDIA::bigint q;
    LiDIA::bigint q_coeff;
    LiDIA::bigint c_coeff;

    maximal_order() {}
    maximal_order(const class step &, const LiDIA::gf_element &);
    void set_c(const LiDIA::bigint &);

    LiDIA::point<LiDIA::gf_element> theta(const LiDIA::point<LiDIA::gf_element> &Q,
                                          const field_extension &ext) const;
    LiDIA::point<LiDIA::gf_element> pi(const LiDIA::point<LiDIA::gf_element> &Q) const;

    LiDIA::bigint norm(const LiDIA::base_vector<LiDIA::bigint> &e) const;
    LiDIA::base_vector<LiDIA::bigint> conjugate(const LiDIA::base_vector<LiDIA::bigint> &e) const;
    LiDIA::bigint_matrix times(const LiDIA::base_vector<LiDIA::bigint> &e) const;

    friend std::istream &operator>>(std::istream &, maximal_order &);
    friend std::ostream &operator<<(std::ostream &, const maximal_order &);
};

std::istream &operator>>(std::istream &in, maximal_order &);
std::ostream &operator<<(std::ostream &out, const maximal_order &);
