#include <array>

#pragma once

#include "basis.h"
#include "step.h"

using namespace std;
using namespace LiDIA;

class maximal_order { // keeps track of info required for the endomorphism theta.
    step step;        // includes psi, phi and omega
    gf_element u;     // Silverman, III, Prop. 1.4

  public:
    bigint c;

    bigint p;
    bigint q;
    bigint q_coeff;
    bigint c_coeff;

    maximal_order() {}
    maximal_order(const class step &, const gf_element &);
    void set_c(const bigint &);

    point<gf_element> theta(const point<gf_element> &Q, const field_extension &ext) const;
    point<gf_element> pi(const point<gf_element> &Q) const;

    bigint norm(const base_vector<bigint> &e) const;
    base_vector<bigint> conjugate(const base_vector<bigint> &e) const;
    bigint_matrix times(const base_vector<bigint> &e) const;

    friend istream &operator>>(istream &, maximal_order &);
    friend ostream &operator<<(ostream &, const maximal_order &);
};

istream &operator>>(istream &in, maximal_order &);
ostream &operator<<(ostream &out, const maximal_order &);
