#include <LiDIA/bigfloat.h>
#include <LiDIA/eco_prime.h>
#include <LiDIA/factorization.h>
#include <LiDIA/finite_fields/sf_gf_polynomial.h>
#include <LiDIA/gf_element.h>
#include <LiDIA/gf_polynomial.h>
#include <LiDIA/random_generator.h>
#include <LiDIA/bigint_lattice.h>
#include <LiDIA/elliptic_curve.h>
#include <LiDIA/point.h>

#include <string>
#include <iostream>
#include <sstream>
#include <array>
#include <vector>

#include <cryptopp/sha.h>
#include <cryptopp/filters.h>
#include <cryptopp/hex.h>

#pragma once

#include "gps.h"
#include "extension.h"
#include "modular.cc"

using namespace std;
using namespace LiDIA;

class torsion_basis;

class step {
    // int le; // degree of isogeny E_{i-1} -> E_i = ord(Q_i)
    // experimenting with storing extremely few shit in here...
    point<gf_element> Q_i; // necessary to store...? the eternal question

    // nota bene... experimentally NOT storing this on-hand. revisit!
    // the below are defined over F_p^2!
    polynomial<gf_element> psi; // division polynomial
    polynomial<gf_element> phi;
    polynomial<gf_element> omega; // omega(x, y) = omega(x) * y

  public:
    elliptic_curve<gf_element> E_i;
    bigint_matrix I_i;

    step(const elliptic_curve<gf_element> &, const bigint_matrix &);
    step(const torsion_basis &, const point<gf_element> &, const array<gf_element, 2> &,
         const step &);
    void isogeny(array<gf_element, 2> &, const field_extension &);
    // it is fed an element of E_{i-1}. inplace

    void generate_ideal(const gf_element &, const torsion_basis &);

    gf_element j_invariant() { return E_i.j_invariant(); } // rename this crap
};

class path : public vector<step> {
  public:
    step dummy;

    path(const elliptic_curve<gf_element> &E_) : dummy(E_, bigint_matrix(4, 4)) {
        // not bothering with n. I_ will only be _read_ from, once, into a lattice...
        bigint_matrix I_(4, 4);
        for (int i = 0; i < 4; i++)
            I_.sto(i, i, 1);
        dummy = step(E_, I_);
    }
    step &back() { // not a reference or const?
        if (empty())
            return dummy;
        else
            return vector<step>::back();
    }
};

class walker {
    int lambda;
    const elliptic_curve<gf_element> &E_;
    // const galois_field &K; // all references?
    const gf_element &iota; // also serve as a ref to K?

    const array<vector<torsion_basis>, 2> &P;

    path random;
    path ideal; // consider tinkering with names here

    string w(bool challenge);
    void random_step(const torsion_basis &);
    void ideal_step(const gf_element &, const torsion_basis &, const bigint_matrix &);
    void reroute();
    void reset();

  public:
    walker(int lambda, const elliptic_curve<gf_element> &, const gf_element &,
           const array<vector<torsion_basis>, 2> &);
    gf_element public_key() { return random.back().j_invariant(); }
    string sign(string message);
    bool verify(string message, string signature);
};
