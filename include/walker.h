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

#include <cryptopp/sha3.h>
#include <cryptopp/filters.h>
#include <cryptopp/hex.h>

#pragma once

#include "extension.h"
#include "modular.h"
#include "order.h"
#include "step.h"

using namespace std;
using namespace LiDIA;

class walker {
    int lambda;
    const elliptic_curve<gf_element> &E_;
    const array<vector<torsion_basis>, 2> &P;
    const maximal_order &O_;

    path random;
    path ideal; // consider tinkering with names here

    string w(bool challenge);

    void update_ideal(const torsion_basis &, const point<gf_element> &);
    void random_step(const torsion_basis &);
    void ideal_step(const torsion_basis &, const bigint_matrix &);

    void reroute();
    void reset();

  public:
    walker(int, const elliptic_curve<gf_element> &, const array<vector<torsion_basis>, 2> &,
           const maximal_order &);
    gf_element public_key() { return random.back().E_i.j_invariant(); }
    // ^^^ requires that we're "between walks"...
    string sign(string message);
    bool verify(string message, string signature);
};
