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

class walker {
    int lambda;
    const LiDIA::elliptic_curve<LiDIA::gf_element> &E_;
    const std::array<std::vector<torsion_basis>, 2> &P;
    const maximal_order &O_;

    path random;
    path ideal; // consider tinkering with names here

    std::string w(bool challenge);

    void update_ideal(const torsion_basis &, const LiDIA::point<LiDIA::gf_element> &);
    void random_step(const torsion_basis &);
    void ideal_step(const torsion_basis &, const LiDIA::bigint_matrix &);

    void reroute();
    void reset();

  public:
    walker(int, const LiDIA::elliptic_curve<LiDIA::gf_element> &,
           const std::array<std::vector<torsion_basis>, 2> &, const maximal_order &);
    LiDIA::gf_element public_key() { return random.back().E_i.j_invariant(); }
    // ^^^ requires that we're "between walks"...
    std::string sign(std::string message);                   // should this return bytes?
    bool verify(std::string message, std::string signature); // should this take bytes?
};
