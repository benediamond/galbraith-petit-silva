#include <LiDIA/eco_prime.h>
#include <LiDIA/gf_element.h>
#include <LiDIA/gf_polynomial.h>
#include <LiDIA/random_generator.h>
#include <LiDIA/finite_fields/sf_gf_polynomial.h>
#include <LiDIA/factorization.h>
#include <LiDIA/bigmod_matrix.h>

#include <vector>

#pragma once

using namespace std;
using namespace LiDIA;

class field_extension {
    bigint p;
    galois_field K;
    galois_field L;
    bigmod_matrix embedding;
    bigmod_matrix restriction;

  public: // warning: naive implementations!
          // can do better: see e.g. De Feo Doliskani Schost
          // https://dl.acm.org/citation.cfm?id=2465956
    field_extension(const galois_field &K, const galois_field &L);
    gf_element embed(const gf_element &a) const;
    polynomial<gf_element> embed(const polynomial<gf_element> &a) const;
    elliptic_curve<gf_element> embed(const elliptic_curve<gf_element> E) const;
    gf_element restrict(const gf_element &b) const;
    polynomial<gf_element> restrict(const polynomial<gf_element> &b) const;

    bigint get_p() const;
    galois_field get_K() const;
    galois_field get_L() const;

    friend istream &operator>>(istream &in, field_extension &ext);
    friend ostream &operator<<(ostream &in, const field_extension &ext);
};
