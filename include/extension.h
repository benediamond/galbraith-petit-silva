#include <LiDIA/eco_prime.h>
#include <LiDIA/gf_element.h>
#include <LiDIA/gf_polynomial.h>
#include <LiDIA/random_generator.h>
#include <LiDIA/finite_fields/sf_gf_polynomial.h>
#include <LiDIA/factorization.h>
#include <LiDIA/bigmod_matrix.h>

#include <vector>

#pragma once

class field_extension {
    LiDIA::bigint p;
    LiDIA::galois_field K;
    LiDIA::galois_field L;
    LiDIA::bigmod_matrix embedding;
    LiDIA::bigmod_matrix restriction;

  public: // warning: naive implementations!
    // can do better: see e.g. De Feo Doliskani Schost
    // https://dl.acm.org/citation.cfm?id=2465956
    field_extension() {}
    field_extension(const LiDIA::galois_field &K, const LiDIA::galois_field &L);
    LiDIA::gf_element embed(const LiDIA::gf_element &a) const;
    LiDIA::polynomial<LiDIA::gf_element> embed(const LiDIA::polynomial<LiDIA::gf_element> &a) const;
    LiDIA::elliptic_curve<LiDIA::gf_element>
    embed(const LiDIA::elliptic_curve<LiDIA::gf_element> E) const;
    LiDIA::gf_element restrict(const LiDIA::gf_element &b) const;
    LiDIA::polynomial<LiDIA::gf_element>
    restrict(const LiDIA::polynomial<LiDIA::gf_element> &b) const;

    LiDIA::bigint get_p() const;
    LiDIA::galois_field get_K() const;
    LiDIA::galois_field get_L() const;

    friend std::istream &operator>>(std::istream &in, field_extension &ext);
    friend std::ostream &operator<<(std::ostream &in, const field_extension &ext);
};
