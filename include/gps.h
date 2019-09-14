#include <LiDIA/eco_prime.h>
#include <LiDIA/gf_element.h>
#include <LiDIA/gf_polynomial.h>
#include <LiDIA/random_generator.h>
#include <LiDIA/finite_fields/sf_gf_polynomial.h>
#include <LiDIA/factorization.h>
#include <LiDIA/gec_complex_multiplication.h>

#include <vector>
#include <array>
#include <fstream>

#pragma once

#include "walker.h"
#include "extension.h"
#include "order.h"
#include "basis.h"

class gps : public LiDIA::eco_prime {
  private:
    int lambda;
    // not storing n!
    LiDIA::galois_field K;

    maximal_order O_;

    LiDIA::bigint Z; // alternative name for B, taken by coeff
    std::array<std::vector<torsion_basis>, 2> P;

  public:
    gps() {}
    gps(int, float);

    void gen_keypair();

    // points on the base curve.

    std::vector<walker> keychain; // making this public... will have "sign" and "verify" methods...

    friend std::istream &operator>>(std::istream &, gps &);
    friend std::ostream &operator<<(std::ostream &, const gps &);
};

std::istream &operator>>(std::istream &, gps &);
std::ostream &operator<<(std::ostream &, const gps &);
