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

using namespace std;
using namespace LiDIA;

class gps : public eco_prime {
  private:
    int lambda;
    // not storing n!
    galois_field K;

    maximal_order O_;

    bigint Z; // alternative name for B, taken by coeff
    array<vector<torsion_basis>, 2> P;

  public:
    gps() {}
    gps(int, float);

    void gen_keypair();

    // points on the base curve.

    vector<walker> keychain; // making this public... will have "sign" and "verify" methods...

    friend istream &operator>>(istream &, gps &);
    friend ostream &operator<<(ostream &, const gps &);
};

istream &operator>>(istream &, gps &);
ostream &operator<<(ostream &, const gps &);
