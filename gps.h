#include <LiDIA/eco_prime.h>
#include <LiDIA/gf_element.h>
#include <LiDIA/gf_polynomial.h>
#include <LiDIA/random_generator.h>
#include <LiDIA/finite_fields/sf_gf_polynomial.h>
#include <LiDIA/factorization.h>

#include <vector>
#include <array>
#include <fstream>

#pragma once

#include "walker.h"
#include "extension.h"

using namespace std;
using namespace LiDIA;

class torsion_basis {
  public:
    int l;
    int le;

    field_extension ext;              // contains K and L
    elliptic_curve<gf_element> E_ext; // weirdly necessary when reading in points... revisit
    array<point<gf_element>, 2> P;    // use the same name again...
    // should any/all of these be references?

    torsion_basis(const int l, const int le, const field_extension &ext,
                  const elliptic_curve<gf_element> E_ext, const array<point<gf_element>, 2> &P)
        : l(l), le(le), ext(ext), E_ext(E_ext), P(P) {}
};

istream &operator>>(istream &, torsion_basis &);
ostream &operator<<(ostream &, const torsion_basis &);

class walker;

class gps : public eco_prime {
  private:
    int lambda;
    // not storing n!

    galois_field K;
    gf_element iota;

    bigint Z; // alternative name for B, taken by coeff
    array<vector<torsion_basis>, 2> P;

  public:
    gps() {}
    gps(int, float);

    void gen_keypair();

    vector<walker> keychain; // making this public... will have "sign" and "verify" methods...

    friend istream &operator>>(istream &, gps &);
    friend ostream &operator<<(ostream &, const gps &);
};

istream &operator>>(istream &, gps &);
ostream &operator<<(ostream &, const gps &);
