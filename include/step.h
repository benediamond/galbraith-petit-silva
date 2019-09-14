#include <LiDIA/eco_prime.h>
#include <LiDIA/gf_element.h>
#include <LiDIA/gf_polynomial.h>
#include <LiDIA/random_generator.h>
#include <LiDIA/finite_fields/sf_gf_polynomial.h>

#include <array>

#pragma once

#include "basis.h"
#include "extension.h"

class step {
    // int le; // degree of isogeny E_{i-1} -> E_i = ord(Q_i)
    // experimenting with storing extremely few shit in here...
    // point<gf_element> Q_i; // necessary to store...? the eternal question
    torsion_basis P; // the basis for the q-torsion.

    // nota bene... experimentally NOT storing this on-hand. revisit!
    // the below are defined over F_p^2!
    LiDIA::polynomial<LiDIA::gf_element> psi; // division polynomial
    LiDIA::polynomial<LiDIA::gf_element> phi;
    LiDIA::polynomial<LiDIA::gf_element> omega; // omega(x, y) = omega(x) * y

  public:
    LiDIA::elliptic_curve<LiDIA::gf_element> E_i;
    LiDIA::bigint_matrix I_i;

    step() {}
    step(const LiDIA::elliptic_curve<LiDIA::gf_element> &, const LiDIA::bigint_matrix &);
    step(const torsion_basis &, const LiDIA::point<LiDIA::gf_element> &, const step &);
    void isogeny(LiDIA::point<LiDIA::gf_element> &, const field_extension &) const;
    // it is fed an element of E_{i-1}. inplace

    // gf_element j_invariant() { return E_i.j_invariant(); } // rename this
    const torsion_basis &get_P() const { return P; }

    friend std::istream &operator>>(std::istream &, step &);
    friend std::ostream &operator<<(std::ostream &, const step &);
};

std::istream &operator>>(std::istream &, step &);
std::ostream &operator<<(std::ostream &, const step &);

class path : public std::vector<step> {
    step dummy;

  public:
    path(const LiDIA::elliptic_curve<LiDIA::gf_element> &E_)
        : dummy(E_, LiDIA::bigint_matrix(4, 4)) {
        for (int i = 0; i < 4; i++)
            dummy.I_i.sto(i, i, 1);
    }

    step &back() { // a _non-const_ reference.
        if (empty())
            return dummy;
        else
            return std::vector<step>::back();
    }
};
