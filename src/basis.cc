#include "basis.h"

torsion_basis::torsion_basis(const int l, const field_extension &ext,
                             const elliptic_curve<gf_element> &E_ext,
                             const array<point<gf_element>, 2> &P)
    : l(l), ext(ext), E_ext(E_ext), P(P) {}

istream &operator>>(istream &in, torsion_basis &P) {
    in >> P.l;
    in >> P.ext;
    gf_element one(P.ext.get_L()), zero(P.ext.get_L());
    one.assign_one();
    elliptic_curve<gf_element> E_ext(one, zero); // necessary so that a4, a6 belong to L.
    in >> E_ext;                                 // redundant unless the base curve changes...?

    for (int j = 0; j < 2; j++) {
        P.P[j].assign_zero(E_ext);
        in >> P.P[j];
    }

    return in;
}

ostream &operator<<(ostream &out, const torsion_basis &P) {
    out << P.l << endl;

    out << P.ext << endl;
    out << P.E_ext;

    for (int j = 0; j < 2; j++)
        out << endl << P.P[j];

    return out;
}
