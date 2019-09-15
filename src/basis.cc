#include "basis.h"

using namespace std;
using namespace LiDIA;

torsion_basis::torsion_basis(const galois_field &K, const int l,
                             const elliptic_curve<gf_element> &E_,
                             const polynomial<gf_element> &curve_eqn,
                             // could just re-determine curve_eqn from E_, but...
                             const array<polynomial<gf_element>, 2> &div_factors)
    : l(l) {
    galois_field L_field((bigint(l)));
    gf_element p_l(L_field);
    p_l.assign(-K.characteristic());
    int degree;
    p_l.order().intify(degree); // note: relies on l being a prime, not a prime power...
    cout << "constructing torsion basis for l = " << l << ", degree = " << degree << endl;
    galois_field L(K.characteristic(), 2 * degree);
    ext = field_extension(K, L); // class member
    polynomial<gf_element> curve_eqn_ext = ext.embed(curve_eqn);
    E_ext = ext.embed(E_);
    factorization<gf_polynomial> div_fact;

    random_generator rg;

    array<polynomial<gf_element>, 2> divisors;
    array<array<gf_element, 2>, 2> P_temp;
    do {
        for (int j = 0; j < 2; j++) {
            factor(div_fact, div_factors[j]); // could use edf
            int r;
            rg >> r;
            divisors[j] = div_fact.prime_base(r % div_fact.no_of_prime_components()).base();

            polynomial<gf_element> lifted_divisor = ext.embed(divisors[j]);
            P_temp[j][0].assign(find_root(lifted_divisor));
            P_temp[j][1].assign(sqrt(curve_eqn_ext(P_temp[j][0])));
            P[j] = point<gf_element>(P_temp[j][0], P_temp[j][1], E_ext);
        }
    } while (bg_algorithm(P[0], P[1], 0, l - 1) != -1);
    // this is even _a question_ only in the l == q case, otherwise will be handled externally
    // in the elkies case, the caller must supply the order-l factors, and can check lin. indce.
    // barring elkies, the caller has no way to assess linear independence of factors, so we must
    // in order to use the same semantics in both cases, we force there to be two poly args
    // this causes extra work in the non-elkies case (ie as currently written, _only_ when l == q)
    // if we _didn't_ care about elkies, we could just pass in the entire divpol, and factor here
    // (this would make sense in the case of prime powers, with vanilla GPS, or maybe a hybrid)
    // until then, this is negligible sacrifice in efficiency for a unified interface / semantics
}

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
