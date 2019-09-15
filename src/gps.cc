#include "gps.h"

using namespace std;
using namespace LiDIA;

gps::gps(int lambda, float epsilon) : lambda(lambda) { // purpose of epsilon?
    bigint p;
    power(p, bigint(2), 4 * lambda);
    do {
        p = next_prime(p);
    } while (p % 24 != 1);
    random_generator rg;

    // i want p to be == 1 modulo 8 (to get the last row of the table)
    // _and_ == 1 (modulo 3), so that (-3 / p) != -1. that way, following
    // Bröker's algorithm, the q chosen won't be 3, and we will use the final case.
    // this avoids a curve with j-inv = 0.
    K.assign(galois_field(p, 2));
    polynomial<gf_element> curve_eqn(K); // will set below
    // see Reinier Bröker, Constructing Supersingular Elliptic Curves
    bigint q(3);
    do {
        q = next_prime(q);
    } while (q % 4 != 3 || jacobi(-q, p) != -1);

    bigfloat temp_min = power(bigfloat(p), 1 + epsilon);
    array<factorization<bigint>, 2> S; // not storing this up top?
    array<bigint, 2> S_temp{1, 1};
    array<bigfloat, 2> estimates{1, 1};

    bigint l(2);
    int parity = 0;
    while (estimates[0] <= temp_min || estimates[1] <= temp_min) {
        l = next_prime(l);
        if (l == q)
            continue; // skip this...?
        S_temp[parity % 2] *= l;
        estimates[parity % 2] *= (l + 1) / (2 * sqrt(bigfloat(l)));
        parity++;
    }
    Z = l;
    S[0] = TrialDiv(S_temp[0]);
    S[1] = TrialDiv(S_temp[1]);
    cout << "debug: p = " << p << ", Z = " << Z << ", S = " << S[0] << ", " << S[1] << endl;

    gec_complex_multiplication gec(bigint(-q));
    gec.compute_class_polynomial(1);
    Fp_polynomial P_K(gec.get_class_polynomial(), p);
    base_vector<bigint> roots(find_roots(P_K)); // guaranteed not to be empty, Lem. 2.3
    cout << "q = " << q << ", P_K: " << P_K << ", roots: " << roots << endl; // debug / temp
    gf_element Ea(K), Eb(K), j(K);           // these will be used in various places
    for (int i = 0; i < roots.size(); i++) { // is there any kind of "foreach" I can use?
        // in practice it seems that there's only 1 root, over F_p, or for that matter at all...
        // where is the guarantee that some j with E(j) == O_0 exists, and further that it is in F_p?
        j.assign(roots[i]);

        gf_element temp;
        temp.assign_one(K);
        temp = 1728 * temp - j;
        A = 3 * j / temp;
        B = 2 * j / temp;
        set_curve(A, B);
        curve_eqn = CurveEqn();
        cout << "curve set: " << E_ << endl;

        polynomial<gf_element> div_pol(K);
        int q_int;
        q.intify(q_int);
        compute_psi(div_pol, q_int);
        divide(div_pol, div_pol, div_pol.lead_coeff()); // necessary?
        torsion_basis Q(K, q_int, E_, curve_eqn, {div_pol, div_pol});

        // "cute" workaround to integrate this and the below. later, we'll have 2 elkies factors
        cout << "q-torsion basis chosen: " << Q.P[0] << ", " << Q.P[1] << endl;

        while (true) {
            // would have been elegant, and actually faster, to not bother with the torsion basis,
            // and just iterate through the factors of the division polynomial above.
            // and yet we will need a full torsion basis soon anyway, to determine if isomorphism
            // and so i guess it's better to just generate the whole thing while we're here,
            // at which point it no longer becomes clean to iterate through distinct subgroups.
            array<int, 2> coeffs;
            do {
                rg >> coeffs[0];
                rg >> coeffs[1];
                coeffs[0] %= q_int; // warning:
                coeffs[1] %= q_int; // modulo bias
            } while (coeffs[0] == 0 && coeffs[1] == 0);

            point<gf_element> Q_i(coeffs[0] * Q.P[0] + coeffs[1] * Q.P[1]);
            step dummy(E_, bigint_matrix(4, 4));
            step theta(Q, Q_i, dummy);

            if (theta.E_i.j_invariant() == j) { // using coefficients alone is not legit?
                gf_element u(K);
                do {
                    u = sqrt(sqrt(A / theta.E_i.get_a4()));
                    power(temp, u, 6);
                } while (temp != B / theta.E_i.get_a6());
                // complicated stuff. see Silverman, III, Prop. 1.4., b.
                // viewing F_q as Z / (q - 1), we're looking for an element in the intersection of
                // the fiber of (A / A') under *4 and the fiber of (B / B') under *6.
                // how do we know that such an element exists? By A^3 B'^2 = A'^3 B^2,
                // these respective fibers in turn live under the same fiber under *12.
                // within this latter fiber, we're looking for elements which in addition
                // have perscribed "residue classes" modulo 4 and modulo 6, respectively.
                // there will be exactly 2 such elements, or half of the fiber (*4)^{-1}(A / A').

                O_ = maximal_order(theta, u); // will fill in c later...!
                break;
            }
        }

        base_vector<gf_element> T_x(find_roots(CurveEqn()));
        // could potentially use extract_sqrt...?
        array<point<gf_element>, 2> T;
        bool success = true;

        field_extension trivial(K, K);
        for (int j = 0; j < 2; j++) { // the array will (?) have length 3, but we only need 2!
            // the third element will be a linear combination of the other two.
            T[j] = point<gf_element>(T_x[j], gf_element(K), E_);
            if (!(T[j] + O_.theta(T[j], trivial)).is_zero() ||
                !(O_.pi(T[j]) - O_.theta(O_.pi(T[j]), trivial)).is_zero()) {
                success = false;
                break;
            }
            galois_field mod_q(q);
            temp.assign_one(mod_q);
            temp.assign(sqrt(-temp / p));
            bigint c(temp.lift_to_Z());
            if (O_.theta(Q.P[j], Q.ext) == c * O_.theta(O_.pi(Q.P[j]), Q.ext))
                O_.set_c(c);
            else if (O_.theta(Q.P[j], Q.ext) == -c * O_.theta(O_.pi(Q.P[j]), Q.ext))
                O_.set_c(-c);
            else {
                success = false;
                break;
            }
        }

        if (success) {
            cout << "Bröker's algorithm success." << endl;
            break;
        } else
            cout << "WARNING: Error during Bröker's alorithm, trying next curve (if one exists?)."
                 << endl;
        // warning: if something is wrong, this could silently run out the whole loop.
    }

    for (int k = 0; k < 2; k++) {                                          // parallelize
        for (lidia_size_t i = 0; i < S[k].no_of_prime_components(); i++) { // parallelize
            l = S[k].prime_base(i).base();
            int l_int;
            l.intify(l_int); // only necessary for very bottom...
            set_prime(l_int);
            build_poly_in_X(meq_pol, jinv);

            array<polynomial<gf_element>, 2> div_factors;
            do {
                for (int j = 0; j < 2; j++) {
                    ftau[0] = find_root(meq_pol);
                    compute_divisor_of_division_polynomial(div_factors[j], Ea, Eb);
                }
            } while (div_factors[0] == div_factors[1]);
            P[k].emplace_back(K, l_int, E_, curve_eqn, div_factors);
        }
    }
}

void gps::gen_keypair() { // kill the return value? revisit this.
    keychain.emplace_back(lambda, E_, P, O_);
}

istream &operator>>(istream &in, gps &gps) {
    in >> gps.lambda;
    in >> gps.K;
    in >> gps.O_;

    in >> gps.Z;

    elliptic_curve<gf_element> E_;
    in >> E_;
    gps.set_curve(E_.get_a4(), E_.get_a6());

    for (int k = 0; k < 2; k++) {
        int r;
        in >> r;
        gps.P[k].resize(r);
        for (int i = 0; i < r; i++) {
            gps.P[k][i].E_ext = E_; // necessary to avoid e = NULL runtime errors
            in >> gps.P[k][i];
        }
    }
    return in;
}

ostream &operator<<(ostream &out, const gps &gps) {
    gf_element::set_output_format(1);

    out << gps.lambda << endl;
    out << gps.K << endl;
    out << gps.O_ << endl;

    out << gps.Z << endl;

    out << gps.E_ << endl;

    for (int k = 0; k < 2; k++) {
        out << gps.P[k].size() << endl;
        for (int i = 0; i < gps.P[k].size(); i++)
            out << gps.P[k][i] << endl;
    }

    return out;
}
