#include "gps.h"

gps::gps(int lambda, float epsilon) : lambda(lambda) { // purpose of epsilon?
    bigint p;
    power(p, bigint(2), 4 * lambda);
    do {
        p = next_prime(p);
    } while (p % 24 != 1);

    // i want p to be == 1 modulo 8 (to get the last row of the table)
    // _and_ == 1 (modulo 3), so that (-3 / p) != -1. that way, following
    // Bröker's algorithm, the q chosen won't be 3, and we will use the final case.
    // this avoids a curve with j-inv = 0.
    K.assign(galois_field(p, 2));

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
        cout << "curve set: " << E_ << endl;

        factorization<gf_polynomial> div_fact;
        polynomial<gf_element> div_pol(K);
        int q_int;
        q.intify(q_int);
        compute_psi(div_pol, q_int);
        divide(div_pol, div_pol, div_pol.lead_coeff());
        factor(div_fact, div_pol); // ONLY time we'll have to factor a divpol. also, could use edf?
        cout << "divpol factored: " << div_fact << endl;

        torsion_basis Q;
        random_generator rg;
        do {
            array<polynomial<gf_element>, 2> divisors;
            for (int j = 0; j < 2; j++) {
                int r;
                rg >> r;
                divisors[j] = div_fact.prime_base(r % div_fact.no_of_prime_components()).base();
            }
            galois_field L(K);
            field_extension ext(K, L);
            elliptic_curve<gf_element> E_ext;
            array<point<gf_element>, 2> P; // shadows class member

            lidia_size_t lcm_degree = lcm(divisors[0].degree(), divisors[1].degree());
            // all factors should be of equal degree.
            if (lcm_degree > 1) {
                L = galois_field(K.characteristic(), 2 * lcm_degree);
                ext = field_extension(K, L);
            }
            // _either way_ i guess i should embed the polys into L, even if it's trivial
            array<array<gf_element, 2>, 2> P_temp;
            for (int j = 0; j < 2; j++) {
                polynomial<gf_element> lifted_divisor = ext.embed(divisors[j]);
                P_temp[j][0].assign(find_root(lifted_divisor));
            }
            polynomial<gf_element> curve_eqn = ext.embed(CurveEqn());
            if (curve_eqn(P_temp[0][0]).is_square() && curve_eqn(P_temp[1][0]).is_square()) {
                for (int j = 0; j < 2; j++)
                    P_temp[j][1].assign(sqrt(curve_eqn(P_temp[j][0])));
            } else {
                L = galois_field(K.characteristic(), 4 * lcm_degree);
                ext = field_extension(K, L);
                curve_eqn = ext.embed(CurveEqn());
                for (int j = 0; j < 2; j++) {
                    polynomial<gf_element> lifted_divisor = ext.embed(divisors[j]);
                    P_temp[j][0].assign(find_root(lifted_divisor));
                    P_temp[j][1].assign(sqrt(curve_eqn(P_temp[j][0])));
                }
            }
            E_ext = ext.embed(E_);
            for (int j = 0; j < 2; j++)
                P[j] = point<gf_element>(P_temp[j][0], P_temp[j][1], E_ext);
            Q = torsion_basis(q_int, ext, E_ext, P);
        } while (bg_algorithm(Q.P[0], Q.P[1], 0, q - 1) != -1);
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

    random_generator rg;
    for (int k = 0; k < 2; k++) {                                          // parallelize
        for (lidia_size_t i = 0; i < S[k].no_of_prime_components(); i++) { // parallelize
            l = S[k].prime_base(i).base();
            int l_int;
            l.intify(l_int); // only necessary for very bottom...
            set_prime(l_int);
            build_poly_in_X(meq_pol, jinv);

            galois_field L(K); // this is disgusting and unnecessary. just to initialize ext
            field_extension ext(K, L);
            elliptic_curve<gf_element> E_ext;
            array<point<gf_element>, 2> P_i;
            do {
                cout << "finding basis for l = " << l << endl;

                array<polynomial<gf_element>, 2> divisors;
                for (int j = 0; j < 2; j++) {
                    ftau[0] = find_root(meq_pol);
                    compute_divisor_of_division_polynomial(divisors[j], Ea, Eb);
                }
                // args Ea and Eb are never used

                lidia_size_t lcm_degree = lcm(divisors[0].degree(), divisors[1].degree());
                // ^^^ max should suffice here, for complicated reasons... (even for e > 1)
                // actually, in the case e == 1, all factors should be of equal degree.
                if (lcm_degree > 1) {
                    L = galois_field(K.characteristic(), 2 * lcm_degree);
                    ext = field_extension(K, L);
                }
                // _either way_ i guess i should embed the polys into L, even if it's trivial
                array<array<gf_element, 2>, 2> P_temp;
                for (int j = 0; j < 2; j++) {
                    polynomial<gf_element> lifted_divisor = ext.embed(divisors[j]);
                    P_temp[j][0] = find_root(lifted_divisor);
                }
                polynomial<gf_element> curve_eqn = ext.embed(CurveEqn());
                if (curve_eqn(P_temp[0][0]).is_square() && curve_eqn(P_temp[1][0]).is_square()) {
                    for (int j = 0; j < 2; j++)
                        P_temp[j][1] = sqrt(curve_eqn(P_temp[j][0]));
                } else {
                    L = galois_field(K.characteristic(), 4 * lcm_degree);
                    ext = field_extension(K, L);
                    curve_eqn = ext.embed(CurveEqn());
                    for (int j = 0; j < 2; j++) {
                        polynomial<gf_element> lifted_divisor = ext.embed(divisors[j]);
                        P_temp[j][0] = find_root(lifted_divisor); // this is extremely slow,
                        // and worse yet, would be completely avoidable if we had a good way to
                        // work with towers of field extensions. we already have the irreducible
                        // polynomial which this thing is going to be a root of.
                        P_temp[j][1] = sqrt(curve_eqn(P_temp[j][0]));
                    }
                }
                E_ext = ext.embed(E_);
                for (int j = 0; j < 2; j++)
                    P_i[j] = point<gf_element>(P_temp[j][0], P_temp[j][1], E_ext);
            } while (bg_algorithm(P_i[0], P_i[1], 0, l - 1) != -1);
            P[k].emplace_back(l_int, ext, E_ext, P_i);
        } // todo: come up with better square root algorithms.
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

    torsion_basis P_(0, field_extension(gps.K, gps.K), E_, {E_, E_});
    // torsion_basis P_(1, field_extension(gps.K, gps.K), gps.E_, {infinity, infinity});
    for (int k = 0; k < 2; k++) {
        int r;
        in >> r;
        // ^^^ need to initialize this with a dummy, to avoid "e == NULL" runtime errors?
        gps.P[k].resize(r, P_);
        for (int i = 0; i < r; i++) {
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
