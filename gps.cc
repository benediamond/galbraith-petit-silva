#include "gps.h"

gps::gps(int lambda, float epsilon) : lambda(lambda) { // purpose of epsilon?
    bigint p;
    power(p, bigint(2), 4 * lambda);
    p.assign(67966);
    do {
        p = next_prime(p);
    } while (p % 4 != 3);
    K.assign(galois_field(p, 2)); // better way to do this...?
    iota.assign_one(K);
    iota = sqrt(-iota);

    bigfloat temp_min = power(bigfloat(p), 1 + epsilon);
    array<factorization<bigint>, 2> S; // not storing this up top?
    floor(Z, 2 * (1 + epsilon) * log(bigfloat(p)));
    while (true) {
        Z++; // = next_prime(params.Z); // inefficient?
        array<bigint, 2> S_temp{1, 1};
        array<bigfloat, 2> estimates{1, 1};
        bigint l(2);
        int parity = 0;
        while (l <= Z) {
            l = next_prime(l);
            bigint e;
            floor(e, log(bigfloat(Z)) / log(bigfloat(l)));
            bigint le;
            power(le, l, e);
            S_temp[parity % 2] *= le;
            estimates[parity % 2] *= power((l + 1) / (2 * sqrt(bigfloat(l))), e);
            parity++;
        }
        if ((estimates[0] > temp_min) && (estimates[1] > temp_min)) {
            S[0] = TrialDiv(S_temp[0]);
            S[1] = TrialDiv(S_temp[1]);
            break;
        }
    }
    // replace this with an operator<< overload
    cout << "debug: p = " << p << ", Z = " << Z << ", S = " << S[0] << ", " << S[1] << endl;

    gf_element Ea(K), Eb(K);
    Ea.assign_one();
    set_curve(Ea, Eb);
    // try to keep agnostic on this curve ahead!
    // in the future, this block can be replaced with something more sophisticated

    bigint l, le; // not my style, but...
    int l_int, le_int, Z_int;
    Z.intify(Z_int);
    lidia_size_t **to_use;
    polynomial<gf_element> ***psi;
    to_use = new lidia_size_t *[Z_int + 1];
    psi = new polynomial<gf_element> **[Z_int + 1];
    for (int i = 0; i <= Z_int; i++) {
        psi[i] = new polynomial<gf_element> *[3];
        to_use[i] = new lidia_size_t[3];
        for (int j = 0; j < 3; j++)
            to_use[i][j] = 0;
    }

    for (int k = 0; k < 2; k++) {
        for (lidia_size_t i = 0; i < S[k].no_of_prime_components(); i++) {
            power(le, S[k].prime_base(i).base(), S[k].prime_exponent(i));
            le.intify(le_int);
            to_use[le_int][0] = 1;
        }
    }
    compute_psi(psi, to_use, Z_int);

    random_generator rg;
    for (int k = 0; k < 2; k++) {                                          // parallelize
        for (lidia_size_t i = 0; i < S[k].no_of_prime_components(); i++) { // parallelize
            l = S[k].prime_base(i).base();
            l.intify(l_int); // only necessary for very bottom...
            power(le, S[k].prime_base(i).base(), S[k].prime_exponent(i));
            le.intify(le_int);

            cout << "l^e = " << le << ", about to factor " << *psi[le_int][0] << endl;
            divide(*psi[le_int][0], *psi[le_int][0], psi[le_int][0]->lead_coeff());
            // ^^^ temp hack for the BUG in factorization.cc, line 646
            factorization<gf_polynomial> div_fact;
            // if (S[k].prime_exponent(i) > 1)
            factor(div_fact, *psi[le_int][0]);
            // else
            // div_fact = edf(*psi[le_int][0], (le_int - 1)/2);

            // why on earth does edf(*psi[le_int][0], (le_int - 1) / 2) hang for l^e = 29?

            galois_field L(K); // this is disgusting and unnecessary. just to initialize ext
            field_extension ext(K, L);
            elliptic_curve<gf_element> E_ext;
            array<point<gf_element>, 2> P_i;
            do {
                cout << "finding basis for le = " << le << endl;
                array<polynomial<gf_element>, 2> divisors;
                for (int j = 0; j < 2; j++) {
                    // generate a random factor of the division polynomial
                    lidia_size_t cur_idx = 0;
                    int random_int;
                    rg >> random_int;
                    random_int %= psi[le_int][0]->degree();
                    while (true) {
                        random_int -= div_fact.prime_base(cur_idx).base().degree();
                        // * div_fact.prime_exponent(cur_idx); should be no multiple roots!
                        // divpol has distinct roots in algebraic closure.
                        if (random_int < 0)
                            break;
                        cur_idx++;
                    }
                    divisors[j] = div_fact.prime_base(cur_idx).base();
                }
                lidia_size_t lcm_degree = lcm(divisors[0].degree(), divisors[1].degree());
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
                        P_temp[j][0] = find_root(lifted_divisor);
                        P_temp[j][1] = sqrt(curve_eqn(P_temp[j][0]));
                    }
                }
                E_ext = ext.embed(E_);
                for (int j = 0; j < 2; j++)
                    P_i[j] = point<gf_element>(P_temp[j][0], P_temp[j][1], E_ext);
                rational_factorization le_fact(S[k].prime_base(i).base(), S[k].prime_exponent(i));
                if (order_point(P_i[0], le_fact) != le || order_point(P_i[1], le_fact) != le)
                    continue;
                bool independent = true;
                for (int d = 0; d < S[k].prime_exponent(i); d++) {
                    bigint ld;
                    power(ld, S[k].prime_base(i).base(), d);
                    if (bg_algorithm(ld * P_i[0], ld * P_i[1], 0, le / ld - 1) != -1) {
                        independent = false;
                        break;
                    }
                }
                if (independent)
                    break;
            } while (true);
            P[k].emplace_back(l_int, le_int, ext, E_ext, P_i);
        }
    }

    for (int i = 0; i <= Z_int; i++)
        for (int j = 0; j < 3; j++) {
            if (to_use[i][j] != 0)
                delete psi[i][j];
        }

    for (int i = 0; i <= Z_int; i++) {
        delete[] psi[i];
        delete[] to_use[i];
    }
    delete[] psi;
    delete[] to_use;
}

void gps::gen_keypair() { // kill the return value? revisit this.
    keychain.emplace_back(lambda, E_, iota, P);
}

istream &operator>>(istream &in, torsion_basis &P) {
    in >> P.l >> P.le;
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
    out << P.l << " " << P.le << endl;

    out << P.ext << endl;
    out << P.E_ext;

    for (int j = 0; j < 2; j++)
        out << endl << P.P[j];

    return out;
}

istream &operator>>(istream &in, gps &gps) {
    in >> gps.lambda;

    in >> gps.K;
    in >> gps.iota; // make sure that verbose output is used!
    in >> gps.Z;

    elliptic_curve<gf_element> E_;
    in >> E_;
    gps.set_curve(E_.get_a4(), E_.get_a6());

    point<gf_element> infinity(gps.E_);
    torsion_basis P_(1, 1, field_extension(gps.K, gps.K), gps.E_, {infinity, infinity});
    for (int k = 0; k < 2; k++) {
        int r;
        in >> r;
        gps.P[k].resize(r, P_);
        for (int i = 0; i < r; i++)
            in >> gps.P[k][i];
    }
    return in;
}

ostream &operator<<(ostream &out, const gps &gps) {
    gf_element::set_output_format(1);

    out << gps.lambda << endl;

    out << gps.K << endl;
    out << gps.iota << endl; // make sure that verbose output is used!
    out << gps.Z << endl;

    out << gps.E_ << endl;

    for (int k = 0; k < 2; k++) {
        out << gps.P[k].size() << endl;
        for (int i = 0; i < gps.P[k].size(); i++)
            out << gps.P[k][i] << endl;
    }

    return out;
}
