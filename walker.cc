#include "walker.h"

using CryptoPP::byte;
using CryptoPP::SHA512;

base_vector<bigint> conjugate(base_vector<bigint> e) { // length 4
    array<bigint, 4> e_bar_arr{e[0] + e[3], -e[1], -e[2], -e[3]};
    base_vector<bigint> e_bar(4);
    for (int i = 0; i < 4; i++)
        e_bar[i].assign(e_bar_arr[i]); // conjugate b_1
    return e_bar;
}

bigint_matrix times(base_vector<bigint> q, bigint p_coeff) {
    array<array<bigint, 4>, 4> O_0_q_arr{
        {{q[0], q[1], q[2], q[3]},
         {-q[1] - q[2], q[0] + q[3], -q[3], q[2]},
         {-p_coeff * q[2], p_coeff * q[3], q[0], -q[1]},
         {-p_coeff * q[3], -p_coeff * q[2], q[1] + q[2], q[0] + q[3]}}};
    bigint_matrix O_0_q(4, 4);
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            O_0_q.sto(i, j, O_0_q_arr[j][i]);
    return O_0_q;
}

step::step(const elliptic_curve<gf_element> &E_, const bigint_matrix &I_)
    : E_i(E_), Q_i(E_), I_i(I_) {}

step::step(const torsion_basis &P, const point<gf_element> &Q_i,
           const array<gf_element, 2> &Q_coords, const step &last)
    : Q_i(Q_i), psi(P.ext.get_L()), phi(P.ext.get_K()), omega(P.ext.get_K()), I_i(last.I_i) {
    // see Galbraith, Thm. 25.1.6., Kohel, Ch. 2.4
    elliptic_curve<gf_element> E_i_minus_1(P.ext.embed(last.E_i));
    point<gf_element> Q(Q_coords[0], Q_coords[1], E_i_minus_1);

    gf_element tG(P.ext.get_L()), wG(P.ext.get_L());
    point<gf_element> Q_accum(Q);
    psi.assign_one();
    polynomial<gf_element> x(P.ext.get_L());
    x.assign_x();
    for (int j = 0; j < (P.le - 1) / 2; j++) { // iterate over G_1
        gf_element F_x(3 * Q_accum.get_x() * Q_accum.get_x() + E_i_minus_1.get_a4());
        gf_element F_y(-2 * Q_accum.get_y());
        gf_element uQ(F_y * F_y);
        gf_element tQ(2 * F_x);
        add(tG, tG, tQ);
        add(wG, wG, uQ + Q_accum.get_x() * tQ);

        multiply(psi, psi, x - Q_accum.get_x());

        add(Q_accum, Q_accum, Q);
    }
    gf_element A_4(P.ext.restrict(E_i_minus_1.get_a4() - 5 * tG));
    gf_element A_6(P.ext.restrict(E_i_minus_1.get_a6() - 7 * wG));
    E_i = elliptic_curve<gf_element>(A_4, A_6);

    psi = P.ext.restrict(psi);

    phi.set_coefficient(4, 3);
    phi.set_coefficient(2 * P.ext.restrict(E_i_minus_1.get_b4()), 1);
    phi.set_coefficient(P.ext.restrict(E_i_minus_1.get_b6()), 0);
    polynomial<gf_element> psi_prime = derivative(psi);

    phi *= psi_prime * psi_prime - derivative(psi_prime) * psi;
    polynomial<gf_element> temp_1(P.ext.get_K()); // stupid, and probably unnecessary
    temp_1.set_coefficient(6, 2);
    temp_1.set_coefficient(P.ext.restrict(E_i_minus_1.get_b4()), 0);
    phi -= temp_1 * psi_prime * psi;
    temp_1.assign_zero(P.ext.get_K()); // arg unnecessary
    temp_1.set_coefficient(P.le, 1);
    temp_1.set_coefficient(2 * psi[psi.degree() - 1], 0); // this won't be used for identity
    phi += temp_1 * psi * psi;

    multiply(omega, derivative(phi), psi);
    subtract(omega, omega, bigint(2) * phi * psi_prime);

    cout << "l^e = " << P.le << ", stepped: " << E_i.j_invariant() << endl;
}

void step::generate_ideal(const gf_element &iota, const torsion_basis &P) {
    bigint_matrix I_i_minus_1(I_i);
    gf_element iota_ext(P.ext.embed(iota));

    point<gf_element> Q_by_2(((P.le + 1) / 2) * Q_i); // only reason class member Q_i necessary
    point<gf_element> i_Q_by_2(-Q_by_2.get_x(), iota_ext * Q_by_2.get_y(), Q_by_2.get_curve());
    point<gf_element> j_Q_by_2(frobenius(Q_by_2));
    point<gf_element> k_Q_by_2(-frobenius(i_Q_by_2)); // k = -ji
    array<point<gf_element>, 4> B_Q_i{Q_i, 2 * i_Q_by_2, i_Q_by_2 + j_Q_by_2, Q_by_2 + k_Q_by_2};

    base_vector<bigint> a_i(4);
    bigint p_coeff((P.ext.get_p() + 1) / 4);
    do {
        while (true) { // { do {
            // coeffs.randomize(P.le); // unfortunately, this function is completely broken
            for (int i = 0; i < 4; i++)
                a_i[i].assign(randomize(bigint(P.le)));
            // could use hensel's lemma, but when f'(r) == 0 (mod l), things get complicated
            // lifting behavior becomes highly nontrivial, and end up needing to tree-search
            // the l-adic integers up to depth e. might as well linearly search l^e
            // this in turn essentially amounts to a fully random trial-and-error search...
            a_i.assign((I_i_minus_1 * bigint_matrix(a_i))(0)); // get 0th column
            bigint norm = a_i[0] * a_i[0] + a_i[0] * a_i[3] + a_i[1] * a_i[1] + a_i[1] * a_i[2] +
                          p_coeff * (a_i[2] * a_i[2] + a_i[3] * a_i[3]);
            if (norm % P.le == 0 && norm % (P.le * P.l) != 0)
                break;
        }
    } while (!((a_i[0] % P.le) * B_Q_i[0] + (a_i[1] % P.le) * B_Q_i[1] +
               (a_i[2] % P.le) * B_Q_i[2] + (a_i[3] % P.le) * B_Q_i[3])
                  .is_zero());
    bigint_lattice I_i_temp(4, 8);
    bigint_matrix O_0_a_i = times(a_i, p_coeff);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++)
            I_i_temp.sto(i, j, P.le * I_i_minus_1[i][j]);
        for (int j = 4; j < 8; j++)
            I_i_temp.sto(i, j, O_0_a_i[i][j - 4]);
    }
    I_i_temp.lll(0.99);
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            I_i.sto(i, j, I_i_temp[i][j]);
    // omitting n once again for now. if you want it back, use det
}

// the "step" has already been constructed (using Velu) and the "coefficients" are ready.
void step::isogeny(array<gf_element, 2> &Q_coords, const field_extension &ext) {
    polynomial<gf_element> lifted_psi(ext.embed(psi));
    polynomial<gf_element> lifted_phi(ext.embed(phi));
    polynomial<gf_element> lifted_omega(ext.embed(omega));

    gf_element psi_x(lifted_psi(Q_coords[0]));
    gf_element psi_x_2, psi_x_3;
    square(psi_x_2, psi_x);
    multiply(psi_x_3, psi_x_2, psi_x);
    gf_element phi_x(lifted_phi(Q_coords[0]));
    gf_element omega_x_y(lifted_omega(Q_coords[0]));
    multiply(omega_x_y, omega_x_y, Q_coords[1]);
    divide(Q_coords[0], phi_x, psi_x_2);
    divide(Q_coords[1], omega_x_y, psi_x_3);
}

void walker::ideal_step(const gf_element &iota, const torsion_basis &P, const bigint_matrix &I) {
    array<point<gf_element>, 2> P_i_by_2{((P.le + 1) / 2) * P.P[0], ((P.le + 1) / 2) * P.P[1]};

    gf_element iota_ext(P.ext.embed(iota));
    array<point<gf_element>, 2> i_P_i_by_2{
        {{-P_i_by_2[0].get_x(), iota_ext * P_i_by_2[0].get_y(), P_i_by_2[0].get_curve()},
         {-P_i_by_2[1].get_x(), iota_ext * P_i_by_2[1].get_y(), P_i_by_2[1].get_curve()}}};
    array<point<gf_element>, 2> j_P_i_by_2{frobenius(P_i_by_2[0]), frobenius(P_i_by_2[1])};
    array<point<gf_element>, 2> k_P_i_by_2{-frobenius(i_P_i_by_2[0]), -frobenius(i_P_i_by_2[1])};
    array<array<point<gf_element>, 4>, 2> B_P_i{
        {{P.P[0], 2 * i_P_i_by_2[0], i_P_i_by_2[0] + j_P_i_by_2[0], P_i_by_2[0] + k_P_i_by_2[0]},
         {P.P[1], 2 * i_P_i_by_2[1], i_P_i_by_2[1] + j_P_i_by_2[1], P_i_by_2[1] + k_P_i_by_2[1]}}};

    rational_factorization le_fact(trialdiv(P.le));
    array<bigint, 2> coeffs;
    for (int j = 0; j < 4; j++) {
        array<point<gf_element>, 2> a_P_i{
            (I[0][j] % P.le) * B_P_i[0][0] + (I[1][j] % P.le) * B_P_i[0][1] +
                (I[2][j] % P.le) * B_P_i[0][2] + (I[3][j] % P.le) * B_P_i[0][3],
            (I[0][j] % P.le) * B_P_i[1][0] + (I[1][j] % P.le) * B_P_i[1][1] +
                (I[2][j] % P.le) * B_P_i[1][2] + (I[3][j] % P.le) * B_P_i[1][3]};
        array<bool, 2> P_i_prim{order_point(a_P_i[0], le_fact) == P.le,
                                order_point(a_P_i[1], le_fact) == P.le};
        if (!P_i_prim[0] && !P_i_prim[1])
            continue; // at least one must annihilate exactly Q...
        // but can we assert that it must be j = 0...? not clear...
        int idx = P_i_prim[0] ? 0 : 1;
        coeffs[idx] = bg_algorithm(a_P_i[idx], a_P_i[1 - idx], 0, P.le - 1, le_fact);
        coeffs[1 - idx] = -1;
        break;
    }

    point<gf_element> Q_i(coeffs[0] * P.P[0] + coeffs[1] * P.P[1]);

    array<gf_element, 2> Q_coords{Q_i.get_x(), Q_i.get_y()};
    for (int j = 0; j < ideal.size(); j++) // carry forward through past walks
        ideal[j].isogeny(Q_coords, P.ext);
    ideal.emplace_back(P, Q_i, Q_coords, ideal.back());
}

void walker::random_step(const torsion_basis &P) {
    array<int, 2> coeffs;
    random_generator rg; // store this somewhere once and for all...?
    do {
        rg >> coeffs[0];
        rg >> coeffs[1];
        coeffs[0] %= P.le; // warning:
        coeffs[1] %= P.le; // modulo bias
    } while (coeffs[0] % P.l == 0 || coeffs[1] % P.l == 0);
    // tricky proof shows this is equivalent to choosing a _subgroup_ uniformly
    point<gf_element> Q_i(coeffs[0] * P.P[0] + coeffs[1] * P.P[1]);
    array<gf_element, 2> Q_coords{Q_i.get_x(), Q_i.get_y()};
    for (int j = 0; j < random.size(); j++) // carry forward through past walks
        random[j].isogeny(Q_coords, P.ext);
    random.emplace_back(P, Q_i, Q_coords, random.back());
    random.back().generate_ideal(iota, P);
}

void walker::reroute() {
    // warning: rewrites steps, but does _not_ regenerate
    bigint_matrix I(random.back().I_i);

    // NOT YET minkowski basis, let alone a random one.
    bigint n(I.det()), m, p(iota.get_field().characteristic()), temp; // silly.
    n.abs();
    sqrt(n, n);
    ceil(m, (log(bigfloat(p))));
    bigint p_coeff((p + 1) / 4);
    array<bigint, 2> S{1, 1};
    for (int k = 0; k < 2; k++) { // different from the usual def. of S...
        S[0] *= P[k][0].le;
        for (int i = 1; i < P[k].size(); i++)
            multiply(S[1], S[1], P[k][i].le);
    }

    base_vector<bigint> del(4);
    while (true) { // if no prime found, then problem...
        for (int i = 0; i < 4; i++)
            del[i].assign(randomize(2 * m) - m); // cleaner than exhaustive, and still fast
        del = (I * bigint_matrix(del))(0);
        bigint norm(del[0] * del[0] + del[0] * del[3] + del[1] * del[1] + del[1] * del[2] +
                    p_coeff * (del[2] * del[2] + del[3] * del[3]));
        if ((norm / n) % 4 == 1 && (norm / n).is_prime()) // so that x_0, y_0 can be found below
            break;
    }

    I = (times(conjugate(del), p_coeff) * I) / n;
    // _not actually necessary_ for J = I' (see paper)!

    temp.assign(I.det());
    temp.abs();
    sqrt(n, temp); // n will temporarily alias n(I'). will be recovered later
    while (true) {
        for (int i = 0; i < 4; i++)
            del[i].assign(randomize(m));
        del = (I * bigint_matrix(del))(0); // reuse del, actually equals "alpha"...
        bigint norm(del[0] * del[0] + del[0] * del[3] + del[1] * del[1] + del[1] * del[2] +
                    p_coeff * (del[2] * del[2] + del[3] * del[3]));
        if (gcd(norm, temp) == n)
            break;
    }

    galois_field F(n, 1);
    gf_element prime(F);
    prime.assign(-1);
    prime.assign(sqrt(prime));

    sqrt(m, n * S[0] / (2 * p)); // includes floor
    base_vector<bigint> b_1(4), b_2_temp(2);
    while (true) { // will keep selecting beta_1 until line 14 is solvable _as is_.
        while (true) {
            b_1[2].assign(randomize(2 * m) - m);
            b_1[3].assign(randomize(2 * m) - m);
            temp.assign(n * S[0] - p * (b_1[2] * b_1[2] + b_1[3] * b_1[3]));
            if (!temp.is_prime())
                continue;
            if (cornacchia(b_1[0], b_1[1], bigint(-4), temp)) {
                divide(b_1[0], b_1[0], 2); // because of LiDIA cornacchia
                break;
            }
        }
        bigmod_matrix b_1_mat(2, 2, n), del_mat(2, 1, n); // embedding into in M_2(Z / NZ)
        // note: using the first row and first col. could fail if identically zero... unlikely
        gf_element x_0(F), y_0(F);
        x_0.assign(1 - p);
        divide(x_0, x_0, 2);
        y_0.assign(1 + p);
        divide(y_0, y_0, 2 * prime);
        bigint x, y;
        x.assign(x_0.lift_to_Z());
        y.assign(y_0.lift_to_Z());
        // x^2 + y^2 should == -p (mod n).

        del[2] += del[2] % 2 == 0 ? 0 : n;
        del[3] += del[3] % 2 == 0 ? 0 : n;
        del_mat.sto(0, 0, del[0] + del[3] / 2 - x * del[2] / 2 - y * del[3] / 2);
        del_mat.sto(1, 0, del[1] + del[2] / 2 - y * del[2] / 2 + x * del[3] / 2);
        // only the first row of the embedding of delta in M_2(Z / NZ)

        array<bigint, 2> b_1_arr{b_1[0] - x * b_1[2] - y * b_1[3],
                                 b_1[1] - y * b_1[2] + x * b_1[3]};
        b_1_mat.sto(0, 0, -b_1_arr[0] * x - b_1_arr[1] * y);
        b_1_mat.sto(0, 1, -b_1_arr[0] * y + b_1_arr[1] * x);
        b_1_mat.sto(1, 0, -b_1_arr[0] * y + b_1_arr[1] * x);
        b_1_mat.sto(1, 1, b_1_arr[0] * x + b_1_arr[1] * y);

        b_2_temp.assign((inv(b_1_mat, temp) * del_mat)(0));

        if (jacobi(p * (b_2_temp[0] * b_2_temp[0] + b_2_temp[1] * b_2_temp[1]) * S[1], n) == 1) {
            subtract(b_1[0], b_1[0], b_1[3]);
            subtract(b_1[1], b_1[1], b_1[2]);
            multiply(b_1[2], b_1[2], 2);
            multiply(b_1[3], b_1[3], 2);
            break;
        }
    }

    prime.assign(S[1]);
    temp.assign(p * (b_2_temp[0] * b_2_temp[0] + b_2_temp[1] * b_2_temp[1]));
    divide(prime, prime, temp);
    prime.assign(sqrt(prime));
    bigint mu(prime.lift_to_Z()); // name lambda already taken

    temp.assign((S[1] - temp * mu * mu) / n);
    base_vector<bigint> b_2(4);
    while (true) {
        b_2[2].randomize(n);
        prime.assign(temp - 2 * p * mu * b_2_temp[0] * b_2[2]);
        b_2[3].assign((prime / (2 * p * mu * b_2_temp[1])).lift_to_Z());
        bigint store_0(mu * b_2_temp[0] + n * b_2[2]);
        bigint store_1(mu * b_2_temp[1] + n * b_2[3]);
        bigint store((S[1] - p * (store_0 * store_0 + store_1 * store_1)) / (n * n));
        // big enough? if not we'll increase S[2] even further.
        if (!store.is_prime())
            continue;
        if (cornacchia(b_2[0], b_2[1], bigint(-4), store)) {
            divide(b_2[0], b_2[0], 2); // because of LiDIA cornacchia
            multiply(b_2[0], b_2[0], n);
            multiply(b_2[1], b_2[1], n);
            add(b_2[2], b_2[2] * n, mu * b_2_temp[0]);
            add(b_2[3], b_2[3] * n, mu * b_2_temp[1]);
            subtract(b_2[0], b_2[0], b_2[3]);
            subtract(b_2[1], b_2[1], b_2[2]);
            multiply(b_2[2], b_2[2], 2);
            multiply(b_2[3], b_2[3], 2);
            break;
        }
    }

    I = (times(conjugate((times(b_2, p_coeff) * bigint_matrix(b_1))(0)), p_coeff) * I) / n;

    // ideal path has already been cleared, at end of previous loop...
    for (int k = 0; k < 2; k++)
        for (int i = 0; i < P[k].size(); i++)
            ideal_step(iota, P[k][i], I); // will push internally
}

void walker::reset() {
    random.erase(random.begin() + P[0].size(), random.end());
    ideal.clear();
}

string walker::w(bool challenge) { // challenge: which path are we dealing with?
    // in the future this can be replaced by some kind of more fancy representation of the isogeny path
    stringstream w;
    if (!challenge) {
        for (int i = 1; i < P[1].size(); i += 2)
            w << random[P[0].size() + i].E_i.j_invariant() << " ";
        if (P[1].size() % 2 == 1)
            w << random[P[0].size() + P[1].size() - 1].E_i.j_invariant() << " ";
    } else {
        for (int i = 1; i < P[0].size() + P[1].size(); i += 2)
            w << ideal[i].E_i.j_invariant() << " ";
        if (ideal.size() % 2 == 1)
            w << ideal[ideal.size() - 1].E_i.j_invariant() << " ";
    }

    return w.str();
}

walker::walker(int lambda, const elliptic_curve<gf_element> &E_, const gf_element &iota,
               const array<vector<torsion_basis>, 2> &P)
    : lambda(lambda), E_(E_), iota(iota), P(P), random(E_), ideal(E_) {
    for (int i = 0; i < P[0].size(); i++)
        random_step(P[0][i]);
}

string g(string w, int bits) {
    byte digest[SHA512::DIGESTSIZE];
    SHA512().CalculateDigest(digest, (const byte *)w.c_str(), w.length());

    CryptoPP::HexEncoder encoder;
    std::string output;
    encoder.Attach(new CryptoPP::StringSink(output));
    encoder.Put(digest, bits / 8 + 1); // number of bytes...?
    encoder.MessageEnd();
    return output;
}

bool decoder(string h, int i) {
    switch (h.at(i / 4)) {
    case '0':
        return false;
    case '1':
        return i % 4 == 0;
    case '2':
        return i % 4 == 1;
    case '3':
        return i % 4 <= 1;
    case '4':
        return i % 4 == 2;
    case '5':
        return i % 2 == 0;
    case '6':
        return i % 4 == 1 || i % 4 == 2;
    case '7':
        return i % 4 != 3;
    case '8':
        return i % 4 == 3;
    case '9':
        return i % 4 == 0 || i % 4 == 3;
    case 'A':
        return i % 2 == 1;
    case 'B':
        return i % 4 != 2;
    case 'C':
        return i % 4 >= 2;
    case 'D':
        return i % 4 != 1;
    case 'E':
        return i % 4 != 0;
    case 'F':
        return true;
    }
    return false; // just to silence the warning
}

string walker::sign(string message) {
    gf_element::set_output_format(0);

    int t = 3 * lambda;
    bigint temp(iota.get_field().characteristic());
    ceil(temp, log(bigfloat(temp)) / log(bigfloat(2)));
    array<int, 2> N;
    (2 * (P[1].size() / 2 + P[1].size() % 2) * temp).intify(N[0]);
    (2 * ((P[0].size() + P[1].size()) / 2 + (P[0].size() + P[1].size()) % 2) * temp).intify(N[1]);

    stringstream ss, gs;
    ss << message << " ";
    ss << random.back().j_invariant() << " ";

    array<vector<string>, 2> rsp_store;
    array<vector<string>, 2> g_store;
    for (int l = 0; l < t; l++) { // parallelize
        cout << "extending the walk..." << endl;
        for (int i = 0; i < P[1].size(); i++)
            random_step(P[1][i]); // extend the walk
        cout << "rerouting..." << endl;
        reroute();
        rsp_store[0].push_back(w(false));
        rsp_store[1].push_back(w(true));

        ss << random.back().j_invariant() << " ";
        g_store[0].push_back(g(rsp_store[0].back(), N[0]));
        g_store[1].push_back(g(rsp_store[1].back(), N[1]));
        gs << g_store[0].back() << " ";
        gs << g_store[1].back() << " ";
        reset();
    }

    CryptoPP::SHA512 H;
    H.Update((const byte *)ss.str().c_str(), ss.str().length());
    H.Update((const byte *)gs.str().c_str(), gs.str().length());
    byte digest[SHA512::DIGESTSIZE];
    H.Final(digest);

    CryptoPP::HexEncoder encoder;
    string h;
    encoder.Attach(new CryptoPP::StringSink(h));
    encoder.Put(digest, t / 8 + 1); // number of bytes...?
    encoder.MessageEnd();

    stringstream sigma;
    sigma << h << " ";
    for (int l = 0; l < t; l++) {
        sigma << rsp_store[decoder(h, l)][l]; // end of w already includes a space
    }
    for (int l = 0; l < t; l++) {
        sigma << g_store[!decoder(h, l)][l] << " ";
    }
    return sigma.str();
}

bool walker::verify(string message, string signature) {
    int t = 3 * lambda;
    galois_field K(iota.get_field());
    modular mod(K);
    bigint temp(K.characteristic());
    ceil(temp, log(bigfloat(temp)) / log(bigfloat(2)));
    array<int, 2> N;
    (2 * (P[1].size() / 2 + P[1].size() % 2) * temp).intify(N[0]);
    (2 * ((P[0].size() + P[1].size()) / 2 + (P[0].size() + P[1].size()) % 2) * temp).intify(N[1]);

    stringstream ss, sigstream(signature);
    ss << message << " ";
    ss << random.back().j_invariant() << " ";
    string h;
    sigstream >> h;
    vector<string> rsp_store;
    gf_element source(K), target(K);
    for (int l = 0; l < t; l++) { // parallelize
        stringstream rsp_stream, temp_stream;
        bool h_i = decoder(h, l);
        if (!h_i) {
            target = public_key();
            for (int i = 1; i < P[1].size(); i += 2) {
                source = target;
                sigstream >> target;
                if (gcd(mod.build(P[1][i - 1].le, source), mod.build(P[1][i].le, target))
                        .degree() == 0)
                    return false;
                rsp_stream << target << " ";
            }
            if (P[1].size() % 2 == 1) {
                source = target;
                sigstream >> target;
                if (!mod.build(P[1][P[1].size() - 1].le, source)(target).is_zero())
                    return false;
                rsp_stream << target << " ";
            }
        } else {
            target = E_.j_invariant();
            for (int i = 1; i < P[0].size(); i += 2) {
                source = target;
                sigstream >> target;
                if (gcd(mod.build(P[0][i - 1].le, source), mod.build(P[0][i].le, target))
                        .degree() == 0)
                    return false;
                rsp_stream << target << " ";
            }
            if (P[0].size() % 2 == 1) {
                source = target;
                sigstream >> target;
                if (gcd(mod.build(P[0][P[0].size() - 1].le, source), mod.build(P[1][0].le, target))
                        .degree() == 0)
                    return false;
                rsp_stream << target << " ";
            }
            for (int i = 1 + P[0].size() % 2; i < P[1].size(); i += 2) {
                source = target;
                sigstream >> target;
                if (gcd(mod.build(P[1][i - 1].le, source), mod.build(P[1][i].le, target))
                        .degree() == 0)
                    return false;
                rsp_stream << target << " ";
            }
            if (P[1].size() % 2 != P[0].size() % 2) {
                source = target;
                sigstream >> target;
                if (!mod.build(P[1][P[1].size() - 1].le, source)(target).is_zero())
                    return false;
                rsp_stream << target << " ";
            }
        }
        rsp_store.push_back(rsp_stream.str());
        ss << target << " ";
    }
    for (int l = 0; l < t; l++) {
        string temp;
        sigstream >> temp;
        if (!decoder(h, l)) {
            ss << g(rsp_store[l], N[0]) << " ";
            ss << temp << " ";
        } else {
            ss << temp << " ";
            ss << g(rsp_store[l], N[1]) << " ";
        }
    }

    byte digest[SHA512::DIGESTSIZE];
    SHA512().CalculateDigest(digest, (const byte *)ss.str().c_str(), ss.str().length());

    CryptoPP::HexEncoder encoder;
    string my_h;
    encoder.Attach(new CryptoPP::StringSink(my_h));
    encoder.Put(digest, t / 8 + 1); // number of bytes...?
    encoder.MessageEnd();
    return h.compare(my_h) == 0;
}
