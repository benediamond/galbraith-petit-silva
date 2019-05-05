#include "walker.h"

using CryptoPP::byte;
using CryptoPP::SHA3_512;

void walker::update_ideal(const torsion_basis &P, const point<gf_element> &Q_i) {
    bigint_matrix I_i_minus_1(random.back().I_i);

    bigint one_over_two((P.l + 1) / 2);
    galois_field mod_l(P.l);
    gf_element temp(mod_l);
    temp.assign(O_.q);
    bigint one_over_q((1 / temp).lift_to_Z());
    point<gf_element> i_Q_i(O_.theta(Q_i, P.ext));
    point<gf_element> j_Q_i(O_.pi(Q_i));
    point<gf_element> k_Q_i(-O_.pi(i_Q_i)); // O_.theta(j_q, P.ext)

    array<point<gf_element>, 4> B_Q_i{one_over_two * (Q_i + i_Q_i), one_over_two * (j_Q_i - k_Q_i),
                                      one_over_q * (i_Q_i - O_.c * k_Q_i), -k_Q_i};
    base_vector<bigint> a_i(4);
    int count = 0;
    do {
        while (true) { // { do {
            // coeffs.randomize(P.le); // unfortunately, this function is completely broken
            for (int i = 0; i < 4; i++)
                a_i[i].assign(randomize(bigint(P.l)));
            a_i.assign((I_i_minus_1 * bigint_matrix(a_i))(0)); // get 0th column
            bigint norm = O_.norm(a_i);
            if (norm % P.l == 0 && norm % (P.l * P.l) != 0)
                break;
        }
    } while (!(a_i[0] % P.l * B_Q_i[0] + a_i[1] % P.l * B_Q_i[1] + a_i[2] % P.l * B_Q_i[2] +
               a_i[3] % P.l * B_Q_i[3])
                  .is_zero());
    bigint_lattice I_i_temp(4, 8);
    bigint_matrix O_0_a_i = O_.times(a_i);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++)
            I_i_temp.sto(i, j, P.l * I_i_minus_1[i][j]);
        for (int j = 4; j < 8; j++)
            I_i_temp.sto(i, j, O_0_a_i[i][j - 4]);
    }
    I_i_temp.lll(0.99);
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            random.back().I_i.sto(i, j, I_i_temp[i][j]);
    // omitting n once again for now. if you want it back, use det
}

void walker::ideal_step(const torsion_basis &P, const bigint_matrix &I_) {
    bigint one_over_two((P.l + 1) / 2);
    galois_field mod_l(P.l);
    gf_element temp(mod_l);
    temp.assign(O_.q);
    bigint one_over_q((1 / temp).lift_to_Z());

    array<point<gf_element>, 2> i_P{O_.theta(P.P[0], P.ext), O_.theta(P.P[1], P.ext)};
    array<point<gf_element>, 2> j_P{O_.pi(P.P[0]), O_.pi(P.P[1])};
    array<point<gf_element>, 2> k_P{-O_.pi(i_P[0]), -O_.pi(i_P[1])};

    array<point<gf_element>, 2> one_i_2{one_over_two * (P.P[0] + i_P[0]),
                                        one_over_two * (P.P[1] + i_P[1])};
    array<point<gf_element>, 2> j_m_k_2{one_over_two * (j_P[0] - k_P[0]),
                                        one_over_two * (j_P[1] - k_P[1])};
    array<point<gf_element>, 2> i_m_ck_q{one_over_q * (i_P[0] - O_.c * k_P[0]),
                                         one_over_q * (i_P[1] - O_.c * k_P[1])};
    array<array<point<gf_element>, 4>, 2> B_P_i{{{one_i_2[0], j_m_k_2[0], i_m_ck_q[0], -k_P[0]},
                                                 {one_i_2[1], j_m_k_2[1], i_m_ck_q[1], -k_P[1]}}};

    rational_factorization l_fact(trialdiv(P.l));
    array<bigint, 2> coeffs;
    for (int j = 0; j < 4; j++) {
        array<point<gf_element>, 2> a_P_i{
            (I_[0][j] % P.l) * B_P_i[0][0] + (I_[1][j] % P.l) * B_P_i[0][1] +
                (I_[2][j] % P.l) * B_P_i[0][2] + (I_[3][j] % P.l) * B_P_i[0][3],
            (I_[0][j] % P.l) * B_P_i[1][0] + (I_[1][j] % P.l) * B_P_i[1][1] +
                (I_[2][j] % P.l) * B_P_i[1][2] + (I_[3][j] % P.l) * B_P_i[1][3]};
        if (a_P_i[0].is_zero() && a_P_i[1].is_zero())
            continue; // at least one must annihilate exactly Q...
        // but can we assert that it must be j = 0...? not clear...
        int idx = a_P_i[0].is_zero() ? 1 : 0;
        coeffs[idx] = bg_algorithm(a_P_i[idx], a_P_i[1 - idx], 0, P.l - 1, l_fact);
        coeffs[1 - idx] = -1;
        break;
    }

    point<gf_element> Q_i(coeffs[0] * P.P[0] + coeffs[1] * P.P[1]);

    // array<gf_element, 2> Q_coords{Q_i.get_x(), Q_i.get_y()};
    for (int j = 0; j < ideal.size(); j++) // carry forward through past walks
        ideal[j].isogeny(Q_i, P.ext);
    ideal.emplace_back(P, Q_i, ideal.back());
}

void walker::random_step(const torsion_basis &P) {
    array<int, 2> coeffs;
    random_generator rg; // store this somewhere once and for all...?
    do {
        rg >> coeffs[0];
        rg >> coeffs[1];
        coeffs[0] %= P.l; // warning:
        coeffs[1] %= P.l; // modulo bias
    } while (coeffs[0] == 0 || coeffs[1] == 0);
    // tricky proof shows this is equivalent to choosing a _subgroup_ uniformly
    point<gf_element> Q_i(coeffs[0] * P.P[0] + coeffs[1] * P.P[1]);
    point<gf_element> Q_isogeny(Q_i);
    for (int j = 0; j < random.size(); j++)  // carry forward through past walks
        random[j].isogeny(Q_isogeny, P.ext); // modifies Q_i in-place
    random.emplace_back(P, Q_isogeny, random.back());
    update_ideal(P, Q_i);
}

void walker::reroute() {
    // warning: rewrites steps, but does _not_ regenerate
    bigint_matrix I_(random.back().I_i);
    // NOT YET minkowski basis, let alone a random one.
    array<bigint, 2> S{1, 1};
    for (int k = 0; k < 2; k++) { // different from the usual def. of S...
        multiply(S[0], S[0], P[k][0].l);
        multiply(S[0], S[0], P[k][1].l);
        for (int i = 2; i < P[k].size(); i++)
            multiply(S[1], S[1], P[k][i].l);
    }

    int tries = 1000; // todo: set this more intelligently.
    while (true) {
        bigint n(I_.det()), m, temp; // silly.
        n.abs();
        sqrt(n, n);
        ceil(m, (log(bigfloat(O_.p))));

        bool success = true;
        base_vector<bigint> del(4);
        while (true) { // if no prime found, then problem...
            for (int i = 0; i < 4; i++)
                del[i].assign(randomize(2 * m) - m); // cleaner than exhaustive, and still fast
            del = (I_ * bigint_matrix(del))(0);
            bigint norm = O_.norm(del); // used to be the simple norm
            if (jacobi(-O_.q, norm / n) == 1 && (norm / n).is_prime())
                // so that x_0, y_0 can be found below
                break;
        }

        bigint_matrix I_n = (O_.times(O_.conjugate(del)) * I_) / n;
        // _not actually necessary_ for J = I' (see paper)!

        temp.assign(I_n.det());
        temp.abs();
        sqrt(n, temp); // n will temporarily alias n(I_'). will be recovered later
        while (true) {
            for (int i = 0; i < 4; i++)
                del[i].assign(randomize(m));
            del = (I_n * bigint_matrix(del))(0); // reuse del, actually equals "alpha"...
            bigint norm = O_.norm(del);
            if (gcd(norm, temp) == n)
                break;
        }

        galois_field mod_n(n, 1);
        gf_element c(mod_n);
        c.assign(-O_.q);
        c.assign(sqrt(c));

        galois_field mod_q(O_.q, 1);
        gf_element k(mod_q);
        k.assign(-del[2]);
        divide(k, k, n);

        // note: using the first row and first col. could fail if identically zero... unlikely
        gf_element x_0(mod_n), y_0(mod_n);
        x_0.assign(1 - O_.p);
        divide(x_0, x_0, 2);
        y_0.assign(1 + O_.p);
        divide(y_0, y_0, 2 * c);
        bigint x, y;
        x.assign(x_0.lift_to_Z());
        y.assign(y_0.lift_to_Z());
        // x^2 + q * y^2 should == -p (mod n).

        del[0] += del[0] % 2 == 0 ? 0 : n;
        del[1] += del[1] % 2 == 0 ? 0 : n;
        del[2] += k.lift_to_Z() * n; // ASSERT: should be divisible by q
        // ^^^ replace delta with a fellow representative (mod N) which also lives in <1, i, j, k>.
        // this is necessarily possible because N and | O_0 / <1, i, j, k> | (4 * q I think) are coprime.

        base_vector<bigint> del_1ijk(4);
        del_1ijk[0] = del[0] / 2;
        del_1ijk[1] = del[0] / 2 + del[2] / O_.q;
        del_1ijk[2] = del[1] / 2;
        del_1ijk[3] = -del[1] / 2 - O_.c * del[2] / O_.q - del[3];
        // re-express this element directly in terms of its coords w.r.t. <1, i, j, k>.

        bigmod_matrix del_mat(1, 2, n); // embedding into in M_2(Z / NZ)
        del_mat.sto(0, 0, del_1ijk[0] + x * del_1ijk[2] + O_.q * y * del_1ijk[3]);
        del_mat.sto(0, 1, del_1ijk[1] + y * del_1ijk[2] - x * del_1ijk[3]);
        // finally, embed del into M_2(Z / NZ) using the matrix on your page...

        sqrt(m, n * S[0] / (2 * O_.p * O_.q)); // this is rather conservative.
        // will work, but could fail to find _anything_ if S[0] isn't big enough.
        // could make things tighter, but could be a pain. revisit
        base_vector<bigint> b_1(4), b_temp(4), b_2_temp(2);
        for (int i = 0; i <= tries; i++) {
            if (i == tries) {
                success = false;
                break; // same as continue
            }
            b_temp[2].assign(randomize(2 * m) - m);
            b_temp[3].assign(randomize(2 * m) - m);
            temp.assign(n * S[0] - O_.p * (b_temp[2] * b_temp[2] + O_.q * b_temp[3] * b_temp[3]));
            if (!temp.is_prime())
                continue;
            if (cornacchia(b_temp[0], b_temp[1], -4 * O_.q, temp)) {
                divide(b_temp[0], b_temp[0], 2); // because of LiDIA cornacchia
                bigmod_matrix b_1_mat(2, 2, n);
                array<bigint, 2> b_1_arr{b_temp[0] + x * b_temp[2] + O_.q * y * b_temp[3],
                                         b_temp[1] + y * b_temp[2] - x * b_temp[3]};
                b_1_mat.sto(0, 0, b_1_arr[0] * x + b_1_arr[1] * O_.q * y);
                b_1_mat.sto(0, 1, b_1_arr[0] * y - b_1_arr[1] * x);
                b_1_mat.sto(1, 0, b_1_arr[0] * O_.q * y - b_1_arr[1] * O_.q * x);
                b_1_mat.sto(1, 1, -b_1_arr[0] * x - b_1_arr[1] * O_.q * y);

                b_2_temp.assign((del_mat * inv(b_1_mat, temp))[0]); // 0th row

                if (jacobi(O_.p * (b_2_temp[0] * b_2_temp[0] + O_.q * b_2_temp[1] * b_2_temp[1]) *
                               S[1],
                           n) == 1) {
                    b_1[0] = 2 * b_temp[0];
                    b_1[1] = 2 * b_temp[2];
                    b_1[2] = O_.q * (b_temp[1] - b_temp[0]);
                    b_1[3] = O_.c * (b_temp[0] - b_temp[1]) - b_temp[2] - b_temp[3];
                    break; // convert b_1 back into coordinates w.r.t. O_0, as opposed to 1 i j and k.
                }
            }
        }
        if (!success) // first cornacchia fail, try again
            continue;

        c.assign(S[1]);
        temp.assign(O_.p * (b_2_temp[0] * b_2_temp[0] + O_.q * b_2_temp[1] * b_2_temp[1]));
        divide(c, c, temp);
        c.assign(sqrt(c));
        bigint mu(c.lift_to_Z()); // name lambda already taken

        temp.assign((S[1] - temp * mu * mu) / n);
        base_vector<bigint> b_2(4);
        for (int i = 0; i <= tries; i++) {
            if (i == tries) {
                success = false;
                break; // same as continue
            }
            b_temp[2].randomize(n); // watch list. try to shrink the sizes of these!
            c.assign(temp - 2 * mu * O_.p * b_2_temp[0] * b_temp[2]);
            b_temp[3].assign((c / (2 * mu * O_.p * O_.q * b_2_temp[1])).lift_to_Z());
            bigint store_0(mu * b_2_temp[0] + n * b_temp[2]);
            bigint store_1(mu * b_2_temp[1] + n * b_temp[3]);
            bigint store((S[1] - O_.p * (store_0 * store_0 + O_.q * store_1 * store_1)) / (n * n));
            // big enough? if not we'll increase S[2] even further.
            if (!store.is_prime())
                continue;
            if (cornacchia(b_temp[0], b_temp[1], -4 * O_.q, store)) {
                divide(b_temp[0], b_temp[0], 2); // because of LiDIA cornacchia
                multiply(b_temp[0], b_temp[0], n);
                multiply(b_temp[1], b_temp[1], n);
                add(b_temp[2], b_temp[2] * n, mu * b_2_temp[0]);
                add(b_temp[3], b_temp[3] * n, mu * b_2_temp[1]);
                b_2[0] = 2 * b_temp[0];
                b_2[1] = 2 * b_temp[2];
                b_2[2] = O_.q * (b_temp[1] - b_temp[0]);
                b_2[3] = O_.c * (b_temp[0] - b_temp[1]) - b_temp[2] - b_temp[3];
                break;
            }
        }
        if (!success) // second cornacchia fail, try again
            continue;

        I_ = (O_.times(O_.conjugate((O_.times(b_2) * bigint_matrix(b_1))(0))) * I_n) / n;
        break;
    }
    // ideal path has already been cleared, at end of previous loop...
    for (int k = 0; k < 2; k++) // key: using 3 here...
        for (int i = 0; i < P[k].size(); i++)
            ideal_step(P[k][i], I_); // will push internally
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

walker::walker(int lambda, const elliptic_curve<gf_element> &E_,
               const array<vector<torsion_basis>, 2> &P, const maximal_order &O_)
    : lambda(lambda), E_(E_), P(P), O_(O_), random(E_), ideal(E_) {
    // not bothering with n. I_ will only be _read_ from, once, into a lattice...
    for (int i = 0; i < P[0].size(); i++)
        random_step(P[0][i]);
}

string g(string w, int bits) {
    byte digest[SHA3_512::DIGESTSIZE];
    SHA3_512().CalculateDigest(digest, (const byte *)w.c_str(), w.length());

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
    bigint temp(O_.p);
    ceil(temp, log(bigfloat(temp)) / log(bigfloat(2)));
    array<int, 2> N;
    (2 * (P[1].size() / 2 + P[1].size() % 2) * temp).intify(N[0]);
    (2 * ((P[0].size() + P[1].size()) / 2 + (P[0].size() + P[1].size()) % 2) * temp).intify(N[1]);

    stringstream ss, gs;
    ss << message << " ";
    ss << public_key() << " ";

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

        ss << public_key() << " ";
        g_store[0].push_back(g(rsp_store[0].back(), N[0]));
        g_store[1].push_back(g(rsp_store[1].back(), N[1]));
        gs << g_store[0].back() << " ";
        gs << g_store[1].back() << " ";
        reset();
    }

    SHA3_512 H;
    H.Update((const byte *)ss.str().c_str(), ss.str().length());
    H.Update((const byte *)gs.str().c_str(), gs.str().length());
    byte digest[SHA3_512::DIGESTSIZE];
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
    galois_field K(public_key().get_field()); // a bit ugly
    modular mod(K);                           // only going to be using primes now. nice
    bigint temp(K.characteristic());
    ceil(temp, log(bigfloat(temp)) / log(bigfloat(2)));
    array<int, 2> N;
    (2 * (P[1].size() / 2 + P[1].size() % 2) * temp).intify(N[0]);
    (2 * ((P[0].size() + P[1].size()) / 2 + (P[0].size() + P[1].size()) % 2) * temp).intify(N[1]);

    stringstream ss, sigstream(signature);
    ss << message << " ";
    ss << public_key() << " ";
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
                if (gcd(mod.build(P[1][i - 1].l, source), mod.build(P[1][i].l, target)).degree() ==
                    0)
                    return false;
                rsp_stream << target << " ";
            }
            if (P[1].size() % 2 == 1) {
                source = target;
                sigstream >> target;
                if (!mod.build(P[1][P[1].size() - 1].l, source)(target).is_zero())
                    return false;
                rsp_stream << target << " ";
            }
        } else {
            target = E_.j_invariant();
            for (int i = 1; i < P[0].size(); i += 2) {
                source = target;
                sigstream >> target;
                if (gcd(mod.build(P[0][i - 1].l, source), mod.build(P[0][i].l, target)).degree() ==
                    0)
                    return false;
                rsp_stream << target << " ";
            }
            if (P[0].size() % 2 == 1) {
                source = target;
                sigstream >> target;
                if (gcd(mod.build(P[0][P[0].size() - 1].l, source), mod.build(P[1][0].l, target))
                        .degree() == 0)
                    return false;
                rsp_stream << target << " ";
            }
            for (int i = 1 + P[0].size() % 2; i < P[1].size(); i += 2) {
                source = target;
                sigstream >> target;
                if (gcd(mod.build(P[1][i - 1].l, source), mod.build(P[1][i].l, target)).degree() ==
                    0)
                    return false;
                rsp_stream << target << " ";
            }
            if ((P[0].size() + P[1].size()) % 2 == 1) {
                source = target;
                sigstream >> target;
                if (!mod.build(P[1][P[1].size() - 1].l, source)(target).is_zero())
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

    byte digest[SHA3_512::DIGESTSIZE];
    SHA3_512().CalculateDigest(digest, (const byte *)ss.str().c_str(), ss.str().length());

    CryptoPP::HexEncoder encoder;
    string my_h;
    encoder.Attach(new CryptoPP::StringSink(my_h));
    encoder.Put(digest, t / 8 + 1); // number of bytes...?
    encoder.MessageEnd();
    return h.compare(my_h) == 0;
}
