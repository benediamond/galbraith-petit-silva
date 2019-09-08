#include <LiDIA/bigint.h>
#include <LiDIA/bigmod.h>

#include <fstream>
#include <iostream>
#include <sstream>

#pragma once

using namespace LiDIA;
using namespace std;

class modular { // made to mimic the interface of meq_prime.
    galois_field K;
    ifstream infile;

  public:
    modular(galois_field K) : K(K) { bigmod::set_modulus(K.characteristic()); }

    polynomial<gf_element> build(bigint l, const gf_element &x) {
        char l_str[l.bit_length() / 3 + 10];
        bigint_to_string(l, l_str);
        infile.close();

        infile.open(string("../src/modular/phi_j_") + l_str + string(".txt"));
        // todo: exception if reading fails! will give cryptic error msg.
        infile.clear();
        infile.seekg(0, ios::beg);

        polynomial<gf_element> f(K);
        rational_factorization l_fact(trialdiv(l));
        int deg;
        (l + 1).intify(deg); // revert: no longer need l^e + l^{e - 1}

        f.set_degree(deg);

        gf_element accumulator;
        accumulator.assign_one(K);
        vector<gf_element> xpow;
        for (int i = 0; i <= deg + 1; i++) {
            xpow.push_back(accumulator);
            multiply(accumulator, accumulator, x);
        }
        string line;
        while (getline(infile, line)) {
            lidia_size_t x_deg, y_deg;
            bigmod coeff;

            replace_if(line.begin(), line.end(),
                       [](char const &c) { return c == '[' || c == ',' || c == ']'; }, ' ');
            stringstream ss(line);
            ss >> x_deg >> y_deg >> coeff;

            f[y_deg] += coeff.mantissa() * xpow[x_deg];
            if (x_deg != y_deg)
                f[x_deg] += coeff.mantissa() * xpow[y_deg];
        }
        return f;
    }
};
