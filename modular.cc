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

    polynomial<gf_element> build(bigint le, const gf_element &x) {
        char le_str[le.bit_length() / 3 + 10];
        bigint_to_string(le, le_str);
        infile.close();
        infile.open(string("modular/phi_j_") + le_str + string(".txt"));
        infile.clear();
        infile.seekg(0, ios::beg);

        polynomial<gf_element> f(K);
        rational_factorization le_fact(trialdiv(le));
        int deg;
        bigint temp;
        power(temp, le_fact.base(0), le_fact.exponent(0) - 1);
        (le + temp).intify(deg);

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
