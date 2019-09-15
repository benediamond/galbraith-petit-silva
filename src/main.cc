#include <LiDIA/random_generator.h>
#include <LiDIA/bigint_matrix.h>
#include <LiDIA/bigmod_matrix.h>
#include <LiDIA/bigint.h>
#include <LiDIA/galois_field.h>
#include <LiDIA/gf_element.h>
#include "gps.h"

#include <sstream>
#include <array>
#include <vector>

using namespace LiDIA;
using namespace std;

int main(int argc, const char *argv[]) {
    // cout << "enter random seed...";
    // unsigned int a;
    // cin >> a;
    // random_generator::seed(a);
    //
    // gps gps(8, 0.001);
    // ofstream out("./src/serializations/l = 4.txt");
    // out << gps;

    gps gps;
    ifstream in("./src/serializations/l = 4.txt");
    in >> gps;

    gps.gen_keypair();
    string signature = gps.keychain[0].sign("Hello World!");
    cout << signature << endl;
    cout << gps.keychain[0].verify("Hello World!", signature) << endl;
    return 0;
}
