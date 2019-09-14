#include "extension.h"

using namespace std;
using namespace LiDIA;

field_extension::field_extension(const galois_field &K, const galois_field &L)
    : K(K), L(L), p(K.characteristic()), embedding(L.degree(), K.degree(), p),
      restriction(K.degree(), L.degree(), p) {
    if (L.characteristic() != p)
        throw std::invalid_argument("characteristics are unequal.");
    if (L.degree() % K.degree() != 0)
        throw std::invalid_argument("deg(K) must divide deg(L).");

    polynomial<gf_element> f_ext(K.irred_polynomial(), L);
    gf_element h = find_root(f_ext);

    gf_element accum(L);
    accum.assign_one();
    for (int j = 0; j < K.degree(); j++) {
        for (int i = 0; i < L.degree(); i++)
            // no good "get_data" routine, inner loop necessary
            embedding.sto(i, j, accum.polynomial_rep()[i]);
        multiply(accum, accum, h);
    }

    cout << "embedding: " << embedding << endl;
    bigint dummy;
    restriction = inv(trans(embedding) * embedding, dummy) * trans(embedding);
    // warning: not always the case that inverse exists! this gets tricky fast
    // see e.g. https://www.sciencedirect.com/science/article/pii/0012365X78901449
    // if fails, try again with different seed?
    cout << "restriction: " << restriction << endl;
}

gf_element field_extension::embed(const gf_element &a) const {
    if (a.get_field() != K)
        throw std::invalid_argument("a doesn't belong to the base field.");
    Fp_polynomial a_poly(a.polynomial_rep());
    bigmod_matrix a_vec(K.degree(), 1, p);
    for (int i = 0; i < K.degree(); i++)
        a_vec.sto(i, 0, a_poly[i]);
    bigmod_matrix b_vec = embedding * a_vec; // or use copy?
    Fp_polynomial b_poly;
    b_poly.set_modulus(p);
    for (int i = 0; i < L.degree(); i++)
        b_poly.set_coefficient(b_vec[i][0], i);

    gf_element b(L);
    b.set_polynomial_rep(b_poly);

    return b;
}

polynomial<gf_element> field_extension::embed(const polynomial<gf_element> &a) const {
    polynomial<gf_element> b(L);
    for (int i = 0; i <= a.degree(); i++)
        b.set_coefficient(embed(a[i]), i);
    return b;
}

elliptic_curve<gf_element> field_extension::embed(const elliptic_curve<gf_element> E) const {
    return elliptic_curve<gf_element>(embed(E.get_a4()), embed(E.get_a6()));
}

gf_element field_extension::restrict(const gf_element &b) const {
    if (b.get_field() != L)
        throw std::invalid_argument("b doesn't belong to the extension field.");
    Fp_polynomial b_poly(b.polynomial_rep());
    bigmod_matrix b_vec(L.degree(), 1, p);
    for (int i = 0; i < L.degree(); i++)
        b_vec.sto(i, 0, b_poly[i]);
    bigmod_matrix a_vec = restriction * b_vec; // or use copy?
    Fp_polynomial a_poly;
    a_poly.set_modulus(p);
    for (int i = 0; i < K.degree(); i++)
        a_poly.set_coefficient(a_vec[i][0], i);

    gf_element a(K);
    a.set_polynomial_rep(a_poly);

    if (embedding * a_vec != b_vec)
        throw std::invalid_argument("b doesn't belong to the image of K.");
    return a;
}

polynomial<gf_element> field_extension::restrict(const polynomial<gf_element> &b) const {
    polynomial<gf_element> a(K);
    for (int i = 0; i <= b.degree(); i++)
        a.set_coefficient(restrict(b[i]), i);
    return a;
}

bigint field_extension::get_p() const { return p; }

galois_field field_extension::get_K() const { return K; }

galois_field field_extension::get_L() const { return L; }

istream &operator>>(istream &in, field_extension &ext) {
    in >> ext.p;
    in >> ext.K;
    in >> ext.L;

    in >> ext.embedding;
    in >> ext.restriction;

    return in;
}

ostream &operator<<(ostream &out, const field_extension &ext) {
    out << ext.p << endl;
    out << ext.K << endl;
    out << ext.L << endl;

    out << ext.embedding.get_no_of_rows() << " " << ext.embedding.get_no_of_columns() << " ";
    out << ext.embedding.get_modulus();

    for (int i = 0; i < ext.embedding.get_no_of_rows(); i++)
        for (int j = 0; j < ext.embedding.get_no_of_columns(); j++)
            out << " " << ext.embedding[i][j];
    out << endl;

    out << ext.restriction.get_no_of_rows() << " " << ext.restriction.get_no_of_columns() << " ";
    out << ext.restriction.get_modulus();

    for (int i = 0; i < ext.restriction.get_no_of_rows(); i++)
        for (int j = 0; j < ext.restriction.get_no_of_columns(); j++)
            out << " " << ext.restriction[i][j];
    // no endl
    return out;
}
