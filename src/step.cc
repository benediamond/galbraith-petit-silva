#include "step.h"

using namespace std;
using namespace LiDIA;

step::step(const elliptic_curve<gf_element> &E_, const bigint_matrix &I_) : E_i(E_), I_i(I_) {}

step::step(const torsion_basis &P, const point<gf_element> &Q, const step &last)
    : P(P), psi(P.ext.get_L()), phi(P.ext.get_K()), omega(P.ext.get_K()), I_i(last.I_i) {
    // see Galbraith, Thm. 25.1.6., Kohel, Ch. 2.4
    elliptic_curve<gf_element> E_i_minus_1(P.ext.embed(last.E_i));
    // Q should _already_ be a point on this curve?

    gf_element tG(P.ext.get_L()), wG(P.ext.get_L());
    point<gf_element> Q_accum(Q);
    psi.assign_one();
    polynomial<gf_element> x(P.ext.get_L());
    x.assign_x();
    for (int j = 0; j < (P.l - 1) / 2; j++) { // iterate over G_1
        gf_element F_x(3 * Q_accum.get_x() * Q_accum.get_x() + E_i_minus_1.get_a4());
        gf_element F_y(-2 * Q_accum.get_y());
        gf_element uQ(F_y * F_y);
        gf_element tQ(2 * F_x);
        add(tG, tG, tQ);
        add(wG, wG, uQ + Q_accum.get_x() * tQ);

        multiply(psi, psi, x - Q_accum.get_x());
        // it looks like this is already sitting around, namely is nothing other than the divpol?

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
    temp_1.set_coefficient(P.l, 1);
    temp_1.set_coefficient(2 * psi[psi.degree() - 1], 0); // this won't be used for identity
    phi += temp_1 * psi * psi;

    multiply(omega, derivative(phi), psi);
    subtract(omega, omega, bigint(2) * phi * psi_prime);

    cout << "l = " << P.l << ", stepped: " << E_i.j_invariant() << endl;
}

// the "step" has already been constructed (using Velu) and the "coefficients" are ready.
void step::isogeny(point<gf_element> &Q, const field_extension &ext) const {
    // opting to do this using an in-place mutation, as opposed to returning a point.
    // slightly nicer to not have to keep an extra variable around i guess,
    // but also prevents the use of "return point<gf_element>(...)".
    if (Q.is_zero()) {
        Q.assign(point<gf_element>(ext.embed(E_i)));
        return;
    }

    polynomial<gf_element> lifted_psi(ext.embed(psi));
    polynomial<gf_element> lifted_phi(ext.embed(phi));
    polynomial<gf_element> lifted_omega(ext.embed(omega));

    gf_element psi_x(lifted_psi(Q.get_x()));
    if (psi_x.is_zero()) {
        Q.assign(point<gf_element>(ext.embed(E_i)));
        return;
    }

    gf_element psi_x_2, psi_x_3;
    square(psi_x_2, psi_x);
    multiply(psi_x_3, psi_x_2, psi_x);
    gf_element phi_x(lifted_phi(Q.get_x()));
    gf_element omega_x_y(lifted_omega(Q.get_x()));
    multiply(omega_x_y, omega_x_y, Q.get_y());

    gf_element temp_x, temp_y; // over L
    divide(temp_x, phi_x, psi_x_2);
    divide(temp_y, omega_x_y, psi_x_3);
    Q.assign(point<gf_element>(temp_x, temp_y, ext.embed(E_i)));
}

istream &operator>>(istream &in, step &step) {
    galois_field K;
    in >> K;

    step.psi.assign_zero(K);
    in >> step.psi;
    step.phi.assign_zero(K);
    in >> step.phi;
    step.omega.assign_zero(K);
    in >> step.omega;

    in >> step.E_i;

    return in;
}

ostream &operator<<(ostream &out, const step &step) {
    out << step.psi.get_field() << endl;
    out << step.psi << endl;
    out << step.phi << endl;
    out << step.omega << endl;

    out << step.E_i << endl;

    return out;
}
