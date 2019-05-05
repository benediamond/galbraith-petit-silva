#include "order.h"

maximal_order::maximal_order(const class step &step, const gf_element &u) : step(step), u(u) {
    p.assign(step.get_P().ext.get_p());
    q.assign(step.get_P().l);
    q_coeff.assign((q + 1) / 4);
    // c_coeff.assign((p * c * c + 1) / q);
}

void maximal_order::set_c(const bigint &new_c) {
    c.assign(new_c);
    c_coeff.assign((p * c * c + 1) / q);
}

point<gf_element> maximal_order::theta(const point<gf_element> &Q,
                                       const field_extension &ext) const {
    // Q can be over any field, which is given by ext.
    // array<gf_element, 2> Q_coords{Q.get_x(), Q.get_y()};
    point<gf_element> Q_theta(Q);
    step.isogeny(Q_theta, ext); // void / in-place
    // todo: what about if this results in the zero point.

    if (Q_theta.is_zero())
        return point<gf_element>(Q.get_curve());
    gf_element lifted_u = ext.embed(u);
    gf_element u_2, u_3;
    square(u_2, lifted_u);
    multiply(u_3, u_2, lifted_u);
    Q_theta.assign(point<gf_element>(u_2 * Q_theta.get_x(), u_3 * Q_theta.get_y(), Q.get_curve()));
    return Q_theta;
}

point<gf_element> maximal_order::pi(const point<gf_element> &Q) const {
    // careful about which degree frobenius this will do.
    // ideally, it will "intelligently" know that the curve's field of definition
    // is actually just 1, even though Q sits on an embedding of the base curve
    // into a field extension. verify this.
    if (Q.is_zero()) // why isn't this baked into Frob?
        return Q;
    return frobenius(Q);
}

bigint maximal_order::norm(const base_vector<bigint> &e) const {
    return q_coeff * (e[0] * e[0] + p * e[1] * e[1]) + c_coeff * e[2] * e[2] +
           q * p * (e[1] * e[3] + e[3] * e[3]) + e[0] * e[2] + p * c * e[1] * e[2] +
           2 * p * c * e[2] * e[3];
}

base_vector<bigint> maximal_order::conjugate(const base_vector<bigint> &e) const { // length 4
    array<bigint, 4> e_bar_arr{e[0], -e[1], -q * e[0] - e[2], c * e[0] - e[3]};
    base_vector<bigint> e_bar(4);
    for (int i = 0; i < 4; i++)
        e_bar[i].assign(e_bar_arr[i]); // conjugate b_1
    return e_bar;
}

bigint_matrix maximal_order::times(const base_vector<bigint> &e) const {
    bigint_matrix O_0_e(4, 4);
    O_0_e.sto(0, 0, -(q - 1) / 2 * e[0] - e[2]);
    O_0_e.sto(1, 0, (q + 1) / 2 * e[1] + c * e[2] + q * e[3]);
    O_0_e.sto(2, 0, q * q_coeff * e[0] + (q + 1) / 2 * e[2]);
    O_0_e.sto(3, 0, -c * q_coeff * e[0] - q_coeff * e[1] - c * e[2] - (q - 1) / 2 * e[3]);
    O_0_e.sto(0, 1, -p * (q + 1) / 2 * e[1] - c * p * e[2] - p * q * e[3]);
    O_0_e.sto(1, 1, -(q - 1) / 2 * e[0] - e[2]);
    O_0_e.sto(2, 1,
              q * p * q_coeff * e[1] + c * p * (q - 1) / 2 * e[2] + p * q * (q - 1) / 2 * e[3]);
    O_0_e.sto(3, 1,
              q_coeff * e[0] - c * p * q_coeff * e[1] + (-(q - 1) / 2 * c_coeff + 1) * e[2] -
                  c * p * (q - 1) / 2 * e[3]);
    O_0_e.sto(0, 2, -e[0] - c * p * e[1] - 2 * c_coeff * e[2] - 2 * c * p * e[3]);
    O_0_e.sto(1, 2, -c * e[0] + e[1] + 2 * e[3]);
    O_0_e.sto(2, 2,
              (q + 1) / 2 * e[0] + c * p * (q + 1) / 2 * e[1] + (c * c * p + 1) * e[2] +
                  c * p * q * e[3]);
    O_0_e.sto(3, 2, -(q + 1) / 2 * c_coeff * e[1] - c * c_coeff * e[2] - (c * c * p + 1) * e[3]);
    O_0_e.sto(0, 3, -p * q * e[1] - 2 * c * p * e[2] - 2 * p * q * e[3]);
    O_0_e.sto(1, 3, -q * e[0] - 2 * e[2]);
    O_0_e.sto(2, 3, p * q * (q + 1) / 2 * e[1] + q * c * p * e[2] + p * q * q * e[3]);
    O_0_e.sto(3, 3,
              (q + 1) / 2 * e[0] - c * p * (q + 1) / 2 * e[1] + (1 - c * c * p) * e[2] -
                  c * p * q * e[3]);
    return O_0_e;
}

istream &operator>>(istream &in, maximal_order &O_) {
    in >> O_.step;
    in >> O_.u;

    in >> O_.c;
    in >> O_.p >> O_.q >> O_.q_coeff >> O_.c_coeff;

    return in;
}

ostream &operator<<(ostream &out, const maximal_order &O_) {
    out << O_.step << endl;
    out << O_.u << endl;
    out << O_.c << endl;

    out << O_.p << " " << O_.q << " " << O_.q_coeff << " " << O_.c_coeff << endl;

    return out;
}
