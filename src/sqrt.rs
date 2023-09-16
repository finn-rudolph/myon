use std::cmp::max;

use log::{info, warn};
use rand::{thread_rng, Rng};
use rug::{
    ops::{NegAssign, Pow},
    Integer,
};

use crate::{
    gfpolynomial::GfPolynomial,
    nt,
    polynomial::{MpPolynomial, Polynomial},
};

// Calculates the algebraic square root of the product of s using q-adic newton iteration.
// Uses divide and conquer to evaluate the product in O(M log n) time, where M is the time needed
// to multiply two numbers in the order of magnitude of the result.
pub fn algebraic_sqrt(s: &MpPolynomial, f: &MpPolynomial) -> Option<MpPolynomial> {
    let mut p: u64 = 3;

    // p must be inert in the number field, which means f must be irreducible mod p.
    while !nt::miller_rabin(p) || !GfPolynomial::from_mp_polynomial(f, p).is_irreducible() {
        p += 1;
    }
    info!("chose the prime for lifting p = {}", p);

    let mut r = MpPolynomial::from(&inv_sqrt_mod_p(
        &GfPolynomial::from_mp_polynomial(&s, p),
        &GfPolynomial::from_mp_polynomial(&f, p),
    ));

    let num_iterations = (s
        .coefficients_ref()
        .iter()
        .fold(Integer::new(), |acc, x| max(acc, x.clone().abs()))
        .significant_bits()
        / p.ilog2())
    .ilog2()
        + 4;

    info!("doing {} iterations of newtons method", num_iterations);

    let mut q = Integer::from(p);

    for _ in 0..num_iterations {
        q.square_mut();
        let mut t = f.mul_mod(s, &f.mul_mod(&r, &r).rem(&q)).rem(&q);

        t[0] -= 3;
        let two_inv = Integer::from(2).invert(&q).unwrap();
        for coefficient in t.coefficients_mut() {
            coefficient.neg_assign();
            *coefficient += &q;
            *coefficient *= &two_inv;
            *coefficient %= &q;
        }

        r = f.mul_mod(&r, &t).rem(&q);

        // Perform a check that r is indeed the inverse square root of s mod q.
        let h = f.mul_mod(s, &f.mul_mod(&r, &r).rem(&q)).rem(&q);
        assert_eq!(h.degree(), 0);
        assert_eq!(h[0], 1);
    }

    let result_mod_q = f.mul_mod(s, &r).rem(&q);
    assert_eq!(
        f.mul_mod(&result_mod_q, &result_mod_q).rem(&q),
        s.clone().rem(&q)
    );

    for i in 0..(1 << f.degree()) {
        let mut result = MpPolynomial::new();
        for j in 0..f.degree() {
            result[j] = if (i >> j) & 1 == 1 {
                result_mod_q[j].clone() - &q
            } else {
                result_mod_q[j].clone()
            };
        }
        if f.mul_mod(&result, &result) == *s {
            return Some(result);
        }
    }

    warn!("newtons method failed");
    return None;
}

pub fn mul_algebraic_integers(integers: &[MpPolynomial], f: &MpPolynomial) -> MpPolynomial {
    if integers.len() == 1 {
        return integers.first().unwrap().clone();
    }
    f.mul_mod(
        &mul_algebraic_integers(&integers[..integers.len() / 2], f),
        &mul_algebraic_integers(&integers[integers.len() / 2..], f),
    )
}

// Compute a square root of s mod p (and, as always, mod f). The algorithm is from Jensen, P. L.
// (2005).
fn inv_sqrt_mod_p(s: &GfPolynomial, f: &GfPolynomial) -> GfPolynomial {
    let p = s.modulus();
    let d = f.degree();
    let mut rng = thread_rng();

    loop {
        let mut u = (GfPolynomial::new(p), GfPolynomial::new(p));
        for i in 0..d {
            u.0[i] = rng.gen_range(0..p);
        }
        while u.0[d - 1] == 0 {
            u.0[d - 1] = rng.gen_range(0..p);
        }
        u.1[0] = p - 1;

        let mut v = (GfPolynomial::new(p), GfPolynomial::new(p));
        v.0[0] = 1;

        let mut e: Integer = (Integer::from(p).pow(d as u32) - 1) / 2;
        while e != 0 {
            if e.is_odd() {
                v = mul_y_polynomials(&u, &v, s, f);
            }
            u = mul_y_polynomials(&u, &u, s, f);
            e >>= 1;
        }

        let g = f.mul_mod(s, &f.mul_mod(&v.1, &v.1));
        if g.degree() == 0 && g[0] == 1 {
            return v.1;
        }
    }
}

// Multiplies u * v modulo y^2 - s, where u, v and s are degree one polynomials in y, whose
// coefficients are polynomials in x. All operations on polynomials in x are done modulo f.
fn mul_y_polynomials(
    u: &(GfPolynomial, GfPolynomial),
    v: &(GfPolynomial, GfPolynomial),
    s: &GfPolynomial,
    f: &GfPolynomial,
) -> (GfPolynomial, GfPolynomial) {
    (
        f.mul_mod(&u.0, &v.0)
            .add(&f.mul_mod(&f.mul_mod(&u.1, &v.1), &s)),
        f.mul_mod(&u.0, &v.1).add(&f.mul_mod(&u.1, &v.0)),
    )
}

pub fn mul_rational_integers(integers: &[Integer]) -> Integer {
    if integers.len() == 1 {
        return integers.first().unwrap().clone();
    }
    mul_rational_integers(&integers[..integers.len() / 2])
        * mul_rational_integers(&integers[integers.len() / 2..])
}

// #[cfg(test)]
// mod tests {
//     use crate::params::MAX_DEGREE;

//     use super::*;

//     // TOOD: move random generation of polynomials (and csc matrix) to test code

//     #[test]
//     fn algebraic_sqrt_random() {
//         for degree in 3..MAX_DEGREE {
//             let f = MpPolynomial::new_random(degree, degree * 20);
//             for bits in (1000..10001).step_by(1000) {
//                 let g = MpPolynomial::new_random(degree - 1, bits);
//                 assert_eq!(g, algebraic_sqrt(&f.mul_mod(&g, &g), &f).unwrap());
//             }
//         }
//     }
// }
