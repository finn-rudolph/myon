use rand::{thread_rng, Rng};
use rug::{ops::Pow, Complete, Integer};

use crate::{
    gfpolynomial::{GfMpPolynomial, GfPolynomial},
    nt,
    polynomial::{MpPolynomial, Polynomial},
};

// Calculates the algebraic square root of the product of 'integers' using q-adic newton iteration.
// Uses divide and conquer to evaluate the product in O(M log n) time, where M is the time needed
// to multiply two numbers in the order of magnitude of the result.
pub fn algebraic_sqrt(integers: &Vec<MpPolynomial>, f: &MpPolynomial) -> MpPolynomial {
    let s = mul_algebraic_integers(integers, f);
    let p = select_p(f);
    let mut r = GfMpPolynomial::from(&inv_sqrt_mod_p(
        &GfPolynomial::from_mp_polynomial(&s, p),
        &GfPolynomial::from_mp_polynomial(&f, p),
    ));

    let max_coeffiecient = s.coefficients_ref().iter().max().unwrap();
    let mut q = Integer::from(p);

    while &q < max_coeffiecient {
        q.square_mut();
        let f_mod_q = GfMpPolynomial::from_mp_polynomial(f, q.clone());
        let mut t = f_mod_q.mul_mod(
            &GfMpPolynomial::from_mp_polynomial(&s, q.clone()),
            &f_mod_q.mul_mod(&r, &r),
        );
        t[0] = ((3u32 - &t[0]).complete() + &q) % &q;
        t = f_mod_q.mul_mod(&r, &t);

        let two_inv = Integer::from(2).invert(&q).unwrap();

        for coefficient in t.coefficients_mut() {
            *coefficient = (coefficient.clone() * &two_inv) % &q;
        }

        r = t;
    }

    MpPolynomial::from(r)
}

fn mul_algebraic_integers(integers: &[MpPolynomial], f: &MpPolynomial) -> MpPolynomial {
    if integers.len() == 1 {
        return integers.first().unwrap().clone();
    }
    f.mul_mod(
        &mul_algebraic_integers(&integers[..integers.len() / 2], f),
        &mul_algebraic_integers(&integers[integers.len() / 2..], f),
    )
}

fn select_p(f: &MpPolynomial) -> u32 {
    let mut p: u32 = 100000007;
    loop {
        // p must be inert in the number field, which means f must be irreducible mod p.
        if nt::miller_rabin(p) && GfPolynomial::from_mp_polynomial(f, p).is_irreducible() {
            return p;
        }
        p += 2;
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
        u.1[0] = 1;

        let mut v = (GfPolynomial::new(p), GfPolynomial::new(p));
        v.0[0] = 1;
        v.1[0] = 1;

        let mut e: Integer = (Integer::from(p).pow(d as u32) - 1) / 2;
        while e != 0 {
            if e.is_odd() {
                v = mul_y_polynomials(&u, &v, s, f);
            }
            u = mul_y_polynomials(&u, &u, s, f);
            e >>= 1;
        }

        let z = f.mul_mod(s, &f.mul_mod(&u.1, &u.1));
        if z.degree() == 0 && z[0] == 1 {
            return u.1;
        }
    }
}

// Caclculates the square root of the product of a set of rational integers.
pub fn rational_sqrt(integers: Vec<Integer>) -> Integer {
    mul_rational_integers(&integers).sqrt()
}

fn mul_rational_integers(integers: &[Integer]) -> Integer {
    if integers.len() == 1 {
        return integers.first().unwrap().clone();
    }
    mul_rational_integers(&integers[..integers.len() / 2])
        * mul_rational_integers(&integers[integers.len() / 2..])
}
