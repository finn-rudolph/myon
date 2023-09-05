use std::cmp::min;

use log::info;
use rug::{
    ops::{NegAssign, Pow},
    Complete, Integer,
};

use crate::{
    lanczos,
    linalg::CscMatrixBuilder,
    nt,
    params::{Params, OVERSQUARENESS},
    polynomial::{self, MpPolynomial, Polynomial},
    sqrt,
};

fn rational_factor_base(m: &Integer, params: &Params) -> Vec<(u32, u32)> {
    let mut base: Vec<(u32, u32)> = Vec::new();

    let mut p: u32 = 2;
    while base.len() < params.rational_base_size {
        if nt::miller_rabin(p) {
            base.push((p, m.mod_u(p)));
        }
        p += 1;
    }

    base
}

fn algebraic_factor_base(f: &MpPolynomial, params: &Params) -> Vec<(u32, u32)> {
    let mut base: Vec<(u32, u32)> = Vec::new();
    let mut p: u32 = 2;

    while base.len() < params.algebraic_base_size {
        if nt::miller_rabin(p) {
            let roots = f.find_roots_mod_p(p);
            base.extend(roots.iter().map(|r| (p, *r)));
        }
        p += 1;
    }

    base.truncate(params.algebraic_base_size);
    base
}

fn quad_char_base(mut p: u32, f: &MpPolynomial, params: &Params) -> Vec<(u32, u32)> {
    let mut base: Vec<(u32, u32)> = Vec::new();
    let f_derivative = f.derivative();

    while base.len() < params.quad_char_base_size {
        if nt::miller_rabin(p) {
            let roots = f.find_roots_mod_p(p);
            for r in roots {
                if f_derivative.evaluate(r) % p != 0 {
                    base.push((p, r));
                }
            }
        }
        p += 1;
    }

    base.truncate(params.quad_char_base_size);
    base
}

fn ilog2_rounded(x: u32) -> u32 {
    ((x * x).ilog2() + 1) >> 1
}

fn line_sieve(b: u32, sieve_array: &mut Vec<i8>, base: &Vec<(u32, u32)>) {
    let a0: i32 = -(sieve_array.len() as i32 / 2);

    for (p, r) in base {
        let log2p = ilog2_rounded(*p) as i8;
        let mut i = ((-((b * r) as i32) % *p as i32 + *p as i32 - a0) % *p as i32) as usize;
        while i < sieve_array.len() {
            sieve_array[i] += log2p;
            i += *p as usize;
        }
    }
}

fn norm(f: &MpPolynomial, a: i32, b: u32) -> Integer {
    let d = f.degree();
    let mut u = Integer::from(1);
    let mut v = Integer::from(b).pow(d as u32);
    let mut result = Integer::new();

    for coefficient in f.coefficients_ref().iter().take(d + 1) {
        result += coefficient * (&u * &v).complete();
        u *= a;
        v /= b;
    }

    result
}

pub fn factorize(n: &Integer) -> Vec<Integer> {
    let params = Params::new(&n);
    let (f, m) = polynomial::select(n, &params);

    info!("set d = {}, m = {}", params.polynomial_degree, &m);
    info!("selected the polynomial {}", &f);

    // Maybe check that the polynomial is irreducible
    let rational_base = rational_factor_base(&m, &params);
    let algebraic_base = algebraic_factor_base(&f, &params);
    let quad_char_base = quad_char_base(algebraic_base.last().unwrap().0 + 1, &f, &params);

    let rational_begin: usize = 1;
    let algebraic_begin = rational_begin + rational_base.len();
    let quad_char_begin = algebraic_begin + algebraic_base.len();
    let base_len = quad_char_begin + quad_char_base.len();

    info!(
        "set up factor base consisting of {} primes on the rational side, {} prime ideals on the \
         algebraic side and {} quadratic characters (total size: {})",
        rational_base.len(),
        algebraic_base.len(),
        quad_char_base.len(),
        base_len
    );

    let mut matrix_builder = CscMatrixBuilder::new();
    matrix_builder.set_num_rows(quad_char_begin + quad_char_base.len());
    let mut relations: Vec<(u32, u32)> = Vec::new();

    let mut rational_sieve_array: Vec<i8> = vec![0; params.sieve_array_size];
    let mut algebraic_sieve_array: Vec<i8> = vec![0; params.sieve_array_size];

    for b in 1.. {
        rational_sieve_array
            .fill(-((ilog2_rounded(b) + m.significant_bits()) as i8) + params.rational_fudge);
        line_sieve(b, &mut rational_sieve_array, &rational_base);

        algebraic_sieve_array.fill(-params.algebraic_threshold);
        line_sieve(b, &mut algebraic_sieve_array, &algebraic_base);

        let a0 = -(params.sieve_array_size as i32 / 2);
        // Consider unsafe access here to avoid bounds checks.
        for i in 0..params.sieve_array_size {
            if rational_sieve_array[i] >= 0 && algebraic_sieve_array[i] >= 0 {
                let a = a0 + i as i32;
                let mut ones_pos: Vec<u32> = Vec::new();

                // Trial divide on the rational side.
                let mut num = a + (b * &m).complete();
                if num < 0 {
                    ones_pos.push(0);
                    num.neg_assign();
                }
                for (i, (p, _)) in rational_base.iter().enumerate() {
                    let e = num.remove_factor_mut(&Integer::from(*p));
                    if e & 1 == 1 {
                        ones_pos.push((rational_begin + i) as u32);
                    }
                }

                // Trial divide on the algebraic side.
                let mut alg_norm = norm(&f, a, b);
                for (i, (p, _)) in algebraic_base.iter().enumerate() {
                    let e = alg_norm.remove_factor_mut(&Integer::from(*p));
                    if e & 1 == 1 {
                        ones_pos.push((algebraic_begin + i) as u32);
                    }
                }

                if num == 1 && alg_norm == 1 {
                    // smooth pair (a, b) found!
                    for (i, (p, s)) in quad_char_base.iter().enumerate() {
                        if nt::legendre(
                            (((a + b as i32 * *s as i32) % *p as i32 + *p as i32) % *p as i32)
                                as u32,
                            *p,
                        ) == p - 1
                        {
                            ones_pos.push((quad_char_begin + i) as u32);
                        }
                    }
                    matrix_builder.add_col(ones_pos);
                    relations.push((a as u32, b));
                }
            }
        }

        info!("collected {} relations", relations.len());
        if relations.len() > base_len + OVERSQUARENESS {
            break;
        }
    }

    info!("collected {} relations", relations.len());

    let (mat, num_dependencies) = lanczos::find_dependencies(&matrix_builder.build());
    let mut factors: Vec<Integer> = Vec::new();

    info!("beginning sqrt phase");

    for i in 0..num_dependencies {
        let mut rational: Vec<Integer> = Vec::new();
        let mut algebraic: Vec<MpPolynomial> = Vec::new();

        for (j, (a, b)) in relations.iter().enumerate() {
            if (mat[j] >> i) & 1 == 1 {
                rational.push(a + (b * &m).complete());
                let mut g = MpPolynomial::new();
                g[0] = Integer::from(*a);
                g[1] = Integer::from(*b);
                algebraic.push(g);
            }
        }

        let a = sqrt::rational_sqrt(&rational);
        let b = sqrt::algebraic_sqrt(&algebraic, &f).evaluate(&m);

        let d = (&a - &b).complete().gcd(&n);
        if d != 1 && d != b {
            factors.push(min(d, (&b - &a).complete().gcd(&n)));
        }
    }

    factors.sort();
    factors.dedup();

    factors
}
