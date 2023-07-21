use std::cmp::max;

use rug::{ops::Pow, Complete, Integer};

use crate::nt;

#[derive(Clone, Copy)]
struct Params {
    rational_base_size: usize,
    algebraic_base_size: usize,
    polynomial_degree: u32,
    sieve_array_size: usize,
    fudge: u8, // different fudge for rational and algebraic side?
}

impl Params {
    const PARAM_TABLE: [(u32, Params); 1] = [(
        128,
        Params {
            rational_base_size: 500,
            algebraic_base_size: 500,
            polynomial_degree: 3,
            sieve_array_size: 10000,
            fudge: 20,
        },
    )];

    fn new(n: &Integer) -> Params {
        let bits = n.significant_bits();

        for (bits_lim, params) in Params::PARAM_TABLE.iter().rev() {
            if bits <= *bits_lim {
                return *params;
            }
        }

        Params::PARAM_TABLE.last().unwrap().1
    }
}

fn rational_factor_base(m: &Integer, params: &Params) -> Vec<(u32, u32)> {
    let mut base: Vec<(u32, u32)> = Vec::new();

    let mut p: u32 = 2;
    while base.len() < params.rational_base_size {
        if nt::is_prime(p) {
            base.push((p, m.mod_u(p)));
        }
        p += 1;
    }

    base
}

// Returns two integers t and m, the selected polynomial is x^d - t. f(m) = 0 mod n.
fn select_polynomial(r: u32, e: u32, s: i32, params: &Params) -> (i32, Integer) {
    let d = params.polynomial_degree;
    let k = (e + d - 1) / d;
    (s * r.pow(k * d - e) as i32, Integer::from(r).pow(k))
}

fn find_polynomial_roots(t: i32, p: u32, params: &Params) -> Vec<u32> {
    let d = params.rational_base_size;
    let mut roots: Vec<u32> = Vec::new();

    for i in 1..p {
        if (i.pow(d as u32) as i32 - t) % p as i32 == 0 {
            roots.push(i);
        }
    }

    roots
}

fn algebraic_factor_base(t: i32, params: &Params) -> Vec<(u32, u32)> {
    let mut base: Vec<(u32, u32)> = Vec::new();
    let mut p: u32 = 2;

    while base.len() < params.algebraic_base_size {
        if nt::is_prime(p) {
            let roots = find_polynomial_roots(t, p, params);
            base.extend(roots.iter().map(|r| (p, *r)));
        }
        p += 1;
    }

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

pub fn factorize(r: u32, e: u32, s: i32) -> Integer {
    let n: Integer = Integer::from(r).pow(e) - s;

    let params = Params::new(&n);
    let d = params.polynomial_degree;
    let (t, m) = select_polynomial(r, e, s, &params);
    let rational_base = rational_factor_base(&m, &params);
    let algebraic_base = algebraic_factor_base(t, &params);

    let mut rational_sieve_array: Vec<i8> = vec![0; params.sieve_array_size];
    let mut algebraic_sieve_array: Vec<i8> = vec![0; params.sieve_array_size];

    for b in 1.. {
        rational_sieve_array.fill(-((ilog2_rounded(b) + m.significant_bits()) as i8));
        line_sieve(b, &mut rational_sieve_array, &rational_base);

        let log_tbd = (d * ilog2_rounded(b) + ilog2_rounded(t.abs() as u32)) as i8;
        let a0 = -(params.sieve_array_size as i32 / 2);
        for i in 0..params.sieve_array_size {
            algebraic_sieve_array[i] = max((d as i32 * (a0 + i as i32)) as i8, log_tbd);
        }
        line_sieve(b, &mut algebraic_sieve_array, &algebraic_base);

        // Consider unsafe access here to avoid bounds checks.
        for i in 0..params.sieve_array_size {
            if rational_sieve_array[i] >= 0 && algebraic_sieve_array[i] >= 0 {
                let a = a0 + i as i32;

                // Trial divide on the rational side.
                let mut x = a + (b * &m).complete();
                for (p, _) in &rational_base {
                    let e = x.remove_factor_mut(&Integer::from(*p));
                }

                // Trial divide on the algrbraic side.
                let mut y = Integer::from(a).pow(d) - t * Integer::from(-(b as i32)).pow(d);
                for (p, _) in &algebraic_base {
                    let e = y.remove_factor_mut(&Integer::from(*p));
                }

                if x == 1 && y == 1 {
                    // smooth pair (a, b) found!
                }
            }
        }
    }

    todo!()
}
