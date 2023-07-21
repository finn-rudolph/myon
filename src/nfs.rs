use rug::{ops::Pow, Integer};

use crate::nt;

#[derive(Clone, Copy)]
struct Params {
    rational_base_size: usize,
    algebraic_base_size: usize,
    polynomial_degree: u32,
}

impl Params {
    const PARAM_TABLE: [(u32, Params); 1] = [(
        128,
        Params {
            rational_base_size: 500,
            algebraic_base_size: 500,
            polynomial_degree: 3,
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

fn rational_factor_base(params: &Params) -> Vec<u32> {
    let mut base: Vec<u32> = Vec::new();

    let mut i: u32 = 2;
    while base.len() < params.rational_base_size {
        if nt::is_prime(i) {
            base.push(i);
        }
        i += 1;
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

pub fn factorize(r: u32, e: u32, s: i32) -> String {
    let n: Integer = Integer::from(r).pow(e) - s;

    let params = Params::new(&n);
    let (t, m) = select_polynomial(r, e, s, &params);
    let rational_base = rational_factor_base(&params);
    let algebraic_base = algebraic_factor_base(t, &params);

    todo!()
}
