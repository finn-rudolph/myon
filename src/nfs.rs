use rug::{ops::Pow, Integer};

use crate::nt;

mod params {
    use rug::Integer;

    static mut RATIONAL_BASE_SIZE: usize = 0;
    static mut ALGEBRAIC_BASE_SIZE: usize = 0;
    static mut POLYNOMIAL_DEGREE: u32 = 0;

    pub fn rational_base_size() -> usize {
        unsafe { RATIONAL_BASE_SIZE }
    }

    pub fn algebraic_base_size() -> usize {
        unsafe { ALGEBRAIC_BASE_SIZE }
    }

    pub fn polynomial_degree() -> u32 {
        unsafe { POLYNOMIAL_DEGREE }
    }

    pub fn set_params(n: &Integer) {}
}

fn rational_factor_base(size: usize) -> Vec<u32> {
    let mut base: Vec<u32> = Vec::new();

    let mut i: u32 = 2;
    while base.len() < size {
        if nt::is_prime(i) {
            base.push(i);
        }
        i += 1;
    }

    base
}

// Returns two integers t and m, the selected polynomial is x^d - t. f(m) = 0 mod n.
fn select_polynomial(r: u32, e: u32, s: i32) -> (i32, Integer) {
    let d = params::polynomial_degree();
    let k = (e + d - 1) / d;
    (s * r.pow(k * d - e) as i32, Integer::from(r).pow(k))
}

fn find_polynomial_roots(t: i32, p: u32) -> Vec<u32> {
    let d = params::rational_base_size();
    let mut roots: Vec<u32> = Vec::new();

    for i in 1..p {
        if i.pow(d as u32) as i32 - t == 0 {
            roots.push(i);
        }
    }

    roots
}

fn algebraic_factor_base(size: usize, t: i32) -> Vec<(u32, u32)> {
    let mut base: Vec<(u32, u32)> = Vec::new();
    let mut p: u32 = 3;

    while base.len() < size {
        if nt::is_prime(p) {
            let roots = find_polynomial_roots(t, p);
            base.extend(roots.iter().map(|r| (p, *r)));
        }
        p += 2;
    }

    base
}

pub fn factorize(r: u32, e: u32, s: i32) -> String {
    const x: u64 = 1 << 43;
    let rational_base = rational_factor_base(500);

    todo!()
}
