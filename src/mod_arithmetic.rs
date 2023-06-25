use rand_xoshiro::{rand_core::RngCore, Xoshiro256PlusPlus};

fn mod_exp(mut a: u32, mut b: u32, m: u32) -> u32 {
    let mut c: u64 = 1;

    while b != 0 {
        if b & 1 == 1 {
            c = (c * a as u64) % m as u64;
        }
        a = (a * a) % m;
        b >>= 1;
    }

    c as u32
}

fn legendre(a: u32, p: u32) -> u32 {
    mod_exp(a, (p - 1) >> 1, p)
}

fn mod_inverse(a: u32, m: u32) -> u32 {
    mod_exp(a, m - 2, m)
}

// Finds a square root of a modulo p using the Tonelli-Shanks algorithm.
pub fn mod_sqrt(mut a: u32, p: u32, xo: &mut Xoshiro256PlusPlus) -> u32 {
    if p & 3 == 3 {
        return mod_exp(a, (p + 1) >> 2, p);
    }

    // About 2 iterations are expected.
    let mut b = xo.next_u32() % p;
    while b == 0 || legendre(b, p) == 1 {
        b = xo.next_u32() % p;
    }

    // Loop invariant: c = b ^ (2 ^ (k - 2)). Before the loop, k = 2, which is possible since p = 1
    // mod 4, so m = (p - 1) / 2^k is an integer.
    let mut m = (p - 1) >> 2;
    let mut correction = 1;
    let mut c = b;
    let mut cinv = mod_inverse(b, p);

    loop {
        if mod_exp(a, m, p) != 1 {
            a = (a * ((c * c) % p)) % p;
            correction = (correction * cinv) % p;
        }
        if m & 1 == 1 {
            break;
        }
        m >>= 1;
        c = (c * c) % p;
        cinv = (cinv * cinv) % p;
    }

    (mod_exp(a, (m + 1) >> 1, p) * correction) % p
}
