use rand_xoshiro::{rand_core::RngCore, Xoshiro256PlusPlus};

// TODO: make this module generic

struct Montgomery {
    n: u32,
    r: u32,
}

impl Montgomery {
    const fn new(n: u32) -> Montgomery {
        let mut mtg = Montgomery { n, r: 1 };
        let mut i = 0;
        while i < 5 {
            mtg.r *= 2 - n * mtg.r;
            i += 1;
        }
        mtg
    }

    const fn reduce(&self, a: u64) -> u32 {
        let q: u64 = (a as u32).wrapping_mul(self.r) as u64;
        let m = q * (self.n as u64) >> 32;
        let x: u32 = ((a >> 32) + self.n as u64 - m) as u32;
        if x >= self.n {
            x - self.n
        } else {
            x
        }
    }

    const fn multiply(&self, a: u32, b: u32) -> u32 {
        self.reduce(a as u64 * b as u64)
    }

    const fn transform(&self, a: u32) -> u32 {
        (((a as u64) << 32) % (self.n as u64)) as u32
    }
}

const fn mod_exp(mut a: u64, mut b: u64, m: u64) -> u64 {
    let mut c: u64 = 1;

    while b != 0 {
        if b & 1 == 1 {
            c = (c * a) % m;
        }
        a = (a * a) % m;
        b >>= 1;
    }

    c
}

const fn legendre(a: u64, p: u64) -> u64 {
    mod_exp(a, (p - 1) >> 1, p)
}

const fn mod_inverse(a: u64, m: u64) -> u64 {
    mod_exp(a, m - 2, m)
}

// Finds a square root of a modulo p using the Tonelli-Shanks algorithm. a and p may not be greater
// than u32::MAX, since multiplication is performed with them.
pub fn tonelli_shanks(mut a: u64, p: u64, xo: &mut Xoshiro256PlusPlus) -> u64 {
    assert!(a <= u32::MAX as u64);
    assert!(p <= u32::MAX as u64);

    if p & 3 == 3 {
        return mod_exp(a, (p + 1) >> 2, p);
    }

    // About 2 iterations are expected.
    let mut b = xo.next_u64() % p;
    while b == 0 || legendre(b, p) == 1 {
        b = xo.next_u64() % p;
    }

    // Loop invariant: c = b ^ (2 ^ (k - 2)). Before the loop, k = 2, which is possible since p = 1
    // mod 4, so m = (p - 1) / 2^k is an integer.
    let mut m = (p - 1) >> 2;
    let mut correction: u64 = 1;
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

// TODO: Add Cipolla's algorithm (it shall be faster sometimes?)
