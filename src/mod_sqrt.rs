use rand_xoshiro::{rand_core::RngCore, Xoshiro256PlusPlus};

// TODO: make this module generic

pub trait ModReduction {
    fn new(n: u32) -> Self;
    fn rem(&self, a: u64) -> u32;
}

struct LemireReduction {
    n: u32, // divisor / modulus
    c: u64, // 2^64 / n + 1
}

// Can quickly compute a / n, a % m and test whether a is divisible by n, for arbitrary a and fixed
// (but not compile-time known) n. Described in Lemire, D, Kaser, O. & Kurz, N. (2019).
impl ModReduction for LemireReduction {
    fn new(n: u32) -> LemireReduction {
        assert!(n != 0);
        LemireReduction {
            n,
            c: u64::MAX / (n as u64) + 1,
        }
    }

    fn rem(&self, a: u64) -> u32 {
        let lowbits = self.c.wrapping_mul(a as u64);
        ((lowbits as u128 * self.n as u128) >> 64) as u32
    }
}

impl LemireReduction {
    fn div(&self, a: u32) -> u32 {
        ((self.c as u128 * a as u128) >> 64) as u32
    }

    fn is_divisible(&self, a: u32) -> bool {
        self.c.wrapping_mul(a as u64) < self.c
    }
}

struct NativeReduction {
    n: u32,
}

impl ModReduction for NativeReduction {
    fn new(n: u32) -> NativeReduction {
        NativeReduction { n }
    }

    fn rem(&self, a: u64) -> u32 {
        (a % self.n as u64) as u32
    }
}

fn mod_exp<T: ModReduction>(mut a: u64, mut b: u64, reduction: &T) -> u64 {
    let mut c: u64 = 1;

    while b != 0 {
        if b & 1 == 1 {
            c = reduction.rem(c * a) as u64;
        }
        a = reduction.rem(a * a) as u64;
        b >>= 1;
    }

    c
}

fn legendre<T: ModReduction>(a: u64, p: u64, lemire: &T) -> u64 {
    mod_exp(a, (p - 1) >> 1, lemire)
}

fn mod_inverse<T: ModReduction>(a: u64, m: u64, lemire: &T) -> u64 {
    mod_exp(a, m - 2, lemire)
}

// Finds a square root of a modulo p using the Tonelli-Shanks algorithm. a and p may not be greater
// than u32::MAX, since multiplication is performed with them.
pub fn mod_sqrt<T: ModReduction>(mut a: u64, p: u64, xo: &mut Xoshiro256PlusPlus) -> u64 {
    assert!(a <= u32::MAX as u64);
    assert!(p <= u32::MAX as u64);

    let reduction = T::new(p as u32);

    if p & 3 == 3 {
        return mod_exp(a, (p + 1) >> 2, &reduction);
    }

    // About 2 iterations are expected.
    let mut b = xo.next_u64() % p;
    while b == 0 || legendre(b, p, &reduction) == 1 {
        b = xo.next_u64() % p;
    }

    // Loop invariant: c = b ^ (2 ^ (k - 2)). Before the loop, k = 2, which is possible since p = 1
    // mod 4, so m = (p - 1) / 2^k is an integer.
    let mut m = (p - 1) >> 2;
    let mut correction: u64 = 1;
    let mut c = b;
    let mut cinv = mod_inverse(b, p, &reduction);

    loop {
        if mod_exp(a, m, &reduction) != 1 {
            a = reduction.rem(a * reduction.rem(c * c) as u64) as u64;
            correction = (correction * cinv) % p;
        }
        if m & 1 == 1 {
            break;
        }
        m >>= 1;
        c = (c * c) % p;
        cinv = (cinv * cinv) % p;
    }

    reduction.rem(mod_exp(a, (m + 1) >> 1, &reduction) * correction) as u64
}

// TODO: Add Cipolla's algorithm (it shall be faster sometimes?)

#[cfg(test)]
mod tests {
    use super::*;
    use rand_xoshiro::{rand_core::SeedableRng, Xoshiro256PlusPlus};

    // Returns true, if (and only if? I'm not sure.) n is a prime. Works for numbers less than
    // u32::MAX.
    fn is_prime(n: u64) -> bool {
        assert!(n <= u32::MAX as u64);

        const MILLER_RABIN_BASES: [u64; 3] = [15, 7363882082, 992620450144556];

        let trailing_zeros = (n - 1).trailing_zeros();
        let u = (n - 1) >> trailing_zeros;
        let lemire = LemireReduction::new(n as u32);

        for mut a in MILLER_RABIN_BASES {
            a = a % n;
            let mut x = mod_exp(a, u, &lemire);
            for _ in 0..trailing_zeros {
                let y = lemire.rem(x * x) as u64;
                if y == 1 && x != 1 && x != n - 1 {
                    return false;
                }
                x = y;
            }
            if x != 1 {
                return false;
            }
        }
        true
    }

    // TODO: find better algorithm
    fn gen_prime(xo: &mut Xoshiro256PlusPlus) -> u32 {
        loop {
            let p = xo.next_u32();
            if is_prime(p as u64) {
                return p;
            }
        }
    }

    #[test]
    fn test_lemire_reduction() {
        let mut xo = Xoshiro256PlusPlus::seed_from_u64(1000000009);

        for _ in 0..1000 {
            let n = xo.next_u32();
            let lemire = LemireReduction::new(n);
            for _ in 0..1000 {
                let a = xo.next_u32();
                assert_eq!(lemire.div(a), a / n);
                assert_eq!(lemire.rem(a as u64), a % n);
                assert_eq!(lemire.is_divisible(a), a % n == 0);
            }
        }
    }

    #[test]
    fn test_tonelli_shanks() {
        let mut xo = Xoshiro256PlusPlus::seed_from_u64(1000000009);

        for _ in 0..1000 {
            let p = gen_prime(&mut xo);
            let lemire = LemireReduction::new(p);
            let mut a = xo.next_u32() % p;
            while legendre(a as u64, p as u64, &lemire) != 1 {
                a = xo.next_u32() % p;
            }
            let x = mod_sqrt::<NativeReduction>(a as u64, p as u64, &mut xo);
            assert_eq!((x * x) % p as u64, a as u64);
        }
    }
}
