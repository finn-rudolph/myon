use core::mem::swap;
use std::ops::{Index, IndexMut};

use rug::{integer::IntegerExt64, ops::Pow, Complete, Integer};

use crate::{
    nt,
    params::MAX_DEGREE,
    polynomial::{MpPolynomial, Polynomial},
};

#[derive(Clone)]
pub struct GfPolynomial {
    coefficients: [u64; MAX_DEGREE + 1],
    modulus: u64,
}

impl GfPolynomial {
    pub fn new(modulus: u64) -> GfPolynomial {
        GfPolynomial {
            coefficients: [0; MAX_DEGREE + 1],
            modulus,
        }
    }

    pub fn from_mp_polynomial(f: &MpPolynomial, modulus: u64) -> GfPolynomial {
        let mut g = GfPolynomial::new(modulus);
        for i in 0..=MAX_DEGREE {
            g[i] = f[i].mod_u64(modulus);
        }
        g
    }

    pub fn modulus(&self) -> u64 {
        self.modulus
    }

    pub fn add(mut self, rhs: &GfPolynomial) -> GfPolynomial {
        for (i, coefficient) in self.coefficients.iter_mut().enumerate() {
            *coefficient = (*coefficient + rhs[i]) % self.modulus;
        }
        self
    }

    // Same routine as in the general polynomial case.
    pub fn mul_mod(&self, f: &GfPolynomial, g: &GfPolynomial) -> GfPolynomial {
        let d = self.degree();
        let p = self.modulus;

        let mut result = GfPolynomial::new(self.modulus);
        for i in 0..d {
            result[i] = (g[i] * f[d - 1]) % p;
        }

        for i in (0..d - 1).rev() {
            let mut leading_coefficient = 0u64;
            for coefficient in result.coefficients_mut().iter_mut().take(d) {
                swap(&mut leading_coefficient, coefficient);
            }

            for j in 0..d {
                result[j] += g[j] * f[i];
                result[j] += p - (leading_coefficient * self[j]) % p;
                result[j] %= p;
            }
        }

        result
    }

    fn rem(mut self, modulus: &GfPolynomial) -> GfPolynomial {
        let (d, e) = (self.degree(), modulus.degree());
        let p = self.modulus;
        let modulus_leading_inv = nt::mod_inv(modulus[e], p);

        for i in (e..=d).rev() {
            let quotient = (self[i] * modulus_leading_inv) % p;
            for j in 0..=e {
                self[j + i - e] = (self[j + i - e] + p - (quotient * modulus[j]) % p) % p;
            }
        }

        self
    }

    fn gcd(self, f: GfPolynomial) -> GfPolynomial {
        if f.degree() == 0 && f[0] == 0 {
            return self;
        }
        let r = self.rem(&f);
        f.gcd(r)
    }

    fn pow_mod(&self, mut f: GfPolynomial, mut e: Integer) -> GfPolynomial {
        let mut g = GfPolynomial::new(self.modulus);
        g[0] = 1;

        while e != 0 {
            if e.is_odd() {
                g = self.mul_mod(&g, &f);
            }
            f = self.mul_mod(&f, &f);
            e >>= 1;
        }

        g
    }

    // Rabin's test of irreducibility for polynomials over finite fields.
    pub fn is_irreducible(&self) -> bool {
        let d = self.degree();
        let p = self.modulus();

        let mut prime_divisors: Vec<u32> = Vec::new();
        for q in 2..=d {
            if d % q == 0 && nt::miller_rabin(q as u64) {
                prime_divisors.push(q as u32);
            }
        }

        let mut x = GfPolynomial::new(p);
        x[1] = 1;

        for q in prime_divisors {
            let mut h = self.pow_mod(x.clone(), Integer::from(p).pow(d as u32 / q));
            h[1] = (h[1] + p - 1) % p; // subtract x
            let g = h.gcd(self.clone());
            if g.degree() != 0 {
                return false;
            }
        }

        let mut g = self.pow_mod(x.clone(), Integer::from(p).pow(d as u32));
        g[1] = (g[1] + p - 1) % p;

        g.degree() == 0 && g[0] == 0
    }
}

impl Polynomial<u64> for GfPolynomial {
    fn degree(&self) -> usize {
        let mut d = MAX_DEGREE;
        while d > 0 && self[d] == 0 {
            d -= 1;
        }
        d
    }

    fn coefficients(self) -> [u64; MAX_DEGREE + 1] {
        self.coefficients
    }

    fn coefficients_ref(&self) -> &[u64; MAX_DEGREE + 1] {
        &self.coefficients
    }

    fn coefficients_mut(&mut self) -> &mut [u64; MAX_DEGREE + 1] {
        &mut self.coefficients
    }
}

impl Index<usize> for GfPolynomial {
    type Output = u64;

    fn index(&self, i: usize) -> &u64 {
        &self.coefficients[i]
    }
}

impl IndexMut<usize> for GfPolynomial {
    fn index_mut(&mut self, i: usize) -> &mut u64 {
        &mut self.coefficients[i]
    }
}

pub struct GfMpPolynomial {
    coefficients: [Integer; MAX_DEGREE + 1],
    modulus: Integer,
}

impl GfMpPolynomial {
    pub fn new(modulus: Integer) -> GfMpPolynomial {
        GfMpPolynomial {
            coefficients: [Integer::ZERO; MAX_DEGREE + 1],
            modulus,
        }
    }

    pub fn from_mp_polynomial(f: &MpPolynomial, modulus: Integer) -> GfMpPolynomial {
        let mut g = GfMpPolynomial::new(modulus);
        for (i, coefficient) in f.coefficients_ref().iter().enumerate() {
            g.coefficients[i] = Integer::from(coefficient) % &g.modulus;
        }
        g
    }

    pub fn modulus(&self) -> &Integer {
        &self.modulus
    }

    // Same routine as in the general polynomial case.
    pub fn mul_mod(&self, f: &GfMpPolynomial, g: &GfMpPolynomial) -> GfMpPolynomial {
        let d = self.degree();
        let p = self.modulus();

        let mut result = GfMpPolynomial::new(self.modulus().clone());
        for i in 0..d {
            result[i] = (&g[i] * &f[d - 1]).complete() % p;
        }

        for i in (0..d - 1).rev() {
            let mut leading_coefficient = Integer::new();
            for coefficient in result.coefficients_mut().iter_mut().take(d) {
                swap(&mut leading_coefficient, coefficient);
            }

            for j in 0..d {
                result[j] += &g[j] * &f[i];
                result[j] += p - (&leading_coefficient * &self[j]).complete() % p;
                result[j] %= p;
            }
        }

        result
    }
}

impl Polynomial<Integer> for GfMpPolynomial {
    fn degree(&self) -> usize {
        let mut d = MAX_DEGREE;
        while d > 0 && self[d] == 0 {
            d -= 1;
        }
        d
    }

    fn coefficients(self) -> [Integer; MAX_DEGREE + 1] {
        self.coefficients
    }

    fn coefficients_ref(&self) -> &[Integer; MAX_DEGREE + 1] {
        &self.coefficients
    }

    fn coefficients_mut(&mut self) -> &mut [Integer; MAX_DEGREE + 1] {
        &mut self.coefficients
    }
}

impl Index<usize> for GfMpPolynomial {
    type Output = Integer;

    fn index(&self, i: usize) -> &Integer {
        &self.coefficients[i]
    }
}

impl IndexMut<usize> for GfMpPolynomial {
    fn index_mut(&mut self, i: usize) -> &mut Integer {
        &mut self.coefficients[i]
    }
}

impl From<&GfPolynomial> for GfMpPolynomial {
    fn from(f: &GfPolynomial) -> GfMpPolynomial {
        let mut g = GfMpPolynomial::new(Integer::from(f.modulus()));
        for (i, coefficient) in f.coefficients_ref().iter().enumerate() {
            g.coefficients[i] = Integer::from(*coefficient);
        }
        g
    }
}
