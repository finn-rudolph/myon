use core::mem::swap;
use std::ops::{Index, IndexMut};

use rug::{ops::Pow, Integer};

use crate::{
    nt,
    polynomial::{self, Polynomial},
};

const MAX_DEGREE: usize = polynomial::MAX_DEGREE;

#[derive(Clone)]
pub struct GfPolynomial {
    coefficients: [u32; MAX_DEGREE + 1],
    modulus: u32,
}

impl GfPolynomial {
    fn new(modulus: u32) -> GfPolynomial {
        GfPolynomial {
            coefficients: [0; MAX_DEGREE + 1],
            modulus,
        }
    }

    pub fn from_polynomial(f: &Polynomial, modulus: u32) -> GfPolynomial {
        let mut g = GfPolynomial::new(modulus);
        for i in 0..=MAX_DEGREE {
            g[i] = f[i].mod_u(modulus);
        }
        g
    }

    pub fn degree(&self) -> usize {
        let mut d = MAX_DEGREE;
        while self[d] == 0 {
            d -= 1;
        }
        d
    }

    // Same routine as in the general polynomial case.
    pub fn mul_mod(&self, f: &GfPolynomial, g: &GfPolynomial) -> GfPolynomial {
        let d = self.degree();
        let p = self.modulus;

        let mut result = GfPolynomial::new(self.modulus);
        for i in 0..d {
            result[i] = ((g[i] as u64 * f[d - 1] as u64) % p as u64) as u32;
        }

        for i in (0..d - 1).rev() {
            let mut leading_coefficient = 0u32;
            for coefficient in &mut result.coefficients {
                swap(&mut leading_coefficient, coefficient);
            }

            for (j, coefficient) in result.coefficients.iter_mut().enumerate() {
                *coefficient +=
                    p - ((leading_coefficient as u64 * self[j] as u64) % p as u64) as u32;
            }

            for j in 0..d {
                result[j] = ((result[j] as u64 + g[j] as u64 * f[i] as u64) % p as u64) as u32;
            }
        }

        result
    }

    fn rem(mut self, modulus: &GfPolynomial) -> GfPolynomial {
        let (d, e) = (self.degree(), modulus.degree());
        let p = self.modulus;

        for i in (e..=d).rev() {
            let quotient = self[i] * nt::mod_inv(self[i], p);
            for j in 0..=e {
                self[j + i - e] = (self[j + i - e] + p
                    - ((quotient as u64 * modulus[j] as u64) % p as u64) as u32)
                    % p;
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
        let p = self.modulus;

        let mut prime_divisors: Vec<u32> = Vec::new();
        for q in 2..=d {
            if d % q == 0 && nt::miller_rabin(q as u32) {
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

        let mut h = self.pow_mod(x.clone(), Integer::from(p).pow(d as u32));
        h[1] = (h[1] + p - 1) % p;

        h.degree() == 0 && h[0] == 0
    }
}

impl Index<usize> for GfPolynomial {
    type Output = u32;

    fn index(&self, i: usize) -> &u32 {
        &self.coefficients[i]
    }
}

impl IndexMut<usize> for GfPolynomial {
    fn index_mut(&mut self, i: usize) -> &mut u32 {
        &mut self.coefficients[i]
    }
}

pub struct GfMpPolynomial {
    coefficients: [Integer; MAX_DEGREE + 1],
    modulus: Integer,
}

impl GfMpPolynomial {
    fn new(modulus: Integer) -> GfMpPolynomial {
        GfMpPolynomial {
            coefficients: [Integer::ZERO; MAX_DEGREE + 1],
            modulus,
        }
    }

    pub fn degree(&self) -> usize {
        let mut d = MAX_DEGREE;
        while self[d] == 0 {
            d -= 1;
        }
        d
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
