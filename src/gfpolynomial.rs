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

    // Multiplies f and g mod self. "mod" means here, whenever we encounter a power of x greater or
    // equal to the degree of self, we subtract some multiple of self to make this power disappear.
    pub fn mul_mod(&self, f: &GfPolynomial, g: &GfPolynomial) -> GfPolynomial {
        let d = self.degree();

        // Initialize result with the leading coefficient of self times rhs. This means, the powers
        // of alpha present in the result are actually shifted up by d - 1.
        let mut result = GfPolynomial::new(self.modulus);
        for i in 0..d {
            result[i] = (&g[i] * &f[d - 1]).into();
        }

        // In each iteration, the leading coffiecient is converted to lower order terms, and the
        // i-th coefficient of self times rhs is added. Thus, the shift in the exponents of powers
        // of alpha is reduced by 1.
        for i in (0..d - 1).rev() {
            let mut leading_coefficient = 0u32;
            for coefficient in &mut result.coefficients {
                swap(&mut leading_coefficient, coefficient);
            }

            for (j, coefficient) in result.coefficients.iter_mut().enumerate() {
                *coefficient -= &leading_coefficient * &self[j];
            }

            for j in 0..d {
                result[j] += &g[j] * &f[i];
            }
        }

        result
    }

    fn rem(&self, modulus: &GfPolynomial) -> GfPolynomial {
        todo!()
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

    pub fn is_irreducible(&self) -> bool {
        let d = self.degree();

        let mut prime_divisors: Vec<u32> = Vec::new();
        for q in 2..=d {
            if d % q == 0 && nt::miller_rabin(q as u32) {
                prime_divisors.push(q as u32);
            }
        }

        let mut x_polynomial = GfPolynomial::new(self.modulus);
        x_polynomial[1] = 1;

        for q in prime_divisors {
            let h = self.pow_mod(
                x_polynomial.clone(),
                Integer::from(self.modulus).pow(d as u32 / q),
            );
        }

        todo!()
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
