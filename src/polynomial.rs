use std::{
    cmp::PartialOrd,
    mem::swap,
    ops::{Index, IndexMut, MulAssign, Rem, Sub},
};

use rug::{ops::Pow, Complete, Integer};

use crate::{nt, params::Params};

pub const MAX_DEGREE: usize = Params::PARAM_TABLE[Params::PARAM_TABLE.len() - 1]
    .1
    .polynomial_degree as usize;

#[derive(Clone)]
pub struct Polynomial([Integer; MAX_DEGREE + 1]);

impl Polynomial {
    fn new() -> Polynomial {
        Polynomial([Integer::ZERO; MAX_DEGREE + 1])
    }

    pub fn degree(&self) -> usize {
        let mut d = MAX_DEGREE;
        while self[d] == 0 {
            d -= 1;
        }
        d
    }

    pub fn evaluate<T: Copy>(&self, x: T) -> Integer
    where
        Integer: MulAssign<T>,
    {
        let mut z = Integer::from(1);
        let mut y = self[0].clone();
        for coefficient in self.0.iter().skip(1) {
            y += coefficient * &z;
            z *= x;
        }
        y
    }

    pub fn derivative(&self) -> Polynomial {
        let mut f = Polynomial::new();
        for (i, coefficient) in self.0.iter().skip(1).enumerate() {
            f[i - 1] = (coefficient * i).complete();
        }
        f
    }

    pub fn find_roots_mod_p<T>(&self, p: T) -> Vec<T>
    where
        Integer: Rem<T> + MulAssign<T>,
        <Integer as Rem<T>>::Output: PartialEq<i32>,
        T: Sub<T, Output = T> + PartialOrd + Copy,
    {
        let mut roots: Vec<T> = Vec::new();

        let i: T = p - p;
        while i < p {
            if self.evaluate(i) % p == 0 {
                roots.push(i);
            }
        }

        roots
    }

    // Multiplies f and g mod self. "mod" means here, whenever we encounter a power of x greater or
    // equal to the degree of self, we subtract some multiple of self to make this power disappear.
    pub fn mul_mod(&self, f: &Polynomial, g: &Polynomial) -> Polynomial {
        let d = self.degree();

        // Initialize result with the leading coefficient of self times rhs. This means, the powers
        // of alpha present in the result are actually shifted up by d - 1.
        let mut result = Polynomial::new();
        for i in 0..d {
            result[i] = (&g[i] * &f[d - 1]).into();
        }

        // In each iteration, the leading coffiecient is converted to lower order terms, and the
        // i-th coefficient of self times rhs is added. Thus, the shift in the exponents of powers
        // of alpha is reduced by 1.
        for i in (0..d - 1).rev() {
            let mut leading_coefficient = Integer::new();
            for coefficient in &mut result.0 {
                swap(&mut leading_coefficient, coefficient);
            }

            for (j, coefficient) in result.0.iter_mut().enumerate() {
                *coefficient -= &leading_coefficient * &self[j];
            }

            for j in 0..d {
                result[j] += &g[j] * &f[i];
            }
        }

        result
    }

    fn pow_mod(&self, mut f: Polynomial, mut e: Integer) -> Polynomial {
        let mut g = Polynomial::new();
        g[0] = Integer::from(1);

        while e != 0 {
            if e.is_odd() {
                g = self.mul_mod(&g, &f);
            }
            f = self.mul_mod(&f, &f);
            e >>= 1;
        }

        g
    }

    fn rem(&self, modulus: &Polynomial) -> Polynomial {
        todo!()
    }

    fn gcd(self, f: Polynomial) -> Polynomial {
        if f.degree() == 0 && f[0] == 0 {
            return self;
        }
        let r = self.rem(&f);
        f.gcd(r)
    }

    pub fn is_irreducible_mod_p(&self, p: u32) -> bool {
        let d = self.degree();

        let mut prime_divisors: Vec<u32> = Vec::new();
        for q in 2..=d {
            if d % q == 0 && nt::miller_rabin(q as u32) {
                prime_divisors.push(q as u32);
            }
        }

        let mut x_polynomial = Polynomial::new();
        x_polynomial[1] = Integer::from(1);

        for q in prime_divisors {
            let h = self.pow_mod(x_polynomial.clone(), Integer::from(p).pow(d as u32 / q));
        }

        todo!()
    }
}

impl Index<usize> for Polynomial {
    type Output = Integer;

    fn index(&self, i: usize) -> &Integer {
        &self.0[i]
    }
}

impl IndexMut<usize> for Polynomial {
    fn index_mut(&mut self, i: usize) -> &mut Integer {
        &mut self.0[i]
    }
}

// Polynomial selection for numbers of the form r^e - s, with small r and |s|. Returns the selected
// polynomial and an integer m, such that f(m) = 0 mod n.
pub fn select_special(r: u32, e: u32, s: i32, params: &Params) -> (Polynomial, Integer) {
    let d = params.polynomial_degree;
    let k = (e + d - 1) / d;

    let mut f = Polynomial::new();
    f[d as usize] = Integer::from(1);
    f[0] = Integer::from(s * r.pow(k * d - e) as i32);
    assert_eq!(f.degree(), params.polynomial_degree as usize);

    (f, Integer::from(r).pow(k))
}
