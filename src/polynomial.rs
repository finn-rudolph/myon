use core::mem::swap;
use std::{
    cmp::PartialOrd,
    ops::{Index, IndexMut, MulAssign, Rem, Sub},
};

use rug::{ops::Pow, Complete, Integer};

use crate::{
    gfpolynomial::GfMpPolynomial,
    params::{Params, MAX_DEGREE},
};

pub trait Polynomial<T>
where
    Self: Index<usize> + IndexMut<usize>,
{
    fn degree(&self) -> usize;

    fn coefficients(self) -> [T; MAX_DEGREE + 1];

    fn coefficients_ref(&self) -> &[T; MAX_DEGREE + 1];

    fn coefficients_mut(&mut self) -> &mut [T; MAX_DEGREE + 1];
}

#[derive(Clone)]
pub struct MpPolynomial([Integer; MAX_DEGREE + 1]);

impl MpPolynomial {
    fn new() -> MpPolynomial {
        MpPolynomial([Integer::ZERO; MAX_DEGREE + 1])
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

    pub fn derivative(&self) -> MpPolynomial {
        let mut f = MpPolynomial::new();
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
    pub fn mul_mod(&self, f: &MpPolynomial, g: &MpPolynomial) -> MpPolynomial {
        let d = self.degree();

        // Initialize result with the leading coefficient of self times rhs. This means, the powers
        // of alpha present in the result are actually shifted up by d - 1.
        let mut result = MpPolynomial::new();
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
}

impl Polynomial<Integer> for MpPolynomial {
    fn degree(&self) -> usize {
        let mut d = MAX_DEGREE;
        while self[d] == 0 {
            d -= 1;
        }
        d
    }

    fn coefficients(self) -> [Integer; MAX_DEGREE + 1] {
        self.0
    }

    fn coefficients_ref(&self) -> &[Integer; MAX_DEGREE + 1] {
        &self.0
    }

    fn coefficients_mut(&mut self) -> &mut [Integer; MAX_DEGREE + 1] {
        &mut self.0
    }
}

impl Index<usize> for MpPolynomial {
    type Output = Integer;

    fn index(&self, i: usize) -> &Integer {
        &self.0[i]
    }
}

impl IndexMut<usize> for MpPolynomial {
    fn index_mut(&mut self, i: usize) -> &mut Integer {
        &mut self.0[i]
    }
}

impl From<GfMpPolynomial> for MpPolynomial {
    fn from(f: GfMpPolynomial) -> MpPolynomial {
        MpPolynomial(f.coefficients())
    }
}

// Polynomial selection for numbers of the form r^e - s, with small r and |s|. Returns the selected
// polynomial and an integer m, such that f(m) = 0 mod n.
pub fn select_special(r: u32, e: u32, s: i32, params: &Params) -> (MpPolynomial, Integer) {
    let d = params.polynomial_degree;
    let k = (e + d as u32 - 1) / d as u32;

    let mut f = MpPolynomial::new();
    f[d as usize] = Integer::from(1);
    f[0] = Integer::from(s * r.pow(k * d as u32 - e) as i32);
    assert_eq!(f.degree(), params.polynomial_degree);

    (f, Integer::from(r).pow(k))
}
