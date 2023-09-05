use core::mem::swap;
use std::{
    fmt::Display,
    ops::{Index, IndexMut, MulAssign},
};

use rug::{Complete, Integer};

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
    pub fn new() -> MpPolynomial {
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
            f[i] = (coefficient * (i + 1)).complete();
        }
        f
    }

    pub fn find_roots_mod_p(&self, p: u64) -> Vec<u64> {
        let mut roots: Vec<u64> = Vec::new();

        let mut i = 0;
        while i < p {
            if self.evaluate(i) % p == 0 {
                roots.push(i);
            }
            i += 1;
        }

        roots
    }

    // Multiplies f and g mod self. "mod" means here, whenever we encounter a power of x greater or
    // equal to the degree of self, we subtract some multiple of self to make this power disappear.
    pub fn mul_mod(&self, f: &MpPolynomial, g: &MpPolynomial) -> MpPolynomial {
        let d = self.degree();

        // Initialize result with the leading coefficient of f times g. This means, the powers
        // of x present in the result are actually shifted up by d - 1.
        let mut result = MpPolynomial::new();
        for i in 0..d {
            result[i] = (&g[i] * &f[d - 1]).into();
        }

        // In each iteration, the leading coffiecient is converted to lower order terms, and the
        // i-th coefficient of self times rhs is added. Thus, the shift in the exponents of powers
        // of x is reduced by 1.
        for i in (0..d - 1).rev() {
            let mut leading_coefficient = Integer::new();
            for coefficient in result.coefficients_mut().iter_mut().take(d) {
                swap(&mut leading_coefficient, coefficient);
            }

            for j in 0..d {
                result[j] += &g[j] * &f[i];
                result[j] -= &leading_coefficient * &self[j];
            }
        }

        result
    }
}

impl Polynomial<Integer> for MpPolynomial {
    fn degree(&self) -> usize {
        let mut d = MAX_DEGREE;
        while d > 0 && self[d] == 0 {
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

impl Display for MpPolynomial {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut s = String::new();
        for i in (0..=self.degree()).rev() {
            s.push_str(&format!("{} * x^{}", self[i], i));
            if i != 0 {
                s.push_str(" + ");
            }
        }
        write!(f, "{}", s)
    }
}

// Naive polynomial selection for general integers. Returns the selected polynomial and an integer
// m, such that f(m) = 0 mod n.
pub fn select(n: &Integer, params: &Params) -> (MpPolynomial, Integer) {
    let d = params.polynomial_degree;
    let m = n.root_ref(d as u32).complete();

    let mut f = MpPolynomial::new();
    let mut x = n.clone();
    for i in 0..=d {
        f[i] = (&x % &m).complete();
        x /= &m;
    }
    assert_eq!(x, 0);

    (f, m)
}
