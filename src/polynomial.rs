use std::{
    cmp::PartialOrd,
    ops::{Index, IndexMut, MulAssign, Rem, Sub},
};

use rug::{ops::Pow, Complete, Integer};

use crate::params::Params;

pub const MAX_DEGREE: usize = Params::PARAM_TABLE[Params::PARAM_TABLE.len() - 1]
    .1
    .polynomial_degree as usize;

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
