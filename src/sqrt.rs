use core::mem::swap;
use std::ops::{Index, IndexMut, Mul};

use rug::{Complete, Integer};

use crate::polynomial::{Polynomial, MAX_DEGREE};

#[derive(Clone)]
struct AlgebraicInt<'a> {
    coefficients: [Integer; MAX_DEGREE + 1],
    min_polynomial: &'a Polynomial,
}

impl AlgebraicInt<'_> {
    pub fn new<'a>(min_polynomial: &'a Polynomial) -> AlgebraicInt<'a> {
        AlgebraicInt {
            coefficients: [Integer::ZERO; MAX_DEGREE + 1],
            min_polynomial,
        }
    }
}

impl Index<usize> for AlgebraicInt<'_> {
    type Output = Integer;

    fn index(&self, i: usize) -> &Integer {
        &self.coefficients[i]
    }
}

impl IndexMut<usize> for AlgebraicInt<'_> {
    fn index_mut(&mut self, i: usize) -> &mut Integer {
        &mut self.coefficients[i]
    }
}

impl<'a> Mul<AlgebraicInt<'a>> for AlgebraicInt<'a> {
    type Output = AlgebraicInt<'a>;

    fn mul(self, rhs: AlgebraicInt<'a>) -> AlgebraicInt<'a> {
        let d = self.min_polynomial.degree();

        // Initialize result with the leading coefficient of self times rhs. This means, the powers
        // of alpha present in the result are actually shifted up by d - 1.
        let mut result = AlgebraicInt::new(self.min_polynomial);
        for i in 0..d {
            result[i] = (&rhs[i] * &self[d - 1]).into();
        }

        // In each iteration, the leading coffiecient is converted to lower order terms, and the
        // i-th coefficient of self times rhs is added. Thus, the shift in the exponents of powers
        // of alpha is reduced by 1.
        for i in (0..d - 1).rev() {
            let mut leading_coefficient = Integer::new();
            for coefficient in &mut result.coefficients {
                swap(&mut leading_coefficient, coefficient);
            }

            for (j, coefficient) in result.coefficients.iter_mut().enumerate() {
                *coefficient -= &leading_coefficient * &self.min_polynomial[j];
            }

            for j in 0..d {
                result[j] += &rhs[j] * &self[i];
            }
        }

        result
    }
}

pub fn algebraic_sqrt(integers: Vec<(u32, u32)>) -> Integer {
    todo!()
}

fn multiply_algebraic_integers<'a>(integers: &[AlgebraicInt<'a>]) -> AlgebraicInt<'a> {
    if integers.len() == 1 {
        return integers.first().unwrap().clone();
    }
    multiply_algebraic_integers(&integers[..integers.len() / 2])
        * multiply_algebraic_integers(&integers[integers.len() / 2..])
}

// Caclculates the square root of the product of a set of rational integers. Multiplies all numbers
// together by divide and conquer in O(M(N) * log n), where M(N) denotes the complexity of
// multiplying two numbers in the order of magnitude of the result, and n the size of the set of
// integers.
pub fn rational_sqrt(integers: Vec<Integer>) -> Integer {
    multiply_rational_integers(integers.as_slice()).sqrt()
}

fn multiply_rational_integers(integers: &[Integer]) -> Integer {
    if integers.len() == 1 {
        return integers.first().unwrap().clone();
    }
    multiply_rational_integers(&integers[..integers.len() / 2])
        * multiply_rational_integers(&integers[integers.len() / 2..])
}
