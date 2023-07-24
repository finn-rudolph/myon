use rug::Integer;

use crate::{gfpolynomial::GfPolynomial, nt, polynomial::Polynomial};

// The algebraic square root by the q-adic newton method. Uses divide and conquer to evaluate the
// product in O(M log n) time, where M is the time needed to multiply two numbers in the order of
// magnitude of the result.
pub fn algebraic_sqrt(integers: Vec<Polynomial>, f: &Polynomial) -> Polynomial {
    let product = mul_algebraic_integers(&integers, f);

    todo!()
}

fn mul_algebraic_integers(integers: &[Polynomial], f: &Polynomial) -> Polynomial {
    if integers.len() == 1 {
        return integers.first().unwrap().clone();
    }
    f.mul_mod(
        &mul_algebraic_integers(&integers[..integers.len() / 2], f),
        &mul_algebraic_integers(&integers[integers.len() / 2..], f),
    )
}

fn select_p(f: &Polynomial) -> u32 {
    let mut p: u32 = 998244353;
    loop {
        // p must be inert in the number field, which means f must be irreducible mod p. Check
        // this by Rabin's test of irreducibility.
        if nt::miller_rabin(p) && GfPolynomial::from_polynomial(f, p).is_irreducible() {
            return p;
        }
        p += 2;
    }
}

fn inverse_sqrt_mod_p(product: &Polynomial, f: &Polynomial) -> Polynomial {
    todo!()
}

// Caclculates the square root of the product of a set of rational integers.
pub fn rational_sqrt(integers: Vec<Integer>) -> Integer {
    mul_rational_integers(&integers).sqrt()
}

fn mul_rational_integers(integers: &[Integer]) -> Integer {
    if integers.len() == 1 {
        return integers.first().unwrap().clone();
    }
    mul_rational_integers(&integers[..integers.len() / 2])
        * mul_rational_integers(&integers[integers.len() / 2..])
}
