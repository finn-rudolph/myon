use rug::Integer;

use crate::polynomial::Polynomial;

// The algebraic square root by the q-adic newton method. Uses divide and conquer to evaluate the
// product in O(M log n) time, where M is the time needed to multiply two numbers in the order of
// magnitude of the result.
pub fn algebraic_sqrt(integers: Vec<Polynomial>, f: &Polynomial) -> Polynomial {
    let product = multiply_algebraic_integers(&integers, f);

    todo!()
}

fn multiply_algebraic_integers(integers: &[Polynomial], f: &Polynomial) -> Polynomial {
    if integers.len() == 1 {
        return integers.first().unwrap().clone();
    }
    f.multiply_mod(
        &multiply_algebraic_integers(&integers[..integers.len() / 2], f),
        &multiply_algebraic_integers(&integers[integers.len() / 2..], f),
    )
}

// Caclculates the square root of the product of a set of rational integers.
pub fn rational_sqrt(integers: Vec<Integer>) -> Integer {
    multiply_rational_integers(&integers).sqrt()
}

fn multiply_rational_integers(integers: &[Integer]) -> Integer {
    if integers.len() == 1 {
        return integers.first().unwrap().clone();
    }
    multiply_rational_integers(&integers[..integers.len() / 2])
        * multiply_rational_integers(&integers[integers.len() / 2..])
}
