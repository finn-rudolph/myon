use rug::Integer;

pub fn algebraic_sqrt(algebraic_base: Vec<(u32, u32)>) -> Integer {
    todo!()
}

// Caclculates the square root of the product of a set of rational integers. Multiplies all numbers
// together by divide and conquer in O(M(N) * log n), where M(N) denotes the complexity of
// multiplying two numbers in the order of magnitude of the result, and n the size of the set of
// integers.
pub fn rational_sqrt(numbers: Vec<Integer>) -> Integer {
    divide_and_conquer_mul(numbers.as_slice()).sqrt()
}

fn divide_and_conquer_mul(numbers: &[Integer]) -> Integer {
    if numbers.len() == 1 {
        return numbers.first().unwrap().clone();
    }
    divide_and_conquer_mul(&numbers[..numbers.len() / 2])
        * divide_and_conquer_mul(&numbers[numbers.len() / 2..])
}
