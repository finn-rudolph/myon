use core::fmt::{self, Formatter};
use core::iter;
use std::error::Error;

use rug::{integer::IsPrime, Complete, Integer};

mod lanczos;
mod linalg;
mod nt;
mod qs;

#[derive(Debug)]
pub struct MyonError {
    message: String,
}

impl MyonError {
    fn new(message: &str) -> MyonError {
        MyonError {
            message: message.to_string(),
        }
    }
}

impl Error for MyonError {}

impl fmt::Display for MyonError {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.message)
    }
}

pub fn factorize(n_str: &str) -> Result<Vec<(String, u32)>, MyonError> {
    // Keep a queue of factors of n that still need to be factored (or recognized as primes).
    let n: Integer = match Integer::parse(n_str) {
        Ok(v) => v.complete(),
        Err(_) => {
            return Err(MyonError::new(
                "The given string could not be parsed as an integer.",
            ))
        }
    };

    if n < 0 {
        return Err(MyonError::new(
            "The given integer is negative. Only positive integers are allowed.",
        ));
    }

    let mut ints_to_factor: Vec<Integer> = vec![n];
    let mut primes: Vec<Integer> = vec![];

    const MILLER_RABIN_ITERATIONS: u32 = 74;

    while !ints_to_factor.is_empty() {
        let mut n = ints_to_factor.pop().unwrap();
        let mut exponent_multiplier: u32 = 1;

        // The quadratic sieve requires that n must not be a power, so extract roots until this is
        // satisfied. (TODO: Maybe optimize this to do it only on primes.)
        let mut i: u32 = 2;
        while n.is_perfect_power() {
            let (root, remainder) = n.root_rem_ref(i).complete();
            if remainder == 0 {
                exponent_multiplier *= i;
                n = root;
            }
            i += 1;
        }

        if n.is_probably_prime(MILLER_RABIN_ITERATIONS) != IsPrime::No {
            primes.extend(iter::repeat(n).take(exponent_multiplier as usize));
        } else {
            let (r, s) = qs::factorize(&n);
            for x in [r, s] {
                ints_to_factor.extend(iter::repeat(x).take(exponent_multiplier as usize));
            }
        }
    }

    // Group equal prime factors together and write each factor and it's exponent in a tuple.
    primes.sort();
    let mut factorization: Vec<(String, u32)> = vec![];

    let mut i = 0;
    while i < primes.len() {
        factorization.push((primes[i].to_string(), 1));
        let mut j = i + 1;
        while j < primes.len() && primes[j] == primes[i] {
            factorization.last_mut().unwrap().1 += 1;
            j += 1;
        }
        i = j;
    }

    Ok(factorization)
}
