use std::collections::VecDeque;

use rug::{integer::IsPrime, Complete, Integer};

mod lanczos;
mod linalg;
mod mod_sqrt;
mod qs;

pub fn factorize(n_str: &str) -> Vec<(String, u32)> {
    let mut primes: Vec<Integer> = vec![];
    let mut composites: VecDeque<Integer> =
        VecDeque::from([Integer::parse(n_str).unwrap().complete()]);

    const MILLER_RABIN_ITERATIONS: u32 = 74;

    while !composites.is_empty() {
        let n = composites.pop_front().unwrap();
        let (r, s) = qs::factorize(&n);

        if r.is_probably_prime(MILLER_RABIN_ITERATIONS) == IsPrime::No {
            composites.push_back(r);
        } else {
            primes.push(r);
        }

        if s.is_probably_prime(MILLER_RABIN_ITERATIONS) == IsPrime::No {
            composites.push_back(s);
        } else {
            primes.push(s);
        }
    }

    primes.sort();
    let mut factorization: Vec<(String, u32)> = vec![];

    let mut i = 0;
    while i < primes.len() {
        factorization.push((primes[i].to_string(), 1));
        let mut j = i + 1;
        while j < primes.len() && primes[j] == primes[i] {
            factorization.last_mut().as_deref_mut().unwrap().1 += 1;
            j += 1;
        }
        i = j;
    }

    factorization
}
