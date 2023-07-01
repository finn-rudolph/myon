use log::{error, info, STATIC_MAX_LEVEL};
use rug::{ops::CompleteRound, Complete, Float, Integer};

use crate::{lanczos, linalg::CscMatrix, nt};

mod params {
    use rug::Integer;

    struct Params {
        bits: usize,
        factor_base_size: usize,
        sieve_array_len: usize,
        polynomial_batch_size: usize,
    }

    impl Params {
        fn new(
            bits: usize,
            factor_base_size: usize,
            sieve_array_len: usize,
            polynomial_batch_size: usize,
        ) -> Params {
            Params {
                bits,
                factor_base_size,
                sieve_array_len,
                polynomial_batch_size,
            }
        }
    }

    const SIQS_PARAMS: [Params; 8] = [
        Params::new(80_, 100_, 5000__, 4),
        Params::new(100, 200_, 25000_, 4),
        Params::new(120, 400_, 25000_, 5),
        Params::new(140, 900_, 50000_, 5),
        Params::new(159, 1200, 100000, 6),
        Params::new(179, 2000, 250000, 6),
        Params::new(199, 3000, 350000, 7),
        Params::new(219, 4500, 500000, 7),
    ];

    fn search_index(n: &Integer) -> usize {
        let bits = n.significant_bits();
        let (mut a, mut b) = (0, SIQS_PARAMS.len() - 1);
        while a < b {
            let mid = (a + b) / 2;
            if SIQS_PARAMS[mid].bits <= bits as usize {
                b = mid;
            } else {
                a = mid + 1;
            }
        }

        a
    }

    pub fn factor_base_size(n: &Integer) -> usize {
        SIQS_PARAMS[search_index(n)].factor_base_size
    }

    pub fn sieve_array_len(n: &Integer) -> usize {
        SIQS_PARAMS[search_index(n)].sieve_array_len
    }

    pub fn polynomial_batch_size(n: &Integer) -> usize {
        SIQS_PARAMS[search_index(n)].polynomial_batch_size
    }
}

struct FactorBaseElem {
    p: u32,
    t: u32, // One solution to t^2 = n mod p (the other one is just -t).
}

impl FactorBaseElem {
    fn new(p: u32, n: &Integer) -> FactorBaseElem {
        FactorBaseElem {
            p,
            t: nt::mod_sqrt((n % p).complete().to_u64().unwrap(), p as u64) as u32,
        }
    }
}

// square_part^2 = smooth_part mod n. A "smooth" relation where smooth_part splits over the factor
// base. ones contains the indices of primes in the factor base, whose exponents in smooth_part's
// factorization are odd.
struct Relation {
    ones: Vec<u32>,
    smooth_part: Integer,
    square_part: Integer,
}

impl Relation {
    fn new(square_part: Integer, smooth_part: Integer, ones: Vec<u32>) -> Relation {
        Relation {
            square_part,
            smooth_part,
            ones,
        }
    }
}

// Represents the polynomial (q^2 x + b)^2 - n.
struct Polynomial {
    q: Integer,
    b: Integer,
    a: Integer,
    c: Integer,
}

impl Polynomial {
    fn new(a: &Integer) -> Polynomial {}

    fn eval(&self, x: usize) -> Integer {
        (&self.a * x).complete() + (x << 1) * &self.b + &self.c
    }
}

fn factor_base(n: &Integer, factor_base_size: usize) -> Vec<FactorBaseElem> {
    let mut factor_base: Vec<FactorBaseElem> = vec![];

    let mut i = 2;
    while factor_base.len() < factor_base_size {
        if nt::is_prime(i) && (i & 1 == 0 || n.legendre(&Integer::from(i)) == 1) {
            factor_base.push(FactorBaseElem::new(i as u32, n));
        }
        i += 1;
    }

    factor_base
}

fn ilog2_rounded(x: u32) -> u32 {
    ((x as u64 * x as u64).ilog2() + 1) >> 1
}

// IDEA: use SIMD gather to increment a bunch of sieve values at once? Probably not helpful
//       since memory is the limit anyway, but try it.
// IDEA: Use a priority queue to eliminate the need for a sieve array. Just reinsert a prime
//       with priority x + p, when it was just incremented position x.
// IDEA: Different method for small numbers (SIMD without gather), especially for 2
// IDEA: somehow do prime powers with the prime itself simultaneously. (e.g. increment a float
//       by 1 / p, convert to integer, add the result * log2(p) and subtract the result from
//       itself => increment by log2(p) every p steps. May be ok to do such expensive operations
//       since memory bandwith is probably the bottleneck.)

// The sieve array contains for each x in [0; m] an approxipation of the sum of logarithms of
// all primes in the factor base that divide (x + sqrt_n)^2 - n.

fn sieve(f: Polynomial, factor_base: &Vec<FactorBaseElem>, m: usize) -> Vec<Relation> {
    let mut sieve_array: Vec<u8> = vec![0; m];

    let a = f.q.clone().square();

    // Don't sieve with 2 for now.
    for &FactorBaseElem { p, t } in factor_base.iter().skip(1) {
        let log2p = ilog2_rounded(p) as u8;
        // Initialization stage (The costly thing we optimize with SIQS)
        let a_inv = nt::mod_inverse((&a % p).complete().to_u64().unwrap(), p as u64);
        let b_mod_p = (&f.b % p).complete().to_u32().unwrap();

        let mut i = (((p + t - b_mod_p) as u64 * a_inv) % p as u64) as usize;
        while i < m {
            sieve_array[i] += log2p;
            i += p as usize;
        }

        i = ((((p << 1) - t - b_mod_p) as u64 * a_inv) % p as u64) as usize;
        while i < m {
            sieve_array[i] += log2p;
            i += p as usize;
        }
    }

    let mut relations: Vec<Relation> = vec![];
    const SIEVE_ERROR_LIMIT: u32 = 20;

    for x in 0..m {
        let mut y = f.eval(x);
        let expected = y.significant_bits(); // log2(f(x))
        if sieve_array[x] as u32 + SIEVE_ERROR_LIMIT >= expected {
            // Candidate for a smooth relation.
            let mut odd_exponent_indices: Vec<u32> = vec![];

            for i in 0..factor_base.len() {
                let mut odd_exponent = false;
                while y.is_divisible_u(factor_base[i].p) {
                    y.div_exact_u_mut(factor_base[i].p);
                    odd_exponent = !odd_exponent;
                }
                if odd_exponent {
                    odd_exponent_indices.push(i as u32);
                }
            }

            if y == 1 {
                relations.push(Relation::new(
                    (&f.a * x).complete() + &f.b,
                    f.eval(x),
                    odd_exponent_indices,
                ));
            }
        }
    }

    relations
}

// TODO: Measure accuracy of sieving heuristic with log.
// TODO: Sieve with prime powers.
pub fn factorize(n: &Integer) -> (Integer, Integer) {
    assert!(!n.is_even());
    assert!(!n.is_perfect_power());

    info!("factor base size set to {}", params::factor_base_size(n));
    let factor_base = factor_base(n, params::factor_base_size(n));
    info!(
        "chose a factor base consisting of {} primes",
        factor_base.len()
    );

    let m = params::sieve_array_len(n);

    let mut relations: Vec<Relation> = vec![];
    let mut q: Integer = ((n << 1u32).complete().sqrt() / m).sqrt().next_prime();
    while relations.len() <= factor_base.len() {
        while n.legendre(&q) != 1 {
            q = q.next_prime();
        }

        q = q.next_prime();
    }

    info!(
        "collected {} smooth relations, starting Block Lanczos",
        relations.len()
    );

    let (x, num_dependencies) = lanczos::find_dependencies(&CscMatrix::new(
        relations.iter().map(|r| r.ones.clone()).collect(),
        factor_base.len(),
    ));
    for i in 0..num_dependencies {
        let (mut a_squared, mut b) = (Integer::from(1), Integer::from(1));
        for j in 0..x.as_ref().len() {
            if (x[j] >> i) & 1 == 1 {
                // The j-th relation is included.
                a_squared *= relations[j].smooth_part.clone().square() - n;
                b *= &relations[j].square_part;
            }
        }

        assert!(a_squared.is_perfect_square());
        assert_eq!(
            (&a_squared % n).complete(),
            (&b.square_ref().complete() % n).complete()
        );
        let c = Integer::from(n.gcd_ref(&(a_squared.sqrt() + b)));
        if c != 1 && &c != n {
            return ((n / &c).complete(), c);
        }
    }

    error!("only trivial divisors found");
    panic!();
}

#[cfg(test)]
mod test {
    use rug::rand::RandState;

    use super::*;

    #[test]
    fn test_quadratic_sieve_general() {}

    #[test]
    fn test_quadratic_sieve_semiprimes() {
        let mut rng = RandState::new();
        rng.seed(&Integer::from(42));
        for _ in 0..10 {
            let p = Integer::from(Integer::random_bits(30, &mut rng)).next_prime();
            let mut q = Integer::from(Integer::random_bits(30, &mut rng)).next_prime();
            while p == q {
                q = Integer::from(Integer::random_bits(30, &mut rng)).next_prime();
            }
            let (s, t) = factorize(&(&p * &q).complete());
            assert!((s == p && t == q) || (s == q && t == p));
        }
    }
}
