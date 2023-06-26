use rand_xoshiro::{rand_core::SeedableRng, Xoshiro256PlusPlus};
use rug::{ops::CompleteRound, Assign, Complete, Float, Integer};

use crate::mod_sqrt;

fn smoothness_bound(n: &Integer) -> usize {
    let l = Float::parse(n.to_string()).unwrap().complete(512).ln();
    return Float::exp(0.5 * l.clone().sqrt() * l.ln().sqrt()).to_f64() as usize;
}

fn factor_base(n: &Integer, smoothness_bound: usize) -> Vec<u32> {
    let mut is_prime: Vec<bool> = vec![1 != 0; smoothness_bound + 1];
    let mut factor_base: Vec<u32> = vec![0; 0];

    for i in 2..smoothness_bound + 1 {
        if is_prime[i] {
            if i & 1 == 0 || n.legendre(&Integer::from(i)) == 1 {
                factor_base.push(i as u32);
            }

            let mut j = 2 * i;
            while j <= smoothness_bound {
                is_prime[j] = false;
                j += i;
                is_prime[j] = false;
                j += i;
            }
        }
    }

    factor_base
}

fn get_sieve_interval_len(smoothness_bound: usize) -> usize {
    100 * smoothness_bound
}

fn ilog2_rounded(x: u32) -> u32 {
    ((x as u64 * x as u64).ilog2() + 1) >> 1
}

// TODO: Measure accuracy of sieving heuristic with log.
// TODO: Sieve with prime powers.
pub fn factorize(n: &Integer) -> (Integer, Integer) {
    let smoothness_bound = smoothness_bound(&n);
    let factor_base = factor_base(&n, smoothness_bound);

    // IDEA: use SIMD gather to increment a bunch of sieve values at once? Probably not helpful
    //       since memory is the limit anyway, but try it.
    // IDEA: Use a priority queue to eliminate the need for a sieve array. Just reinsert a prime
    //       with priority x + p, when it was just incremented position x.
    // IDEA: Different method for small numbers (SIMD without gather)
    // IDEA: somehow do prime powers with the prime itself simultaneously. (e.g. increment a float
    //       by 1 / p, convert to integer, add the result * log2(p) and subtract the result from
    //       itself => increment by log2(p) every p steps. May be ok to do such expensive operations
    //       since memory bandwith is probably the bottleneck.)

    // The sieve array contains for each x in [0; m] an approxipation of the sum of logarithms of
    // all primes in the factor base, that divide (x + sqrt_n)^2 - n.
    let m = get_sieve_interval_len(smoothness_bound);
    let mut sieve_array: Vec<u8> = vec![0; m];

    let mut xo = Xoshiro256PlusPlus::seed_from_u64(42);
    let sqrt_n = n.clone().sqrt();

    const SIEVE_ERROR_LIMIT: u32 = 20;

    for &p in &factor_base {
        let log2p = ilog2_rounded(p);
        let t = mod_sqrt::mod_sqrt((n % p).complete().to_u64().unwrap(), p as u64, &mut xo) as u32;

        let mut i: usize = ((t + p - (&sqrt_n % p).complete().to_u32().unwrap()) % p) as usize;
        while i < m {
            sieve_array[i] += log2p as u8;
            i += p as usize;
        }
        if p != 2 {
            i = (((p - t) + p - (&sqrt_n % p).complete().to_u32().unwrap()) % p) as usize;
            while i < m {
                sieve_array[i] += log2p as u8;
                i += p as usize;
            }
        }
    }

    let log2_sqrt_n: u32 = Float::parse(n.to_string())
        .unwrap()
        .complete(512)
        .log2()
        .to_f32() as u32;

    let mut relations: Vec<Vec<u32>> = vec![];
    for x in 0..m {
        let expected = ilog2_rounded(x as u32) + 1 + log2_sqrt_n;
        if sieve_array[x] as u32 + SIEVE_ERROR_LIMIT >= expected {
            // Candidate for a smooth relation.
            let mut y = Integer::new();
            y.assign(&sqrt_n + x);
            y.square_mut();
            y -= n;

            let mut odd_exponent_indices: Vec<u32> = vec![];

            for i in 0..factor_base.len() {
                let mut odd_exponent = false;
                while y.is_divisible_u(factor_base[i]) {
                    y.div_exact_u_mut(factor_base[i]);
                    odd_exponent = !odd_exponent;
                }
                if y == 1 && odd_exponent {
                    odd_exponent_indices.push(i as u32);
                    break;
                }
            }
            if y == 1 {
                relations.push(odd_exponent_indices);
            }
        }
    }

    (Integer::new(), Integer::new())
}
