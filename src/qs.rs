use rug::{ops::CompleteRound, Float, Integer};

fn smoothness_bound(n: &Integer) -> usize {
    let l = Float::parse(n.to_string()).unwrap().complete(512).ln();
    return Float::exp(0.5 * l.clone().sqrt() * l.ln().sqrt()).to_f64() as usize;
}

fn factor_base(n: &Integer, smoothness_bound: usize) -> Vec<u64> {
    let mut is_prime: Vec<bool> = vec![1 != 0; smoothness_bound + 1];
    let mut factor_base: Vec<u64> = vec![0; 0];

    for i in 2..smoothness_bound + 1 {
        if is_prime[i] {
            if i & 1 == 0 || n.legendre(&Integer::from(i)) == 1 {
                factor_base.push(i as u64);
            }

            let mut j = 2 * i;
            while j <= smoothness_bound {
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

pub fn factorize(n: &Integer) -> (Integer, Integer) {
    let smoothness_bound = smoothness_bound(&n);
    let factor_base = factor_base(&n, smoothness_bound);

    // IDEA: use SIMD gather to increment a bunch of sieve values at once? Probably not helpful
    // since memory is the limit anyway, but try it.
    let sieve_array: Vec<u8> = vec![0; 0];

    (Integer::new(), Integer::new())
}
