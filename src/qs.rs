#[path = "lanczos.rs"]
mod lanczos;

use rug::{ops::CompleteRound, Float, Integer};

fn get_smoothness_bound(n: &Integer) -> u64 {
    let l = Float::parse(n.to_string()).unwrap().complete(512).ln();
    return Float::exp(0.5 * l.clone().sqrt() * l.ln().sqrt()).to_f64() as u64;
}

fn get_factor_base(n: &Integer, smoothness_bound: u64) -> Vec<u64> {
    let mut is_prime: Vec<bool> = vec![1 != 0; smoothness_bound as usize + 1];
    let mut factor_base: Vec<u64> = vec![0; 0];

    for i in 2..smoothness_bound + 1 {
        if is_prime[i as usize] {
            if i & 1 == 0 || n.legendre(&Integer::from(i)) == 1 {
                factor_base.push(i);
            }

            let mut j = 2 * i;
            while j <= smoothness_bound {
                is_prime[j as usize] = false;
                j += i;
            }
        }
    }

    factor_base
}

pub fn factorize(n: &Integer) -> (Integer, Integer) {
    (Integer::new(), Integer::new())
}
