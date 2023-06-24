#[path = "linalg.rs"]
pub mod linalg;

use linalg::blockmatrix;
use linalg::BlockMatrix;
use linalg::CscMatrix;
use linalg::N;

// Finds the largest possible amount of rows / columns, such that the principal
// submatrix of vtav as indicated by d is invertible.
fn max_invertible_submatrix(mut vtav: BlockMatrix, previous_d: u64) -> (u64, BlockMatrix) {
    let mut w_inv = blockmatrix![0; N];
    for i in 0..N {
        w_inv[i] |= 1 << i;
    }
    let mut k: usize = 0;
    for i in 0..N {
        if (previous_d >> i) & 1 == 0 {
            vtav.swap(i, k);
            w_inv.swap(i, k);
            k += 1;
        }
    }
    let mut d = 0;

    // Perform elimination on vtav. All operations are reproduced on w_inv to
    // find an inverse to the submatrix of v, where the i-th row and column is
    // in the submatrix if and only if the i-th bit of d is 1.

    for i in 0..N {
        for j in i..N {
            if (vtav[j] >> i) & 1 == 1 {
                vtav.swap(i, j);
                w_inv.swap(i, j);
                break;
            }
        }

        if (vtav[i] >> i) & 1 == 1 {
            d |= 1 << i;
            for j in 0..N {
                if i != j && (vtav[j] >> i) & 1 == 1 {
                    vtav[j] ^= vtav[i];
                    w_inv[j] ^= w_inv[i];
                }
            }
        } else {
            for j in i..N {
                if (w_inv[j] >> i) & 1 == 1 {
                    vtav.swap(i, j);
                    w_inv.swap(i, j);
                    break;
                }
            }

            assert_eq!((w_inv[i] >> i) & 1, 1u64);

            for j in 0..N {
                if i != j && (w_inv[j] >> i) & 1 == 1 {
                    vtav[j] ^= vtav[i];
                    w_inv[j] ^= w_inv[i];
                }
            }

            vtav[i] = 0;
            w_inv[i] = 0;
        }
    }

    (d, w_inv)
}

fn update_gamma(
    mut gamma: [BlockMatrix; 2],
    c: &BlockMatrix,
    vtav: &BlockMatrix,
    w_inv: &BlockMatrix,
    d: u64,
) -> [BlockMatrix; 2] {
    let mut res: [BlockMatrix; 2] = [blockmatrix![0; N], blockmatrix![0; N]];

    let mut r = blockmatrix![0; N];
    for i in 0..N {
        r[i] = c[i] ^ ((1u64 << i) & !d);
    }

    res[0] = &gamma[0] * &r;
    res[1] = &gamma[0] * w_inv;
    gamma[0] = &gamma[1] * vtav;

    for i in 0..N {
        res[0][i] ^= gamma[0][i] & d;
        res[1][i] ^= gamma[1][i] & !d;
    }

    res
}

// Finds a matrix x, such that a * aT * x = b, using the Block Lanczos
// algorithm.
pub fn lanczos(a: &CscMatrix, b: &BlockMatrix) -> BlockMatrix {
    let n = a.len();

    let mut v = b.clone();
    let mut p = blockmatrix![0; n];
    let mut x = blockmatrix![0; n];
    let mut gamma: [BlockMatrix; 2] = [blockmatrix![0; N], blockmatrix![0; N]];
    gamma[0] = &v.transpose() * &v;

    let mut d: u64 = 64;
    let mut total_d: usize = 0;
    let mut iterations: usize = 0;

    loop {
        let av = &a.transpose() * &(a * &v);
        let vtav = &v.transpose() * &av;
        let vta2v = &av.transpose() * &av;
        let w_inv: BlockMatrix;

        let previous_d = d;
        (d, w_inv) = max_invertible_submatrix(vtav.clone(), d);

        if d == 0 {
            break;
        }
        assert!(w_inv.is_symmetric());

        if total_d + N < a.m {
            assert_eq!(!((1u64 << previous_d.count_zeros()) - 1) | d, !0u64);
        }
        total_d += d.count_ones() as usize;

        let mut tmp = blockmatrix![0; N];
        for i in 0..N {
            tmp[i] = (vta2v[i] & d) ^ (vtav[i] & !d);
        }
        let c = &w_inv * &tmp;

        // Compute new v.

        tmp = &v * &c;
        let pvtav = &p * &vtav;
        let previous_v = v.clone();
        for i in 0..n {
            v[i] = (av[i] & d) ^ (previous_v[i] & !d) ^ tmp[i] ^ (pvtav[i] & d);
        }

        let vw_inv = &previous_v * &w_inv;
        for i in 0..n {
            p[i] = vw_inv[i] ^ (p[i] & !d);
        }

        // Compute new gamma and update x.

        tmp = &vw_inv * &gamma[0].transpose();
        for i in 0..n {
            x[i] ^= tmp[i];
        }
        if iterations != 0 {
            gamma = update_gamma(gamma, &c, &vtav, &w_inv, d);
        } else {
            gamma[0] = &b.transpose() * &v;
            gamma[1] = &b.transpose() * &p;
        }

        iterations += 1;
    }

    x
}

#[cfg(test)]
mod tests {
    use rand_xoshiro::rand_core::RngCore;
    use rand_xoshiro::rand_core::SeedableRng;
    use rand_xoshiro::Xoshiro256StarStar;

    use super::blockmatrix;
    use super::lanczos;
    use super::BlockMatrix;
    use super::CscMatrix;

    fn random_sparse_matrix(xo: &mut Xoshiro256StarStar, n: usize, m: usize) -> CscMatrix {
        let mut end: Vec<u32> = vec![0; 0];
        let mut ones: Vec<u32> = vec![0; 0];
        let avg_ones = m / 5 + (xo.next_u32() as usize % (m / 5 - m / 40));

        for _ in 0..n {
            let weight = xo.next_u32() as usize % std::cmp::min(m, avg_ones << 1);
            for _ in 0..weight {
                ones.push(xo.next_u32() % m as u32);
            }
            end.push(ones.len() as u32);
        }

        CscMatrix::new(m, end, ones)
    }

    fn random_block_matrix(xo: &mut Xoshiro256StarStar, n: usize) -> BlockMatrix {
        let mut a = blockmatrix![0; n];
        for i in 0..n {
            a[i] = xo.next_u64();
        }
        a
    }

    #[test]
    fn test_lanczos() {
        let mut xo = Xoshiro256StarStar::seed_from_u64((1 << 42) - 42);

        for _ in 0..42 {
            let n = (xo.next_u32() as usize % 1900) + 100;
            let m = n - (xo.next_u32() as usize % 20);

            let a = random_sparse_matrix(&mut xo, n, m);
            let b = random_block_matrix(&mut xo, n);
            let x = lanczos(&a, &b);
            let y = &a.transpose() * &(&a * &x);
            for i in 0..n {
                assert_eq!(y[i], b[i]);
            }
        }
    }
}
