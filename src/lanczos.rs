#[path = "linalg.rs"]
pub mod linalg;

use linalg::blockmatrix;
use linalg::BlockMatrix;
use linalg::CscMatrix;
use linalg::N;

fn max_invertible_submatrix(mut vtav: BlockMatrix) -> (u64, BlockMatrix) {
    let mut w_inv = blockmatrix![0; N];
    for i in 0..N {
        w_inv[i] |= 1 << i;
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

    let mut d: u64;

    loop {
        let av = &a.transpose() * &(a * &v);
        let vtav = &v.transpose() * &av;
        let vta2v = &av.transpose() * &av;
        let w_inv: BlockMatrix;

        (d, w_inv) = max_invertible_submatrix(vtav.clone());

        if d == 0 {
            break;
        }
        assert!(w_inv.is_symmetric());

        let mut tmp = blockmatrix![0; N];
        for i in 0..N {
            tmp[i] = (vta2v[i] & d) ^ (vtav[i] & !d);
        }
        let c = &w_inv * &tmp;

        // Compute new gamma and update x.

        tmp = &w_inv * &gamma[0].transpose();
        tmp = &v * &tmp;
        for i in 0..n {
            x[i] ^= tmp[i];
        }

        gamma = update_gamma(gamma, &c, &vtav, &w_inv, d);

        // Compute new v.

        tmp = &v * &c;
        let pvtav = &p * &vtav;
        let ov = v.clone();
        for i in 0..n {
            v[i] = (av[i] & d) ^ (ov[i] & !d) ^ tmp[i] ^ (pvtav[i] & d);
        }

        let vw_inv = &ov * &w_inv;
        for i in 0..n {
            p[i] = vw_inv[i] ^ (p[i] & !d);
        }
    }

    x
}

#[cfg(test)]
mod tests {
    use rand::{thread_rng, Rng};

    use super::blockmatrix;
    use super::lanczos;
    use super::BlockMatrix;
    use super::CscMatrix;

    fn random_sparse_matrix(n: usize, m: usize, avg_ones: usize) -> CscMatrix {
        let mut end: Vec<u32> = vec![0; 0];
        let mut ones: Vec<u32> = vec![0; 0];

        for _ in 0..n {
            let weight = thread_rng().gen_range(0..(avg_ones << 1));
            for _ in 0..weight {
                ones.push(thread_rng().gen_range(0..m) as u32);
            }
            end.push(ones.len() as u32);
        }

        CscMatrix::new(end, ones)
    }

    fn random_block_matrix(n: usize) -> BlockMatrix {
        let mut a = blockmatrix![0; n];
        for i in 0..n {
            a[i] = thread_rng().gen::<u64>();
        }
        a
    }

    #[test]
    fn test_lanczos() {
        for _ in 0..42 {
            let n = thread_rng().gen_range(100..2000);
            let m = thread_rng().gen_range(n - 20..n);
            let avg_ones = thread_rng().gen_range(m / 40..m / 5);

            let a = random_sparse_matrix(n, m, avg_ones);
            let b = random_block_matrix(n);
            let x = lanczos(&a, &b);
            let y = &a * &x;
            for i in 0..n {
                assert_eq!(y[i], b[i]);
            }
        }
    }
}
