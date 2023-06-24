#[path = "linalg.rs"]
pub mod linalg;

use linalg::blockmatrix;
use linalg::BlockMatrix;
use linalg::CscMatrix;
use linalg::N;
use rand_xoshiro::rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;

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

fn update_delta(
    mut delta: [BlockMatrix; 2],
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

    res[0] = &delta[0] * &r;
    res[1] = &delta[0] * w_inv;
    delta[0] = &delta[1] * vtav;

    for i in 0..N {
        res[0][i] ^= delta[0][i] & d;
        res[1][i] ^= delta[1][i] & !d;
    }

    res
}

fn lanczos(a: &CscMatrix, b: &BlockMatrix) -> (BlockMatrix, BlockMatrix) {
    let n = a.len();

    let v0 = &a.transpose() * &(a * b);
    let mut v = v0.clone();
    let mut p = blockmatrix![0; n];
    let mut x = b.clone();
    let mut delta: [BlockMatrix; 2] = [blockmatrix![0; N], blockmatrix![0; N]];
    delta[0] = &v0.transpose() * &v0;

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

        // Update x.

        tmp = &vw_inv * &delta[0].transpose();
        for i in 0..n {
            x[i] ^= tmp[i];
        }
        if iterations != 0 {
            delta = update_delta(delta, &c, &vtav, &w_inv, d);
        } else {
            delta[0] = &v0.transpose() * &v;
            delta[1] = &v0.transpose() * &p;
        }

        iterations += 1;
    }

    (x, v)
}

// Uses x and vm to find vectors in the nullspace of a by elimination. The
// returned BlockMatrix contains the vectors found in the lower order bits, the
// remaining bits are zeroed.
fn combine_columns(a: &CscMatrix, mut x: BlockMatrix, vm: BlockMatrix) -> BlockMatrix {
    let n = a.len();
    let mut r: Vec<Vec<u64>> = (a * &x).explicit_transpose();
    r.append(&mut (a * &vm).explicit_transpose());
    let mut s: Vec<Vec<u64>> = x.explicit_transpose();
    s.append(&mut vm.explicit_transpose());

    let mut i: usize = 0;
    for leading_bit in 0..n {
        if i >= 2 * N {
            break;
        }

        let leading_bit_col: usize = leading_bit / N;
        let leading_bit_mask: u64 = (leading_bit & (N - 1)) as u64;

        for j in i..2 * N {
            if r[i][leading_bit_col] & leading_bit_mask != 0 {
                r.swap(i, j);
                s.swap(i, j);
                break;
            }
        }

        if r[i][leading_bit_col] & leading_bit_mask != 0 {
            for j in i..2 * N {
                if r[j][leading_bit_col] & leading_bit_mask != 0 {
                    for k in 0..r[i].len() {
                        r[j][k] ^= r[i][k];
                        s[j][k] ^= s[i][k];
                    }
                }
            }
            i += 1;
        }
    }

    for j in 0..n {
        x[j] = 0;
        for k in i..N {
            x[j] |= ((s[k][j / N] >> (j & (N - 1))) & 1) << (k - i);
        }
    }

    x
}

// Returns a block of vectors in the null space of a.
pub fn find_dependencies(a: &CscMatrix) -> BlockMatrix {
    let n = a.len();
    let mut xo = Xoshiro256StarStar::seed_from_u64(998244353);
    let y = BlockMatrix::new_random(n, &mut xo);

    let (x, vm) = lanczos(a, &y);
    let null_matrix = combine_columns(a, x, vm);

    let z = a * &null_matrix;
    for i in 0..n {
        assert_eq!(z[i], 0);
    }
    null_matrix
}

#[cfg(test)]
mod tests {
    use rand_xoshiro::rand_core::RngCore;
    use rand_xoshiro::rand_core::SeedableRng;
    use rand_xoshiro::Xoshiro256StarStar;

    use super::lanczos;
    use super::BlockMatrix;
    use super::CscMatrix;

    const NUM_TESTS: usize = 42;
    const MIN_N: usize = 100;
    const MAX_N: usize = 2000;
    const MAX_NM_DIFF: usize = 20;
    const MAX_ONES_FRAC: usize = 5;

    #[test]
    fn test_lanczos() {
        let mut xo = Xoshiro256StarStar::seed_from_u64((1 << 61) - 1);

        for _ in 0..NUM_TESTS {
            let n = (xo.next_u32() as usize % (MAX_N - MIN_N)) + MIN_N;
            let m = n - (xo.next_u32() as usize % MAX_NM_DIFF);

            let a = CscMatrix::new_random(&mut xo, n, m, m / MAX_ONES_FRAC);
            let b = BlockMatrix::new_random(n, &mut xo);
            let (x, vm) = lanczos(&a, &b);
            let y = &a.transpose() * &(&a * &x);
            let z = &a.transpose() * &(&a * &b);
            for i in 0..n {
                assert_eq!(y[i], z[i]);
            }
        }
    }
}
