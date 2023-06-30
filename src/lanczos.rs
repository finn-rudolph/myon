use log::{info, warn};
use rand_xoshiro::{rand_core::SeedableRng, Xoshiro256PlusPlus};

use crate::linalg::{block_matrix, BlockMatrix, CscMatrix, N};

// TODO: add bitmask macro / function / table

// Finds the largest possible amount of rows / columns, such that the principal submatrix of vtav as
// indicated by d is invertible. This function is inspired by the pseudocode in Montgomery, P. L.
// (1995), page 116 and the implementation in msieve.
fn max_invertible_submatrix(mut vtav: BlockMatrix, previous_d: u64) -> (u64, BlockMatrix) {
    let mut w_inv = block_matrix![0; N];
    for i in 0..N {
        w_inv[i] = 1 << i;
    }

    let mut l: usize = 0;
    let mut r: usize = N;
    let mut c = [0; N];
    for i in 0..N {
        if (previous_d >> i) & 1 == 0 {
            c[l] = i;
            l += 1;
        } else {
            r -= 1;
            c[r] = i;
        }
    }
    let mut d = 0;

    // Perform elimination on vtav. All operations are reproduced on w_inv to find an inverse to the
    // submatrix of v, where the i-th row and column is in the submatrix if and only if the i-th bit
    // of d is 1.

    for i in 0..N {
        // Search pivot.
        for j in i..N {
            if (vtav[c[j]] >> c[i]) & 1 == 1 {
                vtav.as_mut().swap(c[i], c[j]);
                w_inv.as_mut().swap(c[i], c[j]);
                break;
            }
        }

        if (vtav[c[i]] >> c[i]) & 1 == 1 {
            // Pivot found.
            d |= 1 << c[i];
            for j in 0..N {
                if i != j && (vtav[c[j]] >> c[i]) & 1 == 1 {
                    vtav[c[j]] ^= vtav[c[i]];
                    w_inv[c[j]] ^= w_inv[c[i]];
                }
            }
        } else {
            for j in i..N {
                if (w_inv[c[j]] >> c[i]) & 1 == 1 {
                    vtav.as_mut().swap(c[i], c[j]);
                    w_inv.as_mut().swap(c[i], c[j]);
                    break;
                }
            }

            assert_eq!((w_inv[c[i]] >> c[i]) & 1, 1u64);

            for j in 0..N {
                if i != j && (w_inv[c[j]] >> c[i]) & 1 == 1 {
                    vtav[c[j]] ^= vtav[c[i]];
                    w_inv[c[j]] ^= w_inv[c[i]];
                }
            }

            vtav[c[i]] = 0;
            w_inv[c[i]] = 0;
        }
    }

    (d, w_inv)
}

// Returns the matrix [v_0T * v_(i + 1) | v_0T * p_(i + 1)] (called delta) for the (i + 1)-th
// iteration using c, vtav w_inv and d from the i-th iteration. delta can be cheaply computed as
// described in Bos, J. W. & Lenstra, A. K. (2017), page 185.
fn update_delta(
    mut delta: [BlockMatrix; 2],
    c: &BlockMatrix,
    vtav: &BlockMatrix,
    w_inv: &BlockMatrix,
    d: u64,
) -> [BlockMatrix; 2] {
    let mut res: [BlockMatrix; 2] = [block_matrix![0; N], block_matrix![0; N]];

    let mut r = block_matrix![0; N];
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

// Finds a matrix x, such that a * x = a * y, but x != y, where A = b * bT.
fn lanczos(b: &CscMatrix, y: &BlockMatrix) -> (BlockMatrix, BlockMatrix) {
    let (n, m) = (b.num_cols(), b.num_rows());

    let v0 = &b.transpose() * &(b * y);
    let mut v = v0.clone();
    let mut p = block_matrix![0; n];
    let mut x = y.clone();
    let mut delta: [BlockMatrix; 2] = [&v0.transpose() * &v0, block_matrix![0; N]];

    let mut d: u64 = !0u64;
    let mut total_d: usize = 0;
    let mut iterations: usize = 0;

    // The recurrence implemented here is based on the decription of Bos, J. W & Lenstra, A. K.
    // (2017), page 184.
    loop {
        let av = &b.transpose() * &(b * &v);
        let vtav = &v.transpose() * &av;
        let vta2v = &av.transpose() * &av;

        let previous_d = d;
        let w_inv: BlockMatrix;
        (d, w_inv) = max_invertible_submatrix(vtav.clone(), d);

        if d == 0 {
            break;
        }
        assert!(w_inv.is_symmetric());

        if total_d + N < m && previous_d | d != !0u64 {
            warn!("Some vectors of v_(i - 1) not included in w_(i - 1) were not included in w_i.");
        }
        total_d += d.count_ones() as usize;

        // Compute c.
        let mut tmp = block_matrix![0; N];
        for i in 0..N {
            tmp[i] = (vta2v[i] & d) ^ (vtav[i] & !d);
        }
        let c = &w_inv * &tmp;

        // Compute v.
        tmp = &v * &c;
        let pvtav = &p * &vtav;
        let vw_inv = &v * &w_inv;
        for i in 0..n {
            v[i] = (av[i] & d) ^ (v[i] & !d) ^ tmp[i] ^ (pvtav[i] & d);
        }

        // Compute p.
        for i in 0..n {
            p[i] = vw_inv[i] ^ (p[i] & !d);
        }

        // Update x.
        tmp = &vw_inv * &delta[0].transpose();
        for i in 0..n {
            x[i] ^= tmp[i];
        }

        // Update delta.
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

// Uses x and vm to find vectors in the nullspace of a by elimination. The returned BlockMatrix
// contains the vectors found in the lower order bits, the remaining bits are zeroed.
fn combine_columns(b: &CscMatrix, mut x: BlockMatrix, vm: BlockMatrix) -> BlockMatrix {
    let (n, m) = (b.num_cols(), b.num_rows());
    let mut r: Vec<Vec<u64>> = (b * &x).explicit_transpose();
    r.append(&mut (b * &vm).explicit_transpose());
    let mut s: Vec<Vec<u64>> = x.explicit_transpose();
    s.append(&mut vm.explicit_transpose());

    let mut i: usize = 0;
    for leading_bit in 0..m {
        if i >= 2 * N {
            break;
        }

        let leading_bit_word: usize = leading_bit / N;
        let leading_bit_mask: u64 = (leading_bit & (N - 1)) as u64;

        // Search pivot.
        let mut j = i;
        while j < 2 * N {
            if r[j][leading_bit_word] & leading_bit_mask != 0 {
                r.swap(i, j);
                s.swap(i, j);
                break;
            }
            j += 1;
        }

        // Eliminate pivot from other columns.
        if r[i][leading_bit_word] & leading_bit_mask != 0 {
            j += 1;
            while j < 2 * N {
                if r[j][leading_bit_word] & leading_bit_mask != 0 {
                    for k in 0..r[i].len() {
                        r[j][k] ^= r[i][k];
                    }
                    for k in 0..s[i].len() {
                        s[j][k] ^= s[i][k];
                    }
                }
                j += 1;
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

// Returns a block of vectors in the nullspace of a.
pub fn find_dependencies(b: &CscMatrix) -> (BlockMatrix, u32) {
    let (n, m) = (b.num_cols(), b.num_rows());

    assert!(
        m < n,
        "matrix has {} rows but only {} columns, so dependencies might not exist.",
        m,
        n
    );

    info!("solving linear system with {} rows and {} columns", m, n);

    let mut xo = Xoshiro256PlusPlus::seed_from_u64(998244353);
    loop {
        let (mut x, vm) = lanczos(b, &BlockMatrix::new_random(n, &mut xo));
        x = combine_columns(b, x, vm);
        let mut u: u64 = 0;
        for i in 0..n {
            u |= x[i];
        }

        if u != 0 {
            // Verify that the vectors of y lie indeed in the nullspace.
            let z = b * &x;
            for i in 0..m {
                assert_eq!(z[i], 0);
            }

            info!(
                "found {} nontrivial vectors in the nullspace",
                u.count_ones()
            );
            return (x, u.count_ones());
        }
        info!("no vectors in the nullspace found, retrying...");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lanczos() {
        let mut xo = Xoshiro256PlusPlus::seed_from_u64((1 << 61) - 1);

        for i in 0..116 {
            let n: usize = 197 + 5 * i;
            let m: usize = n - 19;
            let b = CscMatrix::new_random(n, m, 17, &mut xo);
            let (x, num_dependencies) = find_dependencies(&b);

            assert!(num_dependencies != 0);
            let r = &b * &x;
            for i in 0..m {
                assert_eq!(r[i], 0);
            }
        }
    }
}
