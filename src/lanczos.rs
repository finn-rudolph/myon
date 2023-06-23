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
            d |= 1u64 << i;
            for j in 0..N {
                if i != j && vtav[j] >> i & 1 == 1 {
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
    w_inv: BlockMatrix,
    d: u64,
) -> [BlockMatrix; 2] {
    let mut res: [BlockMatrix; 2] = [blockmatrix![0; N], blockmatrix![0; N]];
    let mut delta: [[BlockMatrix; 2]; 2] = [
        [blockmatrix![0; N], blockmatrix![0; N]],
        [blockmatrix![0; N], blockmatrix![0; N]],
    ];
    for i in 0..N {
        delta[0][0][i] = c[i] ^ ((1u64 << i) & d);
        delta[1][0][i] = vtav[i] & d;
        delta[1][1][i] = d & (1u64 << i);
    }
    delta[0][1] = w_inv;

    res[0] = &gamma[0] * &delta[0][0];
    res[1] = &gamma[0] * &delta[0][1];
    gamma[0] = &gamma[1] * &delta[1][0];
    gamma[1] = &gamma[1] * &delta[1][1];

    for i in 0..N {
        res[0][i] ^= gamma[0][i];
        res[1][i] ^= gamma[1][i];
    }

    res
}

pub fn lanczos(a: &CscMatrix) -> BlockMatrix {
    let n = a.len();

    let mut v = blockmatrix![0; n];
    let mut p = blockmatrix![0; n];
    let mut x = blockmatrix![0; n];
    let mut gamma: [BlockMatrix; 2] = [blockmatrix![0; N], blockmatrix![0; N]];
    gamma[0] = &v.transpose() * &v;

    let mut d: u64;

    x[0] = 998244353; // fill x somewhat randomly (TODO: add RNG later)
    for i in 0..n {
        x[i] = x[i - 1] * x[i - 1];
    }

    loop {
        let av = a * &(a * &v);
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

        gamma = update_gamma(gamma, &c, &vtav, w_inv.clone(), d);

        // Compute new v.

        tmp = &v * &c;
        let pvtav = &p * &vtav;
        for i in 0..n {
            v[i] = (av[i] & d) ^ (v[i] & !d) ^ tmp[i] ^ (pvtav[i] & d);
        }

        let vw_inv = &v * &w_inv;
        for i in 0..n {
            p[i] = vw_inv[i] ^ (p[i] & !d);
        }
    }

    x
}
