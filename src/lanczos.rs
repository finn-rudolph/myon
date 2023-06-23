#[path = "linalg.rs"]
pub mod linalg;

use linalg::blockmatrix;
use linalg::BlockMatrix;
use linalg::CscMatrix;
use linalg::N;

pub fn max_invertible_submatrix(mut vtav: BlockMatrix) -> (u64, BlockMatrix) {
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

pub fn lanczos(a: &CscMatrix) -> BlockMatrix {
    let n = a.len();
    let mut v = blockmatrix![0; n];
    let mut p = blockmatrix![0; n];
    let mut x = blockmatrix![0; n];
    let mut tmp = blockmatrix![0; n];
    let mut d: u64 = 0;

    x[0] = 998244353; // fill x somewhat randomly (TODO: add RNG later)
    for i in 0..n {
        x[i] = x[i - 1] * x[i - 1];
    }

    loop {
        let av = a * &v;
        let vtav = &v.transpose() * &av;
        let vta2v = &v.transpose() * &v;
        let w_inv: BlockMatrix;
        (d, w_inv) = max_invertible_submatrix(vtav);

        if d == 0 {
            break;
        }
        assert!(w_inv.is_symmetric());
    }

    x
}
