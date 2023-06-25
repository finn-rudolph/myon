use std::convert::From;
use std::ops::{Index, IndexMut, Mul};

use rand_xoshiro::{rand_core::RngCore, Xoshiro256PlusPlus};

pub const N: usize = 64;

macro_rules! blockmatrix {
    ( $x:expr; $n:expr ) => {
        BlockMatrix::from(vec![$x; $n])
    };
}

pub(crate) use blockmatrix;

// Column-major sparse matrix storing for each column the ones' positions in a contiguous subsegment
// in 'ones'. The index after the last element of column i is end[i].
pub struct CscMatrix {
    m: usize,      // number of rows
    end: Vec<u32>, // number of columns = end.len()
    ones: Vec<u32>,
}

pub struct CscMatrixTranspose<'a> {
    borrowed: &'a CscMatrix,
}

// A dense binary matrix storing each row as an N-bit integer.
#[derive(Clone)]
pub struct BlockMatrix(Vec<u64>);

pub struct BlockMatrixTranspose<'a> {
    borrowed: &'a BlockMatrix,
}

impl CscMatrix {
    pub fn new(m: usize, end: Vec<u32>, ones: Vec<u32>) -> CscMatrix {
        CscMatrix { m, end, ones }
    }

    pub fn new_random(
        xo: &mut Xoshiro256PlusPlus,
        n: usize,
        m: usize,
        max_ones: usize,
    ) -> CscMatrix {
        let mut end: Vec<u32> = vec![0; 0];
        let mut ones: Vec<u32> = vec![0; 0];

        for _ in 0..n {
            let weight = xo.next_u32() as usize % max_ones;
            for _ in 0..weight {
                ones.push(xo.next_u32() % m as u32);
            }
            end.push(ones.len() as u32);
        }

        CscMatrix::new(m, end, ones)
    }

    pub fn num_cols(&self) -> usize {
        self.end.len()
    }

    pub fn num_rows(&self) -> usize {
        self.m
    }

    // Returns a view on the transposed matrix. The view is tightly bound to the original CscMatrix
    // and is intended to be used only in composition with the '*'-Operator.
    pub fn transpose(&self) -> CscMatrixTranspose {
        CscMatrixTranspose { borrowed: self }
    }
}

impl BlockMatrix {
    pub fn new_random(n: usize, xo: &mut Xoshiro256PlusPlus) -> BlockMatrix {
        let mut a = blockmatrix![0; n];
        for i in 0..n {
            a[i] = xo.next_u64();
        }
        a
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn swap(&mut self, i: usize, j: usize) {
        self.0.swap(i, j);
    }

    // Provides a lightweight view on the transposed matrix, which isn't intendend to be used
    // standalone, but as an argument to the '*'-Operator (on any side).
    pub fn transpose(&self) -> BlockMatrixTranspose {
        BlockMatrixTranspose { borrowed: self }
    }

    // Calculates the transpose explicity as a two-dimensional vector, in row-major format.
    pub fn explicit_transpose(&self) -> Vec<Vec<u64>> {
        let n_words = (self.len() + N - 1) / N;
        let mut res: Vec<Vec<u64>> = vec![vec![0; n_words]; N];

        for i in 0..self.len() {
            for j in 0..N {
                res[j][i / N] |= ((self[i] >> j) & 1) << (i & (N - 1));
            }
        }

        res
    }

    pub fn is_symmetric(&self) -> bool {
        assert_eq!(self.len(), N);
        for i in 0..N {
            for j in 0..N {
                if (self[i] >> j) & 1 != (self[j] >> i) & 1 {
                    return false;
                }
            }
        }
        return true;
    }
}

impl From<Vec<u64>> for BlockMatrix {
    fn from(x: Vec<u64>) -> Self {
        BlockMatrix(x)
    }
}

impl Index<usize> for BlockMatrix {
    type Output = u64;

    fn index(&self, i: usize) -> &u64 {
        &self.0[i]
    }
}

impl IndexMut<usize> for BlockMatrix {
    fn index_mut(&mut self, i: usize) -> &mut u64 {
        &mut self.0[i]
    }
}

impl Mul<&BlockMatrix> for &CscMatrix {
    type Output = BlockMatrix;

    fn mul(self, b: &BlockMatrix) -> BlockMatrix {
        let (n, m) = (self.end.len(), self.m);
        assert_eq!(n, b.len());
        let mut res = blockmatrix![0; m];

        let mut j: usize = 0;
        for i in 0..n {
            while j < self.end[i] as usize {
                res[self.ones[j] as usize] ^= b[i];
                j += 1;
            }
        }

        res
    }
}

impl<'a> Mul<&BlockMatrix> for &CscMatrixTranspose<'a> {
    type Output = BlockMatrix;

    fn mul(self, b: &BlockMatrix) -> BlockMatrix {
        let (n, m) = (self.borrowed.end.len(), self.borrowed.m);
        assert_eq!(m, b.len());
        let mut res = blockmatrix![0; n];

        let mut j: usize = 0;
        for i in 0..n {
            while j < self.borrowed.end[i] as usize {
                res[i] ^= b[self.borrowed.ones[j] as usize];
                j += 1;
            }
        }

        res
    }
}

impl Mul<&BlockMatrix> for &BlockMatrix {
    type Output = BlockMatrix;

    fn mul(self, b: &BlockMatrix) -> BlockMatrix {
        assert_eq!(N, b.len());
        let n = self.len();
        let mut res = blockmatrix![0; n];

        for i in 0..n {
            let mut x = self[i];
            let mut k = 0;
            while x != 0 {
                if (x & 1) != 0 {
                    res[i] ^= b[k];
                }
                x >>= 1;
                k += 1;
            }
        }

        res
    }
}

impl<'a> Mul<&BlockMatrixTranspose<'a>> for &BlockMatrix {
    type Output = BlockMatrix;

    fn mul(self, b: &BlockMatrixTranspose<'a>) -> BlockMatrix {
        let n = self.len();
        assert!(n >= N);
        assert_eq!(N, b.borrowed.len());
        let mut res = blockmatrix![0; n];

        for i in 0..n {
            for j in 0..N {
                res[i] |= (((self[i] & b.borrowed[j]).count_ones() & 1) as u64) << j;
            }
        }

        res
    }
}

impl<'a> Mul<&BlockMatrix> for &BlockMatrixTranspose<'a> {
    type Output = BlockMatrix;

    fn mul(self, b: &BlockMatrix) -> BlockMatrix {
        let n = self.borrowed.len();
        assert_eq!(b.len(), n);
        let mut res = blockmatrix![0; N];

        for i in 0..n {
            let mut x = self.borrowed[i];
            let mut k = 0;
            while x != 0 {
                if (x & 1) != 0 {
                    res[k] ^= b[i];
                }
                x >>= 1;
                k += 1;
            }
        }

        res
    }
}
