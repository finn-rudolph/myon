use core::convert::From;
use core::ops::{Index, IndexMut, Mul};

use rand::{thread_rng, Rng};

pub const N: usize = 64;

// A BlockMatrix of length n filled with x can be created by block_matrix![x; n].
macro_rules! block_matrix {
    ( $x:expr; $n:expr ) => {
        BlockMatrix::from(vec![$x; $n])
    };
}

pub(crate) use block_matrix;

// Column-major sparse matrix storing for each column the ones' positions in a contiguous subsegment
// in 'ones'. The index after the last element of column i is end[i].
pub struct CscMatrix {
    num_rows: usize,
    end: Vec<u32>, // number of columns = end.len()
    ones: Vec<u32>,
}

pub struct CscMatrixTranspose<'a> {
    borrowed: &'a CscMatrix,
}

// A dense binary matrix storing each row as an N-bit integer.
#[repr(transparent)]
#[derive(Clone)]
pub struct BlockMatrix(Vec<u64>);

pub struct BlockMatrixTranspose<'a> {
    borrowed: &'a BlockMatrix,
}

impl CscMatrix {
    // Builds a CscMatrix, given for each column a Vec containing the numbers of rows with a 1.
    pub fn new(ones_positions: Vec<Vec<u32>>, num_rows: usize) -> CscMatrix {
        let mut a = CscMatrix {
            num_rows,
            end: vec![],
            ones: vec![],
        };
        for mut col in ones_positions {
            a.ones.append(&mut col);
            a.end.push(a.ones.len() as u32);
        }
        a
    }

    pub fn new_random(num_cols: usize, num_rows: usize, max_ones: usize) -> CscMatrix {
        let mut end: Vec<u32> = vec![0; 0];
        let mut ones: Vec<u32> = vec![0; 0];
        let mut used: Vec<bool> = vec![false; num_rows];

        let mut rng = thread_rng();

        // Choose the number of nonzero entries for each column at random, then generate the indices
        // of 1s at random, avoiding duplicates in a column.
        for _ in 0..num_cols {
            let weight = rng.gen_range(0..max_ones);
            for _ in 0..weight {
                let mut x = rng.gen_range(0..max_ones) as u32;
                while used[x as usize] {
                    x = rng.gen_range(0..max_ones) as u32;
                }
                ones.push(x);
                used[x as usize] = true;
            }
            end.push(ones.len() as u32);
            for i in 0..weight {
                used[ones[ones.len() - i - 1] as usize] = false;
            }
        }

        CscMatrix {
            num_rows,
            end,
            ones,
        }
    }

    pub fn num_cols(&self) -> usize {
        self.end.len()
    }

    pub fn num_rows(&self) -> usize {
        self.num_rows
    }

    // Returns a view on the transposed matrix. The view is tightly bound to the original CscMatrix
    // and is intended to be used only in composition with the '*'-Operator.
    pub fn transpose(&self) -> CscMatrixTranspose {
        CscMatrixTranspose { borrowed: self }
    }
}

impl BlockMatrix {
    pub fn new_random(n: usize) -> BlockMatrix {
        let mut a = block_matrix![0; n];
        thread_rng().fill(&mut a.as_mut()[..]);
        a
    }

    // Provides a lightweight view on the transposed matrix, which isn't intendend to be used
    // standalone, but as an argument to the '*'-Operator (on any side).
    pub fn transpose(&self) -> BlockMatrixTranspose {
        BlockMatrixTranspose { borrowed: self }
    }

    // Calculates the transpose explicity as a two-dimensional vector, in row-major format.
    pub fn explicit_transpose(&self) -> Vec<Vec<u64>> {
        let n_words = (self.as_ref().len() + N - 1) / N;
        let mut res: Vec<Vec<u64>> = vec![vec![0; n_words]; N];

        for i in 0..self.as_ref().len() {
            for j in 0..N {
                res[j][i / N] |= ((self[i] >> j) & 1) << (i & (N - 1));
            }
        }

        res
    }

    pub fn is_symmetric(&self) -> bool {
        assert_eq!(self.as_ref().len(), N);
        for i in 0..N {
            for j in 0..N {
                if (self[i] >> j) & 1 != (self[j] >> i) & 1 {
                    return false;
                }
            }
        }
        true
    }
}

impl From<Vec<u64>> for BlockMatrix {
    fn from(x: Vec<u64>) -> Self {
        BlockMatrix(x)
    }
}

impl AsRef<Vec<u64>> for BlockMatrix {
    fn as_ref(&self) -> &Vec<u64> {
        &self.0
    }
}

impl AsMut<Vec<u64>> for BlockMatrix {
    fn as_mut(&mut self) -> &mut Vec<u64> {
        &mut self.0
    }
}

impl Index<usize> for BlockMatrix {
    type Output = u64;

    fn index(&self, i: usize) -> &u64 {
        &self.as_ref()[i]
    }
}

impl IndexMut<usize> for BlockMatrix {
    fn index_mut(&mut self, i: usize) -> &mut u64 {
        &mut self.as_mut()[i]
    }
}

// TODO: Optimize all code from here on (unroll, maybe use count trailing zeros to skip zeros in
//       matrix-vector-product, try to remove some branches).

impl Mul<&BlockMatrix> for &CscMatrix {
    type Output = BlockMatrix;

    fn mul(self, b: &BlockMatrix) -> BlockMatrix {
        let (n, m) = (self.num_cols(), self.num_rows());
        assert_eq!(n, b.as_ref().len());
        let mut res = block_matrix![0; m];

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
        let (n, m) = (self.borrowed.num_cols(), self.borrowed.num_rows());
        assert_eq!(m, b.as_ref().len());
        let mut res = block_matrix![0; n];

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
        assert_eq!(N, b.as_ref().len());
        let n = self.as_ref().len();
        let mut res = block_matrix![0; n];

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
        let n = self.as_ref().len();
        assert!(n >= N);
        assert_eq!(N, b.borrowed.as_ref().len());
        let mut res = block_matrix![0; n];

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
        let n = self.borrowed.as_ref().len();
        assert_eq!(b.as_ref().len(), n);
        let mut res = block_matrix![0; N];

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
