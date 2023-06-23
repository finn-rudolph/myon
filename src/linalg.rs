use std::convert::From;
use std::ops::Index;
use std::ops::IndexMut;
use std::ops::Mul;

pub const N: usize = 64;

pub struct CscMatrix {
    end: Vec<u32>,
    ones: Vec<u32>,
}

pub struct CsrMatrix<'a> {
    borrowed: &'a CscMatrix,
}

#[derive(Clone)]
pub struct BlockMatrix(Vec<u64>);
pub struct TransposedBlockMatrix<'a> {
    borrowed: &'a BlockMatrix,
}

impl CscMatrix {
    pub fn new(end: Vec<u32>, ones: Vec<u32>) -> CscMatrix {
        CscMatrix { end, ones }
    }
    pub fn len(&self) -> usize {
        self.end.len()
    }

    pub fn transpose(&self) -> CsrMatrix {
        CsrMatrix { borrowed: self }
    }
}

impl BlockMatrix {
    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn swap(&mut self, i: usize, j: usize) {
        self.0.swap(i, j);
    }

    pub fn transpose(&self) -> TransposedBlockMatrix {
        TransposedBlockMatrix { borrowed: self }
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

macro_rules! blockmatrix {
    ( $x:expr; $n:expr ) => {
        BlockMatrix::from(vec![$x; $n])
    };
}

pub(crate) use blockmatrix;

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
        let n = self.end.len();
        let mut res = BlockMatrix::from(vec![0; n]);

        let mut j = 0;
        for i in 0..n {
            while j < self.end[i] {
                res[i] ^= b[self.ones[j as usize] as usize];
                j += 1;
            }
        }

        res
    }
}

impl<'a> Mul<&BlockMatrix> for &CsrMatrix<'a> {
    type Output = BlockMatrix;

    fn mul(self, b: &BlockMatrix) -> BlockMatrix {
        let n = self.borrowed.end.len();
        let mut res = BlockMatrix::from(vec![0; n]);

        let mut j = 0;
        for i in 0..n {
            while j < self.borrowed.end[i] {
                res[self.borrowed.ones[j as usize] as usize] ^= b[i];
                j += 1;
            }
        }

        res
    }
}

impl Mul<&BlockMatrix> for &BlockMatrix {
    type Output = BlockMatrix;

    fn mul(self, b: &BlockMatrix) -> BlockMatrix {
        let n = self.len();
        let mut res = BlockMatrix::from(vec![0; n]);

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

impl<'a> Mul<&TransposedBlockMatrix<'a>> for &BlockMatrix {
    type Output = BlockMatrix;

    fn mul(self, b: &TransposedBlockMatrix<'a>) -> BlockMatrix {
        let mut res = BlockMatrix::from(vec![0; N]);

        for i in 0..N {
            for j in 0..N {
                res[i] |= (((self[i] & b.borrowed[j]).count_ones() & 1) as u64) << j;
            }
        }

        res
    }
}

impl<'a> Mul<&BlockMatrix> for &TransposedBlockMatrix<'a> {
    type Output = BlockMatrix;

    fn mul(self, b: &BlockMatrix) -> BlockMatrix {
        let n = self.borrowed.len();
        let mut res = BlockMatrix::from(vec![0; N]);

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
