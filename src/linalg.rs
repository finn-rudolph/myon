use std::ops::Mul;

pub struct CscMatrix {
    end: Vec<u32>,
    ones: Vec<u32>,
}

pub struct CsrMatrix<'a> {
    borrowed: &'a CscMatrix,
}

pub struct BlockMatrix(Vec<u64>);
pub struct TransposedBlockMatrix<'a> {
    borrowed: &'a BlockMatrix,
}

impl CscMatrix {
    pub fn transpose(&self) -> CsrMatrix {
        CsrMatrix { borrowed: self }
    }
}

impl BlockMatrix {
    pub fn transpose(&self) -> TransposedBlockMatrix {
        TransposedBlockMatrix { borrowed: self }
    }
}

impl Mul<&BlockMatrix> for &CscMatrix {
    type Output = BlockMatrix;

    fn mul(self, b: &BlockMatrix) -> BlockMatrix {
        let n = self.end.len();
        let mut res = BlockMatrix(vec![0; n]);

        let mut j = 0;
        for i in 0..n {
            while j < self.end[i] {
                res.0[i] ^= b.0[self.ones[j as usize] as usize];
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
        let mut res = BlockMatrix(vec![0; n]);

        let mut j = 0;
        for i in 0..n {
            while j < self.borrowed.end[i] {
                res.0[self.borrowed.ones[j as usize] as usize] ^= b.0[i];
                j += 1;
            }
        }

        res
    }
}

impl Mul<&BlockMatrix> for &BlockMatrix {
    type Output = BlockMatrix;

    fn mul(self, b: &BlockMatrix) -> BlockMatrix {
        let n = self.0.len();
        let mut res = BlockMatrix(vec![0; n]);

        for i in 0..n {
            let mut x = self.0[i];
            let mut k = 0;
            while x != 0 {
                if (x & 1) != 0 {
                    res.0[i] ^= b.0[k];
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
        let mut res = BlockMatrix(vec![0; 64]);

        for i in 0..64 {
            for j in 0..64 {
                res.0[i] |= (((self.0[i] & b.borrowed.0[j]).count_ones() & 1) as u64) << j;
            }
        }

        res
    }
}

impl<'a> Mul<&BlockMatrix> for &TransposedBlockMatrix<'a> {
    type Output = BlockMatrix;

    fn mul(self, b: &BlockMatrix) -> BlockMatrix {
        let n = self.borrowed.0.len();
        let mut res = BlockMatrix(vec![0; 64]);

        for i in 0..n {
            let mut x = self.borrowed.0[i];
            let mut k = 0;
            while x != 0 {
                if (x & 1) != 0 {
                    res.0[k] ^= b.0[i];
                }
                x >>= 1;
                k += 1;
            }
        }

        res
    }
}
