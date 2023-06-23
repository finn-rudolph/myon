use std::ops::Mul;

pub struct CscMatrix {
    start: Vec<u32>,
    ones: Vec<u32>,
}

struct CsrMatrix<'a> {
    borrowed: &'a CscMatrix,
}

pub struct BlockMatrix(Vec<u64>);
struct TransposedBlockMatrix<'a> {
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

impl Mul<BlockMatrix> for CscMatrix {}

impl Mul<BlockMatrix> for CsrMatrix {}

impl Mul<BlockMatrix> for BlockMatrix {}

impl Mul<TransposedBlockMatrix> for BlockMatrix {}

impl Mul<BlockMatrix> for TransposedBlockMatrix {}
