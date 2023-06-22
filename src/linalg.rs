pub struct SparseMatrix {
    end: Vec<u32>,
    ones: Vec<u32>,
}

pub struct BlockMatrix {
    mat: Vec<Vec<u64>>,
}

pub fn mul_ab(a: &BlockMatrix, b: &BlockMatrix) -> BlockMatrix {
    let m = a.mat.len();
    let n = a.mat[0].len();
    let mut res: BlockMatrix = BlockMatrix {
        mat: vec![vec![n as u64]; m],
    };

    for l in 0..m {
        for i in 0..n {
            for j in 0..m {
                let mut x = a.mat[i][j];
                let mut k = 0;
                while x != 0 {
                    if (x & 1) != 0 {
                        res.mat[l][i] ^= b.mat[l][(j << 6) + k];
                    }
                    x >>= 1;
                    k += 1;
                }
            }
        }
    }

    res
}

pub fn mul_atb(a: &BlockMatrix, b: &BlockMatrix) -> BlockMatrix {
    let m = a.mat.len();
    let n = a.mat[0].len();
    let mut res: BlockMatrix = BlockMatrix {
        mat: vec![vec![64 as u64]; m],
    };

    for l in 0..m {
        for i in 0..n {
            for j in 0..m {
                let mut x = a.mat[i][j];
                let mut k = 0;
                while x != 0 {
                    if (x & 1) != 0 {
                        res.mat[l][(j << 6) + k] ^= b.mat[l][i];
                    }
                    x >>= 1;
                    k += 1;
                }
            }
        }
    }

    res
}
