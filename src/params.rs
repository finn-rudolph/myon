use rug::Integer;

pub const MAX_DEGREE: usize = Params::PARAM_TABLE[Params::PARAM_TABLE.len() - 1]
    .1
    .polynomial_degree;

pub const OVERSQUARENESS: usize = 19;

#[derive(Clone, Copy)]
pub struct Params {
    pub rational_base_size: usize,
    pub algebraic_base_size: usize,
    pub quad_char_base_size: usize,
    pub polynomial_degree: usize,
    pub sieve_array_size: usize,
    pub rational_fudge: i8,
    pub algebraic_threshold: i8,
}

impl Params {
    pub const PARAM_TABLE: [(u32, Params); 1] = [(
        64,
        Params {
            rational_base_size: 500,
            algebraic_base_size: 500,
            quad_char_base_size: 50,
            polynomial_degree: 4,
            sieve_array_size: 20000,
            rational_fudge: 7,
            algebraic_threshold: 13,
        },
    )];

    pub fn new(n: &Integer) -> Params {
        let bits = n.significant_bits();

        for (bits_lim, params) in Params::PARAM_TABLE.iter().rev() {
            if bits <= *bits_lim {
                return *params;
            }
        }

        Params::PARAM_TABLE.last().unwrap().1
    }
}
