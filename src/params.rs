use rug::Integer;

pub const MAX_DEGREE: usize = Params::PARAM_TABLE[Params::PARAM_TABLE.len() - 1]
    .1
    .polynomial_degree;

#[derive(Clone, Copy)]
pub struct Params {
    pub rational_base_size: usize,
    pub algebraic_base_size: usize,
    pub quad_char_base_size: usize,
    pub polynomial_degree: usize,
    pub sieve_array_size: usize,
    pub rational_fudge: i8,
    pub algebraic_fudge: i8,
}

impl Params {
    pub const PARAM_TABLE: [(u32, Params); 1] = [(
        128,
        Params {
            rational_base_size: 500,
            algebraic_base_size: 500,
            quad_char_base_size: 100,
            polynomial_degree: 3,
            sieve_array_size: 10000,
            rational_fudge: 20,
            algebraic_fudge: 40,
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
