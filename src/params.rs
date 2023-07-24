use rug::Integer;

pub const MAX_DEGREE: usize = Params::PARAM_TABLE[Params::PARAM_TABLE.len() - 1]
    .1
    .polynomial_degree as usize;

#[derive(Clone, Copy)]
pub struct Params {
    pub rational_base_size: usize,
    pub algebraic_base_size: usize,
    pub quad_char_base_size: usize,
    pub polynomial_degree: u32,
    pub sieve_array_size: usize,
    pub fudge: i8, // different fudge for rational and algebraic side?
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
            fudge: 20,
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
