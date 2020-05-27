use algebra_core::{
    biginteger::{BigInteger768, BigInteger384},
    curves::models::{
        bw6::{BW6Parameters, BW6},
        SWModelParameters,
    },
    field_new,
    fields::FpParameters,
    Fp3,
};

use crate::bw6_761::{Fq, Fq3, Fq3Parameters, Fq6Parameters, FqParameters, Fr, FrParameters};

pub mod g1;
pub mod g2;

#[cfg(test)]
mod tests;

pub use self::{
    g1::{G1Affine, G1Prepared, G1Projective},
    g2::{G2Affine, G2Prepared, G2Projective},
};

pub type BW6_761 = BW6<Parameters>;

pub struct Parameters;

impl BW6Parameters for Parameters {
    const X: BigInteger768 = BigInteger768([0x8508c00000000001, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0]);
    /// `x` is positive.
    const X_IS_NEGATIVE: bool = false;
    const TWIST: Fp3<Self::Fp3Params> = field_new!(Fq3, FQ_ZERO, FQ_ONE, FQ_ZERO);
    // A coefficient of BW6-761 G2 =
    // ```
    // bw6_761_twist_coeff_a = bw6_761_Fq3(bw6_761_Fq::zero(), bw6_761_Fq::zero(),
    //                                  bw6_761_G1::coeff_a);
    //  = (ZERO, ZERO, A_COEFF);
    // ```
    #[rustfmt::skip]
    const TWIST_COEFF_A: Fp3<Self::Fp3Params> = field_new!(Fq3,
        FQ_ZERO,
        FQ_ZERO,
        g1::Parameters::COEFF_A,
    );
    const ATE_LOOP_COUNT: &'static [u64] = &[
        5078512267302010895,
        8121630887260454919,
        6361633214565880833,
        13069376661505132784,
        1923874642746257132,
        1574278065184431084,
    ];
    const ATE_IS_LOOP_COUNT_NEG: bool = false;
    type Fp = Fq;
    type Fp3Params = Fq3Parameters;
    type Fp6Params = Fq6Parameters;
    type G1Parameters = self::g1::Parameters;
    type G2Parameters = self::g2::Parameters;
}

pub const FQ_ZERO: Fq = field_new!(Fq, BigInteger768([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]));
pub const FQ_ONE: Fq = field_new!(Fq, FqParameters::R);
pub const FR_ZERO: Fr = field_new!(Fr, BigInteger384([0, 0, 0, 0, 0, 0]));
pub const FR_ONE: Fr = field_new!(Fr, FrParameters::R);
