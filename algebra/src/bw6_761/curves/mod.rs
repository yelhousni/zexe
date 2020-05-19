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
    /*
    const FINAL_EXPONENT_LAST_CHUNK_1: BigInteger768 = BigInteger768([
        0x3de5800000000089,
        0x832ba4061000003b,
        0xc61c554757551c0c,
        0xc856a0853c9db94c,
        0x2c77d5ac34cb12ef,
        0xad1972339049ce76,
        0x0,
        0x0,
        0x0,
        0x0,
        0x0,
        0x0,
    ]);
    const FINAL_EXPONENT_LAST_CHUNK_W0_IS_NEG: bool = false;
    const FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0: BigInteger768 = BigInteger768([
        0x6f9440000000008c,
        0x1aff40fcf0000082,
        0x9521646d73808c51,
        0x3ba806d298c79fc5,
        0xb521a3d9309c6dd0,
        0x824cd7cfb1e8685a,
        0xa7f6ef02c228c497,
        0xa311dc0a5ef6ff10,
        0x96a147eaf584608d,
        0x828e2c6f9f4f1494,
        0x68f6427062e1b0b,
        0x0,
    ]);
    */
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
