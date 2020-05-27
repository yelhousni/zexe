use crate::{
    bw6_761::*,
    biginteger::BigInteger768 as BigInteger,
};

use algebra_core::curves::{
    bw6,
    bw6::{BW6, BW6Parameters, TwistType},
};

pub mod g1;
pub mod g2;

#[cfg(test)]
mod tests;

pub struct Parameters;

impl BW6Parameters for Parameters {
    const X: BigInteger = BigInteger([0x8508c00000000001, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0]);
    /// `x` is positive.
    const X_IS_NEGATIVE: bool = false;
    const ATE_LOOP_COUNT: &'static [u64] = &[
        0x467a80000000000f,
        0x70b5d44300000007,
        0x58490fb409869401,
        0xb55fc0d440cb48f0,
        0x1ab2f9cb6145aeec,
        0x15d8f58f3501dbec,
    ];
    const ATE_LOOP_COUNT_IS_NEGATIVE: bool = false;
    const TWIST_TYPE: TwistType = TwistType::M;
    type Fp = Fq;
    type Fp3Params = Fq3Parameters;
    type Fp6Params = Fq6Parameters;
    type G1Parameters = g1::Parameters;
    type G2Parameters = g2::Parameters;
}

pub type BW6_761 = BW6<Parameters>;

pub type G1Affine = bw6::G1Affine<Parameters>;
pub type G1Projective = bw6::G1Projective<Parameters>;
pub type G2Affine = bw6::G2Affine<Parameters>;
pub type G2Projective = bw6::G2Projective<Parameters>;
