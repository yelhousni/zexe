use crate::bls12_377::*;
use algebra_core::{
    biginteger::BigInteger384,
    field_new,
    curves::{
    bls12,
    bls12::{Bls12, Bls12Parameters, TwistType}}
};

pub mod g1;
pub mod g2;

#[cfg(test)]
mod tests;

pub struct Parameters;

impl Bls12Parameters for Parameters {
    const X: &'static [u64] = &[0x8508c00000000001];
    /// `x` is positive.
    const X_IS_NEGATIVE: bool = false;
    const TWIST_TYPE: TwistType = TwistType::D;
    const TWIST: Fq2 = field_new!(Fq2,
        field_new!(Fq, BigInteger384([0, 0, 0, 0, 0, 0])),
        field_new!(Fq, BigInteger384([
            202099033278250856u64,
            5854854902718660529u64,
            11492539364873682930u64,
            8885205928937022213u64,
            5545221690922665192u64,
            39800542322357402u64,
        ])),
    );
    type Fp = Fq;
    type Fp2Params = Fq2Parameters;
    type Fp6Params = Fq6Parameters;
    type Fp12Params = Fq12Parameters;
    type G1Parameters = g1::Parameters;
    type G2Parameters = g2::Parameters;
}

pub type Bls12_377 = Bls12<Parameters>;

pub type G1Affine = bls12::G1Affine<Parameters>;
pub type G1Projective = bls12::G1Projective<Parameters>;
pub type G2Affine = bls12::G2Affine<Parameters>;
pub type G2Projective = bls12::G2Projective<Parameters>;
