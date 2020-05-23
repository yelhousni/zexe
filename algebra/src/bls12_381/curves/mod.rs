use crate::bls12_381::*;

use algebra_core::{
    field_new,
    biginteger::BigInteger384,
    curves::bls12::{Bls12, Bls12Parameters, TwistType}
};

pub mod g1;
pub mod g2;

#[cfg(test)]
mod tests;

pub use self::{
    g1::{G1Affine, G1Projective},
    g2::{G2Affine, G2Projective},
};

pub type Bls12_381 = Bls12<Parameters>;

pub struct Parameters;

impl Bls12Parameters for Parameters {
    const X: &'static [u64] = &[0xd201000000010000];
    const X_IS_NEGATIVE: bool = true;
    const TWIST_TYPE: TwistType = TwistType::M;
    // u+1
    const TWIST: Fq2 = field_new!(Fq2,
        field_new!(Fq, BigInteger384([
            8505329371266088957u64,
            17002214543764226050u64,
            6865905132761471162u64,
            8632934651105793861u64,
            6631298214892334189u64,
            1582556514881692819u64,
        ])),
        field_new!(Fq, BigInteger384([
            8505329371266088957u64,
            17002214543764226050u64,
            6865905132761471162u64,
            8632934651105793861u64,
            6631298214892334189u64,
            1582556514881692819u64,
        ])),
    );
    type Fp = Fq;
    type Fp2Params = Fq2Parameters;
    type Fp6Params = Fq6Parameters;
    type Fp12Params = Fq12Parameters;
    type G1Parameters = self::g1::Parameters;
    type G2Parameters = self::g2::Parameters;
}
