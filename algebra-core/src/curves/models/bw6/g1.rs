use crate::{
    bytes::ToBytes,
    curves::{
        bw6::BW6Parameters,
        short_weierstrass_jacobian::{GroupAffine, GroupProjective},
        AffineCurve,
    },
    io::{Result as IoResult, Write},
    Fp3,
};

pub type G1Affine<P> = GroupAffine<<P as BW6Parameters>::G1Parameters>;
pub type G1Projective<P> = GroupProjective<<P as BW6Parameters>::G1Parameters>;

#[derive(Derivative)]
#[derivative(
    Clone(bound = "P: BW6Parameters"),
    Debug(bound = "P: BW6Parameters"),
    PartialEq(bound = "P: BW6Parameters"),
    Eq(bound = "P: BW6Parameters")
)]
pub struct G1Prepared<P: BW6Parameters> {
    pub x: P::Fp,
    pub y: P::Fp,
    pub x_twist: Fp3<P::Fp3Params>,
    pub y_twist: Fp3<P::Fp3Params>,
}

impl<P: BW6Parameters> From<G1Affine<P>> for G1Prepared<P> {
    fn from(g1: G1Affine<P>) -> Self {
        let mut x_twist = P::TWIST.clone();
        x_twist.mul_assign_by_fp(&g1.x);

        let mut y_twist = P::TWIST.clone();
        y_twist.mul_assign_by_fp(&g1.y);

        Self {
            x: g1.x,
            y: g1.y,
            x_twist,
            y_twist,
        }
    }
}

impl<P: BW6Parameters> Default for G1Prepared<P> {
    fn default() -> Self {
        Self::from(G1Affine::<P>::prime_subgroup_generator())
    }
}

impl<P: BW6Parameters> ToBytes for G1Prepared<P> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.x.write(&mut writer)?;
        self.y.write(&mut writer)?;
        self.x_twist.write(&mut writer)?;
        self.y_twist.write(&mut writer)
    }
}
