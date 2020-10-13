use crate::{
    curves::{
        models::{ModelParameters, SWModelParameters},
        short_weierstrass_jacobian::{GroupAffine, GroupProjective},
        PairingEngine,
    },
    fields::{
        fp3::{Fp3, Fp3Parameters},
        fp6_2over3::{Fp6, Fp6Parameters},
        BitIteratorBE, Field, PrimeField, SquareRootField,
    },
    One,
};

use core::marker::PhantomData;

pub type G1Affine<P> = GroupAffine<<P as MNT6Parameters>::G1Parameters>;
pub type G1Projective<P> = GroupProjective<<P as MNT6Parameters>::G1Parameters>;
pub type G2Affine<P> = GroupAffine<<P as MNT6Parameters>::G2Parameters>;
pub type G2Projective<P> = GroupProjective<<P as MNT6Parameters>::G2Parameters>;

pub type GT<P> = Fp6<P>;

pub trait MNT6Parameters: 'static {
    const TWIST: Fp3<Self::Fp3Params>;
    const TWIST_COEFF_A: Fp3<Self::Fp3Params>;
    const ATE_LOOP_COUNT: &'static [u64];
    const FINAL_EXPONENT_LAST_CHUNK_1: <Self::Fp as PrimeField>::BigInt;
    const FINAL_EXPONENT_LAST_CHUNK_W0_IS_NEG: bool;
    const FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0: <Self::Fp as PrimeField>::BigInt;
    type Fp: PrimeField + SquareRootField + Into<<Self::Fp as PrimeField>::BigInt>;
    type Fr: PrimeField + SquareRootField + Into<<Self::Fr as PrimeField>::BigInt>;
    type Fp3Params: Fp3Parameters<Fp = Self::Fp>;
    type Fp6Params: Fp6Parameters<Fp3Params = Self::Fp3Params>;
    type G1Parameters: SWModelParameters<BaseField = Self::Fp, ScalarField = Self::Fr>;
    type G2Parameters: SWModelParameters<
        BaseField = Fp3<Self::Fp3Params>,
        ScalarField = <Self::G1Parameters as ModelParameters>::ScalarField,
    >;
}

#[derive(Derivative)]
#[derivative(Copy, Clone, PartialEq, Eq, Debug, Hash)]
pub struct MNT6<P: MNT6Parameters>(PhantomData<fn() -> P>);

impl<P: MNT6Parameters> MNT6<P> {
    fn ate_miller_loop(p: &G1Affine<P>, q: &G2Affine<P>) -> Fp6<P::Fp6Params> {
        let px = p.x;
        let py = p.y;
        let qx = q.x;
        let qy = q.y;
        let mut py_twist_squared = P::TWIST.square();
        py_twist_squared.mul_assign_by_fp(&py);

        let mut old_rx;
        let mut old_ry;
        let mut rx = qx;
        let mut ry = qy;
        let mut f = Fp6::one();

        // The for loop is executed for all bits (EXCEPT the MSB itself) of
        // mnt6_param_p (skipping leading zeros) in MSB to LSB order
        for bit in BitIteratorBE::without_leading_zeros(P::ATE_LOOP_COUNT).skip(1) {
            old_rx = rx;
            old_ry = ry;

            let old_rx_square = old_rx.square();
            let old_rx_square_3 = old_rx_square.double() + &old_rx_square;
            let old_rx_square_3_a = old_rx_square_3 + &P::TWIST_COEFF_A;
            let old_ry_double_inverse = old_ry.double().inverse().unwrap();

            let gamma = old_rx_square_3_a * &old_ry_double_inverse;
            let gamma_twist = gamma * &P::TWIST;
            let gamma_old_rx = gamma * &old_rx;
            let mut gamma_twist_px = gamma_twist;
            gamma_twist_px.mul_assign_by_fp(&px);

            let x = py_twist_squared;
            let y = gamma_old_rx - &old_ry - &gamma_twist_px;
            let ell_rr_at_p = Fp6::new(x, y);

            rx = gamma.square() - &old_rx.double();
            ry = gamma * &(old_rx - &rx) - &old_ry;
            f = f.square() * &ell_rr_at_p;

            if bit {
                old_rx = rx;
                old_ry = ry;

                let gamma = (old_ry - &qy) * &((old_rx - &qx).inverse().unwrap());
                let gamma_twist = gamma * &P::TWIST;
                let gamma_qx = gamma * &qx;
                let mut gamma_twist_px = gamma_twist;
                gamma_twist_px.mul_assign_by_fp(&px);

                let x = py_twist_squared;
                let y = gamma_qx - &qy - &gamma_twist_px;
                let ell_rq_at_p = Fp6::new(x, y);

                rx = gamma.square() - &old_rx - &qx;
                ry = gamma * &(old_rx - &rx) - &old_ry;
                f = f * &ell_rq_at_p;
            }
        }
        f
    }

    pub fn final_exponentiation(value: &Fp6<P::Fp6Params>) -> GT<P::Fp6Params> {
        let value_inv = value.inverse().unwrap();
        let value_to_first_chunk = Self::final_exponentiation_first_chunk(value, &value_inv);
        let value_inv_to_first_chunk = Self::final_exponentiation_first_chunk(&value_inv, value);
        Self::final_exponentiation_last_chunk(&value_to_first_chunk, &value_inv_to_first_chunk)
    }

    fn final_exponentiation_first_chunk(
        elt: &Fp6<P::Fp6Params>,
        elt_inv: &Fp6<P::Fp6Params>,
    ) -> Fp6<P::Fp6Params> {
        // (q^3-1)*(q+1)

        // elt_q3 = elt^(q^3)
        let mut elt_q3 = elt.clone();
        elt_q3.conjugate();
        // elt_q3_over_elt = elt^(q^3-1)
        let elt_q3_over_elt = elt_q3 * elt_inv;
        // alpha = elt^((q^3-1) * q)
        let mut alpha = elt_q3_over_elt.clone();
        alpha.frobenius_map(1);
        // beta = elt^((q^3-1)*(q+1)
        alpha * &elt_q3_over_elt
    }

    fn final_exponentiation_last_chunk(
        elt: &Fp6<P::Fp6Params>,
        elt_inv: &Fp6<P::Fp6Params>,
    ) -> Fp6<P::Fp6Params> {
        let elt_clone = elt.clone();
        let elt_inv_clone = elt_inv.clone();

        let mut elt_q = elt.clone();
        elt_q.frobenius_map(1);

        let w1_part = elt_q.cyclotomic_exp(&P::FINAL_EXPONENT_LAST_CHUNK_1);
        let w0_part = if P::FINAL_EXPONENT_LAST_CHUNK_W0_IS_NEG {
            elt_inv_clone.cyclotomic_exp(&P::FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0)
        } else {
            elt_clone.cyclotomic_exp(&P::FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0)
        };

        w1_part * &w0_part
    }
}

impl<P: MNT6Parameters> PairingEngine for MNT6<P> {
    type Fr = <P::G1Parameters as ModelParameters>::ScalarField;
    type G1Affine = G1Affine<P>;
    type G1Projective = G1Projective<P>;
    type G1Prepared = G1Affine<P>;
    type G2Affine = G2Affine<P>;
    type G2Projective = G2Projective<P>;
    type G2Prepared = G2Affine<P>;
    type Fq = P::Fp;
    type Fqe = Fp3<P::Fp3Params>;
    type Fqk = Fp6<P::Fp6Params>;

    fn miller_loop<'a, I>(i: I) -> Self::Fqk
    where
        I: IntoIterator<Item = &'a (Self::G1Prepared, Self::G2Prepared)>,
    {
        let mut result = Self::Fqk::one();
        for (p, q) in i {
            result *= &Self::ate_miller_loop(p, q);
        }
        result
    }

    fn final_exponentiation(r: &Self::Fqk) -> Option<Self::Fqk> {
        Some(Self::final_exponentiation(r))
    }
}
