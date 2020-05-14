use crate::{
    biginteger::BigInteger768,
    curves::{models::SWModelParameters, PairingEngine},
    field_new,
    fields::{BitIterator, Field, FpParameters},
    One,
};

use crate::bw6_bis::{Fq, Fq3, Fq6, FqParameters, Fr};

pub mod g1;
pub use self::g1::{G1Affine, G1Projective};

pub mod g2;
pub use self::g2::{G2Affine, G2Projective};

#[cfg(test)]
mod tests;

pub type GT = Fq6;

#[derive(Copy, Clone, Debug)]
pub struct BW6_761;

impl PairingEngine for BW6_761 {
    type Fr = Fr;
    type G1Projective = G1Projective;
    type G1Affine = G1Affine;
    type G1Prepared = G1Affine;
    type G2Projective = G2Projective;
    type G2Affine = G2Affine;
    type G2Prepared = G2Affine;
    type Fq = Fq;
    type Fqe = Fq3;
    type Fqk = Fq6;

    fn miller_loop<'a, I>(i: I) -> Self::Fqk
    where
        I: IntoIterator<Item = &'a (Self::G1Prepared, Self::G2Prepared)>,
    {
        let mut result = Self::Fqk::one();
        for &(ref p, ref q) in i {
            result *= &BW6_761::ate_miller_loop(p, q);
        }
        result
    }

    fn final_exponentiation(r: &Self::Fqk) -> Option<Self::Fqk> {
        Some(BW6_761::final_exponentiation(r))
    }
}

impl BW6_761 {
    pub fn ate_pairing(p: &G1Affine, q: &G2Affine) -> GT {
        BW6_761::final_exponentiation(&BW6_761::ate_miller_loop(p, q))
    }

    fn ate_miller_loop(p: &G1Affine, q: &G2Affine) -> Fq6 {
        let px = p.x;
        let py = p.y;
        let qx = q.x;
        let qy = q.y;
        let mut py_twist_squared = TWIST.square();
        py_twist_squared.mul_assign_by_fp(&py);

        let mut old_rx;
        let mut old_ry;
        let mut rx = qx;
        let mut ry = qy;
        let mut f = Fq6::one();

        // The for loop is executed for all bits (EXCEPT the MSB itself) of
        // bw6_bis_param_p (skipping leading zeros) in MSB to LSB order
        let mut found_one = false;
        for bit in BitIterator::new(ATE_LOOP_COUNT) {
            if !found_one && bit {
                found_one = true;
                continue;
            } else if !found_one {
                continue;
            }

            old_rx = rx;
            old_ry = ry;

            let old_rx_square = old_rx.square();
            let old_rx_square_3 = old_rx_square.double() + &old_rx_square;
            let old_rx_square_3_a = old_rx_square_3 + &g2::Parameters::COEFF_A;
            let old_ry_double_inverse = old_ry.double().inverse().unwrap();

            let gamma = old_rx_square_3_a * &old_ry_double_inverse;
            let gamma_twist = gamma * &TWIST;
            let gamma_old_rx = gamma * &old_rx;
            let mut gamma_twist_px = gamma_twist;
            gamma_twist_px.mul_assign_by_fp(&px);

            let x = py_twist_squared;
            let y = gamma_old_rx - &old_ry - &gamma_twist_px;
            let ell_rr_at_p = Fq6::new(x, y);

            rx = gamma.square() - &old_rx.double();
            ry = gamma * &(old_rx - &rx) - &old_ry;
            f = f.square() * &ell_rr_at_p;

            if bit {
                old_rx = rx;
                old_ry = ry;

                let gamma = (old_ry - &qy) * &((old_rx - &qx).inverse().unwrap());
                let gamma_twist = gamma * &TWIST;
                let gamma_qx = gamma * &qx;
                let mut gamma_twist_px = gamma_twist;
                gamma_twist_px.mul_assign_by_fp(&px);

                let x = py_twist_squared;
                let y = gamma_qx - &qy - &gamma_twist_px;
                let ell_rq_at_p = Fq6::new(x, y);

                rx = gamma.square() - &old_rx - &qx;
                ry = gamma * &(old_rx - &rx) - &old_ry;
                f = f * &ell_rq_at_p;
            }
        }
        f
    }

    fn exp_by_x(elt: &Fq6) -> Fq6 {
        let mut elt_x = elt.clone();
        elt_x = elt_x.cyclotomic_exp(&SEED);
        if SEED_NEG {
            elt_x.conjugate();
        }
        elt_x
    }

    fn final_exponentiation(value: &Fq6) -> GT {
        let value_inv = value.inverse().unwrap();
        let value_to_first_chunk = BW6_761::final_exponentiation_first(value, &value_inv);
        BW6_761::final_exponentiation_last(&value_to_first_chunk)
    }

    fn final_exponentiation_first(elt: &Fq6, elt_inv: &Fq6) -> Fq6 {
        // (q^3-1)*(q+1)

        // elt_q3 = elt^(q^3)
        let mut elt_q3 = elt.clone();
        elt_q3.frobenius_map(3);
        // elt_q3_over_elt = elt^(q^3-1)
        let elt_q3_over_elt = elt_q3 * elt_inv;
        // alpha = elt^((q^3-1) * q)
        let mut alpha = elt_q3_over_elt.clone();
        alpha.frobenius_map(1);
        // beta = elt^((q^3-1)*(q+1)
        alpha * &elt_q3_over_elt
    }

    fn final_exponentiation_last(elt: &Fq6) -> Fq6 {

        /* Algorithm 6 in https://eprint.iacr.org/2020/351.pdf
         *
         * R0(x) := (-103*x^7 + 70*x^6 + 269*x^5 - 197*x^4 - 314*x^3 - 73*x^2 - 263*x - 220)
         * R1(x) := (103*x^9 - 276*x^8 + 77*x^7 + 492*x^6 - 445*x^5 - 65*x^4 + 452*x^3 - 181*x^2 + 34*x + 229)
         *
         * elt ^ R0(SEED) * (elt ^ q) ^ R1(SEED) where elt^q is a Frobenius power in Fq6.
         */

        // steps 1,2,3
        let f0 = elt.clone();
        let mut f0p = f0;
        f0p.frobenius_map(1);
        let f1 = BW6_761::exp_by_x(&f0);
        let mut f1p = f1;
        f1p.frobenius_map(1);
        let f2 = BW6_761::exp_by_x(&f1);
        let mut f2p = f2;
        f2p.frobenius_map(1);
        let f3 = BW6_761::exp_by_x(&f2);
        let mut f3p = f3;
        f3p.frobenius_map(1);
        let f4 = BW6_761::exp_by_x(&f3);
        let mut f4p = f4;
        f4p.frobenius_map(1);
        let f5 = BW6_761::exp_by_x(&f4);
        let mut f5p = f5;
        f5p.frobenius_map(1);
        let f6 = BW6_761::exp_by_x(&f5);
        let mut f6p = f6;
        f6p.frobenius_map(1);
        let f7 = BW6_761::exp_by_x(&f6);
        let mut f7p = f7;
        f7p.frobenius_map(1);

        // step 4
        let f8p = BW6_761::exp_by_x(&f7p);
        let f9p = BW6_761::exp_by_x(&f8p);

        // step 5
        let mut f5p_p3 = f5p;
        f5p_p3.frobenius_map(3);
        let result1 = f3p * &f6p * &f5p_p3;

        // step 6
        let result2 = result1.square();
        let f4_2p = f4 * &f2p;
        let mut tmp1_p3 = f0 * &f1 * &f3 * &f4_2p * &f8p;
        tmp1_p3.frobenius_map(3);
        let result3 = result2 * &f5 * &f0p * &tmp1_p3;

        // step 7
        let result4 = result3.square();
        let mut f7_p3 = f7;
        f7_p3.frobenius_map(3);
        let result5 = result4 * &f9p * &f7_p3;

        // step 8
        let result6 = result5.square();
        let f2_4p = f2 * &f4p;
        let f4_2p_5p = f4_2p * &f5p;
        let mut tmp2_p3 = f2_4p * &f3 * &f3p;
        tmp2_p3.frobenius_map(3);
        let result7 = result6 * &f4_2p_5p * &f6 * &f7p * &tmp2_p3;

        // step 9
        let result8 = result7.square();
        let mut tmp3_p3 = f0p * &f9p;
        tmp3_p3.frobenius_map(3);
        let result9 = result8 * &f0 * &f7 * &f1p * &tmp3_p3;

        // step 10
        let result10 = result9.square();
        let f6p_8p = f6p * &f8p;
        let f5_7p = f5 * &f7p;
        let mut tmp4_p3 = f6p_8p;
        tmp4_p3.frobenius_map(3);
        let result11 = result10 * &f5_7p * &f2p * &tmp4_p3;

        // step 11
        let result12 = result11.square();
        let f3_6 = f3 * &f6;
        let f1_7 = f1 * &f7;
        let mut tmp5_p3 = f1_7 * &f2;
        tmp5_p3.frobenius_map(3);
        let result13 = result12 * &f3_6 * &f9p * &tmp5_p3;

        // step 12
        let result14 = result13.square();
        let mut tmp6_p3 = f4_2p * &f5_7p * &f6p_8p;
        tmp6_p3.frobenius_map(3);
        let result15 = result14 * &f0 * &f0p * &f3p * &f5p * &tmp6_p3;

        // step 13
        let result16 = result15.square();
        let mut tmp7_p3 = f3_6;
        tmp7_p3.frobenius_map(3);
        let result17 = result16 * &f1p * &tmp7_p3;

        // step 14
        let result18 = result17.square();
        let mut tmp8_p3 = f2_4p * &f4_2p_5p * &f9p;
        tmp8_p3.frobenius_map(3);
        let result19 = result18 * &f1_7 * &f5_7p * &f0p * &tmp8_p3;

        result19
    }
}

/// FQ_ZERO = 0
pub const FQ_ZERO: Fq = field_new!(Fq, BigInteger768([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]));

/// FQ_ONE = 1
pub const FQ_ONE: Fq = field_new!(Fq, FqParameters::R);

/// TWIST = (0, 1, 0)
pub const TWIST: Fq3 = field_new!(Fq3, FQ_ZERO, FQ_ONE, FQ_ZERO);

/// ATE_IS_LOOP_COUNT_NEG = false
pub const ATE_IS_LOOP_COUNT_NEG: bool = false;

/// ATE_LOOP_COUNT =
/// 3362637538168598222219435186298528655381674028954528064283340709388076588006567983337308081752755143497537638367247
pub const ATE_LOOP_COUNT: &'static [u64] = &[
    5078512267302010895,
    8121630887260454919,
    6361633214565880833,
    13069376661505132784,
    1923874642746257132,
    1574278065184431084,
];

/// SEED_NEG = false
pub const SEED_NEG: bool = false;

/// SEED =
/// 9586122913090633729
pub const SEED: BigInteger768 = BigInteger768([0x8508c00000000001, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0]);

/// FINAL_EXPONENT_LAST_CHUNK_W0_IS_NEG = false
pub const FINAL_EXPONENT_LAST_CHUNK_W0_IS_NEG: bool = false;

/// FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0 =
/// 2156695813352724974824326851054479880127610960548355747044807332080688727374737671308314095389122345740953981240668571898337613282699493372314698360451061276517306188376803619985090458895588556562724088277106828
pub const FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0: BigInteger768 = BigInteger768([
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

/// FINAL_EXPONENT_LAST_CHUNK_W1 =
/// 26642435879335816683987677701488073867751118270052650655942102502312977592501693353047140953112195348280268661194889
pub const FINAL_EXPONENT_LAST_CHUNK_W1: BigInteger768 = BigInteger768([
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
