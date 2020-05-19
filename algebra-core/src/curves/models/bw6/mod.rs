use crate::{
    curves::{
        models::{ModelParameters, SWModelParameters},
        PairingEngine,
    },
    fields::{
        fp3::{Fp3, Fp3Parameters},
        fp6_2over3::{Fp6, Fp6Parameters},
        BitIterator, Field, PrimeField, SquareRootField,
    },
    One, Zero,
};

use core::marker::PhantomData;

pub mod g1;
pub mod g2;

use self::g2::{AteAdditionCoefficients, AteDoubleCoefficients, G2ProjectiveExtended};
pub use self::{
    g1::{G1Affine, G1Prepared, G1Projective},
    g2::{G2Affine, G2Prepared, G2Projective},
};

pub type GT<P> = Fp6<P>;

pub trait BW6Parameters: 'static {
    const X: <Self::Fp as PrimeField>::BigInt;
    const X_IS_NEGATIVE: bool;
    const TWIST: Fp3<Self::Fp3Params>;
    const TWIST_COEFF_A: Fp3<Self::Fp3Params>;
    const ATE_LOOP_COUNT: &'static [u64];
    const ATE_IS_LOOP_COUNT_NEG: bool;
    // const FINAL_EXPONENT_LAST_CHUNK_1: <Self::Fp as PrimeField>::BigInt;
    // const FINAL_EXPONENT_LAST_CHUNK_W0_IS_NEG: bool;
    // const FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0: <Self::Fp as PrimeField>::BigInt;
    type Fp: PrimeField + SquareRootField + Into<<Self::Fp as PrimeField>::BigInt>;
    type Fp3Params: Fp3Parameters<Fp = Self::Fp>;
    type Fp6Params: Fp6Parameters<Fp3Params = Self::Fp3Params>;
    type G1Parameters: SWModelParameters<BaseField = Self::Fp>;
    type G2Parameters: SWModelParameters<
        BaseField = Fp3<Self::Fp3Params>,
        ScalarField = <Self::G1Parameters as ModelParameters>::ScalarField,
    >;
}

#[derive(Derivative)]
#[derivative(Copy, Clone, PartialEq, Eq, Debug, Hash)]
pub struct BW6<P: BW6Parameters>(PhantomData<fn() -> P>);

impl<P: BW6Parameters> BW6<P> {
    fn doubling_step_for_flipped_miller_loop(
        r: &G2ProjectiveExtended<P>,
    ) -> (G2ProjectiveExtended<P>, AteDoubleCoefficients<P>) {
        let a = r.t.square();
        let b = r.x.square();
        let c = r.y.square();
        let d = c.square();
        let e = (r.x + &c).square() - &b - &d;
        let f = (b + &b + &b) + &(P::TWIST_COEFF_A * &a);
        let g = f.square();

        let d_eight = d.double().double().double();

        let x = -(e + &e + &e + &e) + &g;
        let y = -d_eight + &(f * &(e + &e - &x));
        let z = (r.y + &r.z).square() - &c - &r.z.square();
        let t = z.square();

        let r2 = G2ProjectiveExtended { x, y, z, t };
        let coeff = AteDoubleCoefficients {
            c_h: (r2.z + &r.t).square() - &r2.t - &a,
            c_4c: c + &c + &c + &c,
            c_j: (f + &r.t).square() - &g - &a,
            c_l: (f + &r.x).square() - &g - &b,
        };

        (r2, coeff)
    }

    fn mixed_addition_step_for_flipped_miller_loop(
        x: &Fp3<P::Fp3Params>,
        y: &Fp3<P::Fp3Params>,
        r: &G2ProjectiveExtended<P>,
    ) -> (G2ProjectiveExtended<P>, AteAdditionCoefficients<P>) {
        let a = y.square();
        let b = r.t * x;
        let d = ((r.z + y).square() - &a - &r.t) * &r.t;
        let h = b - &r.x;
        let i = h.square();
        let e = i + &i + &i + &i;
        let j = h * &e;
        let v = r.x * &e;
        let l1 = d - &(r.y + &r.y);

        let x = l1.square() - &j - &(v + &v);
        let y = l1 * &(v - &x) - &(j * &(r.y + &r.y));
        let z = (r.z + &h).square() - &r.t - &i;
        let t = z.square();

        let r2 = G2ProjectiveExtended { x, y, z, t };
        let coeff = AteAdditionCoefficients { c_l1: l1, c_rz: z };

        (r2, coeff)
    }

    pub fn ate_miller_loop(p: &G1Prepared<P>, q: &G2Prepared<P>) -> Fp6<P::Fp6Params> {
        let l1_coeff = Fp3::new(p.x, P::Fp::zero(), P::Fp::zero()) - &q.x_over_twist;

        let mut f = <Fp6<P::Fp6Params>>::one();

        let mut dbl_idx: usize = 0;
        let mut add_idx: usize = 0;

        let mut found_one = false;

        for bit in BitIterator::new(P::ATE_LOOP_COUNT) {
            // code below gets executed for all bits (EXCEPT the MSB itself) of
            // bw6_param_p (skipping leading zeros) in MSB to LSB order
            if !found_one && bit {
                found_one = true;
                continue;
            } else if !found_one {
                continue;
            }

            let dc = &q.double_coefficients[dbl_idx];
            dbl_idx += 1;

            let g_rr_at_p = Fp6::new(
                -dc.c_4c - &(dc.c_j * &p.x_twist) + &dc.c_l,
                dc.c_h * &p.y_twist,
            );

            f = f.square() * &g_rr_at_p;

            if bit {
                let ac = &q.addition_coefficients[add_idx];
                add_idx += 1;

                let g_rq_at_p = Fp6::new(
                    ac.c_rz * &p.y_twist,
                    -(q.y_over_twist * &ac.c_rz + &(l1_coeff * &ac.c_l1)),
                );
                f *= &g_rq_at_p;
            }
        }

        if P::ATE_IS_LOOP_COUNT_NEG {
            let ac = &q.addition_coefficients[add_idx];

            let g_rnegr_at_p = Fp6::new(
                ac.c_rz * &p.y_twist,
                -(q.y_over_twist * &ac.c_rz + &(l1_coeff * &ac.c_l1)),
            );
            f = (f * &g_rnegr_at_p).inverse().unwrap();
        }

        f
    }

    fn exp_by_x(mut f: Fp6<P::Fp6Params>) -> Fp6<P::Fp6Params> {
        f = f.cyclotomic_exp(&P::X);
        if P::X_IS_NEGATIVE {
            f.conjugate();
        }
        f
    }

    pub fn final_exponentiation(f: &Fp6<P::Fp6Params>) -> Fp6<P::Fp6Params> {
        // final_exponent = (q^6-1)/r
        //                = (q^3-1)*(q+1) * (q^2-q+1)/r
        //                =    easy_part  *  hard_part

        let f_inv = f.inverse().unwrap();
        let f_to_first_chunk = Self::final_exponentiation_first_chunk(f, &f_inv);
        Self::final_exponentiation_last_chunk(&f_to_first_chunk)
    }

    fn final_exponentiation_first_chunk(
        f: &Fp6<P::Fp6Params>,
        f_inv: &Fp6<P::Fp6Params>,
    ) -> Fp6<P::Fp6Params> {
        // (q^3-1)*(q+1)

        // f_q3 = f^(q^3)
        let mut f_q3 = f.clone();
        f_q3.frobenius_map(3);
        // f_q3_over_f = f^(q^3-1)
        let f_q3_over_f = f_q3 * f_inv;
        // alpha = f^((q^3-1) * q)
        let mut alpha = f_q3_over_f.clone();
        alpha.frobenius_map(1);
        // beta = f^((q^3-1)*(q+1)
        alpha * &f_q3_over_f
    }

    fn final_exponentiation_last_chunk(f: &Fp6<P::Fp6Params>) -> Fp6<P::Fp6Params> {
        // hard_part
        // From https://eprint.iacr.org/2020/351.pdf, Alg.6
        // R0(x) := (-103*x^7 + 70*x^6 + 269*x^5 - 197*x^4 - 314*x^3 - 73*x^2 - 263*x - 220)
        // R1(x) := (103*x^9 - 276*x^8 + 77*x^7 + 492*x^6 - 445*x^5 - 65*x^4 + 452*x^3 - 181*x^2 + 34*x + 229)
        // f ^ R0(u) * (f ^ q) ^ R1(u) in a 2-NAF multi-exp fashion.

        // steps 1,2,3
        let f0 = f.clone();
        let mut f0p = f0;
        f0p.frobenius_map(1);
        let f1 = Self::exp_by_x(f0);
        let mut f1p = f1;
        f1p.frobenius_map(1);
        let f2 = Self::exp_by_x(f1);
        let mut f2p = f2;
        f2p.frobenius_map(1);
        let f3 = Self::exp_by_x(f2);
        let mut f3p = f3;
        f3p.frobenius_map(1);
        let f4 = Self::exp_by_x(f3);
        let mut f4p = f4;
        f4p.frobenius_map(1);
        let f5 = Self::exp_by_x(f4);
        let mut f5p = f5;
        f5p.frobenius_map(1);
        let f6 = Self::exp_by_x(f5);
        let mut f6p = f6;
        f6p.frobenius_map(1);
        let f7 = Self::exp_by_x(f6);
        let mut f7p = f7;
        f7p.frobenius_map(1);

        // step 4
        let f8p = Self::exp_by_x(f7p);
        let f9p = Self::exp_by_x(f8p);

        // step 5
        let mut f5p_p3 = f5p;
        f5p_p3.frobenius_map(3);
        let result1 = f3p * &f6p * &f5p_p3;

        // step 6
        // TODO: cyclotomic_square?
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

impl<P: BW6Parameters> PairingEngine for BW6<P> {
    type Fr = <P::G1Parameters as ModelParameters>::ScalarField;
    type G1Projective = G1Projective<P>;
    type G1Affine = G1Affine<P>;
    type G1Prepared = G1Prepared<P>;
    type G2Projective = G2Projective<P>;
    type G2Affine = G2Affine<P>;
    type G2Prepared = G2Prepared<P>;
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
