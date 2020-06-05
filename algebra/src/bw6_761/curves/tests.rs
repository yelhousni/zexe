use std::time::Instant;
use algebra_core::{
    test_rng, AffineCurve, Field, One, PairingEngine, PrimeField, ProjectiveCurve,
    msm::VariableBaseMSM, UniformRand, Zero,
};
use rand::{Rng, SeedableRng};
use rand_xorshift::XorShiftRng;

use crate::bw6_761::*;

use crate::tests::{curves::*, groups::*};

fn naive_var_base_msm<G: AffineCurve>(
    bases: &[G],
    scalars: &[<G::ScalarField as PrimeField>::BigInt],
) -> G::Projective {
    let mut acc = G::Projective::zero();

    for (base, scalar) in bases.iter().zip(scalars.iter()) {
        acc += &base.mul(*scalar);
    }
    acc
}

#[test]
fn test_g1_projective_curve() {
    curve_tests::<G1Projective>();

    sw_tests::<g1::Parameters>();
}

#[test]
fn test_g1_projective_group() {
    let mut rng = test_rng();
    let a: G1Projective = rng.gen();
    let b: G1Projective = rng.gen();
    group_test(a, b);
}

#[test]
fn test_g1_generator() {
    let generator = G1Affine::prime_subgroup_generator();
    assert!(generator.is_on_curve());
    assert!(generator.is_in_correct_subgroup_assuming_on_curve());
}

#[test]
fn test_g2_projective_curve() {
    curve_tests::<G2Projective>();

    sw_tests::<g2::Parameters>();
}

#[test]
fn test_g2_projective_group() {
    let mut rng = test_rng();
    let a: G2Projective = rng.gen();
    let b: G2Projective = rng.gen();
    group_test(a, b);
}

#[test]
fn test_g2_generator() {
    let generator = G2Affine::prime_subgroup_generator();
    assert!(generator.is_on_curve());
    assert!(generator.is_in_correct_subgroup_assuming_on_curve());
}

#[test]
fn test_bilinearity() {
    let mut rng = test_rng();
    let a: G1Projective = rng.gen();
    let b: G2Projective = rng.gen();
    let s: Fr = rng.gen();

    let sa = a.mul(s);
    let sb = b.mul(s);

    let ans1 = BW6_761::pairing(sa, b);
    let ans2 = BW6_761::pairing(a, sb);
    let ans3 = BW6_761::pairing(a, b).pow(s.into_repr());

    assert_eq!(ans1, ans2);
    assert_eq!(ans2, ans3);

    assert_ne!(ans1, Fq6::one());
    assert_ne!(ans2, Fq6::one());
    assert_ne!(ans3, Fq6::one());

    assert_eq!(ans1.pow(Fr::characteristic()), Fq6::one());
    assert_eq!(ans2.pow(Fr::characteristic()), Fq6::one());
    assert_eq!(ans3.pow(Fr::characteristic()), Fq6::one());
}

#[test]
fn test_g1_msm() {
    const SAMPLES: usize = 1 << 10;

    let mut rng = XorShiftRng::seed_from_u64(234872845u64);

    let v = (0..SAMPLES)
        .map(|_| Fr::rand(&mut rng).into_repr())
        .collect::<Vec<_>>();
    let g = (0..SAMPLES)
        .map(|_| G1Projective::rand(&mut rng).into_affine())
        .collect::<Vec<_>>();

    let naive = naive_var_base_msm(g.as_slice(), v.as_slice());
    let fast = VariableBaseMSM::multi_scalar_mul(g.as_slice(), v.as_slice());

    assert_eq!(naive.into_affine(), fast.into_affine());
}

#[test]
fn test_g1_msm_unequal_numbers() {
    const SAMPLES: usize = 1 << 10;

    let mut rng = XorShiftRng::seed_from_u64(234872845u64);

    let v = (0..SAMPLES - 1)
        .map(|_| Fr::rand(&mut rng).into_repr())
        .collect::<Vec<_>>();
    let g = (0..SAMPLES)
        .map(|_| G1Projective::rand(&mut rng).into_affine())
        .collect::<Vec<_>>();

    let naive = naive_var_base_msm(g.as_slice(), v.as_slice());
    let fast = VariableBaseMSM::multi_scalar_mul(g.as_slice(), v.as_slice());

    assert_eq!(naive.into_affine(), fast.into_affine());
}

#[test]
fn test_g1_msm_single() {
    let mut rng = XorShiftRng::seed_from_u64(234872845u64);

    let scalar = Fr::rand(&mut rng).into_repr();
    let base = G1Projective::rand(&mut rng).into_affine();

    let start = Instant::now();
    let mul = base.mul(scalar);
    let duration = start.elapsed();
    println!("G1MUL: {:?}", duration);

    let start = Instant::now();
    let msm = VariableBaseMSM::multi_scalar_mul(&[base], &[scalar]);
    let duration = start.elapsed();
    println!("G1MSM-1: {:?}", duration);

    assert_eq!(mul.into_affine(), msm.into_affine());
}

#[test]
fn test_g2_msm() {
    const SAMPLES: usize = 1 << 10;

    let mut rng = XorShiftRng::seed_from_u64(234872845u64);

    let v = (0..SAMPLES)
        .map(|_| Fr::rand(&mut rng).into_repr())
        .collect::<Vec<_>>();
    let g = (0..SAMPLES)
        .map(|_| G2Projective::rand(&mut rng).into_affine())
        .collect::<Vec<_>>();

    let naive = naive_var_base_msm(g.as_slice(), v.as_slice());
    let fast = VariableBaseMSM::multi_scalar_mul(g.as_slice(), v.as_slice());

    assert_eq!(naive.into_affine(), fast.into_affine());
}

#[test]
fn test_g2_msm_unequal_numbers() {
    const SAMPLES: usize = 1 << 10;

    let mut rng = XorShiftRng::seed_from_u64(234872845u64);

    let v = (0..SAMPLES - 1)
        .map(|_| Fr::rand(&mut rng).into_repr())
        .collect::<Vec<_>>();
    let g = (0..SAMPLES)
        .map(|_| G2Projective::rand(&mut rng).into_affine())
        .collect::<Vec<_>>();

    let naive = naive_var_base_msm(g.as_slice(), v.as_slice());
    let fast = VariableBaseMSM::multi_scalar_mul(g.as_slice(), v.as_slice());

    assert_eq!(naive.into_affine(), fast.into_affine());
}

#[test]
fn test_g2_msm_single() {
    let mut rng = XorShiftRng::seed_from_u64(234872845u64);

    let scalar = Fr::rand(&mut rng).into_repr();
    let base = G2Projective::rand(&mut rng).into_affine();

    let mul = base.mul(scalar);
    let msm = VariableBaseMSM::multi_scalar_mul(&[base], &[scalar]);

    assert_eq!(mul.into_affine(), msm.into_affine());
}
