use std::{
    char,
    ops::{Neg, Sub},
    str::FromStr,
};

use num_bigint::BigInt;
use num_traits::{FromPrimitive, One, ToPrimitive};
use pairings_derive::polynomial_inverse;

use crate::{
    field_element::FieldElement,
    finite_field::{FiniteField, NonExtendedField},
    polynomial::{AbstractCoeff, Polynomial},
};

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub struct Mod13;

impl FiniteField for Mod13 {
    type T = i64;

    fn modulus() -> Self::T {
        13
    }

    fn one() -> Self::T {
        1
    }

    fn zero() -> Self::T {
        0
    }
}

impl NonExtendedField for Mod13 {
    fn to_bits_be(s: i64) -> Vec<u8> {
        let mut max = i64::BITS - s.leading_zeros();

        let mut res = vec![0u8; max.try_into().unwrap()];
        for i in (0..max).rev() {
            res[i as usize] = (s >> i & 1) as u8;
        }

        res
    }

    fn to_uint(s: i64) -> Option<usize> {
        Some(s.try_into().unwrap())
    }

    fn from_uint(s: usize) -> Option<i64> {
        Some(s.try_into().unwrap())
    }
}

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub struct Mod13_2;

impl FiniteField for Mod13_2 {
    type T = Polynomial<F13>;

    fn modulus() -> Self::T {
        Polynomial::new(vec![F13::new(2), F13::new(0), F13::new(1)])
    }

    fn zero() -> Self::T {
        Polynomial::new(vec![F13::zero()])
    }

    fn one() -> Self::T {
        Polynomial::new(vec![F13::one()])
    }

    #[polynomial_inverse]
    fn inverse(value: &Self::T) -> Self::T;
}

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub struct Mod13_4;

impl FiniteField for Mod13_4 {
    type T = Polynomial<F13>;

    fn modulus() -> Self::T {
        Polynomial::new(vec![
            F13::new(2),
            F13::zero(),
            F13::zero(),
            F13::zero(),
            F13::one(),
        ])
    }

    fn zero() -> Self::T {
        Polynomial::new(vec![F13::zero()])
    }

    fn one() -> Self::T {
        Polynomial::new(vec![F13::one()])
    }

    #[polynomial_inverse]
    fn inverse(value: &Self::T) -> Self::T;
}

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub struct BN254;

impl FiniteField for BN254 {
    type T = BigInt;

    fn modulus() -> Self::T {
        BigInt::from_str(
            "21888242871839275222246405745257275088696311157297823662689037894645226208583",
        )
        .unwrap()
    }

    fn one() -> Self::T {
        BigInt::from_i64(1).unwrap()
    }

    fn zero() -> Self::T {
        BigInt::from_i64(0).unwrap()
    }
}

impl NonExtendedField for BN254 {
    fn to_bits_be(s: BigInt) -> Vec<u8> {
        let max = s.bits();

        let mut res = vec![0u8; max.try_into().unwrap()];
        for i in (0..max).rev() {
            res[i as usize] = (s.clone() >> i & BigInt::one()).to_u8().unwrap();
        }

        res
    }

    fn to_uint(s: BigInt) -> Option<usize> {
        s.to_usize()
    }

    fn from_uint(s: usize) -> Option<BigInt> {
        BigInt::from_usize(s)
    }
}

pub type F13 = FieldElement<Mod13>;
pub type F13_2 = FieldElement<Mod13_2>;
pub type F13_4 = FieldElement<Mod13_4>;
pub type FBN254 = FieldElement<BN254>;
