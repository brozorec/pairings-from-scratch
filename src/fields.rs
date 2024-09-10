use std::str::FromStr;

use derive_lib::polynomial_inverse;
use num_bigint::BigInt;
use num_traits::{FromPrimitive, One, ToPrimitive};

use crate::{
    field_element::FieldElement,
    finite_field::{FiniteField, NonExtendedField},
    polynomial::Polynomial,
};

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub struct Ff13;

impl FiniteField for Ff13 {
    type T = i16;

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

impl NonExtendedField for Ff13 {
    fn to_bits_be(s: i16) -> Vec<u8> {
        let max = i16::BITS - s.leading_zeros();

        let mut res = vec![0u8; max.try_into().unwrap()];
        for i in (0..max).rev() {
            res[i as usize] = (s >> i & 1) as u8;
        }

        res
    }

    fn to_uint(s: Self::T) -> Option<usize> {
        Some(s.try_into().unwrap())
    }

    fn from_uint(s: usize) -> Option<Self::T> {
        Some(s.try_into().unwrap())
    }
}

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub struct Ff13_2;

impl FiniteField for Ff13_2 {
    type T = Polynomial<Fe13>;

    fn modulus() -> Self::T {
        Polynomial::from(vec![2, 0, 1])
    }

    fn zero() -> Self::T {
        Polynomial::new(vec![Fe13::zero()])
    }

    fn one() -> Self::T {
        Polynomial::new(vec![Fe13::one()])
    }

    #[polynomial_inverse]
    fn inverse(value: &Self::T) -> Self::T;
}

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub struct Ff13_4;

impl FiniteField for Ff13_4 {
    type T = Polynomial<Fe13>;

    fn modulus() -> Self::T {
        Polynomial::from(vec![2, 0, 0, 0, 1])
    }

    fn zero() -> Self::T {
        Polynomial::new(vec![Fe13::zero()])
    }

    fn one() -> Self::T {
        Polynomial::new(vec![Fe13::one()])
    }

    #[polynomial_inverse]
    fn inverse(value: &Self::T) -> Self::T;
}

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub struct FfBn254;

impl FiniteField for FfBn254 {
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

impl NonExtendedField for FfBn254 {
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

pub type Fe13 = FieldElement<Ff13>;
pub type Fe13_2 = FieldElement<Ff13_2>;
pub type Fe13_4 = FieldElement<Ff13_4>;
pub type FeBn254 = FieldElement<FfBn254>;
