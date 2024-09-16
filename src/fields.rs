use std::str::FromStr;

use derive_lib::polynomial_inverse;
use num_bigint::BigInt;
use num_traits::{FromPrimitive, One, ToPrimitive, Zero};

use crate::{
    field_element::FieldElement,
    finite_field::{FiniteField, NonExtendedField},
    polynomial::Polynomial,
};

// ---------------- Ff 13 ---------------------

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
    fn to_bits(s: i16) -> Vec<bool> {
        let max = i16::BITS - s.leading_zeros();
        let mut res = vec![false; max as usize];
        for i in 0..max {
            res[i as usize] = (s >> (max - 1 - i) & 1) != 0;
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

// ---------------- Ff 43 ---------------------

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub struct Ff43;

impl FiniteField for Ff43 {
    type T = i64;

    fn modulus() -> Self::T {
        43
    }

    fn one() -> Self::T {
        1
    }

    fn zero() -> Self::T {
        0
    }
}

impl NonExtendedField for Ff43 {
    fn to_bits(s: i64) -> Vec<bool> {
        let max = i64::BITS - s.leading_zeros();
        let mut res = vec![false; max as usize];
        for i in 0..max {
            res[i as usize] = (s >> (max - 1 - i) & 1) != 0;
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
pub struct Ff43_6;

impl FiniteField for Ff43_6 {
    type T = Polynomial<Fe43>;

    fn modulus() -> Self::T {
        // 
        Polynomial::from(vec![6, 0, 0, 0, 0, 0, 1])
    }

    fn zero() -> Self::T {
        Polynomial::new(vec![Fe43::zero()])
    }

    fn one() -> Self::T {
        Polynomial::new(vec![Fe43::one()])
    }

    #[polynomial_inverse]
    fn inverse(value: &Self::T) -> Self::T;
}



// ---------------- Ff Bn254 ---------------------

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
    fn to_bits(s: BigInt) -> Vec<bool> {
        let max = s.bits();
        let mut res = vec![false; max as usize];
        for i in 0..max {
            res[i as usize] = (s.clone() >> (max - 1 - i) & BigInt::one()) != BigInt::zero();
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
pub type Fe43 = FieldElement<Ff43>;
pub type Fe43_6 = FieldElement<Ff43_6>;
pub type FeBn254 = FieldElement<FfBn254>;
