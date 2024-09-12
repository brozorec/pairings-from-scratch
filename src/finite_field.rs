use std::{
    fmt::Display,
    ops::{Add, BitAnd, Div, Mul, Neg, Rem, Shr, Sub},
    usize,
};

use num_traits::Pow;

pub trait FiniteField: Copy + Eq {
    type T: Clone
        + PartialEq
        + Default
        + Display
        + Neg
        + Add<Output = Self::T>
        + Sub<Output = Self::T>
        + Mul<Output = Self::T>
        + Div<Output = Self::T>
        + Rem<Output = Self::T>;

    /// Returns the modulus of the field.
    fn modulus() -> Self::T;

    /// Returns the additive identity (zero element).
    fn zero() -> Self::T;

    /// Returns the multiplicative identity (one element).
    fn one() -> Self::T;

    /// Reduces a given element under the modulus of the field.
    fn reduce(value: Self::T) -> Self::T {
        (Self::modulus() + value) % Self::modulus()
    }

    /// Computes the multiplicative inverse of an element within the field.
    /// Adding a default implementation, but when T is a polynomial,
    /// it must be re-implemented.
    fn inverse(value: &Self::T) -> Self::T {
        let zero = Self::zero();
        // xgcd
        let mut r0 = Self::modulus();
        let mut r1 = value.clone();
        let mut t0 = zero.clone();
        let mut t1 = Self::one();

        while r1 != zero {
            let quotient = r0.clone() / r1.clone();

            let r2 = r0 - quotient.clone() * r1.clone();
            r0 = r1;
            r1 = r2;

            let t2 = t0 - quotient * t1.clone();
            t0 = t1;
            t1 = t2;
        }

        Self::reduce(t0)
    }
}

pub trait NonExtendedField:
    FiniteField<T: Shr<usize, Output = Self::T> + BitAnd + Pow<usize, Output = Self::T>>
{
    fn to_bits(s: Self::T) -> Vec<bool>;

    fn to_uint(s: Self::T) -> Option<usize>;

    fn from_uint(s: usize) -> Option<Self::T>;
}
