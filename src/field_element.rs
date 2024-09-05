use core::fmt;
use std::{
    fmt::{Debug, Display},
    ops::{Add, Div, Mul, Neg, Rem, Sub, Shr, BitAnd},
    usize, mem,
};

use crate::finite_field::{FiniteField, NonExtendedField};

#[derive(Clone, Eq, PartialEq)]
pub struct FieldElement<M: FiniteField> {
    pub value: M::T,
}

impl<M: FiniteField> FieldElement<M> {
    pub fn new(value: M::T) -> Self {
        FieldElement {
            value: M::reduce(value),
        }
    }

    pub fn zero() -> Self {
        Self::new(M::zero())
    }

    pub fn is_zero(&self) -> bool {
        M::zero() == self.value
    }

    pub fn one() -> Self {
        Self::new(M::one())
    }

    pub fn inverse(&self) -> Self {
        Self::new(M::inverse(&self.value))
    }

    // Montgomery ladder
    pub fn pow<S: NonExtendedField>(&self, exp: S::T) -> Self {
        if exp == S::zero() {
            return Self::one();
        }

        let mut r0 = Self::one();
        let mut r1 = self.clone();

        for bit in S::to_bits_be(exp).iter() {
            if *bit == 1 {
                r0 = r1.clone() * r0;
            }

            r1 = r1.clone() * r1.clone();
        }

        r0
    }
}

impl<M: FiniteField> Display for FieldElement<M> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        //write!(f, "F{}({})", M::modulus(), self.value)
        write!(f, "{}", self.value)
    }
}

impl<M: FiniteField> Debug for FieldElement<M> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.value)
    }
}

// Implement arithmetic operations for FieldElement
impl<M: FiniteField> Add for FieldElement<M> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self::new(self.value + rhs.value)
    }
}

impl<M: FiniteField> Sub for FieldElement<M> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Self::new(M::modulus() + self.value - rhs.value)
    }
}

impl<M: FiniteField> Mul for FieldElement<M> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self::new(self.value * rhs.value)
    }
}

impl<M: FiniteField> Div for FieldElement<M> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        self * rhs.inverse()
    }
}

impl<M: FiniteField> Neg for FieldElement<M> {
    type Output = Self;

    fn neg(self) -> Self {
        Self::new(M::modulus() - self.value)
    }
}

impl<M: FiniteField> Default for FieldElement<M> {
    fn default() -> Self {
        Self::zero()
    }
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use num_bigint::BigInt;
    use num_traits::FromPrimitive;

    use crate::{
        fields::{F13, F13_2, FBN254, Mod13, BN254},
        polynomial::Polynomial,
    };

    use super::*;

    #[test]
    fn test_finite_field_pow() {
        assert_eq!(F13::new(8).pow::<Mod13>(169), F13::new(8));
        assert_eq!(F13::new(8).pow::<Mod13>(13), F13::new(8));
        assert_eq!(F13::new(3).pow::<Mod13>(5), F13::new(9));
        assert_eq!(F13::new(3).pow::<Mod13>(6), F13::one());
        assert_eq!(F13::new(3).pow::<Mod13>(0), F13::one());
        assert_eq!(F13::new(4).pow::<Mod13>(7), F13::new(4));
        assert_eq!(F13::new(0).pow::<Mod13>(4), F13::zero());
    }

    #[test]
    fn test_ext_finite_field_pow() {
        let element: F13_2 = Polynomial::from(vec![3]).into();
        let res: F13_2 = Polynomial::from(vec![9]).into();
        assert_eq!(element.pow::<Mod13>(5), res);
    }

    #[test]
    fn test_ext_finite_field_inverse() {
        let ext_field_element: F13_2 = Polynomial::from(vec![3, 5]).into();
        let inverse = ext_field_element.inverse();
        assert_eq!(
            inverse.value,
            Polynomial::from(vec![6, 3])
        );

        // Multiply the element by its inverse and check if we get the identity element
        let product = inverse * ext_field_element;
        let identity = Polynomial::new(vec![F13::one()]);
        assert_eq!(product.value, identity);
    }

    #[test]
    fn test_ext_finite_field_mul() {
        let element1: F13_2 = Polynomial::from(vec![7, 3]).into();
        let element2: F13_2 = Polynomial::from(vec![5, 6]).into();

        let product = element1 * element2;
        assert_eq!(
            product.value,
            Polynomial::new(vec![F13::new(12), F13::new(5)])
        )
    }

    #[test]
    fn test_finite_field_big_add() {
        let a = FBN254::new(BigInt::from_i64(-1).unwrap());
        let b = FBN254::new(BigInt::from_i64(-1).unwrap());
        let c = FBN254::new(
            BigInt::from_str(
                "21888242871839275222246405745257275088696311157297823662689037894645226208581",
            )
            .unwrap(),
        );
        assert_eq!(a + b, c);
    }

    #[test]
    fn test_finite_field_big_inverse() {
        let num = BigInt::from_str("21888242871839275222246405745").unwrap();
        let a = FBN254::new(num);
        let b = a.inverse();
        assert_eq!(a * b, FBN254::one());
    }

    #[test]
    fn test_finite_field_add() {
        let a = F13::new(7);
        let b = F13::new(10);
        let c = F13::new(4); // (7 + 10) % 13 = 17 % 13 = 4
        assert_eq!(a + b, c);
    }

    #[test]
    fn test_finite_field_sub() {
        let a = F13::new(7);
        let b = F13::new(10);
        let c = F13::new(10); // (7 - 10) % 13 = -3 % 13 = 10
        assert_eq!(a - b, c);
    }

    #[test]
    fn test_finite_field_mul() {
        let a = F13::new(7);
        let b = F13::new(10);
        let c = F13::new(5); // (7 * 10) % 13 = 70 % 13 = 5
        assert_eq!(a * b, c);
    }

    #[test]
    fn test_finite_field_div() {
        let a = F13::new(7);
        let b = F13::new(3);
        let c = F13::new(11); // (7 / 3) % 13 = 7 * 9 % 13 = 63 % 13 = 11
        assert_eq!(a / b, c);
    }

    #[test]
    fn test_finite_field_neg() {
        let a = F13::new(7);
        let b = F13::new(6); // -7 % 13 = 6
        assert_eq!(-a, b);
    }

    #[test]
    fn test_finite_field_inverse() {
        let a = F13::new(7);
        let inv_a = F13::new(2); // The inverse of 7 mod 13 is 2
        assert_eq!(a.inverse(), inv_a);
    }

    #[test]
    fn test_finite_field_identity_addition() {
        let zero = F13::zero();
        let a = F13::new(5);
        assert_eq!(a.clone() + zero.clone(), a.clone());
        assert_eq!(zero + a.clone(), a);
    }

    #[test]
    fn test_finite_field_identity_multiplication() {
        let one = F13::one();
        let a = F13::new(5);
        assert_eq!(a.clone() * one.clone(), a.clone());
        assert_eq!(one * a.clone(), a);
    }

    #[test]
    fn test_finite_field_associative_addition() {
        let a = F13::new(3);
        let b = F13::new(4);
        let c = F13::new(5);
        assert_eq!((a.clone() + b.clone()) + c.clone(), a + (b + c));
    }

    #[test]
    fn test_finite_field_associative_multiplication() {
        let a = F13::new(3);
        let b = F13::new(4);
        let c = F13::new(5);
        assert_eq!((a.clone() * b.clone()) * c.clone(), a * (b * c));
    }

    #[test]
    fn test_finite_field_distributive_property() {
        let a = F13::new(3);
        let b = F13::new(4);
        let c = F13::new(5);
        assert_eq!(
            a.clone() * (b.clone() + c.clone()),
            (a.clone() * b) + (a * c)
        );
    }
}
