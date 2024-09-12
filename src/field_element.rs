use core::fmt;
use std::{
    fmt::{Debug, Display},
    ops::{Add, Div, Mul, Neg, Sub},
};

use crate::finite_field::{FiniteField, NonExtendedField};

#[derive(Clone, Eq, PartialEq)]
pub struct FieldElement<M: FiniteField>(M::T);

impl<M: FiniteField> FieldElement<M> {
    pub fn new(value: M::T) -> Self {
        FieldElement(M::reduce(value))
    }

    pub fn value(&self) -> &M::T {
        &self.0
    }

    pub fn zero() -> Self {
        Self::new(M::zero())
    }

    pub fn is_zero(&self) -> bool {
        M::zero() == self.0
    }

    pub fn one() -> Self {
        Self::new(M::one())
    }

    pub fn inverse(&self) -> Self {
        Self::new(M::inverse(&self.0))
    }

    // Exponentiation is computed with the square-and-multiply algorithm,
    // which iterates over the bits in the expansion of the exponent,
    // squares an accumulator variable in each iteration,
    // and additionally multiplies it by the base element if the bit is set.
    pub fn pow<S: NonExtendedField>(&self, exp: S::T) -> Self {
        if exp == S::zero() {
            return Self::one();
        }

        let mut acc = Self::one();
        let base = self.clone();

        for bit in S::to_bits(exp).iter() {
            acc = acc.clone() * acc.clone();

            if *bit {
                acc = acc * base.clone();
            }
        }

        acc
    }
}

impl<M: FiniteField> Display for FieldElement<M> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl<M: FiniteField> Debug for FieldElement<M> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

// Implement arithmetic operations for FieldElement

/// Addition of two field elements is simply the sum of their values,
/// followed by reduction modulo the field’s modulus.
impl<M: FiniteField> Add for FieldElement<M> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self::new(self.0 + rhs.0)
    }
}

/// Negation in a finite field is defined as the
/// additive inverse of an element where a + -a ≡ p mod p
impl<M: FiniteField> Neg for FieldElement<M> {
    type Output = Self;

    fn neg(self) -> Self {
        Self::new(M::modulus() - self.0)
    }
}

/// Subtraction is equivalent to adding the negation of
/// the second element to the first.
impl<M: FiniteField> Sub for FieldElement<M> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        self + -rhs
    }
}

/// Multiplication of two field elements follows the usual
/// rules for integers, with the result reduced modulo the
/// field’s modulus.
impl<M: FiniteField> Mul for FieldElement<M> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self::new(self.0 * rhs.0)
    }
}

/// Division is equivalent to multiplying the left-hand side
/// by the multiplicative inverse of the right-hand side.
impl<M: FiniteField> Div for FieldElement<M> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        self * rhs.inverse()
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
        fields::{Fe13, Fe13_2, FeBn254, Ff13},
        polynomial::Polynomial,
    };

    #[test]
    fn test_finite_field_add() {
        let a = Fe13::new(7);
        let b = Fe13::new(10);
        let c = Fe13::new(4); // (7 + 10) % 13 = 17 % 13 = 4
        assert_eq!(a + b, c);
    }

    #[test]
    fn test_finite_field_big_add() {
        let a = FeBn254::new(BigInt::from_i64(-1).unwrap());
        let b = FeBn254::new(BigInt::from_i64(-1).unwrap());
        let c = FeBn254::new(
            BigInt::from_str(
                "21888242871839275222246405745257275088696311157297823662689037894645226208581",
            )
            .unwrap(),
        );
        assert_eq!(a + b, c);
    }

    #[test]
    fn test_finite_field_sub() {
        let a = Fe13::new(7);
        let b = Fe13::new(10);
        let c = Fe13::new(10); // (7 - 10) % 13 = -3 % 13 = 10
        assert_eq!(a - b, c);
    }

    #[test]
    fn test_finite_field_mul() {
        let a = Fe13::new(7);
        let b = Fe13::new(10);
        let c = Fe13::new(5); // (7 * 10) % 13 = 70 % 13 = 5
        assert_eq!(a * b, c);
    }

    #[test]
    fn test_ext_finite_field_mul() {
        let element1: Fe13_2 = Polynomial::from(vec![7, 3]).into();
        let element2: Fe13_2 = Polynomial::from(vec![5, 6]).into();

        let product = element1 * element2;
        assert_eq!(*product.value(), Polynomial::from(vec![12, 5]).into())
    }

    #[test]
    fn test_finite_field_div() {
        let a = Fe13::new(7);
        let b = Fe13::new(3);
        let c = Fe13::new(11); // (7 / 3) % 13 = 7 * 9 % 13 = 63 % 13 = 11
        assert_eq!(a / b, c);
    }

    #[test]
    fn test_finite_field_neg() {
        let a = Fe13::new(7);
        let b = Fe13::new(6); // -7 % 13 = 6
        assert_eq!(-a, b);
    }

    #[test]
    fn test_finite_field_inverse() {
        let a = Fe13::new(7);
        let inv_a = Fe13::new(2); // The inverse of 7 mod 13 is 2
        assert_eq!(a.inverse(), inv_a);
    }

    #[test]
    fn test_finite_field_big_inverse() {
        let num = BigInt::from_str("21888242871839275222246405745").unwrap();
        let a = FeBn254::new(num);
        let b = a.inverse();
        assert_eq!(a * b, FeBn254::one());
    }

    #[test]
    fn test_ext_finite_field_inverse() {
        let ext_field_element: Fe13_2 = Polynomial::from(vec![3, 5]).into();
        let inverse = ext_field_element.inverse();
        assert_eq!(*inverse.value(), Polynomial::from(vec![6, 3]));

        // Multiply the element by its inverse and check if we get the identity element
        let product = inverse * ext_field_element;
        let identity = Polynomial::new(vec![Fe13::one()]);
        assert_eq!(*product.value(), identity);
    }

    #[test]
    fn test_finite_field_identity_add() {
        let zero = Fe13::zero();
        let a = Fe13::new(5);
        assert_eq!(a.clone() + zero.clone(), a.clone());
        assert_eq!(zero + a.clone(), a);
    }

    #[test]
    fn test_finite_field_identity_mul() {
        let one = Fe13::one();
        let a = Fe13::new(5);
        assert_eq!(a.clone() * one.clone(), a.clone());
        assert_eq!(one * a.clone(), a);
    }

    #[test]
    fn test_finite_field_associative_add() {
        let a = Fe13::new(3);
        let b = Fe13::new(4);
        let c = Fe13::new(5);
        assert_eq!((a.clone() + b.clone()) + c.clone(), a + (b + c));
    }

    #[test]
    fn test_finite_field_associative_mul() {
        let a = Fe13::new(3);
        let b = Fe13::new(4);
        let c = Fe13::new(5);
        assert_eq!((a.clone() * b.clone()) * c.clone(), a * (b * c));
    }

    #[test]
    fn test_finite_field_distributive_property() {
        let a = Fe13::new(3);
        let b = Fe13::new(4);
        let c = Fe13::new(5);
        assert_eq!(
            a.clone() * (b.clone() + c.clone()),
            (a.clone() * b) + (a * c)
        );
    }

    #[test]
    fn test_finite_field_pow() {
        assert_eq!(Fe13::new(8).pow::<Ff13>(169), Fe13::new(8));
        assert_eq!(Fe13::new(8).pow::<Ff13>(13), Fe13::new(8));
        assert_eq!(Fe13::new(3).pow::<Ff13>(5), Fe13::new(9));
        assert_eq!(Fe13::new(3).pow::<Ff13>(6), Fe13::one());
        assert_eq!(Fe13::new(3).pow::<Ff13>(0), Fe13::one());
        assert_eq!(Fe13::new(4).pow::<Ff13>(7), Fe13::new(4));
        assert_eq!(Fe13::new(0).pow::<Ff13>(4), Fe13::zero());
    }

    #[test]
    fn test_ext_finite_field_pow() {
        let element: Fe13_2 = Polynomial::from(vec![3]).into();
        let res: Fe13_2 = Polynomial::from(vec![9]).into();
        assert_eq!(element.pow::<Ff13>(5), res);
    }
}
