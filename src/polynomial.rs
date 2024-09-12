use std::fmt::{self, Display};
use std::ops::{Add, Div, Mul, Neg, Rem, Sub};

use crate::field_element::FieldElement;
use crate::finite_field::FiniteField;

pub trait Coefficient:
    Clone
    + Default
    + PartialEq
    + Display
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
{
}

impl<T> Coefficient for T where
    T: Clone
        + Default
        + PartialEq
        + Display
        + Add<Output = T>
        + Sub<Output = T>
        + Mul<Output = T>
        + Div<Output = T>
        + Neg<Output = T>
{
}

#[derive(Clone, Default, Debug)]
pub struct Polynomial<C: Coefficient>(Vec<C>);

impl<C: Coefficient> Polynomial<C> {
    /// Creates a new polynomial, trimming trailing zeros.
    pub fn new(coeffs: Vec<C>) -> Self {
        let leading = coeffs.iter().rposition(|c| c != &C::default()).unwrap_or(0);
        let coefficients = &coeffs[..=leading];
        Self(coefficients.to_vec())
    }

    pub fn from_coefficients(coeffs: &[C]) -> Self {
        Self::new(coeffs.to_vec())
    }

    /// Returns a slice of the polynomial's coefficients.
    #[inline]
    pub fn coefficients(&self) -> &[C] {
        &self.0
    }

    /// Returns the coefficient of the highest-degree term.
    #[inline]
    pub fn leading_coefficient(&self) -> &C {
        &self.coefficients()[self.degree()]
    }

    /// Calculates the degree of the polynomial.
    #[inline]
    pub fn degree(&self) -> usize {
        self.coefficients()
            .iter()
            .rposition(|coeff| coeff != &C::default())
            .unwrap_or(0)
    }

    /// Checks if the polynomial is zero.
    pub fn is_zero(&self) -> bool {
        self.coefficients()
            .iter()
            .all(|coeff| *coeff == C::default())
    }

    /// Performs polynomial long division.
    fn div_mod(&self, divisor: &Self) -> (Self, Self) {
        let mut quotient = Polynomial::new(vec![C::default()]);
        let mut remainder = self.clone();

        while remainder.degree() >= divisor.degree()
            && *remainder.leading_coefficient() != C::default()
        {
            let degree = remainder.degree() - divisor.degree();
            let mut coeffs = vec![C::default(); degree + 1];

            let rem_coeff = remainder.leading_coefficient().clone();
            let div_coeff = divisor.leading_coefficient().clone();
            let leading_coeff = rem_coeff / div_coeff;

            coeffs[degree] = leading_coeff;
            let term = Polynomial::new(coeffs);

            quotient = quotient + term.clone();
            remainder = remainder - term * divisor.clone();
        }
        (quotient, remainder)
    }
}

impl<T: Coefficient> PartialEq for Polynomial<T> {
    fn eq(&self, other: &Self) -> bool {
        if self.degree() != other.degree() {
            return false;
        }

        self.0.iter().zip(other.0.iter()).all(|(a, b)| *a == *b)
    }
}

impl<T: Coefficient> Display for Polynomial<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for (i, coeff) in self.0.iter().enumerate() {
            if *coeff != T::default() {
                if i != 0 {
                    write!(f, " + ")?;
                }

                if i == 0 {
                    write!(f, "{}", coeff)?;
                } else if i == 1 {
                    write!(f, "{}*x", coeff)?;
                } else {
                    write!(f, "{}*x^{}", coeff, i)?;
                }
            }
        }
        Ok(())
    }
}

impl<T: Coefficient> FromIterator<T> for Polynomial<T> {
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        Self(iter.into_iter().collect())
    }
}

// instead of Polynomial::new(vec![Fe13::new(12), Fe13::new(4)])
// just do Polynomial::from(vec![12, 4])
impl<M: FiniteField> From<Vec<M::T>> for Polynomial<FieldElement<M>> {
    fn from(value: Vec<M::T>) -> Self {
        value.into_iter().map(|v| FieldElement::new(v)).collect()
    }
}

impl<A: Coefficient, M: FiniteField<T = Self>> Into<FieldElement<M>> for Polynomial<A> {
    fn into(self) -> FieldElement<M> {
        FieldElement::<M>::new(self)
    }
}

// Implement operator overloading

impl<T: Coefficient> Add for Polynomial<T> {
    type Output = Polynomial<T>;

    fn add(self, other: Polynomial<T>) -> Polynomial<T> {
        let len = self.degree().max(other.degree()) + 1;
        let mut coeffs = vec![T::default(); len];

        for (i, c) in self.coefficients().iter().enumerate() {
            coeffs[i] = c.clone();
        }
        for (i, c) in other.coefficients().iter().enumerate() {
            coeffs[i] = coeffs[i].clone() + c.clone();
        }

        Polynomial::new(coeffs)
    }
}

impl<T: Coefficient> Sub for Polynomial<T> {
    type Output = Polynomial<T>;

    fn sub(self, other: Polynomial<T>) -> Polynomial<T> {
        self + -other
    }
}

impl<T: Coefficient> Mul for Polynomial<T> {
    type Output = Polynomial<T>;

    fn mul(self, other: Polynomial<T>) -> Polynomial<T> {
        let n = self.degree() + other.degree() + 1;
        let mut coeffs = vec![T::default(); n];

        for (i, a) in self.coefficients().iter().enumerate() {
            if *a == T::default() {
                continue;
            }
            for (j, b) in other.coefficients().iter().enumerate() {
                coeffs[i + j] = coeffs[i + j].clone() + a.clone() * b.clone();
            }
        }

        Polynomial::new(coeffs)
    }
}

impl<T: Coefficient> Mul<T> for Polynomial<T> {
    type Output = Polynomial<T>;

    fn mul(self, scalar: T) -> Polynomial<T> {
        let coeffs = self.0.iter().map(|v| v.clone() * scalar.clone()).collect();

        Polynomial::new(coeffs)
    }
}

impl<T: Coefficient> Neg for Polynomial<T> {
    type Output = Polynomial<T>;

    fn neg(self) -> Polynomial<T> {
        let neg_coeffs = self.0.iter().map(|v| -v.clone()).collect();
        Polynomial::new(neg_coeffs)
    }
}

impl<T: Coefficient> Div for Polynomial<T> {
    type Output = Polynomial<T>;

    fn div(self, other: Polynomial<T>) -> Polynomial<T> {
        self.div_mod(&other).0
    }
}

impl<T: Coefficient> Rem for Polynomial<T> {
    type Output = Polynomial<T>;

    fn rem(self, other: Polynomial<T>) -> Polynomial<T> {
        self.div_mod(&other).1
    }
}

#[cfg(test)]
mod tests {
    use crate::{fields::Fe13, finite_field::FiniteField};

    use super::*;

    #[test]
    fn test_polynomial_creation() {
        let p = Polynomial::new(vec![1, 2, 3]);
        assert_eq!(p.coefficients(), &[1, 2, 3]);
    }

    #[test]
    fn test_polynomial_degree() {
        let p = Polynomial::new(vec![0, 0, 3, 0, 5]);
        assert_eq!(p.degree(), 4);
    }

    #[test]
    fn test_polynomial_leading_coefficient() {
        let p = Polynomial::new(vec![0, 0, 3, 0, 5]);
        assert_eq!(*p.leading_coefficient(), 5);
    }

    #[test]
    fn test_polynomial_is_zero() {
        let p = Polynomial::new(vec![0, 0, 3, 0, 5]);
        assert_eq!(p.is_zero(), false);
        let p = Polynomial::new(vec![Fe13::zero()]);
        assert_eq!(p.is_zero(), true);
    }

    #[test]
    fn test_polynomial_addition() {
        let p1 = Polynomial::new(vec![1, 2, 3]); // Represents 1 + 2x + 3x^2
        let p2 = Polynomial::new(vec![3, 4, 5]); // Represents 3 + 4x + 5x^2
        let sum = p1.add(p2); // Should represent 4 + 6x + 8x^2
        assert_eq!(sum.coefficients(), &[4, 6, 8]);
    }

    #[test]
    fn test_polynomial_subtraction() {
        let p1 = Polynomial::new(vec![5, 7, 9]); // Represents 5 + 7x + 9x^2
        let p2 = Polynomial::new(vec![1, 2, 3]); // Represents 1 + 2x + 3x^2
        let difference = p1.sub(p2); // Should represent 4 + 5x + 6x^2
        assert_eq!(difference.coefficients(), &[4, 5, 6]);
    }

    #[test]
    fn test_polynomial_multiplication() {
        let p1 = Polynomial::new(vec![1, 2]); // Represents 1 + 2x
        let p2 = Polynomial::new(vec![3, 4]); // Represents 3 + 4x
        let product = p1.mul(p2); // Should represent 3 + 10x + 8x^2
        assert_eq!(product.coefficients(), &[3, 10, 8]);
    }

    #[test]
    fn test_polynomial_division() {
        let dividend = Polynomial::new(vec![1, -3, 2]); // Represents 1 - 3x + 2x^2
        let divisor = Polynomial::new(vec![1, -1]); // Represents 1 - x
        let (quotient, remainder) = dividend.div_mod(&divisor); // Should represent quotient = 2x - 1, remainder = 1
        assert_eq!(quotient.coefficients(), &[1, -2]);
        assert_eq!(remainder.coefficients(), &[0]);
    }

    #[test]
    fn test_polynomial_zero_division() {
        let p1 = Polynomial::new(vec![1, 2, 3]);
        let zero_poly = Polynomial::new(vec![0]);
        let result = std::panic::catch_unwind(|| {
            let _ = p1 / zero_poly;
        });
        assert!(result.is_err()); // Division by zero should panic or handle accordingly
    }

    #[test]
    fn test_polynomial_zero_addition() {
        let p1 = Polynomial::new(vec![1, 2, 3]);
        let zero_poly = Polynomial::new(vec![0, 0, 0]);
        let sum = p1.clone() + zero_poly; // Should result in p1 itself
        assert_eq!(sum.coefficients(), p1.coefficients());
    }

    #[test]
    fn test_polynomial_zero_subtraction() {
        let p1 = Polynomial::new(vec![1, 2, 3]);
        let zero_poly = Polynomial::new(vec![0, 0, 0]);
        let difference = p1.clone() - zero_poly; // Should result in p1 itself
        assert_eq!(difference.coefficients(), p1.coefficients());
    }

    #[test]
    fn test_polynomial_addition_ff() {
        let p1 = Polynomial::new(vec![Fe13::new(1), Fe13::new(2), Fe13::new(9)]); // Represents 1 + 2x + 3x^2
        let p2 = Polynomial::new(vec![Fe13::new(3), Fe13::new(4), Fe13::new(5)]); // Represents 3 + 4x + 5x^2
        let sum = p1 + p2; // Should represent 4 + 6x + 8x^2
        assert_eq!(
            sum.coefficients(),
            &[Fe13::new(4), Fe13::new(6), Fe13::new(1)]
        );
    }

    #[test]
    fn test_polynomial_multiplication_ff() {
        let p1 = Polynomial::new(vec![Fe13::new(7), Fe13::new(3)]);
        let p2 = Polynomial::new(vec![Fe13::new(5), Fe13::new(6)]);
        let product = p1 * p2;
        assert_eq!(
            product.coefficients(),
            &[Fe13::new(9), Fe13::new(5), Fe13::new(5)]
        );
    }

    #[test]
    fn test_polynomial_division_ff() {
        let p1 = Polynomial::new(vec![Fe13::new(2), Fe13::new(5), Fe13::new(7), Fe13::new(3)]);
        let p2 = Polynomial::new(vec![Fe13::new(1), Fe13::new(4), Fe13::new(1)]);
        let quotient = p1 / p2;
        assert_eq!(quotient.coefficients(), &[Fe13::new(8), Fe13::new(3)]);
    }

    #[test]
    fn test_polynomial_division_ff_case_1() {
        let p1 = Polynomial::new(vec![
            Fe13::new(6),
            Fe13::new(3),
            Fe13::new(10),
            Fe13::new(7),
        ]);
        let p2 = Polynomial::new(vec![Fe13::new(3), Fe13::new(2)]);
        let (quotient, remainder) = p1.div_mod(&p2);
        assert_eq!(
            quotient.coefficients(),
            &[Fe13::new(10), Fe13::new(3), Fe13::new(10)]
        );
        assert_eq!(remainder.coefficients(), &[Fe13::new(2)]);
    }

    #[test]
    fn test_polynomial_division_ff_case_2() {
        let p1 = Polynomial::new(vec![Fe13::new(8), Fe13::new(10), Fe13::new(12)]);
        let p2 = Polynomial::new(vec![Fe13::new(4), Fe13::new(1)]);
        let (quotient, remainder) = p1.div_mod(&p2);
        assert_eq!(quotient.coefficients(), &[Fe13::new(1), Fe13::new(12)]);
        assert_eq!(remainder.coefficients(), &[Fe13::new(4)]);
    }

    #[test]
    fn test_polynomial_division_ff_case_3() {
        let p1 = Polynomial::new(vec![Fe13::new(3), Fe13::new(5)]);
        let p2 = Polynomial::new(vec![Fe13::new(6)]);
        let (quotient, remainder) = p1.div_mod(&p2);
        assert_eq!(quotient.coefficients(), &[Fe13::new(7), Fe13::new(3)]);
        assert_eq!(remainder.coefficients(), &[Fe13::new(0)]);
    }
}
