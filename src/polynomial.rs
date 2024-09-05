use std::default;
use std::fmt::{self, Display};
use std::ops::{Add, Div, Mul, Neg, Rem, Sub};

use crate::field_element::FieldElement;
use crate::finite_field::FiniteField;

pub trait AbstractCoeff:
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

impl<T> AbstractCoeff for T where
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
pub struct Polynomial<T>(Vec<T>);

impl<T: AbstractCoeff> Polynomial<T> {
    pub fn new(coeffs: Vec<T>) -> Self {
        let leading = coeffs.iter().rposition(|c| c != &T::default()).unwrap_or(0);
        let coefficients = &coeffs[..=leading];
        Self(coefficients.to_vec())
    }

    pub fn from_coefficients(coeffs: &[T]) -> Self {
        Self::new(coeffs.to_vec())
    }

    pub fn coefficients(&self) -> &[T] {
        &self.0
    }

    pub fn leading_coefficient(&self) -> &T {
        &self.coefficients()[self.degree()]
    }

    pub fn degree(&self) -> usize {
        self.coefficients()
            .iter()
            .rposition(|coeff| coeff != &T::default())
            .unwrap_or(0)
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_zero(&self) -> bool {
        self.coefficients()
            .iter()
            .all(|coeff| *coeff == T::default())
    }

    fn div_mod(&self, divisor: &Self) -> (Self, Self) {
        let dividend = self.clone();
        let n = dividend.degree() + divisor.degree() + 1;

        let mut quotient = Polynomial::new(vec![T::default(); n]);
        let mut remainder = self.clone();

        while remainder.degree() >= divisor.degree()
            && *remainder.leading_coefficient() != T::default()
        {
            let leading_coeff =
                remainder.leading_coefficient().clone() / divisor.leading_coefficient().clone();
            let degree = remainder.degree() - divisor.degree();
            let mut coeffs = vec![T::default(); degree + 1];
            coeffs[degree] = leading_coeff;
            let term = Polynomial::new(coeffs);

            quotient = quotient + term.clone();
            remainder = remainder - (term * divisor.clone());
        }
        (quotient, remainder)
    }
}

impl<T: AbstractCoeff> PartialEq for Polynomial<T> {
    fn eq(&self, other: &Self) -> bool {
        if self.degree() != other.degree() { return false; }

        self.0.iter().zip(other.0.iter()).all(|(a, b)| *a == *b)
    }
}

impl<T> FromIterator<T> for Polynomial<T> {
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        Self(iter.into_iter().collect())
    }
}

impl<T: AbstractCoeff> Display for Polynomial<T> {
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

impl<A: AbstractCoeff, M: FiniteField<T = Self>> Into<FieldElement<M>> for Polynomial<A> {
    fn into(self) -> FieldElement<M> {
        FieldElement::<M>::new(self)
    }
}

// Implement operator overloading

impl<T: AbstractCoeff> Add for Polynomial<T> {
    type Output = Polynomial<T>;

    fn add(self, other: Polynomial<T>) -> Polynomial<T> {
        let len = self.len().max(other.len());
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

impl<T: AbstractCoeff> Sub for Polynomial<T> {
    type Output = Polynomial<T>;

    fn sub(self, other: Polynomial<T>) -> Polynomial<T> {
        self + (-other)
    }
}

impl<T: AbstractCoeff> Mul for Polynomial<T> {
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

impl<T: AbstractCoeff> Mul<T> for Polynomial<T> {
    type Output = Polynomial<T>;

    fn mul(self, scalar: T) -> Polynomial<T> {
        let coeffs = self.0.iter().map(|v| v.clone() * scalar.clone()).collect();

        Polynomial::new(coeffs)
    }
}

impl<T: AbstractCoeff> Neg for Polynomial<T> {
    type Output = Polynomial<T>;

    fn neg(self) -> Polynomial<T> {
        let neg_coeffs = self.0.iter().map(|v| -v.clone()).collect();
        Polynomial::new(neg_coeffs)
    }
}

impl<T: AbstractCoeff> Div for Polynomial<T> {
    type Output = Polynomial<T>;

    fn div(self, other: Polynomial<T>) -> Polynomial<T> {
        self.div_mod(&other).0
    }
}

impl<T: AbstractCoeff> Rem for Polynomial<T> {
    type Output = Polynomial<T>;

    fn rem(self, other: Polynomial<T>) -> Polynomial<T> {
        self.div_mod(&other).1
    }
}
