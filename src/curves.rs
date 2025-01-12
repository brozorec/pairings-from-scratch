use crate::{
    elliptic_curve::{AffinePoint, EllipticCurve},
    field_element::FieldElement,
    fields::{Ff13, Ff13_4, Ff43_6, Ff43},
    finite_field::FiniteField,
    pairing::Pairing,
    polynomial::Polynomial,
};

#[derive(Debug, Clone, PartialEq)]
pub struct TinyJJ;

impl EllipticCurve for TinyJJ {
    type BaseField = Ff13_4;
    type ScalarField = Ff13;

    fn a() -> FieldElement<Self::BaseField> {
        Polynomial::from(vec![8]).into()
    }

    fn b() -> FieldElement<Self::BaseField> {
        Polynomial::from(vec![8]).into()
    }

    fn generator() -> AffinePoint<Self> {
        // 5*x^3 + 12*x^2 + 2*x + 8 and 5*x^2 + x
        AffinePoint::new_xy(
            Polynomial::from(vec![8, 2, 12, 5]).into(),
            Polynomial::from(vec![0, 1, 5]).into(),
        )
    }

    fn embedding_degree() -> usize {
        4
    }

    fn order() -> <Self::ScalarField as FiniteField>::T {
        28800
    }

    fn r() -> <Self::ScalarField as FiniteField>::T {
        5
    }
}

// BLS6_6
#[derive(Debug, Clone, PartialEq)]
pub struct MoonMath;

impl EllipticCurve for MoonMath {
    type BaseField = Ff43_6;
    type ScalarField = Ff43;

    fn a() -> FieldElement<Self::BaseField> {
        Polynomial::from(vec![0]).into()
    }

    fn b() -> FieldElement<Self::BaseField> {
        Polynomial::from(vec![6]).into()
    }

    fn generator() -> AffinePoint<Self> {
        // 7*v^2, 16*v^3
        AffinePoint::new_xy(
            Polynomial::from(vec![0, 0, 7]).into(),
            Polynomial::from(vec![0, 0, 0, 16]).into(),
        )
    }

    fn embedding_degree() -> usize {
        6
    }

    fn order() -> <Self::ScalarField as FiniteField>::T {
        6321251664
    }

    fn r() -> <Self::ScalarField as FiniteField>::T {
        13
    }
}
