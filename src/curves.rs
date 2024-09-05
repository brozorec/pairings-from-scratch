use crate::{
    elliptic_curve::{AffinePoint, EllipticCurve},
    field_element::FieldElement,
    fields::{Mod13, Mod13_4, F13},
    finite_field::FiniteField,
    polynomial::Polynomial, pairing::Pairing,
};

#[derive(Debug, Clone, PartialEq)]
pub struct TinyJJ;

impl EllipticCurve for TinyJJ {
    type BaseField = Mod13_4;
    type ScalarField = Mod13;

    fn a() -> FieldElement<Self::BaseField> {
        Polynomial::new(vec![F13::new(8)]).into()
    }

    fn b() -> FieldElement<Self::BaseField> {
        Polynomial::new(vec![F13::new(8)]).into()
    }

    fn generator() -> AffinePoint<Self> {
        // 5*x^3 + 12*x^2 + 2*x + 8 and 5*x^2 + x
        AffinePoint::new_xy(
            Polynomial::new(vec![F13::new(8), F13::new(2), F13::new(12), F13::new(5)]).into(),
            Polynomial::new(vec![F13::zero(), F13::one(), F13::new(5)]).into(),
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

impl Pairing for TinyJJ {}
