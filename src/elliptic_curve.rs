use std::{
    ops::{Add, Mul, Neg},
    usize,
};

use num_traits::Pow;

use crate::{
    field_element::FieldElement,
    finite_field::{FiniteField, NonExtendedField},
    pairing::Pairing,
};

pub trait EllipticCurve: Clone + PartialEq {
    type BaseField: FiniteField;
    type ScalarField: FiniteField;

    fn a() -> FieldElement<Self::BaseField>;
    fn b() -> FieldElement<Self::BaseField>;

    fn embedding_degree() -> usize;
    fn generator() -> AffinePoint<Self>;

    // numer of points in the extened curve
    fn order() -> <Self::ScalarField as FiniteField>::T;

    // this biggest cofactor of the order of the NON-extended curve
    fn r() -> <Self::ScalarField as FiniteField>::T;

    fn pairing(p: &AffinePoint<Self>, q: &AffinePoint<Self>) -> FieldElement<Self::BaseField>
    where
        Self: Pairing,
    {
        Pairing::pairing(p, q)
    }
}

#[derive(Debug, PartialEq, Clone)]
pub enum AffinePoint<E: EllipticCurve> {
    Infinity,
    XY(FieldElement<E::BaseField>, FieldElement<E::BaseField>),
}

impl<E: EllipticCurve> AffinePoint<E> {
    pub fn new_xy(x: FieldElement<E::BaseField>, y: FieldElement<E::BaseField>) -> Self {
        // check coordinates are valid, if not return Infinity
        if Self::is_on_curve(&x, &y) {
            AffinePoint::XY(x, y)
        } else {
            Self::new_inf()
        }
    }

    pub fn new_inf() -> Self {
        AffinePoint::Infinity
    }

    pub fn is_on_curve(x: &FieldElement<E::BaseField>, y: &FieldElement<E::BaseField>) -> bool {
        y.clone() * y.clone() == x.clone() * x.clone() * x.clone() + E::a() * x.clone() + E::b()
    }

    /// Returns the x and y coordinates of this affine point.
    pub fn xy(&self) -> Option<(FieldElement<E::BaseField>, FieldElement<E::BaseField>)> {
        match self {
            AffinePoint::XY(x, y) => Some((x.clone(), y.clone())),
            _ => None,
        }
    }

    /// Returns the x coordinate of this affine point.
    pub fn x(&self) -> Option<FieldElement<E::BaseField>> {
        self.xy().map(|(x, _)| x)
    }

    /// Returns the y coordinate of this affine point.
    pub fn y(&self) -> Option<FieldElement<E::BaseField>> {
        self.xy().map(|(_, y)| y)
    }

    /// Is `self` the point at infinity?
    pub fn is_inf(&self) -> bool {
        self.xy().is_none()
    }

    pub fn double(&self) -> Self {
        match self {
            AffinePoint::XY(x, y) => {
                // explain why from claude
                if y.is_zero() {
                    return Self::Infinity;
                }

                let x_pow_2 = x.clone() * x.clone();

                // (3*x^2 + a) / 2*y
                let m = (x_pow_2.clone() + x_pow_2.clone() + x_pow_2 + E::a())
                    / (y.clone() + y.clone());

                let new_x = m.clone() * m.clone() - x.clone() - x.clone();
                let new_y = m.clone() * (x.clone() - new_x.clone()) - y.clone();

                AffinePoint::new_xy(new_x, new_y)
            }
            _ => AffinePoint::Infinity,
        }
    }

    pub fn trace_map(&self) -> Self
    where
        E::ScalarField: NonExtendedField,
    {
        match self {
            AffinePoint::XY(x, y) => {
                let mut point = AffinePoint::XY(x.clone(), y.clone());
                for i in 1..E::embedding_degree() {
                    let power = E::ScalarField::modulus().pow(i);
                    let new_x = x.pow::<E::ScalarField>(power.clone());
                    let new_y = y.pow::<E::ScalarField>(power);
                    point = point + AffinePoint::new_xy(new_x, new_y);
                }
                point
            }
            _ => AffinePoint::Infinity,
        }
    }
}

impl<E: EllipticCurve> Add for AffinePoint<E> {
    type Output = AffinePoint<E>;

    fn add(self, other: Self) -> Self::Output {
        if self == other {
            return self.double();
        }

        if self == -other.clone() {
            return Self::Infinity;
        }

        match self {
            AffinePoint::XY(x1, y1) => {
                if let AffinePoint::XY(x2, y2) = other {
                    let m = (y2 - y1.clone()) / (x2.clone() - x1.clone());
                    let x = m.clone() * m.clone() - x1.clone() - x2.clone();
                    let y = y1 + m * (x.clone() - x1);
                    Self::new_xy(x, -y)
                } else {
                    // other == Infinity
                    Self::new_xy(x1, y1)
                }
            }
            // self == Infinity
            _ => other,
        }
    }
}

impl<E: EllipticCurve<ScalarField: NonExtendedField>> Mul<<E::ScalarField as FiniteField>::T>
    for AffinePoint<E>
{
    type Output = AffinePoint<E>;

    fn mul(self, scalar: <E::ScalarField as FiniteField>::T) -> Self::Output {
        match self {
            AffinePoint::XY(x, y) => {
                let mut point = Self::XY(x.clone(), y.clone());
                for bit in E::ScalarField::to_bits_be(scalar).iter().rev().skip(1) {
                    point = point.double();

                    if *bit == 1 {
                        point = point + Self::XY(x.clone(), y.clone());
                    }
                }
                point
            }
            _ => AffinePoint::Infinity,
        }
    }
}

impl<E: EllipticCurve> Neg for AffinePoint<E> {
    type Output = AffinePoint<E>;

    fn neg(self) -> Self::Output {
        match self {
            AffinePoint::XY(x, y) => AffinePoint::XY(x, -y),
            _ => AffinePoint::Infinity,
        }
    }
}

pub fn get_all_points<E: EllipticCurve<ScalarField: NonExtendedField>>() -> Vec<AffinePoint<E>> {
    let max = E::ScalarField::to_uint(E::order()).unwrap();
    let mut result = vec![E::generator(); max];
    let mut acc = E::generator();
    for i in (1..max).into_iter() {
        acc = acc.clone() + E::generator();
        result[i] = acc.clone();
    }
    result
}

#[cfg(test)]
mod tests {
    use crate::{curves::TinyJJ, fields::Fe13_4, polynomial::Polynomial};

    use super::*;

    //#[test]
    //fn test_ec_all_points() {
    //let points = get_all_points::<TinyJJ>();
    //assert!(points.len() == TinyJJ::order().try_into().unwrap());
    //}

    #[test]
    fn test_ec_generator() {
        let prod = TinyJJ::generator() * TinyJJ::order();
        assert!(prod.is_inf());
    }

    #[test]
    fn test_ec_scalar_mul() {
        let p1 = AffinePoint::<TinyJJ>::new_xy(
            Polynomial::from(vec![8]).into(),
            Polynomial::from(vec![8]).into(),
        );
        let x: Fe13_4 = Polynomial::from(vec![7]).into();
        let y: Fe13_4 = Polynomial::from(vec![2]).into();
        let prod = p1.clone() * 123;
        assert!(prod.x().unwrap() == x);
        assert!(prod.y().unwrap() == y);

        let prod2 = p1 * 5;
        assert!(prod2.is_inf());
    }

    #[test]
    fn test_ec_add() {
        let p1 = AffinePoint::<TinyJJ>::new_xy(
            Polynomial::from(vec![7]).into(),
            Polynomial::from(vec![11]).into(),
        );
        let p2 = AffinePoint::<TinyJJ>::new_xy(
            Polynomial::from(vec![12]).into(),
            Polynomial::from(vec![8]).into(),
        );
        let x: Fe13_4 = Polynomial::from(vec![11]).into();
        let y: Fe13_4 = Polynomial::from(vec![7]).into();
        let sum = p1 + p2;
        assert!(sum.x().unwrap() == x);
        assert!(sum.y().unwrap() == y);
    }

    #[test]
    fn test_ec_neg() {
        let p1 = AffinePoint::<TinyJJ>::new_xy(
            Polynomial::from(vec![8]).into(),
            Polynomial::from(vec![8]).into(),
        );
        let neg = -p1;
        let x: Fe13_4 = Polynomial::from(vec![8]).into();
        let y: Fe13_4 = Polynomial::from(vec![5]).into();
        assert!(neg.x().unwrap() == x);
        assert!(neg.y().unwrap() == y);
    }

    #[test]
    fn test_ec_add_to_inf() {
        let p1 = AffinePoint::<TinyJJ>::new_xy(
            Polynomial::from(vec![8]).into(),
            Polynomial::from(vec![8]).into(),
        );
        let sum = p1.clone() + -p1;
        assert!(sum.is_inf());
    }

    #[test]
    fn test_ec_double() {
        let p1 = AffinePoint::<TinyJJ>::new_xy(
            Polynomial::from(vec![7]).into(),
            Polynomial::from(vec![11]).into(),
        );
        let x: Fe13_4 = Polynomial::from(vec![8]).into();
        let y: Fe13_4 = Polynomial::from(vec![5]).into();
        let sum = p1.double();
        assert!(sum.x().unwrap() == x);
        assert!(sum.y().unwrap() == y);
    }

    #[test]
    fn test_ec_double_to_inf() {
        let p = TinyJJ::generator() * 14400;
        assert!(p.double().is_inf());
    }

    #[test]
    fn test_ec_trace_map() {
        let p1 = AffinePoint::<TinyJJ>::new_xy(
            Polynomial::from(vec![8]).into(),
            Polynomial::from(vec![8]).into(),
        );
        let x: Fe13_4 = Polynomial::from(vec![8]).into();
        let y: Fe13_4 = Polynomial::from(vec![5]).into();
        let tm = p1.trace_map();
        assert!(tm.x().unwrap() == x);
        assert!(tm.y().unwrap() == y);
    }
}
