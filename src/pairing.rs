use num_traits::Pow;

use crate::{
    curves::TinyJJ,
    elliptic_curve::{AffinePoint, EllipticCurve},
    field_element::FieldElement,
    finite_field::{FiniteField, NonExtendedField},
};

pub fn linefunc<E: EllipticCurve>(
    p: AffinePoint<E>,
    q: AffinePoint<E>,
    t: AffinePoint<E>,
) -> FieldElement<E::BaseField> {
    // Infinity not allowed
    let (x1, y1) = p.xy().unwrap();
    let (x2, y2) = q.xy().unwrap();
    let (xt, yt) = t.xy().unwrap();

    if x1 != x2 {
        let m = (y2 - y1.clone()) / (x2 - x1.clone());
        m * (xt - x1) - (yt - y1)
    } else if y1 == y2 {
        // (3 * x1**2 + a) / (2 * y1)
        let x_pow_2 = x1.clone() * x1.clone();
        let m = (x_pow_2.clone() + x_pow_2.clone() + x_pow_2 + E::a()) / (y1.clone() + y1.clone());
        m * (xt - x1) - (yt - y1)
    } else {
        xt - x1
    }
}

pub trait Pairing: EllipticCurve<ScalarField: NonExtendedField> {
    fn pairing(p: AffinePoint<Self>, q: AffinePoint<Self>) -> FieldElement<Self::BaseField> {
        // TODO: check points are valid
        let mut point = p.clone();
        let mut f = FieldElement::<Self::BaseField>::one();
        let bits = Self::ScalarField::to_bits_be(Self::r());
        for bit in bits.iter().rev().skip(1) {
            let f_new = linefunc(point.clone(), point.clone(), q.clone());
            f = f.clone() * f * f_new;
            point = point.clone().double();

            if *bit == 1 {
                let f_new = linefunc(point.clone(), p.clone(), q.clone());
                f = f * f_new;
                point = point + p.clone()
            }
        }

        assert!(point == p.clone() * Self::r());

        let k = Self::embedding_degree();
        let q = <Self::ScalarField as FiniteField>::modulus();
        let one = <Self::ScalarField as FiniteField>::one();
        let f_exp = (q.pow(k) - one) / Self::r();

        f.pow::<Self::ScalarField>(f_exp)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        curves::TinyJJ, elliptic_curve::EllipticCurve, fields::F13_4, polynomial::Polynomial,
    };

    #[test]
    fn test_pairing() {
        let p = AffinePoint::<TinyJJ>::new_xy(
            Polynomial::from(vec![8]).into(),
            Polynomial::from(vec![8]).into(),
        );
        let q = AffinePoint::<TinyJJ>::new_xy(
            Polynomial::from(vec![7, 0, 4]).into(),
            Polynomial::from(vec![0, 10, 0, 5]).into(),
        );
        let result: F13_4 = Polynomial::from(vec![3, 7, 7, 6]).into();
        assert!(result == Pairing::pairing(p, q));
    }

    #[test]
    fn test_pairing_bilinearity() {
        let p = AffinePoint::<TinyJJ>::new_xy(
            Polynomial::from(vec![8]).into(),
            Polynomial::from(vec![8]).into(),
        );
        let q = AffinePoint::<TinyJJ>::new_xy(
            Polynomial::from(vec![7, 0, 4]).into(),
            Polynomial::from(vec![0, 10, 0, 5]).into(),
        );
        assert!(Pairing::pairing(p.clone(), q.double()) == Pairing::pairing(p.double(), q));
    }

    #[test]
    fn test_pairing_linefunc() {
        let one = TinyJJ::generator();
        let two = TinyJJ::generator() * 2;
        let three = TinyJJ::generator() * 3;
        let negone = TinyJJ::generator() * (TinyJJ::order() - 1);
        let negtwo = TinyJJ::generator() * (TinyJJ::order() - 2);
        let negthree = TinyJJ::generator() * (TinyJJ::order() - 3);

        assert_eq!(
            linefunc(one.clone(), two.clone(), one.clone()),
            F13_4::zero()
        );
        assert_eq!(
            linefunc(one.clone(), two.clone(), two.clone()),
            F13_4::zero()
        );
        assert_ne!(
            linefunc(one.clone(), two.clone(), three.clone()),
            F13_4::zero()
        );
        assert_eq!(linefunc(one.clone(), two.clone(), negthree), F13_4::zero());
        assert_eq!(
            linefunc(one.clone(), negone.clone(), one.clone()),
            F13_4::zero()
        );
        assert_eq!(
            linefunc(one.clone(), negone.clone(), negone.clone()),
            F13_4::zero()
        );
        assert_ne!(linefunc(one.clone(), negone, two.clone()), F13_4::zero());
        assert_eq!(
            linefunc(one.clone(), one.clone(), one.clone()),
            F13_4::zero()
        );
        assert_ne!(linefunc(one.clone(), one.clone(), two), F13_4::zero());
        assert_eq!(linefunc(one.clone(), one, negtwo), F13_4::zero());
    }
}
