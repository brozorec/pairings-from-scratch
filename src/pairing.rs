use num_traits::Pow;

use crate::{
    elliptic_curve::{AffinePoint, EllipticCurve},
    field_element::FieldElement,
    finite_field::{FiniteField, NonExtendedField},
    logger::{log_table_row, log_table_titles},
};

pub fn dist_relationship<E: EllipticCurve>(
    p: &AffinePoint<E>,
    q: &AffinePoint<E>,
    t: &AffinePoint<E>,
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

pub trait Pairing: EllipticCurve {
    fn is_valid_g1(p: &AffinePoint<Self>) -> bool {
        let k = Self::ScalarField::from_uint(Self::embedding_degree()).unwrap();
        p.trace_map() == p.clone() * k
    }

    fn is_valid_g2(q: &AffinePoint<Self>) -> bool {
        q.trace_map().is_inf()
    }

    fn miller_loop(p: &AffinePoint<Self>, q: &AffinePoint<Self>) -> FieldElement<Self::BaseField> {
        // 1. It first carries out an implicit multiplication of 𝑃 by 𝑟, using
        // the standard double-and-add algorithm for point multiplication.
        // The result of this multiplication is the point at infinity,
        // as 𝑃 is of order 𝑟.
        //
        // 2. At each step in the process, a value from the base field 𝑓 is
        // calculated from a distance relationship between the current line
        // and the point 𝑄.
        //
        // 3. This value is multiplicatively accumulated, and its ﬁnal value
        // is the output of Miller’s loop.
        let mut point = p.clone();
        let mut f = FieldElement::<Self::BaseField>::one();

        let bits = Self::ScalarField::to_bits(Self::r());
        log_table_titles();
        for bit in bits.iter().skip(1) {
            let f_new = dist_relationship(&point, &point, &q);
            f = f.clone() * f * f_new.clone();
            point = point.clone().double();

            log_table_row(*bit, &f_new, &f, &point);

            if *bit {
                let f_new = dist_relationship(&point, &p, &q);
                f = f * f_new.clone();
                point = point + p.clone();

                log_table_row(*bit, &f_new, &f, &point);
            }
        }

        assert!(point.is_inf());
        // same as:
        //assert!(point == p.clone() * Self::r());
        f
    }

    fn final_exponentiation(f: FieldElement<Self::BaseField>) -> FieldElement<Self::BaseField> {
        // The output value 𝑓 from the Miller loop must of the same order
        // as the points 𝑃 and 𝑄 which are of order 𝑟. That’s why we must
        // raise 𝑓 to some power to bring down its order to 𝑟.
        //
        // It turns out that if we raise 𝑓 by the power of (𝑞^𝑘 − 1) / 𝑟,
        // it will eliminate all multiples of order 𝑟 and we'll get a
        // value from the field of order 𝑟.
        let k = Self::embedding_degree();
        let q = Self::ScalarField::modulus();
        let one = Self::ScalarField::one();
        let f_exp = (q.pow(k) - one) / Self::r();

        f.pow::<Self::ScalarField>(f_exp)
    }

    fn tate_pairing(p: &AffinePoint<Self>, q: &AffinePoint<Self>) -> FieldElement<Self::BaseField> {
        assert!(Self::is_valid_g1(p), "p is not a G1 point");
        assert!(Self::is_valid_g2(q), "q is not a G2 point");

        let f = Self::miller_loop(p, q);
        Self::final_exponentiation(f)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        curves::TinyJJ, elliptic_curve::EllipticCurve, fields::Fe13_4, polynomial::Polynomial,
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
        let result: Fe13_4 = Polynomial::from(vec![3, 7, 7, 6]).into();
        assert!(result == Pairing::tate_pairing(&p, &q));
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
        assert!(Pairing::tate_pairing(&p, &q.double()) == Pairing::tate_pairing(&p.double(), &q));
    }

    #[test]
    fn test_pairing_dist_relationship() {
        let one = TinyJJ::generator();
        let two = TinyJJ::generator() * 2;
        let three = TinyJJ::generator() * 3;
        let negone = TinyJJ::generator() * (TinyJJ::order() - 1);
        let negtwo = TinyJJ::generator() * (TinyJJ::order() - 2);
        let negthree = TinyJJ::generator() * (TinyJJ::order() - 3);

        assert_eq!(dist_relationship(&one, &two, &one), Fe13_4::zero());
        assert_eq!(dist_relationship(&one, &two, &two), Fe13_4::zero());
        assert_ne!(dist_relationship(&one, &two, &three), Fe13_4::zero());
        assert_eq!(dist_relationship(&one, &two, &negthree), Fe13_4::zero());
        assert_eq!(dist_relationship(&one, &negone, &one), Fe13_4::zero());
        assert_eq!(dist_relationship(&one, &negone, &negone), Fe13_4::zero());
        assert_ne!(dist_relationship(&one, &negone, &two), Fe13_4::zero());
        assert_eq!(dist_relationship(&one, &one, &one), Fe13_4::zero());
        assert_ne!(dist_relationship(&one, &one, &two), Fe13_4::zero());
        assert_eq!(dist_relationship(&one, &one, &negtwo), Fe13_4::zero());
    }
}
