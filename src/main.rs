#![allow(unused)]

use crate::{
    curves::TinyJJ, elliptic_curve::AffinePoint, fields::Fe13_4, pairing::Pairing,
    polynomial::Polynomial,
};

mod curves;
mod elliptic_curve;
mod field_element;
mod fields;
mod finite_field;
mod pairing;
mod polynomial;
mod logger;

fn main() {
    let p = AffinePoint::<TinyJJ>::new_xy(
        Polynomial::from(vec![8]).into(),
        Polynomial::from(vec![8]).into(),
    );
    let q = AffinePoint::<TinyJJ>::new_xy(
        Polynomial::from(vec![7, 0, 4]).into(),
        Polynomial::from(vec![0, 10, 0, 5]).into(),
    );
    let result: Fe13_4 = p.pairing(&q);
    println!("Output from pairing p and q: {}", result);
}
