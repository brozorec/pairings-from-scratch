use pairings_from_scratch::{curves::TinyJJ, elliptic_curve::AffinePoint, polynomial::Polynomial};

fn main() {
    let p = AffinePoint::<TinyJJ>::new_xy(
        Polynomial::from(vec![8]).into(),
        Polynomial::from(vec![8]).into(),
    );
    let q = AffinePoint::<TinyJJ>::new_xy(
        Polynomial::from(vec![7, 0, 4]).into(),
        Polynomial::from(vec![0, 10, 0, 5]).into(),
    );
    let result = p.pairing(&q);
    println!("Output from pairing p and q: {}", result);
}
