use std::env;

use crate::{
    elliptic_curve::{AffinePoint, EllipticCurve},
    field_element::FieldElement,
};

pub fn is_log_mode() -> bool {
    env::var("LOG_MODE").is_ok()
}

pub fn log_table_titles() {
    if is_log_mode() {
        println!(
            "{0: <5} | {1: <10} | {2: <25} | {3: <25} | {4: <15}",
            "bit",
            "operation",
            "distance relationship",
            "accumulated distance",
            "accumulator point"
        );
        println!(
            "{0:-<5} | {1:-<10} | {2:-<25} | {3:-<25} | {4:-<15}",
            "", "", "", "", ""
        );
    }
}

pub fn log_table_row<E: EllipticCurve>(
    bit: &str,
    op: &str,
    f_new: &FieldElement<E::BaseField>,
    f: &FieldElement<E::BaseField>,
    point: &AffinePoint<E>,
) {
    if is_log_mode() {
        match point {
            AffinePoint::XY(x, y) => {
                println!(
                    "{: <5} | {: <10} | {: <25} | {: <25} | ({:?}, {:?})",
                    bit,
                    op,
                    f_new.to_string(),
                    f.to_string(),
                    x,
                    y
                );
            }
            AffinePoint::Infinity => {
                println!(
                    "{: <5} | {: <10} | {: <25} | {: <25} | Infinity ",
                    bit,
                    op,
                    f_new.to_string(),
                    f.to_string()
                );
                println!(
                    "{0:-<5} | {1:-<10} | {2:-<25} | {3:-<25} | {4:-<15}",
                    "", "", "", "", ""
                );
                println!("Output from Miller's loop: {:?}", f);
            }
        }
    }
}
