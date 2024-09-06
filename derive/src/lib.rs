use proc_macro::TokenStream;
use quote::quote;

#[proc_macro_attribute]
pub fn polynomial_inverse(_metadata: TokenStream, _input: TokenStream) -> TokenStream {
    let method = quote! {
        fn inverse(value: &Self::T) -> Self::T {
            let mut r0 = Self::modulus();
            let mut r1 = value.clone();

            let mut t0 = Self::zero();
            let mut t1 = Self::one();

            while !r1.is_zero() {
                let quotient = r0.clone() / r1.clone();

                let r2 = r0 - quotient.clone() * r1.clone();
                r0 = r1;
                r1 = r2;

                let t2 = t0 - quotient * t1.clone();
                t0 = t1;
                t1 = t2;
            }

            let lead_inverse = r0.leading_coefficient().inverse();
            t0 * lead_inverse
        }
    };
    TokenStream::from(method)
}
