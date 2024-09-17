# Elliptic Curve Pairings from Scratch in Rust

This repository is a hands-on implementation of elliptic curve pairings from scratch in Rust, inspired by the article series:

- [Pairings for the Rest of Us, Part 1: Finding G1 and G2](https://hackmd.io/@brozorec/pairings-for-the-rest-of-us-1)
- [Pairings for the Rest of Us, Part 2: The Tate Pairing](https://hackmd.io/@brozorec/pairings-for-the-rest-of-us-2)
- [Pairings for the Rest of Us, Part 3: From Scratch in Rust](https://hackmd.io/@brozorec/pairings-for-the-rest-of-us-3)

This implementation starts with basic finite field arithmetic and builds up to implementing elliptic curve pairings, including the Miller loop and final exponentiation. The code is designed for learning and is **not optimized** for performance or real-world usage with large curves.

## Features

- **Finite Field Arithmetic**: Basic operations such as addition, subtraction, multiplication, division, and inversion over finite fields.
- **Elliptic Curve Operations**: Point addition, doubling, and scalar multiplication on elliptic curves.
- **Pairings**: Implementation of Miller's algorithm and final exponentiation for computing elliptic curve pairings.
- **Curves**: Supports the **TinyJubJub** curve for simplicity, as well as larger curves like **BLS6_6** (MoonMath).  

## Usage

### Running the Example

To run the pairing example with logging enabled:

```bash
LOG_MODE=true cargo run --release
```

This will show detailed logs of the pairing process, including the step-by-step Miller loop operations.

Sample output:
```
bit   | operation  | distance relationship     | accumulated distance      | accumulator point
----- | ---------- | ------------------------- | ------------------------- | ---------------
0     | double     | 2 + 3*x + 11*x^2 + 8*x^3  | 2 + 3*x + 11*x^2 + 8*x^3  | (7, 11)
1     | add        | 11 + 3*x + 1*x^2 + 8*x^3  | 2 + 6*x + 10*x^2 + 12*x^3 | (8, 5)
1     | add        | 12 + 4*x^2                | 9 + 2*x + 11*x^2 + 12*x^3 | Infinity 
----- | ---------- | ------------------------- | ------------------------- | ---------------
Output from Miller's loop: 9 + 2*x + 11*x^2 + 12*x^3
Output from pairing p and q: 3 + 7*x + 7*x^2 + 6*x^3
```

### Experiment with Larger Curves

You can experiment with the **MoonMath** (BLS6_6) curve by changing the curve in the code. For example, replace **TinyJubJub** with **MoonMath** in the following lines:

```rust
let p = AffinePoint::<MoonMath>::new_xy(
    Polynomial::from(vec![...]).into(),
    Polynomial::from(vec![...]).into(),
);
let q = AffinePoint::<MoonMath>::new_xy(
    Polynomial::from(vec![...]).into(),
    Polynomial::from(vec![...]).into(),
);
let result: Fe64_6 = p.pairing(&q);
```

Then run the pairing as before:

```bash
LOG_MODE=true cargo run --release
```

## Important Notes

- This code is meant for learning and is **not optimized**. It works well for small curves like **TinyJubJub**, but larger curves such as **BN254** may cause performance issues or errors.
- Future improvements will address these limitations, with a focus on handling larger curves and optimizing performance.

## Next Steps

In the next part of this series, we’ll work on resolving some of the issues encountered when using larger curves like **BN254**, as well as improving the overall performance of the pairing operations.

## Contributions

I welcome feedback, corrections, and suggestions! Feel free to submit pull requests or open issues if you’d like to contribute or improve the project.

## Connect with Me

If you have any questions or feedback, you can reach me on X (Twitter) at [@BoyanBarakov](https://twitter.com/BoyanBarakov).
