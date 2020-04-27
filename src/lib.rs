//! `haven` is a modern lattice cryptography library.
//!
//! It aims to provide safe, ergonomic and fast cryptographic building blocks.
//!
//! The focus right now is on (Fully) Homomorphic Encryption schemes, however the overall goal is
//! to support various lattice-based cryptographic capabilities.
//!
//! **NOTE:** It's early days for `haven`. Everything you currently get is highly experimental.
//! Please don't use it in production! Given its state I'd love to see all the different ways it
//! might break!
//!
//! Feel free to reach out if you run into any problems or just want to provide some feedback.

#[allow(non_snake_case)]
/// Implementation of the [GSW](https://eprint.iacr.org/2013/340.pdf) Homomorphic Encryption scheme.
pub struct GSW {
    /// Modulus
    q: usize,
    /// Lattice dimension
    n: usize,
    /// Error distribution
    x: usize,
    /// Automatically computed via `\mathcal{O}(n log q)`
    m: usize,
    /// Automatically computed via `\lfloor log_2 q \rfloor + 1`
    l: usize,
    /// Automatically computed via `(n + 1) \times l`
    N: usize,
}

#[allow(non_snake_case)]
impl GSW {
    /// Create a new `GSW` instance with the given security parameters.
    pub fn new(q: usize, n: usize, x: usize) -> Self {
        let m: usize = n * (q as f64).log2() as usize; // TODO: Check if `log2` is sufficient here
        let l: usize = (q as f64).log2().floor() as usize + 1;
        let N: usize = (n + 1) * l;
        GSW { q, n, x, m, l, N }
    }
}

#[cfg(test)]
mod tests {
    use crate::GSW;

    #[test]
    fn gsw_instance_creation() {
        let q = 65536;
        let n = 3;
        let x = 8;
        let gsw = GSW::new(q, n, x);

        assert_eq!(gsw.q, q);
        assert_eq!(gsw.n, n);
        assert_eq!(gsw.x, x);
        assert_eq!(gsw.m, 48);
        assert_eq!(gsw.l, 17);
        assert_eq!(gsw.N, 68);
    }
}
