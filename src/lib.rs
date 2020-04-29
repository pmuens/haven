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

#![allow(non_snake_case)]

#[macro_use]
extern crate ndarray;
extern crate rand;

use ndarray::{arr1, stack, Array1, Array2, Axis};
use rand::distributions::Standard;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

/// Creates a new vector of length `(vector.len() * range)` and computes the powers of
/// 2 in `modulus` from range `0` to `range` for each value in `vector`.
fn powers_of_2(modulus: usize, range: usize, vector: &Array1<usize>) -> Array1<usize> {
    let vector_vec: Vec<usize> = vector.to_vec();
    let mut result_vec: Vec<usize> = Vec::with_capacity(vector_vec.len() * range);

    for elem in vector_vec.iter() {
        for i in 0..range {
            result_vec.push((elem * (2usize.pow(i as u32))) % modulus);
        }
    }

    arr1(&result_vec)
}

/// Secret Key information (includes the secret key vector).
pub struct SecretKey {
    t: Array1<usize>,
    s: Array1<usize>,
    v: Array1<usize>,
}

impl SecretKey {
    /// Generate a new SecretKey.
    pub fn new(q: usize, n: usize, l: usize) -> Self {
        // Init our Rng
        let mut rng: StdRng = StdRng::from_entropy();

        // t: Randomly generate a vector of length `n` with values in `(mod q)`
        let t_values: Vec<usize> = (0..n)
            .map(|_| rng.sample::<usize, Standard>(Standard) % q)
            .collect();
        let t = arr1(&t_values);

        // s: Create a vector of length `n + 1`, negate `t` values and add `q (mod q)`, prepend a `1`
        let mut s_values = Vec::with_capacity(t_values.len() + 1);
        s_values.push(1);
        let neg_t_values = t_values.into_iter().map(|t_val| (q - t_val) % q);
        s_values.extend(neg_t_values);
        let s = arr1(&s_values);

        // v: Transform `s` via `powers_of_2` function
        let v = powers_of_2(q, l, &s);

        SecretKey { t, s, v }
    }
}

/// Public Key information (includes the public key matrix).
pub struct PublicKey {
    B: Array2<usize>,
    e: Array1<usize>,
    b: Array1<usize>,
    A: Array2<usize>,
}

impl PublicKey {
    /// Generate a new PublicKey.
    pub fn new(q: usize, n: usize, _x: usize, m: usize, t: Array1<usize>) -> Self {
        // Init our Rng
        let mut rng: StdRng = StdRng::from_entropy();

        // B: Generate an `n x m` matrix of random values in `(mod q)`
        let B_values: Vec<usize> = (0..n * m)
            .map(|_| rng.sample::<usize, Standard>(Standard) % q)
            .collect();
        let B: Array2<usize> = Array2::from_shape_vec((n, m), B_values).unwrap();

        // e: Randomly generate a vector of length `m` with values in `(mod q)`
        // TODO: Use `x` to draw random variables from distribution
        let e_values: Vec<usize> = (0..m)
            .map(|_| rng.sample::<usize, Standard>(Standard) % q)
            .collect();
        let e = arr1(&e_values);

        // b: `b = t dot B + e (mod q)`
        let b: Array1<usize> = ((t.dot(&B) % q) + &e) % q;

        // A: `A` has `m + 1` rows, the first row is `b`, the rest is `B`
        // TODO: We need to turn the vector into a matrix so that we can use stacking. Not sure if
        //  this is the best solution...
        let b_2d: Array2<usize> = Array2::from_shape_vec((1, m), b.to_vec()).unwrap();
        let A = stack(Axis(0), &[b_2d.view(), B.view()]).unwrap();

        PublicKey { B, e, b, A }
    }
}

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
    /// Automatically generated when instantiated via `new`
    secret_key: SecretKey,
    /// Automatically generated when instantiated via `new`
    public_key: PublicKey,
}

impl GSW {
    /// Create a new `GSW` instance with the given security parameters.
    /// **NOTE:** This will automatically generate a `PublicKey` and `SecretKey` key pair.
    pub fn new(q: usize, n: usize, x: usize) -> Self {
        let m: usize = n * (q as f64).log2() as usize; // TODO: Check if `log2` is sufficient here
        let l: usize = (q as f64).log2().floor() as usize + 1;
        let N: usize = (n + 1) * l;
        let secret_key = SecretKey::new(q, n, l);
        let public_key = PublicKey::new(q, n, x, m, secret_key.t.clone());
        GSW {
            q,
            n,
            x,
            m,
            l,
            N,
            secret_key,
            public_key,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{powers_of_2, GSW};
    use ndarray::arr1;

    fn create_gsw() -> GSW {
        let q = 65536;
        let n = 3;
        let x = 8;
        GSW::new(q, n, x)
    }

    #[test]
    fn gsw_security_params() {
        let gsw = create_gsw();

        assert_eq!(gsw.q, 65536);
        assert_eq!(gsw.n, 3);
        assert_eq!(gsw.x, 8);
        assert_eq!(gsw.m, 48);
        assert_eq!(gsw.l, 17);
        assert_eq!(gsw.N, 68);
    }

    #[test]
    fn gsw_keygen_sk() {
        let gsw = create_gsw();
        let sk = &gsw.secret_key;

        let t = sk.t.to_owned();
        let t_vec = t.to_vec();
        let s = sk.s.to_owned();
        let s_vec = s.to_vec();
        let v = sk.v.to_owned();
        let v_vec = v.to_vec();

        // Check the vector shapes
        assert_eq!(t.shape(), [gsw.n]);
        assert_eq!(s.shape(), [gsw.n + 1]);
        assert_eq!(v.shape(), [s.shape()[0] * gsw.l]);

        // `t`
        // All values in `t` must be positive and in `(mod q)`
        assert!(t_vec.iter().all(|val| val >= &0 && val <= &gsw.q));

        // `s`
        // The first value of `s` must be 1
        assert_eq!(s_vec[0], 1);
        // All elements of `s` need to be `-t[i] + q (mod q)` (except the first one which is the `1`)
        assert!(s_vec
            .iter()
            .skip(1)
            .enumerate()
            .all(|(i, val)| val >= &0 && val == &(gsw.q - t_vec[i])));

        // `v`
        // The first element must 1 because the first element of `s` is 1 and `1 * (2 ** 0) % q = 1`
        assert_eq!(v_vec[0], 1);
        // The "last" computed element for the first value in `s` must be 0 because `1 * (2 ** l - 1) % q = 0`
        assert_eq!(v_vec[gsw.l - 1], 0);
        // All values of `v` should be positive and in `(mod q)`
        assert!(v_vec.iter().all(|val| val >= &0 && val <= &gsw.q));
    }

    #[test]
    fn gsw_keygen_pk() {
        let gsw = create_gsw();
        let pk = &gsw.public_key;

        let B = pk.B.to_owned();
        let B_vec = B.clone().into_raw_vec();
        let e = pk.e.to_owned();
        let e_vec = e.to_vec();
        let b = pk.b.to_owned();
        let b_vec = b.to_vec();
        let A = pk.A.to_owned();

        // Check the shapes
        assert_eq!(B.shape(), [gsw.n, gsw.m]);
        assert_eq!(e.shape(), [gsw.m]);
        assert_eq!(b.shape(), [gsw.m]);
        assert_eq!(A.shape(), [gsw.n + 1, gsw.m]);

        // `B`
        // All values in `B `should be positive and in `(mod q)`
        assert!(B_vec.iter().all(|val| val >= &0 && val <= &gsw.q));

        // `e`
        // All values in `e` should be positive and in `(mod q)`
        assert!(e_vec.iter().all(|val| val >= &0 && val <= &gsw.q));

        // `b`
        // All values in `b` should be positive and in `(mod q)`
        assert!(b_vec.iter().all(|val| val >= &0 && val <= &gsw.q));

        // `A`
        // The first row should be `b`
        assert_eq!(A.slice(s![0..1, ..]).to_owned().into_raw_vec(), b_vec);
        // The last rows should be `B`
        assert_eq!(A.slice(s![1..gsw.n + 1, ..]).view(), B.view());
    }

    #[test]
    fn powers_of_two() {
        let vector_vec = vec![1, 60782, 26640, 17805];
        let vector = arr1(&vector_vec);

        let result = powers_of_2(65536, 17, &vector);
        assert_eq!(result.shape(), [68]);
        assert_eq!(
            result.to_vec(),
            vec![
                1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 0,
                60782, 56028, 46520, 27504, 55008, 44480, 23424, 46848, 28160, 56320, 47104, 28672,
                57344, 49152, 32768, 0, 0, 26640, 53280, 41024, 16512, 33024, 512, 1024, 2048,
                4096, 8192, 16384, 32768, 0, 0, 0, 0, 0, 17805, 35610, 5684, 11368, 22736, 45472,
                25408, 50816, 36096, 6656, 13312, 26624, 53248, 40960, 16384, 32768, 0
            ]
        );
    }
}
