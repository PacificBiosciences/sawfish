#![allow(dead_code)]

use std::iter::Sum;

use num::{Float, NumCast};

pub fn error_prob_to_phred(prob: f64) -> f64 {
    -10f64 * prob.log10().max(f64::MIN_10_EXP as f64)
}

pub fn ln_error_prob_to_phred(ln_prob: f64) -> f64 {
    -10f64 * (ln_prob / std::f64::consts::LN_10).max(f64::MIN_10_EXP as f64)
}

pub fn error_prob_to_qphred(prob: f64) -> i32 {
    error_prob_to_phred(prob).round() as i32
}

pub fn ln_error_prob_to_qphred(ln_prob: f64) -> i32 {
    ln_error_prob_to_phred(ln_prob).round() as i32
}

/// Standardize ln-transformed unnormalized prob distro input
///
/// Returns the index of the most probable component
///
pub fn normalize_ln_distro<F: Float>(x: &mut [F]) -> Option<usize> {
    if x.is_empty() {
        return None;
    }

    let mut max_index = 0;
    let mut max_p = *x.first().unwrap();
    for (index, p) in x.iter().skip(1).enumerate() {
        if *p > max_p {
            max_p = *p;
            max_index = index + 1;
        }
    }

    let mut sum = NumCast::from(0).unwrap();
    for p in x.iter_mut() {
        *p = (*p - max_p).exp();
        sum = sum + *p;
    }

    for p in x.iter_mut() {
        *p = *p / sum;
    }

    Some(max_index)
}

/// Get the complement of pdf[index] from a normalized prob distro
///
/// As pdf[index] approaches 1, computing the complement as 1 - pdf[index] starts to significantly
/// degrade precision. Instead the value is found by summing the rest of the pdf.
///
pub fn get_complement_prob<F: Float + Sum<F>>(pdf: &[F], index: usize) -> F {
    pdf.iter()
        .enumerate()
        .filter(|(i, _)| *i != index)
        .map(|(_, p)| *p)
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_prob_to_qphred() {
        let q = error_prob_to_qphred(0.001);
        assert_eq!(q, 30);
    }

    #[test]
    fn test_ln_error_prob_to_qphred() {
        let q = ln_error_prob_to_qphred(0.001.ln());
        assert_eq!(q, 30);
    }

    #[test]
    fn test_normalize_ln_distro() {
        let x = [0.001, 0.001, 0.002, 0.001];
        let mut x = x.into_iter().map(|x: f64| x.ln()).collect::<Vec<_>>();

        let max_index = normalize_ln_distro(&mut x);
        assert_eq!(max_index, Some(2));
        approx::assert_ulps_eq!(x[0], 0.2, max_ulps = 4);
        approx::assert_ulps_eq!(x[2], 0.4, max_ulps = 4);
    }

    #[test]
    fn test_get_complement_prob() {
        let x: [f32; 3] = [0.9999999, 0.00000005, 0.00000005];
        let x0c = get_complement_prob(&x, 0);
        approx::assert_ulps_eq!(x0c, 0.0000001, max_ulps = 4);
    }
}
