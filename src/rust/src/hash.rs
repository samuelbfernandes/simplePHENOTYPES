//! Stable content hash for pedigree keys (DECISION-024).
//!
//! FNV-1a with the 128-bit parameters of the FNV specification, over the UTF-8
//! bytes of a canonical text encoding built in R. Deterministic and independent
//! of R, package and dependency versions, so a pedigree key computed today is
//! the key the same individual gets after an upgrade.

use extendr_api::prelude::*;

const FNV128_OFFSET: u128 = 0x6c62272e07bb014262b821756295c58d;
const FNV128_PRIME: u128 = 0x0000000001000000000000000000013b;

/// FNV-1a-128 of a byte slice.
pub fn fnv1a_128(bytes: &[u8]) -> u128 {
    let mut h = FNV128_OFFSET;
    for &b in bytes {
        h ^= u128::from(b);
        h = h.wrapping_mul(FNV128_PRIME);
    }
    h
}

/// Hex (32 lower-case digits) FNV-1a-128 of each input string.
/// @noRd
#[extendr]
pub fn stable_hash_core(x: Vec<String>) -> Vec<String> {
    x.iter()
        .map(|s| format!("{:032x}", fnv1a_128(s.as_bytes())))
        .collect()
}

extendr_module! {
    mod hash;
    fn stable_hash_core;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_input_is_the_offset_basis() {
        assert_eq!(fnv1a_128(b""), FNV128_OFFSET);
    }

    #[test]
    fn distinct_inputs_differ_and_repeat_equally() {
        assert_ne!(fnv1a_128(b"a"), fnv1a_128(b"b"));
        assert_eq!(fnv1a_128(b"founder|A|P1"), fnv1a_128(b"founder|A|P1"));
        assert_ne!(fnv1a_128(b"founder|A|P1"), fnv1a_128(b"founder|B|P1"));
    }
}
