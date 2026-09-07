//! Bit-packed haplotypes and genome layout for the isqg meiosis port.
//!
//! Bits are stored in **ascending map-position order**: bit `j` is the locus
//! with rank `j` in the sorted map. isqg itself stores the reverse
//! (`bit = n-1-rank`) and calls `std::reverse` on export; storing ascending
//! yields the same observable output with no reversal, and sidesteps isqg's
//! mis-assignment of bit indices for tied map positions.
//!
//! Nothing in this module calls an RNG (DECISION-012).

/// Number of 64-bit words needed to hold `n` bits.
///
/// `usize::div_ceil` is 1.73; this crate's MSRV is 1.65.
#[inline]
pub const fn n_words(n: usize) -> usize {
    (n + 63) >> 6
}

/// A fixed-length bit vector.
///
/// INVARIANT: `words.len() == n_words(n)`, and every bit at index `>= n` in the
/// final word is zero. Complementing (`flip_all`) is what breaks this if the
/// tail is not re-cleared, and a derived `PartialEq` would then report two
/// logically equal values as unequal.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Bits {
    words: Vec<u64>,
    n: usize,
}

impl Bits {
    pub fn zeros(n: usize) -> Self {
        Bits {
            words: vec![0u64; n_words(n)],
            n,
        }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.n
    }

    #[inline]
    pub fn get(&self, i: usize) -> bool {
        debug_assert!(i < self.n);
        debug_assert!(self.tail_is_clean());
        self.words[i >> 6] >> (i & 63) & 1 == 1
    }

    #[inline]
    pub fn set(&mut self, i: usize, v: bool) {
        debug_assert!(i < self.n);
        let (w, bit) = (i >> 6, 1u64 << (i & 63));
        if v {
            self.words[w] |= bit;
        } else {
            self.words[w] &= !bit;
        }
    }

    /// XOR every bit in `[from, n)`.
    ///
    /// This is the one bulk operation meiosis needs. `from >= n` is a no-op,
    /// which is what makes `breaks == n_loci` fall out naturally instead of
    /// needing a special case — and avoids the `u64 << 64` that a shift-based
    /// implementation would hit there.
    pub fn toggle_from(&mut self, from: usize) {
        if from >= self.n {
            return;
        }
        let first = from >> 6;
        // `from & 63 <= 63`, so this shift is always well-defined.
        self.words[first] ^= !0u64 << (from & 63);
        for w in self.words.iter_mut().skip(first + 1) {
            *w = !*w;
        }
        self.clear_tail();
    }

    pub fn flip_all(&mut self) {
        for w in self.words.iter_mut() {
            *w = !*w;
        }
        self.clear_tail();
    }

    /// Parse a string of `one`/other characters, ascending: `code[j]` -> bit `j`.
    pub fn from_code(code: &str, one: char) -> Self {
        let mut b = Bits::zeros(code.chars().count());
        for (j, c) in code.chars().enumerate() {
            if c == one {
                b.set(j, true);
            }
        }
        b
    }

    /// Render ascending as a string, for direct comparison with isqg's
    /// `spc$gamete()` output.
    pub fn to_code(&self, one: char, zero: char) -> String {
        debug_assert!(self.tail_is_clean());
        (0..self.n)
            .map(|j| if self.get(j) { one } else { zero })
            .collect()
    }

    /// Re-zero bits beyond `n` in the final word. Must be called by every
    /// operation that can dirty them (`ones`, `toggle_from`, `flip_all`).
    #[inline]
    fn clear_tail(&mut self) {
        let rem = self.n & 63;
        // Guard is mandatory: `1u64 << 64` panics in debug and wraps in release.
        if rem != 0 {
            let last = self.words.len() - 1;
            self.words[last] &= (1u64 << rem) - 1;
        }
    }

    #[cfg(test)]
    fn tail_is_clean(&self) -> bool {
        let rem = self.n & 63;
        rem == 0 || self.words[self.words.len() - 1] >> rem == 0
    }

    #[cfg(not(test))]
    #[inline]
    fn tail_is_clean(&self) -> bool {
        true
    }
}

/// Chromosome boundaries and map positions for one species.
///
/// `positions` holds every locus, concatenated chromosome by chromosome in
/// ascending order; `chr_start` has `n_chr + 1` entries delimiting them.
pub struct GenomeLayout {
    chr_start: Vec<usize>,
    positions: Vec<f64>,
}

impl GenomeLayout {
    /// # Panics
    /// If the per-chromosome counts do not describe `positions`, or any
    /// chromosome is empty. isqg takes a chromosome's length from its last map
    /// position, which is undefined for an empty map, so R must never send one.
    pub fn new(loci_per_chr: &[i32], positions: &[f64]) -> Self {
        let mut chr_start = Vec::with_capacity(loci_per_chr.len() + 1);
        let mut acc = 0usize;
        chr_start.push(0);
        for &c in loci_per_chr {
            assert!(c > 0, "each chromosome needs at least one locus, got {}", c);
            acc += c as usize;
            chr_start.push(acc);
        }
        assert_eq!(
            acc,
            positions.len(),
            "loci_per_chr sums to {} but {} positions were supplied",
            acc,
            positions.len()
        );
        GenomeLayout {
            chr_start,
            positions: positions.to_vec(),
        }
    }

    #[inline]
    pub fn n_chr(&self) -> usize {
        self.chr_start.len() - 1
    }

    #[inline]
    pub fn n_loci(&self) -> usize {
        self.positions.len()
    }

    #[inline]
    pub fn chr_positions(&self, c: usize) -> &[f64] {
        &self.positions[self.chr_start[c]..self.chr_start[c + 1]]
    }

    #[inline]
    pub fn chr_offset(&self, c: usize) -> usize {
        self.chr_start[c]
    }
}

/// Write one individual's `-1/0/1` genotype into a loci-major buffer.
///
/// Coding matches isqg: `1` = both strands carry allele A, `0` = heterozygous,
/// `-1` = neither. Layout is `out[locus * n_prog + prog]`, the same row-major
/// convention as `numeric.rs` (`out[snp * n_samp + samp]`).
pub fn write_genotype(cis: &Bits, trans: &Bits, out: &mut [i32], prog: usize, n_prog: usize) {
    for j in 0..cis.len() {
        let (c, t) = (cis.get(j), trans.get(j));
        out[j * n_prog + prog] = if c && t {
            1
        } else if c || t {
            0
        } else {
            -1
        };
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // The whole trailing-bit bug class lives at these boundaries.
    const SIZES: [usize; 8] = [1, 37, 63, 64, 65, 127, 128, 129];

    #[test]
    fn all_ones_via_flip_all_keeps_the_tail_clean() {
        // The terminal Bernoulli flip is the production path that dirties the
        // trailing word, so it is the one that must stay clean at every width.
        for &n in &SIZES {
            let mut b = Bits::zeros(n);
            b.flip_all();
            assert!(b.tail_is_clean(), "dirty tail at n={}", n);
            assert!((0..n).all(|j| b.get(j)), "missing bit at n={}", n);
            assert_eq!(b.to_code('1', '0'), "1".repeat(n));
        }
    }

    #[test]
    fn flip_all_is_an_involution_and_keeps_tail_clean() {
        for &n in &SIZES {
            let mut b = Bits::zeros(n);
            b.set(0, true);
            if n > 1 {
                b.set(n - 1, true);
            }
            let original = b.clone();
            b.flip_all();
            assert!(b.tail_is_clean(), "dirty tail at n={}", n);
            b.flip_all();
            // Equality here is exactly what a dirty tail would break.
            assert_eq!(b, original, "not an involution at n={}", n);
        }
    }

    #[test]
    fn toggle_from_covers_the_suffix_only() {
        for &n in &SIZES {
            for from in 0..=n {
                let mut b = Bits::zeros(n);
                b.toggle_from(from);
                assert!(b.tail_is_clean(), "dirty tail at n={} from={}", n, from);
                for j in 0..n {
                    assert_eq!(b.get(j), j >= from, "n={} from={} j={}", n, from, j);
                }
            }
        }
    }

    #[test]
    fn toggle_from_at_or_past_n_is_a_noop() {
        // `breaks == n_loci` is reachable whenever a chiasma lands exactly on
        // the last marker, because isqg counts positions <= the chiasma.
        for &n in &SIZES {
            let mut b = Bits::zeros(n);
            b.toggle_from(n);
            b.toggle_from(n + 1);
            b.toggle_from(n + 1000);
            assert_eq!(b, Bits::zeros(n), "n={}", n);
        }
    }

    #[test]
    fn toggle_from_twice_cancels() {
        for &n in &SIZES {
            for from in [0, n / 2, n] {
                let mut b = Bits::zeros(n);
                b.toggle_from(from);
                b.toggle_from(from);
                assert_eq!(b, Bits::zeros(n), "n={} from={}", n, from);
            }
        }
    }

    #[test]
    fn code_round_trips_ascending() {
        let code = "1001110100001111011";
        let b = Bits::from_code(code, '1');
        assert_eq!(b.len(), code.len());
        assert_eq!(b.to_code('1', '0'), code);
        // Ascending: character 0 of the string is bit 0.
        assert!(b.get(0));
        assert!(!b.get(1));
    }

    #[test]
    fn genotype_uses_loci_major_layout() {
        // One AA locus in a field of aa, at an asymmetric index, pins both the
        // coding and the flat-index formula.
        let n_loci = 37;
        let n_prog = 5;
        let mut cis = Bits::zeros(n_loci);
        let mut trans = Bits::zeros(n_loci);
        cis.set(11, true);
        trans.set(11, true);
        cis.set(4, true); // het at locus 4

        let mut out = vec![0i32; n_loci * n_prog];
        write_genotype(&cis, &trans, &mut out, 3, n_prog);

        assert_eq!(out[11 * n_prog + 3], 1, "AA locus in the wrong cell");
        assert_eq!(out[4 * n_prog + 3], 0, "het locus in the wrong cell");
        assert_eq!(out[3], -1, "locus 0 of progeny 3");
        // A transposed write would have landed at prog * n_loci + locus.
        assert_eq!(out[3 * n_loci + 11], 0, "looks transposed");
    }

    #[test]
    fn genotype_is_symmetric_in_the_two_strands_only_for_homozygotes() {
        let mut cis = Bits::zeros(4);
        let trans = Bits::zeros(4);
        cis.set(1, true);
        let mut a = vec![0i32; 4];
        let mut b = vec![0i32; 4];
        write_genotype(&cis, &trans, &mut a, 0, 1);
        write_genotype(&trans, &cis, &mut b, 0, 1);
        // -1/0/1 cannot distinguish phase; that is why the parity harness also
        // compares raw gamete masks.
        assert_eq!(a, b);
    }

    #[test]
    #[should_panic(expected = "positions")]
    fn layout_rejects_mismatched_counts() {
        GenomeLayout::new(&[3, 2], &[0.0, 0.1, 0.2, 0.3]);
    }

    #[test]
    #[should_panic(expected = "at least one locus")]
    fn layout_rejects_empty_chromosome() {
        GenomeLayout::new(&[2, 0], &[0.0, 0.1]);
    }

    #[test]
    fn layout_slices_chromosomes_in_order() {
        let g = GenomeLayout::new(&[11, 1, 25], &vec![0.5; 37]);
        assert_eq!(g.n_chr(), 3);
        assert_eq!(g.n_loci(), 37);
        assert_eq!(g.chr_positions(0).len(), 11);
        assert_eq!(g.chr_positions(1).len(), 1);
        assert_eq!(g.chr_positions(2).len(), 25);
        assert_eq!(g.chr_offset(2), 12);
    }
}
