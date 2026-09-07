//! Meiosis, crossing and double-haploid production — the isqg port.
//!
//! **This module never calls an RNG (DECISION-012).** R draws every random
//! quantity, in isqg's exact order, and passes the results in:
//!
//! ```text
//! per meiosis event, per chromosome ascending:
//!     n_x       ~ rpois(1, L)              L = LAST map position, in Morgans
//!     chiasmata ~ sort(runif(n_x, 0, L))   not drawn when n_x == 0
//!     flip      ~ rbinom(1, 1, 0.5)        ALWAYS drawn, even when n_x == 0
//! ```
//!
//! Everything here is the deterministic remainder: the XOR chain that turns
//! chiasmata into an ancestry mask, and the assembly of gametes into progeny.

use crate::genome::{write_genotype, Bits, GenomeLayout};
use extendr_api::prelude::*;

// R represents NA_integer_ as i32::MIN.
const R_NA_INT: i32 = i32::MIN;

/// Number of map positions `<= chiasma` — isqg's `std::upper_bound`.
///
/// No tolerance is applied, and none may be: an epsilon here is a parity break
/// by construction.
#[inline]
pub fn breaks_at(positions: &[f64], chiasma: f64) -> usize {
    positions.partition_point(|&p| p <= chiasma)
}

/// Ancestry mask for ONE chromosome. Bit `j` set => take locus `j` from `cis`.
///
/// Each chiasma toggles every locus at or downstream of `breaks_at`. Note the
/// bound is `j >= breaks`, not `j > breaks`: isqg toggles bit indices
/// `[0, n-1-breaks]` in its reversed layout, which maps to ascending ranks
/// `[breaks, n-1]`.
///
/// The chiasmata are consumed exactly as supplied — never sorted, deduplicated
/// or filtered. XOR is commutative so order is irrelevant, and duplicate draws
/// are meant to cancel.
pub fn chromosome_mask(positions: &[f64], chiasmata: &[f64], flip: bool) -> Bits {
    let mut mask = Bits::zeros(positions.len());
    for &x in chiasmata {
        mask.toggle_from(breaks_at(positions, x));
    }
    if flip {
        mask.flip_all();
    }
    mask
}

/// The pre-drawn randomness for a run, decoded from R's flat vectors.
///
/// Ragged data (chiasmata counts vary per chromosome per event) is passed as
/// one concatenated `chiasmata` slice plus a `counts` slice, ordered
/// `counts[event * n_chr + chr]` — the same nesting as the R draw loop, so the
/// two cannot drift apart. `offsets` is the prefix sum, computed once; nothing
/// else in the crate indexes `counts` directly.
pub struct EventStream<'a> {
    chiasmata: &'a [f64],
    flips: &'a [i32],
    offsets: Vec<usize>,
    n_chr: usize,
    n_events: usize,
}

impl<'a> EventStream<'a> {
    /// # Panics
    /// On any inconsistency between `chiasmata`, `counts` and `flips`. The
    /// crate convention is to validate in R, but a `counts`/`chiasmata`
    /// mismatch silently shifts every subsequent block and yields
    /// plausible-looking output — exactly the failure mode of commit 4167402.
    /// These are the only assertions in the module.
    pub fn new(chiasmata: &'a [f64], counts: &'a [i32], flips: &'a [i32], n_chr: usize) -> Self {
        assert!(n_chr > 0, "n_chr must be positive");
        assert_eq!(
            counts.len(),
            flips.len(),
            "counts has {} entries but flips has {}",
            counts.len(),
            flips.len()
        );
        assert_eq!(
            counts.len() % n_chr,
            0,
            "counts length {} is not a multiple of n_chr {}",
            counts.len(),
            n_chr
        );

        let mut offsets = Vec::with_capacity(counts.len() + 1);
        let mut acc = 0usize;
        offsets.push(0);
        for &c in counts {
            assert!(
                c >= 0 && c != R_NA_INT,
                "crossover counts must be non-negative and non-NA, got {}",
                c
            );
            acc += c as usize;
            offsets.push(acc);
        }
        assert_eq!(
            acc,
            chiasmata.len(),
            "counts sum to {} but {} chiasmata were supplied",
            acc,
            chiasmata.len()
        );
        for &f in flips {
            assert!(f != R_NA_INT, "flips must not contain NA");
        }

        EventStream {
            chiasmata,
            flips,
            offsets,
            n_chr,
            n_events: counts.len() / n_chr,
        }
    }

    #[inline]
    pub fn n_events(&self) -> usize {
        self.n_events
    }

    #[inline]
    fn slot(&self, event: usize, chr: usize) -> usize {
        event * self.n_chr + chr
    }

    #[inline]
    pub fn block(&self, event: usize, chr: usize) -> &'a [f64] {
        let k = self.slot(event, chr);
        &self.chiasmata[self.offsets[k]..self.offsets[k + 1]]
    }

    #[inline]
    pub fn flip(&self, event: usize, chr: usize) -> bool {
        self.flips[self.slot(event, chr)] != 0
    }
}

/// Whole-genome ancestry mask for one meiosis event, chromosomes ascending.
pub fn genome_mask(layout: &GenomeLayout, events: &EventStream, event: usize) -> Bits {
    let mut mask = Bits::zeros(layout.n_loci());
    for c in 0..layout.n_chr() {
        let chr_mask = chromosome_mask(
            layout.chr_positions(c),
            events.block(event, c),
            events.flip(event, c),
        );
        let base = layout.chr_offset(c);
        for j in 0..chr_mask.len() {
            mask.set(base + j, chr_mask.get(j));
        }
    }
    mask
}

/// Draw a gamete: `(mask & cis) | (!mask & trans)`.
pub fn recombine(cis: &Bits, trans: &Bits, mask: &Bits) -> Bits {
    let mut out = Bits::zeros(cis.len());
    for j in 0..cis.len() {
        out.set(
            j,
            if mask.get(j) {
                cis.get(j)
            } else {
                trans.get(j)
            },
        );
    }
    out
}

/// The mating designs. They differ only in how many meiosis events each
/// progeny consumes and how the resulting gametes are paired.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Design {
    Cross,
    SelfCross,
    Dh,
}

impl Design {
    pub fn from_str(s: &str) -> Design {
        match s {
            "dh" => Design::Dh,
            "selfcross" => Design::SelfCross,
            _ => Design::Cross,
        }
    }

    /// `dh` makes one gamete and doubles it; the others make two.
    pub fn events_per_progeny(self) -> usize {
        match self {
            Design::Dh => 1,
            _ => 2,
        }
    }
}

/// Produce `n_prog` progeny, returning a loci-major `-1/0/1` genotype buffer
/// (`out[locus * n_prog + prog]`).
///
/// Events are consumed progeny-major, parent-minor — isqg generates progeny one
/// at a time, drawing the first parent's meiosis then the second's, rather than
/// all of one parent's then all of the other's. For `SelfCross`, pass the same
/// parent twice (isqg selfs against a deep copy, which is still two independent
/// meioses).
#[allow(clippy::too_many_arguments)]
pub fn mate_haplotypes(
    layout: &GenomeLayout,
    p1_cis: &Bits,
    p1_trans: &Bits,
    p2_cis: &Bits,
    p2_trans: &Bits,
    events: &EventStream,
    design: Design,
    n_prog: usize,
) -> Vec<(Bits, Bits)> {
    let per = design.events_per_progeny();
    let mut progeny = Vec::with_capacity(n_prog);

    for i in 0..n_prog {
        let pair = if design == Design::Dh {
            let g = recombine(p1_cis, p1_trans, &genome_mask(layout, events, i));
            (g.clone(), g)
        } else {
            let a = recombine(p1_cis, p1_trans, &genome_mask(layout, events, per * i));
            let b = recombine(p2_cis, p2_trans, &genome_mask(layout, events, per * i + 1));
            (a, b)
        };
        progeny.push(pair);
    }

    progeny
}

/// As [`mate_haplotypes`], projected onto a loci-major `-1/0/1` genotype buffer.
///
/// The projection is lossy: it cannot distinguish the two phases of a
/// heterozygote. Callers that will breed from the progeny must keep the
/// haplotypes instead.
#[allow(clippy::too_many_arguments)]
pub fn mate_core(
    layout: &GenomeLayout,
    p1_cis: &Bits,
    p1_trans: &Bits,
    p2_cis: &Bits,
    p2_trans: &Bits,
    events: &EventStream,
    design: Design,
    n_prog: usize,
) -> Vec<i32> {
    let progeny = mate_haplotypes(
        layout, p1_cis, p1_trans, p2_cis, p2_trans, events, design, n_prog,
    );
    let mut out = vec![0i32; layout.n_loci() * n_prog];
    for (i, (cis, trans)) in progeny.iter().enumerate() {
        write_genotype(cis, trans, &mut out, i, n_prog);
    }
    out
}

// ---------------------------------------------------------------------------
// extendr boundary — no logic beyond argument marshalling.
// ---------------------------------------------------------------------------

/// Produce progeny genotypes from pre-drawn meiosis randomness.
///
/// Chiasmata are concatenated and delimited by `counts`, ordered
/// `counts[event * n_chr + chr]`. Events run progeny-major: for "cross" and
/// "selfcross", event `2i` is parent 1's gamete for progeny `i` and `2i+1` is
/// parent 2's; for "dh" there is one event per progeny.
///
/// @param loci_per_chr Integer vector of loci counts, chromosomes ascending.
/// @param positions    Map positions in Morgans, concatenated, ascending within chromosome.
/// @param p1_cis       Parent 1 cis strand as a '0'/'1' string, ascending.
/// @param p1_trans     Parent 1 trans strand.
/// @param p2_cis       Parent 2 cis strand (same as p1 for selfcross/dh).
/// @param p2_trans     Parent 2 trans strand.
/// @param chiasmata    Concatenated crossover positions.
/// @param counts       Crossovers per (event, chromosome).
/// @param flips        0/1 strand-choice per (event, chromosome).
/// @param design       "cross", "selfcross", or "dh".
/// @param n_prog       Number of progeny.
/// @return Integer vector length n_loci * n_prog, loci-major (-1/0/1).
/// @noRd
#[extendr]
#[allow(clippy::too_many_arguments)]
pub fn meiosis_core(
    loci_per_chr: &[i32],
    positions: &[f64],
    p1_cis: &str,
    p1_trans: &str,
    p2_cis: &str,
    p2_trans: &str,
    chiasmata: &[f64],
    counts: &[i32],
    flips: &[i32],
    design: &str,
    n_prog: i32,
) -> Vec<i32> {
    let layout = GenomeLayout::new(loci_per_chr, positions);
    let events = EventStream::new(chiasmata, counts, flips, layout.n_chr());
    mate_core(
        &layout,
        &Bits::from_code(p1_cis, '1'),
        &Bits::from_code(p1_trans, '1'),
        &Bits::from_code(p2_cis, '1'),
        &Bits::from_code(p2_trans, '1'),
        &events,
        Design::from_str(design),
        n_prog as usize,
    )
}

/// Raw ancestry masks, one '0'/'1' string per meiosis event.
///
/// Parity hook mirroring isqg's `spc$gamete()`. Unlike the -1/0/1 genotype —
/// a lossy 3-valued projection that cannot distinguish a cis/trans swap — the
/// mask pins the recombination algorithm on its own.
///
/// @param loci_per_chr Integer vector of loci counts, chromosomes ascending.
/// @param positions    Map positions in Morgans, concatenated.
/// @param chiasmata    Concatenated crossover positions.
/// @param counts       Crossovers per (event, chromosome).
/// @param flips        0/1 strand-choice per (event, chromosome).
/// @return Character vector, one mask per event, ascending map order.
/// @noRd
#[extendr]
pub fn gamete_masks_core(
    loci_per_chr: &[i32],
    positions: &[f64],
    chiasmata: &[f64],
    counts: &[i32],
    flips: &[i32],
) -> Vec<String> {
    let layout = GenomeLayout::new(loci_per_chr, positions);
    let events = EventStream::new(chiasmata, counts, flips, layout.n_chr());
    (0..events.n_events())
        .map(|e| genome_mask(&layout, &events, e).to_code('1', '0'))
        .collect()
}

/// Progeny haplotypes from pre-drawn meiosis randomness.
///
/// Same arguments as [`meiosis_core`], but returns the phased strands rather
/// than the `-1/0/1` genotype: element `2i` is progeny `i`'s first strand and
/// `2i + 1` its second, each a '0'/'1' string in ascending map order.
///
/// This is what breeding programmes need. A genotype loses the phase of every
/// heterozygote, so reconstructing strands from it would make an F1 -- which is
/// heterozygous everywhere -- come out with one all-allele-1 strand and one
/// all-allele-2 strand, and every later generation would then recombine
/// haplotypes the founders never had.
///
/// @param loci_per_chr Integer vector of loci counts, chromosomes ascending.
/// @param positions    Map positions in Morgans, concatenated.
/// @param p1_cis       Parent 1 first strand as a '0'/'1' string, ascending.
/// @param p1_trans     Parent 1 second strand.
/// @param p2_cis       Parent 2 first strand (same as p1 for selfcross/dh).
/// @param p2_trans     Parent 2 second strand.
/// @param chiasmata    Concatenated crossover positions.
/// @param counts       Crossovers per (event, chromosome).
/// @param flips        0/1 strand-choice per (event, chromosome).
/// @param design       "cross", "selfcross", or "dh".
/// @param n_prog       Number of progeny.
/// @return Character vector of length 2 * n_prog.
/// @noRd
#[extendr]
#[allow(clippy::too_many_arguments)]
pub fn mate_haplotypes_core(
    loci_per_chr: &[i32],
    positions: &[f64],
    p1_cis: &str,
    p1_trans: &str,
    p2_cis: &str,
    p2_trans: &str,
    chiasmata: &[f64],
    counts: &[i32],
    flips: &[i32],
    design: &str,
    n_prog: i32,
) -> Vec<String> {
    let layout = GenomeLayout::new(loci_per_chr, positions);
    let events = EventStream::new(chiasmata, counts, flips, layout.n_chr());
    let progeny = mate_haplotypes(
        &layout,
        &Bits::from_code(p1_cis, '1'),
        &Bits::from_code(p1_trans, '1'),
        &Bits::from_code(p2_cis, '1'),
        &Bits::from_code(p2_trans, '1'),
        &events,
        Design::from_str(design),
        n_prog as usize,
    );
    let mut out = Vec::with_capacity(progeny.len() * 2);
    for (cis, trans) in progeny {
        out.push(cis.to_code('1', '0'));
        out.push(trans.to_code('1', '0'));
    }
    out
}

extendr_module! {
    mod meiosis;
    fn meiosis_core;
    fn mate_haplotypes_core;
    fn gamete_masks_core;
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pos() -> Vec<f64> {
        vec![0.0, 0.5, 1.0, 1.5, 2.0]
    }

    #[test]
    fn breaks_counts_positions_at_or_below_the_chiasma() {
        let p = pos();
        assert_eq!(breaks_at(&p, -0.1), 0);
        assert_eq!(breaks_at(&p, 0.25), 1);
        assert_eq!(breaks_at(&p, 0.75), 2);
        assert_eq!(breaks_at(&p, 9.0), 5);
    }

    #[test]
    fn a_chiasma_exactly_on_a_marker_leaves_it_upstream() {
        // isqg uses upper_bound (<=), so a marker at the breakpoint stays with
        // the upstream segment. Measure-zero under runif, but pinned here
        // because any grid-snapped breakpoint hits it every time.
        let p = pos();
        assert_eq!(breaks_at(&p, 0.5), 2);
        assert_eq!(breaks_at(&p, 2.0), 5); // == n_loci: toggle must be a no-op
    }

    #[test]
    fn zero_chiasmata_gives_all_zero_or_all_one_by_the_flip() {
        let p = pos();
        assert_eq!(chromosome_mask(&p, &[], false).to_code('1', '0'), "00000");
        // The Bernoulli is drawn even with no crossovers; dropping it would
        // desynchronise every later draw.
        assert_eq!(chromosome_mask(&p, &[], true).to_code('1', '0'), "11111");
    }

    #[test]
    fn one_chiasma_toggles_the_downstream_segment_inclusively() {
        let p = pos();
        // breaks_at(0.75) == 2, so loci 2,3,4 flip — NOT 3,4.
        assert_eq!(
            chromosome_mask(&p, &[0.75], false).to_code('1', '0'),
            "00111"
        );
    }

    #[test]
    fn a_chiasma_upstream_of_the_first_marker_toggles_everything() {
        // Reachable whenever a chromosome's first marker is not at 0, since
        // chiasmata are drawn on (0, L).
        let p = vec![0.15, 0.5, 1.0];
        assert_eq!(chromosome_mask(&p, &[0.05], false).to_code('1', '0'), "111");
    }

    #[test]
    fn two_chiasmata_restore_the_intervening_segment() {
        let p = pos();
        let m = chromosome_mask(&p, &[0.75, 1.75], false);
        assert_eq!(m.to_code('1', '0'), "00110");
    }

    #[test]
    fn duplicate_chiasmata_cancel_and_order_is_irrelevant() {
        let p = pos();
        assert_eq!(
            chromosome_mask(&p, &[0.75, 0.75], false).to_code('1', '0'),
            "00000"
        );
        assert_eq!(
            chromosome_mask(&p, &[1.75, 0.75], false),
            chromosome_mask(&p, &[0.75, 1.75], false)
        );
    }

    #[test]
    fn single_locus_chromosome() {
        let p = vec![0.8];
        assert_eq!(chromosome_mask(&p, &[], false).to_code('1', '0'), "0");
        assert_eq!(chromosome_mask(&p, &[0.4], false).to_code('1', '0'), "1");
        // A chiasma at or past the only marker gives breaks == n: no toggle.
        assert_eq!(chromosome_mask(&p, &[0.8], false).to_code('1', '0'), "0");
    }

    #[test]
    fn event_stream_slices_ragged_blocks_in_order() {
        // 2 events x 2 chromosomes; counts are (e,c)-ordered.
        let chiasmata = [0.1, 0.2, 0.3, 0.4];
        let counts = [1i32, 0, 2, 1];
        let flips = [0i32, 1, 1, 0];
        let ev = EventStream::new(&chiasmata, &counts, &flips, 2);

        assert_eq!(ev.n_events(), 2);
        assert_eq!(ev.block(0, 0), &[0.1]);
        assert_eq!(ev.block(0, 1), &[] as &[f64]);
        assert_eq!(ev.block(1, 0), &[0.2, 0.3]);
        assert_eq!(ev.block(1, 1), &[0.4]);
        assert!(!ev.flip(0, 0));
        assert!(ev.flip(0, 1));
        assert!(!ev.flip(1, 1));
    }

    #[test]
    #[should_panic(expected = "counts sum to")]
    fn event_stream_rejects_a_count_chiasmata_mismatch() {
        EventStream::new(&[0.1, 0.2], &[1i32, 0], &[0i32, 0], 2);
    }

    #[test]
    #[should_panic(expected = "non-negative")]
    fn event_stream_rejects_negative_counts() {
        EventStream::new(&[0.1], &[-1i32, 2], &[0i32, 0], 2);
    }

    #[test]
    fn genome_mask_concatenates_chromosomes_ascending() {
        let layout = GenomeLayout::new(&[3, 2], &[0.0, 0.5, 1.0, 0.0, 1.0]);
        // chr0: no crossover, no flip -> 000. chr1: no crossover, flip -> 11.
        let ev = EventStream::new(&[], &[0i32, 0], &[0i32, 1], 2);
        assert_eq!(genome_mask(&layout, &ev, 0).to_code('1', '0'), "00011");
    }

    #[test]
    fn recombine_picks_cis_where_the_mask_is_set() {
        let cis = Bits::from_code("1111", '1');
        let trans = Bits::from_code("0000", '1');
        let mask = Bits::from_code("1010", '1');
        assert_eq!(recombine(&cis, &trans, &mask).to_code('1', '0'), "1010");
    }

    #[test]
    fn dh_progeny_are_fully_homozygous() {
        let layout = GenomeLayout::new(&[4], &[0.0, 0.3, 0.6, 0.9]);
        let ev = EventStream::new(&[0.45], &[1i32], &[0i32], 1);
        let p_cis = Bits::from_code("1100", '1');
        let p_trans = Bits::from_code("0011", '1');
        let out = mate_core(
            &layout,
            &p_cis,
            &p_trans,
            &p_cis,
            &p_trans,
            &ev,
            Design::Dh,
            1,
        );
        assert!(out.iter().all(|&g| g != 0), "DH produced a heterozygote");
    }

    #[test]
    fn cross_takes_cis_from_the_first_parent() {
        let layout = GenomeLayout::new(&[2], &[0.0, 1.0]);
        // No crossovers, no flips: each gamete is the parent's trans strand
        // (mask all zero), so the phase is fully determined by parent order.
        let ev = EventStream::new(&[], &[0i32, 0], &[0i32, 0], 1);
        let p1_cis = Bits::from_code("11", '1');
        let p1_trans = Bits::from_code("11", '1');
        let p2_cis = Bits::from_code("00", '1');
        let p2_trans = Bits::from_code("00", '1');

        let ab = mate_core(
            &layout,
            &p1_cis,
            &p1_trans,
            &p2_cis,
            &p2_trans,
            &ev,
            Design::Cross,
            1,
        );
        // p1 contributes A, p2 contributes a => heterozygous everywhere.
        assert_eq!(ab, vec![0, 0]);
    }

    #[test]
    fn events_are_consumed_progeny_major() {
        // Two progeny, one chromosome, no crossovers. Flips per event:
        // [p1 of prog0, p2 of prog0, p1 of prog1, p2 of prog1] = [0,0,1,1].
        // Progeny 0 takes both trans strands; progeny 1 takes both cis.
        let layout = GenomeLayout::new(&[2], &[0.0, 1.0]);
        let ev = EventStream::new(&[], &[0i32; 4], &[0i32, 0, 1, 1], 1);
        let cis = Bits::from_code("11", '1');
        let trans = Bits::from_code("00", '1');
        let out = mate_core(&layout, &cis, &trans, &cis, &trans, &ev, Design::Cross, 2);
        // Loci-major, n_prog = 2: [locus0_prog0, locus0_prog1, locus1_prog0, ...]
        assert_eq!(out, vec![-1, 1, -1, 1]);
    }

    #[test]
    fn dh_haplotypes_are_segments_of_the_parent_strands_not_the_mask() {
        // A parent heterozygous at every locus, with an ALTERNATING allele
        // pattern so the strands cannot be confused with a mask. Each doubled
        // haploid must be a segmentwise copy of one parental strand.
        //
        // This is the shape of a real bug: reconstructing progeny strands from
        // a -1/0/1 genotype loses the phase of every heterozygote and yields
        // progeny that carry the ancestry mask instead of parental alleles.
        let layout = GenomeLayout::new(&[8], &[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]);
        let cis = Bits::from_code("10101010", '1');
        let trans = Bits::from_code("01010101", '1');

        // One crossover at 0.35 => breaks at index 4; no terminal flip.
        let ev = EventStream::new(&[0.35], &[1i32], &[0i32], 1);
        let prog = mate_haplotypes(&layout, &cis, &trans, &cis, &trans, &ev, Design::Dh, 1);

        let (a, b) = &prog[0];
        assert_eq!(a, b, "a doubled haploid must have identical strands");
        // Loci 0-3 from trans, loci 4-7 from cis.
        assert_eq!(a.to_code('1', '0'), "01011010");
    }

    #[test]
    fn mate_core_agrees_with_mate_haplotypes() {
        let layout = GenomeLayout::new(&[4], &[0.0, 0.3, 0.6, 0.9]);
        let ev = EventStream::new(&[0.45, 0.2], &[1i32, 1], &[0i32, 1], 1);
        let p_cis = Bits::from_code("1100", '1');
        let p_trans = Bits::from_code("0110", '1');

        let hap = mate_haplotypes(
            &layout,
            &p_cis,
            &p_trans,
            &p_cis,
            &p_trans,
            &ev,
            Design::Cross,
            1,
        );
        let geno = mate_core(
            &layout,
            &p_cis,
            &p_trans,
            &p_cis,
            &p_trans,
            &ev,
            Design::Cross,
            1,
        );
        let mut expected = vec![0i32; 4];
        write_genotype(&hap[0].0, &hap[0].1, &mut expected, 0, 1);
        assert_eq!(geno, expected);
    }

    #[test]
    fn design_parsing_and_event_budget() {
        assert_eq!(Design::from_str("dh"), Design::Dh);
        assert_eq!(Design::from_str("selfcross"), Design::SelfCross);
        assert_eq!(Design::from_str("cross"), Design::Cross);
        assert_eq!(Design::from_str("anything else"), Design::Cross);
        assert_eq!(Design::Dh.events_per_progeny(), 1);
        assert_eq!(Design::Cross.events_per_progeny(), 2);
        assert_eq!(Design::SelfCross.events_per_progeny(), 2);
    }
}
