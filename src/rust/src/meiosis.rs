//! Meiosis, crossing and double-haploid production — the isqg port.
//!
//! **This module never calls an RNG (DECISION-012).** R draws every random
//! quantity, in isqg's exact order, and passes the results in:
//!
//! ```text
//! per meiosis event, per chromosome ascending:
//!     n_x       ~ rpois(1, L)              L = LAST map position, in Morgans
//!                                          (not the span last - first)
//!     chiasmata ~ sort(runif(n_x, 0, L))   not drawn when n_x == 0
//!     flip      ~ rbinom(1, 1, 0.5)        ALWAYS drawn, even when n_x == 0
//! ```
//!
//! Everything here is the deterministic remainder: the XOR chain that turns
//! chiasmata into an ancestry mask, and the assembly of gametes into progeny.
//!
//! The kernel does not care how the chiasmata were distributed: with the
//! default Poisson draws above, or with the optional two-pathway gamma
//! interference model (R's `.draw_meiosis_interference()`, SPEC-0020 item 3), R
//! hands it sorted positions in `[0, L]` and the same `counts`/`flips`
//! contract. `mate_many_core()` is the batched, integer-in/integer-out entry
//! point (SPEC-0020 item 2): many matings, one shared event stream, no strings.
//!
//! # Error handling
//!
//! **Nothing in this module panics on bad input.** On the toolchains where the
//! Rust panic runtime cannot unwind through R's frames (for example a
//! gcc-linked macOS build) a panic aborts the whole R process, and `extendr`
//! turns a returned `Err` into a panic unless the `result_list` feature is on.
//! So the crate enables `result_list`, every `#[extendr]` entry point returns
//! `Result<_, String>`, and the R wrappers (`R/extendr-wrappers.R`) re-raise the
//! error as an ordinary R error. Every inconsistency between the arguments —
//! strand length or alphabet against the layout, the event budget against
//! `n_prog * events_per_progeny`, non-finite or out-of-range chiasmata,
//! non-binary flips — is reported here, where the layout is known, as well as
//! in the R callers.

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
    /// Validate `chiasmata`, `counts` and `flips` against the layout.
    ///
    /// A `counts`/`chiasmata` mismatch silently shifts every subsequent block
    /// and yields plausible-looking output — exactly the failure mode of commit
    /// 4167402 — so it is an error. So is anything a valid R draw cannot
    /// produce: a negative or NA count, a flip other than 0 or 1, a
    /// non-finite chiasma, or a chiasma outside `[0, L]` of its chromosome
    /// (R draws `runif(k, 0, L)`; a value below 0 used to toggle the whole
    /// chromosome and one above `L` was silently a no-op).
    pub fn try_new(
        layout: &GenomeLayout,
        chiasmata: &'a [f64],
        counts: &'a [i32],
        flips: &'a [i32],
    ) -> Result<Self, String> {
        let n_chr = layout.n_chr();
        if n_chr == 0 {
            return Err("n_chr must be positive".to_string());
        }
        if counts.len() != flips.len() {
            return Err(format!(
                "counts has {} entries but flips has {}",
                counts.len(),
                flips.len()
            ));
        }
        if counts.len() % n_chr != 0 {
            return Err(format!(
                "counts length {} is not a multiple of n_chr {}",
                counts.len(),
                n_chr
            ));
        }

        let mut offsets = Vec::with_capacity(counts.len() + 1);
        let mut acc = 0usize;
        offsets.push(0);
        for &c in counts {
            if c < 0 {
                return Err(format!(
                    "crossover counts must be non-negative and non-NA, got {}",
                    if c == R_NA_INT {
                        "NA".to_string()
                    } else {
                        c.to_string()
                    }
                ));
            }
            acc = acc
                .checked_add(c as usize)
                .ok_or_else(|| "crossover counts overflow".to_string())?;
            offsets.push(acc);
        }
        if acc != chiasmata.len() {
            return Err(format!(
                "counts sum to {} but {} chiasmata were supplied",
                acc,
                chiasmata.len()
            ));
        }
        for &f in flips {
            if f != 0 && f != 1 {
                return Err(format!(
                    "flips must be 0 or 1, got {}",
                    if f == R_NA_INT {
                        "NA".to_string()
                    } else {
                        f.to_string()
                    }
                ));
            }
        }
        for slot in 0..counts.len() {
            let len = layout.chr_len(slot % n_chr);
            for &x in &chiasmata[offsets[slot]..offsets[slot + 1]] {
                // NaN fails both comparisons, so this also rejects it.
                if !(0.0..=len).contains(&x) {
                    return Err(format!(
                        "chiasma positions must be finite and within [0, L] of their chromosome (L = {}), got {}",
                        len, x
                    ));
                }
            }
        }

        Ok(EventStream {
            chiasmata,
            flips,
            offsets,
            n_chr,
            n_events: counts.len() / n_chr,
        })
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
    /// `"cross"`, `"selfcross"` or `"dh"`; anything else is an error (it used
    /// to fall back to `Cross` silently).
    pub fn parse(s: &str) -> Result<Design, String> {
        match s {
            "cross" => Ok(Design::Cross),
            "selfcross" => Ok(Design::SelfCross),
            "dh" => Ok(Design::Dh),
            other => Err(format!(
                "design must be \"cross\", \"selfcross\" or \"dh\", got {:?}",
                other
            )),
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

/// The four parental strands, validated against the layout.
pub struct Parents {
    p1_cis: Bits,
    p1_trans: Bits,
    p2_cis: Bits,
    p2_trans: Bits,
}

impl Parents {
    pub fn parse(
        layout: &GenomeLayout,
        p1_cis: &str,
        p1_trans: &str,
        p2_cis: &str,
        p2_trans: &str,
    ) -> Result<Self, String> {
        let n = layout.n_loci();
        Ok(Parents {
            p1_cis: Bits::parse_strand(p1_cis, n, "p1_cis")?,
            p1_trans: Bits::parse_strand(p1_trans, n, "p1_trans")?,
            p2_cis: Bits::parse_strand(p2_cis, n, "p2_cis")?,
            p2_trans: Bits::parse_strand(p2_trans, n, "p2_trans")?,
        })
    }
}

/// Convert R's `n_prog` to a count (NA and negatives are errors).
fn progeny_count(n_prog: i32) -> Result<usize, String> {
    usize::try_from(n_prog).map_err(|_| {
        format!(
            "n_prog must be a non-negative whole number, got {}",
            if n_prog == R_NA_INT {
                "NA".to_string()
            } else {
                n_prog.to_string()
            }
        )
    })
}

/// Produce `n_prog` progeny.
///
/// Events are consumed progeny-major, parent-minor — isqg generates progeny one
/// at a time, drawing the first parent's meiosis then the second's, rather than
/// all of one parent's then all of the other's. For `SelfCross`, pass the same
/// parent twice (isqg selfs against a deep copy, which is still two independent
/// meioses).
///
/// The event stream must hold **exactly** `n_prog * events_per_progeny` events:
/// surplus events used to be ignored silently (so a stale draw looked like a
/// valid run) and a deficit indexed out of bounds.
#[allow(clippy::too_many_arguments)]
pub fn mate_haplotypes(
    layout: &GenomeLayout,
    parents: &Parents,
    events: &EventStream,
    design: Design,
    n_prog: usize,
) -> Result<Vec<(Bits, Bits)>, String> {
    let per = design.events_per_progeny();
    let needed = n_prog
        .checked_mul(per)
        .ok_or_else(|| "n_prog is too large".to_string())?;
    if events.n_events() != needed {
        return Err(format!(
            "{} progeny need {} meiosis events ({} per progeny) but {} were supplied",
            n_prog,
            needed,
            per,
            events.n_events()
        ));
    }
    layout
        .n_loci()
        .checked_mul(n_prog)
        .ok_or_else(|| "n_prog is too large".to_string())?;

    let mut progeny = Vec::with_capacity(n_prog);
    for i in 0..n_prog {
        let pair = if design == Design::Dh {
            let g = recombine(
                &parents.p1_cis,
                &parents.p1_trans,
                &genome_mask(layout, events, i),
            );
            (g.clone(), g)
        } else {
            let a = recombine(
                &parents.p1_cis,
                &parents.p1_trans,
                &genome_mask(layout, events, per * i),
            );
            let b = recombine(
                &parents.p2_cis,
                &parents.p2_trans,
                &genome_mask(layout, events, per * i + 1),
            );
            (a, b)
        };
        progeny.push(pair);
    }

    Ok(progeny)
}

/// As [`mate_haplotypes`], projected onto a loci-major `-1/0/1` genotype buffer
/// (`out[locus * n_prog + prog]`).
///
/// The projection is lossy: it cannot distinguish the two phases of a
/// heterozygote. Callers that will breed from the progeny must keep the
/// haplotypes instead.
pub fn mate_core(
    layout: &GenomeLayout,
    parents: &Parents,
    events: &EventStream,
    design: Design,
    n_prog: usize,
) -> Result<Vec<i32>, String> {
    let progeny = mate_haplotypes(layout, parents, events, design, n_prog)?;
    let mut out = vec![0i32; layout.n_loci() * n_prog];
    for (i, (cis, trans)) in progeny.iter().enumerate() {
        write_genotype(cis, trans, &mut out, i, n_prog);
    }
    Ok(out)
}

/// One mating of a batch (see [`mate_many`]): which parental strands it uses,
/// its design and its number of progeny.
#[derive(Clone, Copy, Debug)]
pub struct Mating {
    pub p1_cis: usize,
    pub p1_trans: usize,
    pub p2_cis: usize,
    pub p2_trans: usize,
    pub design: Design,
    pub n_prog: usize,
}

/// Produce the progeny of a whole batch of matings in one call.
///
/// The batched, integer-I/O form of [`mate_haplotypes`]: every mating draws
/// its meiosis events from one shared [`EventStream`], consecutively and in
/// plan order (a mating of design `d` with `n` progeny consumes
/// `n * d.events_per_progeny()` events, progeny-major, parent-minor, exactly as
/// [`mate_haplotypes`]), so a batch is the same computation as running its
/// matings one after the other with their draws concatenated.
///
/// `strands` holds `n_strands` parental strands of `layout.n_loci()` 0/1
/// entries each, strand after strand, in the CALLER's marker order; `order[j]`
/// is the caller's index of the locus with ascending map rank `j`, so the
/// output is in the caller's order too and no permutation is needed on the R
/// side. Returns `(cis, trans)`, each `n_loci * total_progeny` integers, one
/// progeny per column (progeny of the matings concatenated in order).
pub fn mate_many(
    layout: &GenomeLayout,
    order: &[usize],
    strands: &[i32],
    n_strands: usize,
    matings: &[Mating],
    events: &EventStream,
) -> Result<(Vec<i32>, Vec<i32>), String> {
    let n = layout.n_loci();
    if order.len() != n {
        return Err(format!(
            "order has {} entries but the layout has {} loci",
            order.len(),
            n
        ));
    }
    let mut seen = vec![false; n];
    for &o in order {
        if o >= n || seen[o] {
            return Err("order must be a permutation of the loci".to_string());
        }
        seen[o] = true;
    }
    if n_strands.checked_mul(n) != Some(strands.len()) {
        return Err(format!(
            "strands has {} entries but {} strands of {} loci need {}",
            strands.len(),
            n_strands,
            n,
            n_strands.saturating_mul(n)
        ));
    }
    if let Some(k) = strands.iter().position(|&v| v != 0 && v != 1) {
        return Err(format!(
            "parental strands must hold only 0 and 1 (strand {} locus {} does not)",
            k / n.max(1) + 1,
            k % n.max(1) + 1
        ));
    }
    let mut needed = 0usize;
    let mut total = 0usize;
    for (k, m) in matings.iter().enumerate() {
        for &s in &[m.p1_cis, m.p1_trans, m.p2_cis, m.p2_trans] {
            if s >= n_strands {
                return Err(format!(
                    "mating {} refers to strand {} but only {} strands were supplied",
                    k + 1,
                    s + 1,
                    n_strands
                ));
            }
        }
        needed = m
            .n_prog
            .checked_mul(m.design.events_per_progeny())
            .and_then(|e| needed.checked_add(e))
            .ok_or_else(|| "n_prog is too large".to_string())?;
        total = total
            .checked_add(m.n_prog)
            .ok_or_else(|| "n_prog is too large".to_string())?;
    }
    if events.n_events() != needed {
        return Err(format!(
            "the matings need {} meiosis events but {} were supplied",
            needed,
            events.n_events()
        ));
    }
    let cells = n
        .checked_mul(total)
        .ok_or_else(|| "n_prog is too large".to_string())?;

    // Gamete from a parent's strand pair under an ancestry mask, written into
    // `out` (one progeny's strand, caller's marker order).
    let gamete = |mask: &Bits, cis: usize, trans: usize, out: &mut [i32]| {
        let c = &strands[cis * n..(cis + 1) * n];
        let t = &strands[trans * n..(trans + 1) * n];
        for (j, &o) in order.iter().enumerate() {
            out[o] = if mask.get(j) { c[o] } else { t[o] };
        }
    };

    let mut cis_out = vec![0i32; cells];
    let mut trans_out = vec![0i32; cells];
    let mut ev = 0usize;
    let mut col = 0usize;
    for m in matings {
        for _ in 0..m.n_prog {
            let range = col * n..(col + 1) * n;
            if m.design == Design::Dh {
                let mask = genome_mask(layout, events, ev);
                ev += 1;
                gamete(&mask, m.p1_cis, m.p1_trans, &mut cis_out[range.clone()]);
                trans_out[range.clone()].copy_from_slice(&cis_out[range]);
            } else {
                let mask_a = genome_mask(layout, events, ev);
                let mask_b = genome_mask(layout, events, ev + 1);
                ev += 2;
                gamete(&mask_a, m.p1_cis, m.p1_trans, &mut cis_out[range.clone()]);
                gamete(&mask_b, m.p2_cis, m.p2_trans, &mut trans_out[range]);
            }
            col += 1;
        }
    }
    Ok((cis_out, trans_out))
}

// ---------------------------------------------------------------------------
// extendr boundary — argument marshalling and validation only.
// ---------------------------------------------------------------------------

/// Produce progeny genotypes from pre-drawn meiosis randomness.
///
/// Chiasmata are concatenated and delimited by `counts`, ordered
/// `counts[event * n_chr + chr]`. Events run progeny-major: for "cross" and
/// "selfcross", event `2i` is parent 1's gamete for progeny `i` and `2i+1` is
/// parent 2's; for "dh" there is one event per progeny. The number of events
/// must be exactly `n_prog` times that.
///
/// Returns `list(ok, err)` (feature `result_list`); the R wrapper raises `err`
/// as an R error.
///
/// @param loci_per_chr Integer vector of loci counts, chromosomes ascending.
/// @param positions    Map positions in Morgans, concatenated, non-decreasing within chromosome.
/// @param p1_cis       Parent 1 cis strand as a '0'/'1' string of length n_loci, ascending.
/// @param p1_trans     Parent 1 trans strand.
/// @param p2_cis       Parent 2 cis strand (same as p1 for selfcross/dh).
/// @param p2_trans     Parent 2 trans strand.
/// @param chiasmata    Concatenated crossover positions, each within [0, L] of its chromosome.
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
) -> std::result::Result<Vec<i32>, String> {
    let layout = GenomeLayout::try_new(loci_per_chr, positions)?;
    let events = EventStream::try_new(&layout, chiasmata, counts, flips)?;
    let parents = Parents::parse(&layout, p1_cis, p1_trans, p2_cis, p2_trans)?;
    mate_core(
        &layout,
        &parents,
        &events,
        Design::parse(design)?,
        progeny_count(n_prog)?,
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
/// @param chiasmata    Concatenated crossover positions, each within [0, L] of its chromosome.
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
) -> std::result::Result<Vec<String>, String> {
    let layout = GenomeLayout::try_new(loci_per_chr, positions)?;
    let events = EventStream::try_new(&layout, chiasmata, counts, flips)?;
    Ok((0..events.n_events())
        .map(|e| genome_mask(&layout, &events, e).to_code('1', '0'))
        .collect())
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
/// @param p1_cis       Parent 1 first strand as a '0'/'1' string of length n_loci, ascending.
/// @param p1_trans     Parent 1 second strand.
/// @param p2_cis       Parent 2 first strand (same as p1 for selfcross/dh).
/// @param p2_trans     Parent 2 second strand.
/// @param chiasmata    Concatenated crossover positions, each within [0, L] of its chromosome.
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
) -> std::result::Result<Vec<String>, String> {
    let layout = GenomeLayout::try_new(loci_per_chr, positions)?;
    let events = EventStream::try_new(&layout, chiasmata, counts, flips)?;
    let parents = Parents::parse(&layout, p1_cis, p1_trans, p2_cis, p2_trans)?;
    let progeny = mate_haplotypes(
        &layout,
        &parents,
        &events,
        Design::parse(design)?,
        progeny_count(n_prog)?,
    )?;
    let mut out = Vec::with_capacity(progeny.len() * 2);
    for (cis, trans) in progeny {
        out.push(cis.to_code('1', '0'));
        out.push(trans.to_code('1', '0'));
    }
    Ok(out)
}

/// Progeny of a whole batch of matings in one call, integer in and out.
///
/// The vectorised, string-free form of `mate_haplotypes_core()` (SPEC-0020 item 2):
/// `strands` is an integer 0/1 vector of `n_strands` parental strands of
/// `length(positions)` loci each (strand after strand), in the caller's marker
/// order; `order` (1-based) maps ascending map rank to the caller's index, so
/// the progeny come back in the caller's order. `mating` holds, per mating, the
/// 1-based strand indices `p1_cis, p1_trans, p2_cis, p2_trans` (4 consecutive
/// entries); `design` and `n_prog` have one entry per mating. The meiosis events
/// are shared and consumed in mating order (`n_prog * (1 for "dh", else 2)` per
/// mating, progeny-major), so the result equals running the matings one by one
/// with their draws concatenated. The kernel draws nothing (DECISION-012).
///
/// @param loci_per_chr Integer vector of loci counts, chromosomes ascending.
/// @param positions    Map positions in Morgans, ascending map order, concatenated.
/// @param order        1-based permutation: caller index of the locus with each map rank.
/// @param strands      Integer 0/1 vector, `n_strands` strands of n_loci entries.
/// @param n_strands    Number of parental strands in `strands`.
/// @param mating       4 * n_matings 1-based strand indices (p1_cis, p1_trans, p2_cis, p2_trans).
/// @param design       Character vector, one of "cross", "selfcross", "dh" per mating.
/// @param n_prog       Progeny per mating.
/// @param chiasmata    Concatenated crossover positions, each within [0, L] of its chromosome.
/// @param counts       Crossovers per (event, chromosome).
/// @param flips        0/1 strand-choice per (event, chromosome).
/// @return `list(cis, trans)`: integer vectors of n_loci * sum(n_prog) 0/1 entries, one progeny per column.
/// @noRd
#[extendr]
#[allow(clippy::too_many_arguments)]
pub fn mate_many_core(
    loci_per_chr: &[i32],
    positions: &[f64],
    order: &[i32],
    strands: &[i32],
    n_strands: i32,
    mating: &[i32],
    design: Vec<String>,
    n_prog: &[i32],
    chiasmata: &[f64],
    counts: &[i32],
    flips: &[i32],
) -> std::result::Result<List, String> {
    let layout = GenomeLayout::try_new(loci_per_chr, positions)?;
    let events = EventStream::try_new(&layout, chiasmata, counts, flips)?;
    let n_matings = n_prog.len();
    if design.len() != n_matings || mating.len() != 4 * n_matings {
        return Err(format!(
            "{} matings need {} designs and {} strand indices, got {} and {}",
            n_matings,
            n_matings,
            4 * n_matings,
            design.len(),
            mating.len()
        ));
    }
    let to_index = |v: i32, what: &str| -> Result<usize, String> {
        if v >= 1 {
            Ok((v - 1) as usize)
        } else {
            Err(format!(
                "{} must be a positive 1-based index, got {}",
                what, v
            ))
        }
    };
    let order = order
        .iter()
        .map(|&v| to_index(v, "order"))
        .collect::<Result<Vec<_>, _>>()?;
    let n_strands = usize::try_from(n_strands)
        .map_err(|_| format!("n_strands must be non-negative, got {}", n_strands))?;
    let mut matings = Vec::with_capacity(n_matings);
    for k in 0..n_matings {
        matings.push(Mating {
            p1_cis: to_index(mating[4 * k], "mating")?,
            p1_trans: to_index(mating[4 * k + 1], "mating")?,
            p2_cis: to_index(mating[4 * k + 2], "mating")?,
            p2_trans: to_index(mating[4 * k + 3], "mating")?,
            design: Design::parse(&design[k])?,
            n_prog: progeny_count(n_prog[k])?,
        });
    }
    let (cis, trans) = mate_many(&layout, &order, strands, n_strands, &matings, &events)?;
    Ok(list!(cis = cis, trans = trans))
}

extendr_module! {
    mod meiosis;
    fn meiosis_core;
    fn mate_many_core;
    fn mate_haplotypes_core;
    fn gamete_masks_core;
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pos() -> Vec<f64> {
        vec![0.0, 0.5, 1.0, 1.5, 2.0]
    }

    fn layout(loci: &[i32], positions: &[f64]) -> GenomeLayout {
        GenomeLayout::try_new(loci, positions).unwrap()
    }

    fn stream<'a>(
        l: &GenomeLayout,
        chiasmata: &'a [f64],
        counts: &'a [i32],
        flips: &'a [i32],
    ) -> EventStream<'a> {
        EventStream::try_new(l, chiasmata, counts, flips).unwrap()
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
        let l = layout(&[2, 2], &[0.0, 1.0, 0.0, 1.0]);
        let chiasmata = [0.1, 0.2, 0.3, 0.4];
        let counts = [1i32, 0, 2, 1];
        let flips = [0i32, 1, 1, 0];
        let ev = stream(&l, &chiasmata, &counts, &flips);

        assert_eq!(ev.n_events(), 2);
        assert_eq!(ev.block(0, 0), &[0.1]);
        assert_eq!(ev.block(0, 1), &[] as &[f64]);
        assert_eq!(ev.block(1, 0), &[0.2, 0.3]);
        assert_eq!(ev.block(1, 1), &[0.4]);
        assert!(!ev.flip(0, 0));
        assert!(ev.flip(0, 1));
        assert!(!ev.flip(1, 1));
    }

    fn err_of(r: Result<EventStream, String>) -> String {
        match r {
            Ok(_) => panic!("expected an error"),
            Err(e) => e,
        }
    }

    #[test]
    fn event_stream_rejects_a_count_chiasmata_mismatch() {
        let l = layout(&[2, 2], &[0.0, 1.0, 0.0, 1.0]);
        let e = err_of(EventStream::try_new(
            &l,
            &[0.1, 0.2],
            &[1i32, 0],
            &[0i32, 0],
        ));
        assert!(e.contains("counts sum to"), "{}", e);
    }

    #[test]
    fn event_stream_rejects_negative_and_na_counts() {
        let l = layout(&[2, 2], &[0.0, 1.0, 0.0, 1.0]);
        let e = err_of(EventStream::try_new(&l, &[0.1], &[-1i32, 2], &[0i32, 0]));
        assert!(e.contains("non-negative"), "{}", e);
        let e = err_of(EventStream::try_new(&l, &[], &[R_NA_INT, 0], &[0i32, 0]));
        assert!(e.contains("non-negative"), "{}", e);
    }

    #[test]
    fn event_stream_rejects_non_binary_flips() {
        let l = layout(&[2], &[0.0, 1.0]);
        for bad in [2i32, -1, R_NA_INT] {
            let e = err_of(EventStream::try_new(&l, &[], &[0i32], &[bad]));
            assert!(e.contains("flips must be 0 or 1"), "{}", e);
        }
    }

    #[test]
    fn event_stream_rejects_ragged_shapes() {
        let l = layout(&[2, 2], &[0.0, 1.0, 0.0, 1.0]);
        assert!(EventStream::try_new(&l, &[], &[0i32, 0, 0], &[0i32, 0, 0]).is_err());
        assert!(EventStream::try_new(&l, &[], &[0i32, 0], &[0i32]).is_err());
    }

    #[test]
    fn event_stream_rejects_non_finite_and_out_of_range_chiasmata() {
        let l = layout(&[3], &[0.0, 0.5, 1.0]);
        for bad in [
            f64::NAN,
            f64::INFINITY,
            f64::NEG_INFINITY,
            -1.0,
            5.0,
            1.0 + 1e-12,
        ] {
            let e = err_of(EventStream::try_new(&l, &[bad], &[1i32], &[0i32]));
            assert!(e.contains("chiasma"), "{}: {}", bad, e);
        }
        // the boundaries themselves are legal: 0 and L = the last position
        assert!(EventStream::try_new(&l, &[0.0, 1.0], &[2i32], &[0i32]).is_ok());
    }

    #[test]
    fn chiasma_range_is_per_chromosome() {
        let l = layout(&[2, 2], &[0.0, 0.2, 0.0, 2.0]);
        // 1.0 is inside chr 2 (L = 2) but outside chr 1 (L = 0.2)
        assert!(EventStream::try_new(&l, &[1.0], &[0i32, 1], &[0i32, 0]).is_ok());
        assert!(EventStream::try_new(&l, &[1.0], &[1i32, 0], &[0i32, 0]).is_err());
    }

    #[test]
    fn genome_mask_concatenates_chromosomes_ascending() {
        let l = layout(&[3, 2], &[0.0, 0.5, 1.0, 0.0, 1.0]);
        // chr0: no crossover, no flip -> 000. chr1: no crossover, flip -> 11.
        let ev = stream(&l, &[], &[0i32, 0], &[0i32, 1]);
        assert_eq!(genome_mask(&l, &ev, 0).to_code('1', '0'), "00011");
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
        let l = layout(&[4], &[0.0, 0.3, 0.6, 0.9]);
        let ev = stream(&l, &[0.45], &[1i32], &[0i32]);
        let p = Parents::parse(&l, "1100", "0011", "1100", "0011").unwrap();
        let out = mate_core(&l, &p, &ev, Design::Dh, 1).unwrap();
        assert!(out.iter().all(|&g| g != 0), "DH produced a heterozygote");
    }

    #[test]
    fn cross_takes_cis_from_the_first_parent() {
        let l = layout(&[2], &[0.0, 1.0]);
        // No crossovers, no flips: each gamete is the parent's trans strand
        // (mask all zero), so the phase is fully determined by parent order.
        let ev = stream(&l, &[], &[0i32, 0], &[0i32, 0]);
        let p = Parents::parse(&l, "11", "11", "00", "00").unwrap();
        let ab = mate_core(&l, &p, &ev, Design::Cross, 1).unwrap();
        // p1 contributes A, p2 contributes a => heterozygous everywhere.
        assert_eq!(ab, vec![0, 0]);
    }

    #[test]
    fn events_are_consumed_progeny_major() {
        // Two progeny, one chromosome, no crossovers. Flips per event:
        // [p1 of prog0, p2 of prog0, p1 of prog1, p2 of prog1] = [0,0,1,1].
        // Progeny 0 takes both trans strands; progeny 1 takes both cis.
        let l = layout(&[2], &[0.0, 1.0]);
        let ev = stream(&l, &[], &[0i32; 4], &[0i32, 0, 1, 1]);
        let p = Parents::parse(&l, "11", "00", "11", "00").unwrap();
        let out = mate_core(&l, &p, &ev, Design::Cross, 2).unwrap();
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
        let l = layout(&[8], &[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]);
        let p = Parents::parse(&l, "10101010", "01010101", "10101010", "01010101").unwrap();

        // One crossover at 0.35 => breaks at index 4; no terminal flip.
        let ev = stream(&l, &[0.35], &[1i32], &[0i32]);
        let prog = mate_haplotypes(&l, &p, &ev, Design::Dh, 1).unwrap();

        let (a, b) = &prog[0];
        assert_eq!(a, b, "a doubled haploid must have identical strands");
        // Loci 0-3 from trans, loci 4-7 from cis.
        assert_eq!(a.to_code('1', '0'), "01011010");
    }

    #[test]
    fn mate_core_agrees_with_mate_haplotypes() {
        let l = layout(&[4], &[0.0, 0.3, 0.6, 0.9]);
        let ev = stream(&l, &[0.45, 0.2], &[1i32, 1], &[0i32, 1]);
        let p = Parents::parse(&l, "1100", "0110", "1100", "0110").unwrap();

        let hap = mate_haplotypes(&l, &p, &ev, Design::Cross, 1).unwrap();
        let geno = mate_core(&l, &p, &ev, Design::Cross, 1).unwrap();
        let mut expected = vec![0i32; 4];
        write_genotype(&hap[0].0, &hap[0].1, &mut expected, 0, 1);
        assert_eq!(geno, expected);
    }

    #[test]
    fn the_event_budget_must_match_n_prog_exactly() {
        let l = layout(&[3], &[0.0, 0.1, 0.2]);
        let p = Parents::parse(&l, "111", "000", "111", "000").unwrap();
        // 4 events for 1 cross progeny (needs 2): surplus used to be ignored
        let ev = stream(&l, &[], &[0i32; 4], &[0i32; 4]);
        assert!(mate_haplotypes(&l, &p, &ev, Design::Cross, 1).is_err());
        // 1 event for 1 cross progeny: deficit used to index out of bounds
        let ev = stream(&l, &[], &[0i32], &[0i32]);
        assert!(mate_haplotypes(&l, &p, &ev, Design::Cross, 1).is_err());
        // zero events, one progeny
        let ev = stream(&l, &[], &[], &[]);
        assert!(mate_haplotypes(&l, &p, &ev, Design::Dh, 1).is_err());
        // exact budgets pass: dh needs 1 per progeny, self/cross 2
        let ev = stream(&l, &[], &[0i32; 3], &[0i32; 3]);
        assert!(mate_haplotypes(&l, &p, &ev, Design::Dh, 3).is_ok());
        let ev = stream(&l, &[], &[0i32; 4], &[0i32; 4]);
        assert!(mate_haplotypes(&l, &p, &ev, Design::SelfCross, 2).is_ok());
    }

    #[test]
    fn progeny_count_rejects_na_and_negative() {
        assert_eq!(progeny_count(0), Ok(0));
        assert_eq!(progeny_count(7), Ok(7));
        assert!(progeny_count(-1).is_err());
        assert!(progeny_count(R_NA_INT).is_err());
    }

    #[test]
    fn parents_must_match_the_layout() {
        let l = layout(&[3], &[0.0, 0.1, 0.2]);
        assert!(Parents::parse(&l, "11", "000", "111", "000").is_err());
        assert!(Parents::parse(&l, "1111", "000", "111", "000").is_err());
        assert!(Parents::parse(&l, "121", "000", "111", "000").is_err());
        assert!(Parents::parse(&l, "111", "000", "111", "000").is_ok());
    }

    #[test]
    fn design_parsing_and_event_budget() {
        assert_eq!(Design::parse("dh"), Ok(Design::Dh));
        assert_eq!(Design::parse("selfcross"), Ok(Design::SelfCross));
        assert_eq!(Design::parse("cross"), Ok(Design::Cross));
        assert!(Design::parse("anything else").is_err());
        assert_eq!(Design::Dh.events_per_progeny(), 1);
        assert_eq!(Design::Cross.events_per_progeny(), 2);
        assert_eq!(Design::SelfCross.events_per_progeny(), 2);
    }

    #[test]
    fn entry_points_return_errors_not_panics() {
        // the six reconciliation abort probes, at the Rust level
        let r = mate_haplotypes_core(&[1], &[0.0], "1", "1", "1", "1", &[], &[], &[], "cross", 1);
        assert!(r.is_err());
        let r = gamete_masks_core(&[3], &[0.0, 0.5, 1.0], &[0.1, 0.2], &[1], &[0]);
        assert!(r.is_err());
        let r = mate_haplotypes_core(
            &[3],
            &[0.0, 0.1, 0.2],
            "111",
            "000",
            "111",
            "000",
            &[],
            &[0, 0],
            &[0, 0],
            "cross",
            2,
        );
        assert!(r.is_err());
        let r = mate_haplotypes_core(
            &[3],
            &[0.0, 0.1, 0.2],
            "11111111",
            "00000000",
            "11111111",
            "00000000",
            &[],
            &[0, 0],
            &[0, 0],
            "cross",
            1,
        );
        assert!(r.is_err());
        let r = mate_haplotypes_core(
            &[3],
            &[0.0, 0.1, 0.2],
            "111",
            "000",
            "111",
            "000",
            &[],
            &[0, 0],
            &[0, 0],
            "cross",
            -1,
        );
        assert!(r.is_err());
        let r = meiosis_core(
            &[3],
            &[0.0, 0.1, 0.2],
            "111",
            "000",
            "111",
            "000",
            &[],
            &[0, 0],
            &[0, 0],
            "bogus",
            1,
        );
        assert!(r.is_err());
    }

    // ---- mate_many: the batched integer-I/O path ---------------------------

    fn bits_to_i32(b: &Bits) -> Vec<i32> {
        (0..b.len()).map(|j| b.get(j) as i32).collect()
    }

    /// Two chromosomes, 7 loci, an UNSORTED caller order, two matings (a cross
    /// of strands (0,1)x(2,3) with 2 progeny, then a dh of strands (4,5) with
    /// 2 progeny): a batch must equal the matings run one at a time.
    #[test]
    fn mate_many_equals_sequential_mate_haplotypes() {
        let loci = [4i32, 3];
        let positions = [0.0, 0.2, 0.5, 0.9, 0.1, 0.4, 1.0];
        let l = layout(&loci, &positions);
        // the caller stores the loci in a scrambled order
        let order = [3usize, 0, 5, 1, 6, 2, 4]; // order[rank] = caller index
        let strand_codes = [
            "1010110", "0110011", "1111000", "0001111", "1001001", "0110110",
        ];
        let n = 7;
        // caller-order strands: caller index order[rank] holds the code char at rank
        let mut strands = vec![0i32; 6 * n];
        for (s, code) in strand_codes.iter().enumerate() {
            for (rank, ch) in code.chars().enumerate() {
                strands[s * n + order[rank]] = (ch == '1') as i32;
            }
        }
        // events: mating 1 = 2 progeny x 2 events, mating 2 = 2 progeny x 1 event
        let counts = [
            1i32, 0, 2, 1, 0, 0, 1, 1, // mating 1, events 0..3 (2 chr each)
            0, 2, 1, 0, // mating 2, events 0..1
        ];
        let chiasmata = [0.3, 0.2, 0.6, 0.5, 0.35, 0.8, 0.2, 0.7, 0.9];
        let flips = [0i32, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0];
        let total: i32 = counts.iter().sum();
        assert_eq!(total as usize, chiasmata.len());
        let ev = stream(&l, &chiasmata, &counts, &flips);
        let matings = [
            Mating {
                p1_cis: 0,
                p1_trans: 1,
                p2_cis: 2,
                p2_trans: 3,
                design: Design::Cross,
                n_prog: 2,
            },
            Mating {
                p1_cis: 4,
                p1_trans: 5,
                p2_cis: 4,
                p2_trans: 5,
                design: Design::Dh,
                n_prog: 2,
            },
        ];
        let (cis, trans) = mate_many(&l, &order, &strands, 6, &matings, &ev).unwrap();
        assert_eq!(cis.len(), n * 4);

        // reference: the same events through mate_haplotypes, mating by mating
        let ev1 = stream(&l, &chiasmata[..6], &counts[..8], &flips[..8]);
        let p1 = Parents::parse(&l, "1010110", "0110011", "1111000", "0001111").unwrap();
        let r1 = mate_haplotypes(&l, &p1, &ev1, Design::Cross, 2).unwrap();
        let ev2 = stream(&l, &chiasmata[6..], &counts[8..], &flips[8..]);
        let p2 = Parents::parse(&l, "1001001", "0110110", "1001001", "0110110").unwrap();
        let r2 = mate_haplotypes(&l, &p2, &ev2, Design::Dh, 2).unwrap();
        let expect: Vec<(Bits, Bits)> = r1.into_iter().chain(r2).collect();
        for (col, (c, t)) in expect.iter().enumerate() {
            let (c, t) = (bits_to_i32(c), bits_to_i32(t));
            for rank in 0..n {
                assert_eq!(cis[col * n + order[rank]], c[rank], "cis col {col}");
                assert_eq!(trans[col * n + order[rank]], t[rank], "trans col {col}");
            }
        }
    }

    #[test]
    fn mate_many_rejects_inconsistent_inputs() {
        let l = layout(&[3], &[0.0, 0.5, 1.0]);
        let ev = stream(&l, &[], &[0i32, 0], &[0i32, 0]);
        let m = |design, n_prog| Mating {
            p1_cis: 0,
            p1_trans: 1,
            p2_cis: 0,
            p2_trans: 1,
            design,
            n_prog,
        };
        let strands = [1i32, 0, 1, 0, 1, 0];
        let id = [0usize, 1, 2];
        assert!(mate_many(&l, &id, &strands, 2, &[m(Design::Cross, 1)], &ev).is_ok());
        // event budget
        assert!(mate_many(&l, &id, &strands, 2, &[m(Design::Cross, 2)], &ev).is_err());
        assert!(mate_many(&l, &id, &strands, 2, &[m(Design::Dh, 1)], &ev).is_err());
        // order is not a permutation / wrong length
        assert!(mate_many(&l, &[0, 0, 2], &strands, 2, &[m(Design::Cross, 1)], &ev).is_err());
        assert!(mate_many(&l, &[0, 1], &strands, 2, &[m(Design::Cross, 1)], &ev).is_err());
        // strand bookkeeping and alphabet
        assert!(mate_many(&l, &id, &strands, 3, &[m(Design::Cross, 1)], &ev).is_err());
        assert!(mate_many(&l, &id, &[1, 2, 1, 0, 1, 0], 2, &[m(Design::Cross, 1)], &ev).is_err());
        let mut bad = m(Design::Cross, 1);
        bad.p2_trans = 2;
        assert!(mate_many(&l, &id, &strands, 2, &[bad], &ev).is_err());
    }

    #[test]
    fn mate_many_core_returns_errors_not_panics() {
        let r = mate_many_core(
            &[3],
            &[0.0, 0.1, 0.2],
            &[1, 2, 3],
            &[1, 0, 1, 0, 1, 0],
            2,
            &[1, 2, 1, 2],
            vec!["bogus".to_string()],
            &[1],
            &[],
            &[0, 0],
            &[0, 0],
        );
        assert!(r.is_err());
        let r = mate_many_core(
            &[3],
            &[0.0, 0.1, 0.2],
            &[0, 2, 3],
            &[1, 0, 1, 0, 1, 0],
            2,
            &[1, 2, 1, 2],
            vec!["cross".to_string()],
            &[1],
            &[],
            &[0, 0],
            &[0, 0],
        );
        assert!(r.is_err());
    }
}
