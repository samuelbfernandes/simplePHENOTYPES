//! Coalescent founder generator (DECISION-050, docs/SPEC-coalescent.md).
//!
//! A sequential Markov coalescent (SMC', Marjoram & Wall 2006) walked along one
//! chromosome scaled to [0, 1). Time is in units of 4 N0 generations (the ms /
//! MaCS convention): a pair of lineages coalesces at rate 2 / lambda(u), where
//! lambda(u) N0 is the population size at time u (piecewise constant, ms `-eN`).
//! With a marginal tree of total branch length L, mutations (infinite sites) and
//! recombination points arrive along the sequence at rates theta L and rho L, so
//! the next event is an exponential race between the two.
//!
//! At a recombination point a uniform point (branch b, time t) is chosen on the
//! tree, the lineage above t floats and re-coalesces at rate 2 k(u) / lambda(u)
//! with the k(u) lineages of the current tree -- its own branch included (the
//! SMC' rule; joining it leaves the tree unchanged). Otherwise b is pruned at t
//! and regrafted (a subtree prune and regraft that reuses b's old parent node).
//!
//! This is the only place in the package where Rust draws random numbers: a
//! xoshiro256++ generator seeded through SplitMix64 from one seed per chromosome
//! that R draws from its own stream (DECISION-050), so a run is reproducible.

use extendr_api::prelude::*;

// ---------------------------------------------------------------------------
// PRNG: xoshiro256++ (Blackman & Vigna), seeded by SplitMix64
// ---------------------------------------------------------------------------

pub struct Rng {
    s: [u64; 4],
}

fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

impl Rng {
    pub fn new(seed: u64) -> Self {
        let mut st = seed;
        let s = [
            splitmix64(&mut st),
            splitmix64(&mut st),
            splitmix64(&mut st),
            splitmix64(&mut st),
        ];
        Rng { s }
    }

    pub fn next_u64(&mut self) -> u64 {
        let result = self.s[0]
            .wrapping_add(self.s[3])
            .rotate_left(23)
            .wrapping_add(self.s[0]);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        result
    }

    /// Uniform on [0, 1) with 53 random bits.
    pub fn unif(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / 9_007_199_254_740_992.0)
    }

    /// Standard exponential.
    pub fn exp1(&mut self) -> f64 {
        -(1.0 - self.unif()).ln()
    }

    /// Uniform integer in 0..n (n > 0), by rejection (no modulo bias).
    pub fn below(&mut self, n: usize) -> usize {
        let n64 = n as u64;
        let zone = u64::MAX - (u64::MAX % n64);
        loop {
            let x = self.next_u64();
            if x < zone {
                return (x % n64) as usize;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Demography: piecewise-constant relative size lambda(u)
// ---------------------------------------------------------------------------

/// Range of a relative population size and of a change time (4 N0 units).
/// Outside it the coalescence rates or waiting times leave double precision
/// (a rate of Inf collapses a waiting time to 0, a time of Inf breaks the walk).
pub const MIN_SIZE: f64 = 1e-9;
pub const MAX_SIZE: f64 = 1e9;
pub const MAX_TIME: f64 = 1e12;

pub struct Demography {
    /// Change times, ascending, times[0] = 0.
    times: Vec<f64>,
    /// Relative sizes, sizes[i] on [times[i], times[i + 1]).
    sizes: Vec<f64>,
}

impl Demography {
    pub fn new(change_times: &[f64], change_sizes: &[f64]) -> std::result::Result<Self, String> {
        if change_times.len() != change_sizes.len() {
            return Err("history times and sizes differ in length".into());
        }
        let mut times = vec![0.0];
        let mut sizes = vec![1.0];
        for (i, (&t, &x)) in change_times.iter().zip(change_sizes).enumerate() {
            if !(t > 0.0 && t <= MAX_TIME && (MIN_SIZE..=MAX_SIZE).contains(&x)) {
                return Err(format!(
                    "history times must be in (0, {MAX_TIME:e}] and sizes in \
                     [{MIN_SIZE:e}, {MAX_SIZE:e}] (relative to N0, time in 4 N0 generations)"
                ));
            }
            if i > 0 && t <= change_times[i - 1] {
                return Err("history times must be strictly increasing".into());
            }
            times.push(t);
            sizes.push(x);
        }
        Ok(Demography { times, sizes })
    }

    /// Epoch index of time u.
    fn epoch(&self, u: f64) -> usize {
        match self.times.partition_point(|&t| t <= u) {
            0 => 0,
            i => i - 1,
        }
    }

    fn next_change(&self, epoch: usize) -> f64 {
        self.times.get(epoch + 1).copied().unwrap_or(f64::INFINITY)
    }
}

// ---------------------------------------------------------------------------
// Fenwick tree over branch lengths (sampling a point uniformly on the tree)
// ---------------------------------------------------------------------------

struct Fenwick {
    tree: Vec<f64>,
    vals: Vec<f64>,
}

impl Fenwick {
    fn new(n: usize) -> Self {
        Fenwick {
            tree: vec![0.0; n + 1],
            vals: vec![0.0; n],
        }
    }

    fn set(&mut self, i: usize, v: f64) {
        let delta = v - self.vals[i];
        self.vals[i] = v;
        let mut j = i + 1;
        while j < self.tree.len() {
            self.tree[j] += delta;
            j += j & j.wrapping_neg();
        }
    }

    /// Sum of all values (prefix sum over the whole tree; `rebuild` limits drift).
    fn total(&self) -> f64 {
        let mut j = self.vals.len();
        let mut acc = 0.0;
        while j > 0 {
            acc += self.tree[j];
            j &= j - 1;
        }
        acc
    }

    /// Index whose cumulative interval contains `target` (0 <= target < total).
    fn find(&self, mut target: f64) -> usize {
        let n = self.vals.len();
        let mut pos = 0usize;
        let mut step = n.next_power_of_two();
        while step > 0 {
            let nxt = pos + step;
            if nxt <= n && self.tree[nxt] <= target {
                pos = nxt;
                target -= self.tree[nxt];
            }
            step >>= 1;
        }
        // guard against floating-point drift: return a node with positive length
        let mut i = pos.min(n - 1);
        if self.vals[i] <= 0.0 {
            i = (0..n).rev().find(|&j| self.vals[j] > 0.0).unwrap_or(i);
        }
        i
    }

    fn rebuild(&mut self) {
        let vals = std::mem::take(&mut self.vals);
        let n = vals.len();
        *self = Fenwick::new(n);
        for (i, v) in vals.into_iter().enumerate() {
            self.set(i, v);
        }
    }
}

// ---------------------------------------------------------------------------
// The marginal tree
// ---------------------------------------------------------------------------

const NONE: usize = usize::MAX;

/// Error for coalescence times that overflow f64 (an extreme `history`).
const TIME_OVERFLOW: &str = "coalescence rates or times left double precision \
                             (non-finite): the history is too extreme";

/// Deme label of a node whose time is at or above the join (one population).
const MERGED: u8 = 2;

/// Two isolated subpopulations that merge into one at `join` (ms `-I 2 n1 n2`
/// with no migration, then `-ej join 2 1`): haplotypes `0..n_first` sample
/// deme 0 and the rest deme 1. Both demes take the population size of the
/// history (ms `-eN` applies to every subpopulation).
#[derive(Clone, Copy)]
pub struct Split {
    pub n_first: usize,
    pub join: f64,
}

pub struct Tree {
    n: usize,
    parent: Vec<usize>,
    children: Vec<[usize; 2]>,
    time: Vec<f64>,
    root: usize,
    /// (time, node) of the internal nodes, ascending in time.
    internal: Vec<(f64, usize)>,
    fen: Fenwick,
    events: u64,
    split: Option<Split>,
    /// Deme of each node: 0 / 1 below the join (all its leaves sample that
    /// deme, as there is no migration), MERGED at or above it.
    deme: Vec<u8>,
    /// Times of the internal nodes below the join, per deme, ascending.
    deme_times: [Vec<f64>; 2],
    /// Number of sampled haplotypes per deme.
    deme_n: [usize; 2],
}

impl Tree {
    /// Kingman coalescent tree of n >= 2 samples under the demography, with an
    /// optional split into two isolated demes below `split.join`.
    pub fn kingman(
        n: usize,
        dem: &Demography,
        split: Option<Split>,
        rng: &mut Rng,
    ) -> std::result::Result<Self, String> {
        if let Some(sp) = split {
            if sp.n_first == 0 || sp.n_first >= n {
                return Err("a split needs at least one haplotype in each subpopulation".into());
            }
            if !(sp.join > 0.0 && sp.join <= MAX_TIME) {
                return Err(format!("the split time must be in (0, {MAX_TIME:e}]"));
            }
        }
        let m = 2 * n - 1;
        let mut deme = vec![MERGED; m];
        let mut deme_n = [0, 0];
        if let Some(sp) = split {
            for (i, d) in deme.iter_mut().enumerate().take(n) {
                *d = u8::from(i >= sp.n_first);
            }
            deme_n = [sp.n_first, n - sp.n_first];
        }
        let mut tree = Tree {
            n,
            parent: vec![NONE; m],
            children: vec![[NONE, NONE]; m],
            time: vec![0.0; m],
            root: NONE,
            internal: Vec::with_capacity(n - 1),
            fen: Fenwick::new(m),
            events: 0,
            split,
            deme,
            deme_times: [Vec::new(), Vec::new()],
            deme_n,
        };
        let mut u = 0.0;
        let mut ep = dem.epoch(0.0);
        let mut next_node = n;
        let mut join_at = |tree: &mut Tree, a: usize, b: usize, u: f64, label: u8| -> usize {
            let p = next_node;
            next_node += 1;
            tree.time[p] = u;
            tree.children[p] = [a, b];
            tree.parent[a] = p;
            tree.parent[b] = p;
            tree.internal.push((u, p));
            tree.deme[p] = label;
            if label != MERGED {
                tree.deme_times[label as usize].push(u);
            }
            p
        };
        let mut active: Vec<usize> = match split {
            None => (0..n).collect(),
            Some(sp) => {
                // isolated demes: each coalesces at k_d (k_d - 1) / lambda
                let mut act: [Vec<usize>; 2] =
                    [(0..sp.n_first).collect(), (sp.n_first..n).collect()];
                loop {
                    let k0 = act[0].len() as f64;
                    let k1 = act[1].len() as f64;
                    let r0 = k0 * (k0 - 1.0) / dem.sizes[ep];
                    let r1 = k1 * (k1 - 1.0) / dem.sizes[ep];
                    let rate = r0 + r1;
                    let change = dem.next_change(ep);
                    let bound = change.min(sp.join);
                    let dt = if rate > 0.0 {
                        if !rate.is_finite() {
                            return Err(TIME_OVERFLOW.into());
                        }
                        rng.exp1() / rate
                    } else {
                        f64::INFINITY
                    };
                    if u + dt >= bound {
                        u = bound;
                        if change <= sp.join {
                            ep += 1;
                        }
                        if sp.join <= change {
                            break;
                        }
                        continue;
                    }
                    u += dt;
                    let d = usize::from(rng.unif() * rate >= r0);
                    let i = rng.below(act[d].len());
                    let a = act[d].swap_remove(i);
                    let j = rng.below(act[d].len());
                    let b = act[d].swap_remove(j);
                    let p = join_at(&mut tree, a, b, u, d as u8);
                    act[d].push(p);
                }
                let [a0, a1] = act;
                a0.into_iter().chain(a1).collect()
            }
        };
        while active.len() > 1 {
            let k = active.len() as f64;
            let rate = k * (k - 1.0) / dem.sizes[ep];
            if !(rate.is_finite() && rate > 0.0) {
                return Err(TIME_OVERFLOW.into());
            }
            let dt = rng.exp1() / rate;
            if !(u + dt).is_finite() {
                return Err(TIME_OVERFLOW.into());
            }
            let change = dem.next_change(ep);
            if u + dt >= change {
                u = change;
                ep += 1;
                continue;
            }
            u += dt;
            let i = rng.below(active.len());
            let a = active.swap_remove(i);
            let j = rng.below(active.len());
            let b = active.swap_remove(j);
            let p = join_at(&mut tree, a, b, u, MERGED);
            active.push(p);
        }
        tree.root = active[0];
        for c in 0..m {
            tree.refresh_len(c);
        }
        if !tree.total_length().is_finite() {
            return Err(TIME_OVERFLOW.into());
        }
        Ok(tree)
    }

    fn branch_len(&self, c: usize) -> f64 {
        match self.parent[c] {
            NONE => 0.0,
            p => self.time[p] - self.time[c],
        }
    }

    fn refresh_len(&mut self, c: usize) {
        let l = self.branch_len(c);
        self.fen.set(c, l);
    }

    pub fn total_length(&self) -> f64 {
        self.fen.total()
    }

    pub fn tmrca(&self) -> f64 {
        self.time[self.root]
    }

    /// Number of lineages of the tree at time u > 0.
    #[cfg(test)]
    fn lineages_at(&self, u: f64) -> usize {
        let above = self.internal.len() - self.internal.partition_point(|&(t, _)| t <= u);
        1 + above
    }

    /// A uniform point on the tree: (branch child node, time).
    fn sample_point(&self, rng: &mut Rng) -> (usize, f64) {
        let total = self.total_length();
        let c = self.fen.find(rng.unif() * total);
        let lo = self.time[c];
        let hi = self.time[self.parent[c]];
        (c, lo + rng.unif() * (hi - lo))
    }

    /// Re-coalescence time of a lineage floating from time t on branch b
    /// (SMC'). Below a split's join it can only meet lineages of its own deme.
    fn float_time(
        &self,
        t: f64,
        b: usize,
        dem: &Demography,
        rng: &mut Rng,
    ) -> std::result::Result<f64, String> {
        let mut u = t;
        let mut ep = dem.epoch(u);
        if let Some(sp) = self.split {
            if u < sp.join {
                let d = self.deme[b] as usize;
                let list = &self.deme_times[d];
                let mut idx = list.partition_point(|&x| x <= u);
                loop {
                    // lineages of deme d at u: its samples minus its coalescences so far
                    let k = (self.deme_n[d] - idx) as f64;
                    let rate = 2.0 * k / dem.sizes[ep];
                    if !(rate.is_finite() && rate > 0.0) {
                        return Err(TIME_OVERFLOW.into());
                    }
                    let next_node = list.get(idx).copied().unwrap_or(f64::INFINITY);
                    let change = dem.next_change(ep);
                    let bound = next_node.min(change).min(sp.join);
                    let dt = rng.exp1() / rate;
                    if u + dt < bound {
                        return Ok(u + dt);
                    }
                    u = bound;
                    if change <= bound {
                        ep += 1;
                    }
                    if next_node <= bound {
                        idx += 1;
                    }
                    if sp.join <= bound {
                        break;
                    }
                }
            }
        }
        let mut idx = self.internal.partition_point(|&(x, _)| x <= u);
        loop {
            let k = (1 + self.internal.len() - idx) as f64;
            let rate = 2.0 * k / dem.sizes[ep];
            if !(rate.is_finite() && rate > 0.0) {
                return Err(TIME_OVERFLOW.into());
            }
            let next_node = self.internal.get(idx).map_or(f64::INFINITY, |&(t, _)| t);
            let change = dem.next_change(ep);
            let bound = next_node.min(change);
            let dt = rng.exp1() / rate;
            if !(u + dt).is_finite() {
                // only an overflow can leave the last epoch above the root
                return Err(TIME_OVERFLOW.into());
            }
            if u + dt < bound {
                return Ok(u + dt);
            }
            u = bound;
            if change <= next_node {
                ep += 1;
            }
            if next_node <= change {
                idx += 1;
            }
        }
    }

    /// The lineages present at time tau: child nodes whose branch spans tau
    /// (the root counts above its own time). O(n); the tests check the O(1)
    /// sampler against it.
    #[cfg(test)]
    fn lineages_crossing(&self, tau: f64) -> Vec<usize> {
        (0..self.parent.len())
            .filter(|&c| {
                let below = self.time[c] <= tau;
                match self.parent[c] {
                    NONE => below,
                    p => below && tau < self.time[p],
                }
            })
            .collect()
    }

    fn remove_internal(&mut self, t: f64, node: usize) {
        let mut i = self.internal.partition_point(|&(x, _)| x < t);
        while i < self.internal.len() && self.internal[i].0 == t {
            if self.internal[i].1 == node {
                self.internal.remove(i);
                return;
            }
            i += 1;
        }
    }

    fn insert_internal(&mut self, t: f64, node: usize) {
        let i = self.internal.partition_point(|&(x, _)| x < t);
        self.internal.insert(i, (t, node));
    }

    /// A uniform lineage among those present at time tau: an O(log n) search
    /// plus O(1) expected rejection draws (keeping `internal` sorted still
    /// shifts O(n) entries per regraft, a fast memmove). Candidates are the root and both children of every internal node
    /// older than tau (2m + 1 of them for m such nodes; k = m + 1 lineages);
    /// a candidate is a lineage at tau iff its branch starts at or below tau,
    /// so rejection keeps a uniform draw with acceptance k / (2k - 1) >= 1/2.
    fn sample_lineage_at(&self, tau: f64, b: usize, rng: &mut Rng) -> usize {
        if let Some(sp) = self.split {
            if tau < sp.join {
                // below the join: a uniform lineage of b's deme (O(n) scan; only
                // split runs re-coalesce below the join)
                let d = self.deme[b];
                let cands: Vec<usize> = (0..self.parent.len())
                    .filter(|&c| {
                        self.deme[c] == d
                            && self.time[c] <= tau
                            && self.parent[c] != NONE
                            && tau < self.time[self.parent[c]]
                    })
                    .collect();
                return cands[rng.below(cands.len())];
            }
        }
        let idx = self.internal.partition_point(|&(x, _)| x <= tau);
        let m = self.internal.len() - idx;
        loop {
            let r = rng.below(2 * m + 1);
            let c = if r == 2 * m {
                self.root
            } else {
                self.children[self.internal[idx + r / 2].1][r % 2]
            };
            if self.time[c] <= tau {
                return c;
            }
        }
    }

    /// One SMC' recombination event.
    pub fn recombine(
        &mut self,
        dem: &Demography,
        rng: &mut Rng,
    ) -> std::result::Result<(), String> {
        self.events += 1;
        let (b, t) = self.sample_point(rng);
        let tau = self.float_time(t, b, dem, rng)?;
        let mut c = self.sample_lineage_at(tau, b, rng);
        let p = self.parent[b];
        if c == b {
            return Ok(()); // rejoined its own branch: the tree is unchanged
        }
        // prune: remove p, connecting b's sibling s to p's parent gp
        let s = if self.children[p][0] == b {
            self.children[p][1]
        } else {
            self.children[p][0]
        };
        let gp = self.parent[p];
        self.parent[s] = gp;
        if gp == NONE {
            self.root = s;
        } else {
            let slot = if self.children[gp][0] == p { 0 } else { 1 };
            self.children[gp][slot] = s;
        }
        if c == p {
            c = s; // p's lineage above its time is s's lineage once p is gone
        }
        self.remove_internal(self.time[p], p);
        if self.deme[p] != MERGED {
            let d = self.deme[p] as usize;
            let old = self.time[p];
            let i = self.deme_times[d].partition_point(|&x| x < old);
            self.deme_times[d].remove(i);
        }
        // regraft: p becomes the parent of b and c at tau
        let cp = self.parent[c];
        self.parent[p] = cp;
        if cp == NONE {
            self.root = p;
        } else {
            let slot = if self.children[cp][0] == c { 0 } else { 1 };
            self.children[cp][slot] = p;
        }
        self.children[p] = [b, c];
        self.parent[b] = p;
        self.parent[c] = p;
        self.time[p] = tau;
        self.insert_internal(tau, p);
        self.deme[p] = match self.split {
            Some(sp) if tau < sp.join => self.deme[b],
            _ => MERGED,
        };
        if self.deme[p] != MERGED {
            let d = self.deme[p] as usize;
            let i = self.deme_times[d].partition_point(|&x| x < tau);
            self.deme_times[d].insert(i, tau);
        }
        for node in [b, c, s, p] {
            self.refresh_len(node);
        }
        if self.events % 65_536 == 0 {
            self.fen.rebuild();
        }
        if !self.total_length().is_finite() {
            return Err(TIME_OVERFLOW.into());
        }
        Ok(())
    }

    /// Leaves below node c, as a bit set of n bits.
    fn leaf_bits(&self, c: usize) -> Vec<u64> {
        let mut bits = vec![0u64; (self.n + 63) / 64];
        let mut stack = vec![c];
        while let Some(x) = stack.pop() {
            if x < self.n {
                bits[x / 64] |= 1u64 << (x % 64);
            } else {
                stack.extend_from_slice(&self.children[x]);
            }
        }
        bits
    }

    /// A mutation on a uniform point of the tree: its carrier bit set.
    fn mutate(&self, rng: &mut Rng) -> Vec<u64> {
        let total = self.total_length();
        let c = self.fen.find(rng.unif() * total);
        self.leaf_bits(c)
    }
}

// ---------------------------------------------------------------------------
// Walking one chromosome
// ---------------------------------------------------------------------------

/// Upper bound on mutation + recombination events per chromosome; a run that
/// would exceed it stops with an error instead of running (practically) forever.
pub const MAX_EVENTS: u64 = 2_000_000_000;

pub struct Sites {
    pub positions: Vec<f64>,
    pub carriers: Vec<Vec<u64>>,
    pub n_total: u64,
    pub first_tmrca: f64,
    pub first_length: f64,
    /// T_MRCA of the tree at the end of the chromosome (position 1).
    pub last_tmrca: f64,
}

/// Simulate one chromosome of n haplotypes. `keep` = 0 keeps every
/// segregating site; otherwise a uniform subset of `keep` sites (reservoir
/// sampling), returned in position order.
pub fn simulate_chromosome(
    n: usize,
    theta: f64,
    rho: f64,
    dem: &Demography,
    split: Option<Split>,
    keep: usize,
    seed: u64,
) -> std::result::Result<Sites, String> {
    let mut rng = Rng::new(seed);
    let mut tree = Tree::kingman(n, dem, split, &mut rng)?;
    let first_tmrca = tree.tmrca();
    let first_length = tree.total_length();
    let mut positions = Vec::new();
    let mut carriers = Vec::new();
    let mut n_total: u64 = 0;
    let mut x = 0.0;
    let both = theta + rho;
    if !both.is_finite() {
        return Err("theta + rho must be finite".into());
    }
    // expected number of events on the first tree's scale: refuse at once when
    // it already passes the budget (the loop budget catches later growth)
    if both * first_length > MAX_EVENTS as f64 {
        return Err(format!(
            "about {:.3e} mutation and recombination events expected on one chromosome \
             (more than {MAX_EVENTS}): theta and rho are too large for this sample",
            both * first_length
        ));
    }
    let mut n_events: u64 = 0;
    if both > 0.0 {
        loop {
            n_events += 1;
            if n_events > MAX_EVENTS {
                return Err(format!(
                    "more than {MAX_EVENTS} mutation and recombination events on one \
                     chromosome: theta and rho are too large for this sample"
                ));
            }
            let l = tree.total_length();
            x += rng.exp1() / (both * l);
            if x >= 1.0 {
                break;
            }
            if rng.unif() * both < theta {
                n_total += 1;
                if keep == 0 || positions.len() < keep {
                    positions.push(x);
                    carriers.push(tree.mutate(&mut rng));
                } else {
                    let j = (rng.unif() * n_total as f64) as u64;
                    if (j as usize) < keep {
                        positions[j as usize] = x;
                        carriers[j as usize] = tree.mutate(&mut rng);
                    }
                }
            } else {
                tree.recombine(dem, &mut rng)?;
            }
        }
    }
    let mut order: Vec<usize> = (0..positions.len()).collect();
    order.sort_by(|&a, &b| positions[a].total_cmp(&positions[b]));
    Ok(Sites {
        positions: order.iter().map(|&i| positions[i]).collect(),
        carriers: order.iter().map(|&i| carriers[i].clone()).collect(),
        n_total,
        first_tmrca,
        first_length,
        last_tmrca: tree.tmrca(),
    })
}

/// Coalescent haplotypes of one chromosome.
///
/// Returns list(pos, hap, n_total, tmrca, length): `pos` in [0, 1); `hap` an
/// integer vector of length n_sites * n_hap, sites x haplotypes column-major
/// (1 = derived allele); `n_total` the number of segregating sites before
/// subsetting; `tmrca` and `length` of the tree at position 0, `tmrca_end`
/// of the tree at position 1. `split_n_first` > 0 splits the sample: the
/// first `split_n_first` haplotypes in one subpopulation, the rest in another,
/// isolated until they merge at `split_time` (4 N0 units; ms -I 2 / -ej).
/// @noRd
#[allow(clippy::too_many_arguments)]
#[extendr]
pub fn coalescent_chromosome_core(
    n_hap: i32,
    theta: f64,
    rho: f64,
    history_times: Vec<f64>,
    history_sizes: Vec<f64>,
    seg_sites: i32,
    seed: f64,
    split_n_first: i32,
    split_time: f64,
) -> std::result::Result<List, String> {
    if n_hap < 2 {
        return Err("n_hap must be at least 2".into());
    }
    if !(theta.is_finite() && theta >= 0.0 && rho.is_finite() && rho >= 0.0) {
        return Err("theta and rho must be finite and non-negative".into());
    }
    if seg_sites < 0 {
        return Err("seg_sites must be non-negative".into());
    }
    if !((0.0..=4_294_967_295.0).contains(&seed) && seed.fract() == 0.0) {
        return Err("seed must be a whole number in [0, 2^32)".into());
    }
    let dem = Demography::new(&history_times, &history_sizes)?;
    let n = n_hap as usize;
    let split = if split_n_first > 0 {
        Some(Split {
            n_first: split_n_first as usize,
            join: split_time,
        })
    } else {
        None
    };
    let sites = simulate_chromosome(n, theta, rho, &dem, split, seg_sites as usize, seed as u64)?;
    let s = sites.positions.len();
    let mut hap = vec![0i32; s * n];
    for (i, bits) in sites.carriers.iter().enumerate() {
        for h in 0..n {
            if bits[h / 64] >> (h % 64) & 1 == 1 {
                hap[i + h * s] = 1;
            }
        }
    }
    Ok(list!(
        pos = sites.positions,
        hap = hap,
        n_total = sites.n_total as f64,
        tmrca = sites.first_tmrca,
        length = sites.first_length,
        tmrca_end = sites.last_tmrca
    ))
}

extendr_module! {
    mod coalescent;
    fn coalescent_chromosome_core;
}

#[cfg(test)]
mod tests {
    use super::*;

    fn harmonic(n: usize) -> f64 {
        (1..n).map(|i| 1.0 / i as f64).sum()
    }

    #[test]
    fn xoshiro_reference_stream() {
        // xoshiro256++ from state {1, 2, 3, 4}; the first two outputs were checked
        // by hand from the recurrence (rotl(s0 + s3, 23) + s0, then the update).
        let mut r = Rng { s: [1, 2, 3, 4] };
        assert_eq!(r.next_u64(), 41_943_041);
        assert_eq!(r.next_u64(), 58_720_359);
        assert_eq!(r.next_u64(), 3_588_806_011_781_223);
    }

    #[test]
    fn unif_and_below_are_in_range() {
        let mut r = Rng::new(7);
        for _ in 0..10_000 {
            let u = r.unif();
            assert!((0.0..1.0).contains(&u));
            assert!(r.below(3) < 3);
        }
    }

    #[test]
    fn kingman_tmrca_mean_constant_size() {
        // E[T_MRCA] = 1 - 1/n in units of 4N (= 2(1 - 1/n) in 2N units)
        let dem = Demography::new(&[], &[]).unwrap();
        let n = 10;
        let reps = 20_000;
        let mut rng = Rng::new(1);
        let mean: f64 = (0..reps)
            .map(|_| Tree::kingman(n, &dem, None, &mut rng).unwrap().tmrca())
            .sum::<f64>()
            / reps as f64;
        let expected = 1.0 - 1.0 / n as f64;
        // sd of T_MRCA < 1 here; 5 standard errors
        assert!(
            (mean - expected).abs() < 5.0 / (reps as f64).sqrt(),
            "mean {mean}"
        );
    }

    #[test]
    fn kingman_length_mean_is_harmonic() {
        // E[L] = sum_{i<n} 1/i in units of 4N, so E[S] = theta * that
        let dem = Demography::new(&[], &[]).unwrap();
        let n = 8;
        let reps = 20_000;
        let mut rng = Rng::new(2);
        let mean: f64 = (0..reps)
            .map(|_| {
                Tree::kingman(n, &dem, None, &mut rng)
                    .unwrap()
                    .total_length()
            })
            .sum::<f64>()
            / reps as f64;
        assert!((mean - harmonic(n)).abs() < 0.03, "mean {mean}");
    }

    #[test]
    fn pair_tmrca_under_a_size_change() {
        // size 1 on [0, t1), l2 after; pair rate 2 / lambda, so the survival is
        // e^{-2u} before t1 and e^{-2 t1} e^{-2 (u - t1) / l2} after, and
        // E[T] = integral of the survival = (1 - e^{-2 t1}) / 2 + e^{-2 t1} l2 / 2.
        let (t1, l2) = (0.3, 5.0);
        let dem = Demography::new(&[t1], &[l2]).unwrap();
        let surv_t1 = (-2.0 * t1).exp();
        let expected = 0.5 * (1.0 - surv_t1) + surv_t1 * l2 / 2.0;
        let reps = 40_000;
        let mut rng = Rng::new(3);
        let mean: f64 = (0..reps)
            .map(|_| Tree::kingman(2, &dem, None, &mut rng).unwrap().tmrca())
            .sum::<f64>()
            / reps as f64;
        assert!(
            (mean - expected).abs() < 0.03,
            "mean {mean} expected {expected}"
        );
    }

    /// Walk recombinations only (theta = 0) to the end of the sequence and
    /// return the tree there.
    fn tree_at_end(n: usize, rho: f64, dem: &Demography, rng: &mut Rng) -> Tree {
        let mut tree = Tree::kingman(n, dem, None, rng).unwrap();
        let mut x = 0.0;
        loop {
            x += rng.exp1() / (rho * tree.total_length());
            if x >= 1.0 {
                return tree;
            }
            tree.recombine(dem, rng).unwrap();
        }
    }

    #[test]
    fn smc_prime_keeps_the_marginal_tree_kingman() {
        // The tree at a fixed *position* is Kingman whatever rho: E[L] = sum 1/i.
        // (After a fixed number of *events* it is length-biased instead,
        // E[L^2] / E[L], because events arrive at rate rho L.)
        let dem = Demography::new(&[], &[]).unwrap();
        let n = 6;
        let reps = 20_000;
        let mut rng = Rng::new(4);
        let mut acc = 0.0;
        let mut acc_t = 0.0;
        for _ in 0..reps {
            let tree = tree_at_end(n, 30.0, &dem, &mut rng);
            acc += tree.total_length();
            acc_t += tree.tmrca();
        }
        let mean = acc / reps as f64;
        // sd(L) = sqrt(sum 1/i^2) ~ 1.21; 5 standard errors
        assert!(
            (mean - harmonic(n)).abs() < 5.0 * 1.21 / (reps as f64).sqrt(),
            "mean {mean}"
        );
        let mean_t = acc_t / reps as f64;
        assert!(
            (mean_t - (1.0 - 1.0 / n as f64)).abs() < 0.03,
            "tmrca {mean_t}"
        );
    }

    #[test]
    fn marginal_tree_follows_the_demography_after_recombination() {
        // pair TMRCA at the end of a recombining sequence, size change at t1
        let (t1, l2) = (0.3, 5.0);
        let dem = Demography::new(&[t1], &[l2]).unwrap();
        let surv_t1 = (-2.0 * t1).exp();
        let expected = 0.5 * (1.0 - surv_t1) + surv_t1 * l2 / 2.0;
        let reps = 40_000;
        let mut rng = Rng::new(6);
        let mean: f64 = (0..reps)
            .map(|_| tree_at_end(2, 10.0, &dem, &mut rng).tmrca())
            .sum::<f64>()
            / reps as f64;
        assert!(
            (mean - expected).abs() < 0.05,
            "mean {mean} expected {expected}"
        );
    }

    #[test]
    fn tree_stays_consistent_through_recombination() {
        let dem = Demography::new(&[0.5], &[3.0]).unwrap();
        let mut rng = Rng::new(5);
        let mut tree = Tree::kingman(12, &dem, None, &mut rng).unwrap();
        for _ in 0..5_000 {
            tree.recombine(&dem, &mut rng).unwrap();
            // every leaf reaches the root; parents are older than children
            for leaf in 0..tree.n {
                let mut x = leaf;
                let mut steps = 0;
                while tree.parent[x] != NONE {
                    assert!(tree.time[tree.parent[x]] > tree.time[x]);
                    x = tree.parent[x];
                    steps += 1;
                    assert!(steps < 2 * tree.n);
                }
                assert_eq!(x, tree.root);
            }
            assert_eq!(tree.internal.len(), tree.n - 1);
            assert!(tree.internal.windows(2).all(|w| w[0].0 <= w[1].0));
            assert!(tree.internal.iter().all(|&(t, v)| tree.time[v] == t));
            let len: f64 = (0..tree.parent.len()).map(|c| tree.branch_len(c)).sum();
            assert!((len - tree.total_length()).abs() < 1e-9);
            assert_eq!(tree.lineages_at(tree.tmrca() + 1.0), 1);
        }
    }

    #[test]
    fn seeded_runs_reproduce() {
        let dem = Demography::new(&[], &[]).unwrap();
        let a = simulate_chromosome(10, 20.0, 15.0, &dem, None, 0, 99).unwrap();
        let b = simulate_chromosome(10, 20.0, 15.0, &dem, None, 0, 99).unwrap();
        assert_eq!(a.positions, b.positions);
        assert_eq!(a.carriers, b.carriers);
    }

    #[test]
    fn reservoir_keeps_the_requested_number_in_order() {
        let dem = Demography::new(&[], &[]).unwrap();
        let s = simulate_chromosome(20, 200.0, 50.0, &dem, None, 30, 11).unwrap();
        assert!(s.n_total > 30);
        assert_eq!(s.positions.len(), 30);
        assert!(s.positions.windows(2).all(|w| w[0] <= w[1]));
        // every kept site segregates
        for bits in &s.carriers {
            let ones: u32 = bits.iter().map(|w| w.count_ones()).sum();
            assert!(ones >= 1 && (ones as usize) < 20);
        }
    }

    #[test]
    fn overflowing_rates_are_errors_not_hangs() {
        let dem = Demography::new(&[], &[]).unwrap();
        assert!(simulate_chromosome(2, 1e308, 1e308, &dem, None, 0, 1).is_err());
    }

    #[test]
    fn extreme_history_sizes_are_errors_not_panics() {
        // Codex probe: time 5e-324, size 1e308; seed 41 overflowed to Inf and
        // stepped past the last epoch (index panic). Now an error for every seed.
        // and time 5e-324, size 1e-308 collapsed waiting times to 0. Both are
        // now refused when the history is built.
        assert!(Demography::new(&[5e-324], &[1e308]).is_err());
        assert!(Demography::new(&[5e-324], &[1e-308]).is_err());
        assert!(Demography::new(&[1e13], &[1.0]).is_err());
        assert!(Demography::new(&[1e12], &[1e9]).is_ok());
        assert!(Demography::new(&[1e-300], &[1e-9]).is_ok());
        // the edges of the range still give finite trees for every seed
        let dem = Demography::new(&[1e-300], &[1e9]).unwrap();
        for seed in 0..200 {
            let s = simulate_chromosome(4, 0.0, 0.0, &dem, None, 0, seed).unwrap();
            assert!(s.first_tmrca.is_finite() && s.first_tmrca > 0.0);
        }
        let dem = Demography::new(&[1e-300], &[1e-9]).unwrap();
        for seed in 0..200 {
            let s = simulate_chromosome(4, 0.0, 0.0, &dem, None, 0, seed).unwrap();
            assert!(s.first_tmrca > 1e-300);
        }
    }

    #[test]
    fn lineage_sampler_is_uniform_over_the_lineages_at_tau() {
        let dem = Demography::new(&[], &[]).unwrap();
        let mut rng = Rng::new(21);
        let tree = Tree::kingman(9, &dem, None, &mut rng).unwrap();
        for &tau in &[
            1e-6,
            0.05,
            0.2,
            0.5,
            tree.tmrca() * 0.99,
            tree.tmrca() + 1.0,
        ] {
            let crossing = tree.lineages_crossing(tau);
            let draws = 60_000;
            let mut counts = vec![0usize; tree.parent.len()];
            for _ in 0..draws {
                let c = tree.sample_lineage_at(tau, 0, &mut rng);
                assert!(crossing.contains(&c));
                counts[c] += 1;
            }
            let expected = draws as f64 / crossing.len() as f64;
            for &c in &crossing {
                // 6 binomial standard deviations
                let sd = (expected * (1.0 - 1.0 / crossing.len() as f64)).sqrt();
                assert!((counts[c] as f64 - expected).abs() < 6.0 * sd + 1.0);
            }
        }
    }

    #[test]
    fn two_locus_correlation_is_smc_prime_not_smc() {
        // Pair TMRCAs at the two ends of a sequence of scaled length rho. The
        // exact ARG gives Corr = (rho + 18) / (rho^2 + 13 rho + 18) and the SMC
        // 1 / (1 + rho); SMC' lies just below the ARG value (Wilton et al. 2015).
        let dem = Demography::new(&[], &[]).unwrap();
        let rho = 1.0;
        let reps = 40_000;
        let mut rng = Rng::new(31);
        let (mut sx, mut sy, mut sxx, mut syy, mut sxy) = (0.0, 0.0, 0.0, 0.0, 0.0);
        for _ in 0..reps {
            let mut tree = Tree::kingman(2, &dem, None, &mut rng).unwrap();
            let a = tree.tmrca();
            let mut x = 0.0;
            loop {
                x += rng.exp1() / (rho * tree.total_length());
                if x >= 1.0 {
                    break;
                }
                tree.recombine(&dem, &mut rng).unwrap();
            }
            let b = tree.tmrca();
            sx += a;
            sy += b;
            sxx += a * a;
            syy += b * b;
            sxy += a * b;
        }
        let r = reps as f64;
        let cov = sxy / r - (sx / r) * (sy / r);
        let corr = cov / ((sxx / r - (sx / r).powi(2)) * (syy / r - (sy / r).powi(2))).sqrt();
        let arg = (rho + 18.0) / (rho * rho + 13.0 * rho + 18.0);
        let smc = 1.0 / (1.0 + rho);
        let se = (1.0 - corr * corr) / r.sqrt();
        assert!(corr > smc + 4.0 * se, "corr {corr} vs SMC {smc}");
        assert!(corr < arg + 3.0 * se, "corr {corr} vs ARG {arg}");
    }

    fn check_split_invariants(tree: &Tree) {
        let sp = tree.split.unwrap();
        for v in 0..tree.parent.len() {
            if v < tree.n {
                assert_eq!(tree.deme[v], u8::from(v >= sp.n_first));
            } else if tree.time[v] < sp.join {
                let [a, b] = tree.children[v];
                assert!(tree.deme[v] < MERGED);
                assert_eq!(tree.deme[a], tree.deme[v]);
                assert_eq!(tree.deme[b], tree.deme[v]);
            } else {
                assert_eq!(tree.deme[v], MERGED);
            }
        }
        for d in 0..2 {
            let mut want: Vec<f64> = (tree.n..tree.parent.len())
                .filter(|&v| tree.deme[v] == d as u8)
                .map(|v| tree.time[v])
                .collect();
            want.sort_by(f64::total_cmp);
            assert_eq!(tree.deme_times[d], want);
        }
        assert!(tree.tmrca() >= sp.join);
    }

    #[test]
    fn split_keeps_demes_isolated_until_the_join() {
        let dem = Demography::new(&[0.4], &[2.0]).unwrap();
        let split = Some(Split {
            n_first: 5,
            join: 0.3,
        });
        let mut rng = Rng::new(41);
        let mut tree = Tree::kingman(11, &dem, split, &mut rng).unwrap();
        check_split_invariants(&tree);
        for _ in 0..5_000 {
            tree.recombine(&dem, &mut rng).unwrap();
            check_split_invariants(&tree);
        }
    }

    #[test]
    fn between_deme_pair_tmrca_is_join_plus_half() {
        // one haplotype per deme: no coalescence before the join, then rate 2
        // at constant size, so E[T] = join + 1/2 -- at position 0 and, after
        // recombination along the sequence, at position 1 too.
        let dem = Demography::new(&[], &[]).unwrap();
        let join = 0.7;
        let split = Some(Split { n_first: 1, join });
        let reps = 40_000;
        let mut rng = Rng::new(51);
        let (mut a, mut b) = (0.0, 0.0);
        for _ in 0..reps {
            let mut tree = Tree::kingman(2, &dem, split, &mut rng).unwrap();
            a += tree.tmrca();
            let mut x = 0.0;
            loop {
                x += rng.exp1() / (5.0 * tree.total_length());
                if x >= 1.0 {
                    break;
                }
                tree.recombine(&dem, &mut rng).unwrap();
            }
            b += tree.tmrca();
        }
        // sd of T = 1/2; 5 standard errors
        let tol = 5.0 * 0.5 / (reps as f64).sqrt();
        assert!(
            (a / reps as f64 - (join + 0.5)).abs() < tol,
            "start {}",
            a / reps as f64
        );
        assert!(
            (b / reps as f64 - (join + 0.5)).abs() < tol,
            "end {}",
            b / reps as f64
        );
    }

    #[test]
    fn within_deme_sample_before_a_late_join_is_kingman() {
        // all but one haplotype in deme 0 and a join far in the past: the deme-0
        // subtree is Kingman, E[T_MRCA of deme 0] = 1 - 1/n0, and the tree root
        // sits above the join
        let dem = Demography::new(&[], &[]).unwrap();
        let n0 = 6;
        let split = Some(Split {
            n_first: n0,
            join: 50.0,
        });
        let reps = 20_000;
        let mut rng = Rng::new(61);
        let mut acc = 0.0;
        for _ in 0..reps {
            let tree = Tree::kingman(n0 + 1, &dem, split, &mut rng).unwrap();
            acc += tree.deme_times[0].last().copied().unwrap();
        }
        let mean = acc / reps as f64;
        assert!((mean - (1.0 - 1.0 / n0 as f64)).abs() < 0.02, "mean {mean}");
    }

    #[test]
    fn split_inputs_are_validated() {
        let dem = Demography::new(&[], &[]).unwrap();
        let mut rng = Rng::new(1);
        assert!(Tree::kingman(
            4,
            &dem,
            Some(Split {
                n_first: 0,
                join: 1.0
            }),
            &mut rng
        )
        .is_err());
        assert!(Tree::kingman(
            4,
            &dem,
            Some(Split {
                n_first: 4,
                join: 1.0
            }),
            &mut rng
        )
        .is_err());
        assert!(Tree::kingman(
            4,
            &dem,
            Some(Split {
                n_first: 2,
                join: 0.0
            }),
            &mut rng
        )
        .is_err());
        assert!(Tree::kingman(
            4,
            &dem,
            Some(Split {
                n_first: 2,
                join: f64::NAN
            }),
            &mut rng
        )
        .is_err());
    }
}
