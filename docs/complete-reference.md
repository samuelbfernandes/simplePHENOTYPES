Complete reference: every option in one place
================

This vignette exercises **every option added in version 2**, end to end,
in one runnable document. It is deliberately exhaustive rather than
gentle: if you want a guided introduction, start with
`vignette("simplePHENOTYPES-v2")`, then
`vignette("breeding-populations")` and `vignette("genetic-maps")`.

Everything below actually runs when this vignette is built, so the
printed output is real, not illustrative.

``` r
library(simplePHENOTYPES)
data("SNP55K_maize282_maf04")
geno <- SNP55K_maize282_maf04
dim(geno)
#> [1] 10650   285
```

The first five columns are marker metadata (`snp`, `allele`, `chr`,
`pos`, `cm`); the rest are individuals.

``` r
geno[1:5, 1:8]
#>           snp allele chr     pos        cm 4226 4722 33-16
#> 1 ss196422159    A/G   1  379844 0.0000000   -1   -1    -1
#> 2 ss196422171    G/A   1  613257 0.2482394   -1   -1    -1
#> 3 ss196422173    G/A   1  659354 0.2972609    1    1    -1
#> 4 ss196422186    G/C   1  992572 0.6515820    1   -1    -1
#> 5 ss196500940    G/A   1 2044264 1.7694445    1    1     1
```

# 1. Genotype input

## Converting from other formats

`as_numeric()` reads HapMap, VCF, GDS and PLINK bed/ped and returns the
numeric format used throughout.

``` r
hmp <- data.frame(
  `rs#`   = c("snp1", "snp2", "snp3"),
  alleles = c("C/T", "A/G", "G/T"),
  chrom   = c(1, 1, 2),
  pos     = c(100, 500, 200),
  strand  = "+", `assembly#` = NA, center = NA,
  protLSID = NA, assayLSID = NA, panelLSID = NA, QCcode = NA,
  L1 = c("C", "A", "G"), L2 = c("T", "G", "T"),
  L3 = c("C", "A", "T"), L4 = c("T", "G", "G"),
  check.names = FALSE, stringsAsFactors = FALSE
)

as_numeric(hmp, to_r = TRUE, verbose = FALSE)
#>    snp allele chr pos cm L1 L2 L3 L4
#> 1 snp1    C/T   1 100 NA  1 -1  1 -1
#> 2 snp2    A/G   1 500 NA  1 -1  1 -1
#> 3 snp3    G/T   2 200 NA  1 -1 -1  1
```

Coding schemes and allele orientation:

``` r
# 0 / 1 / 2 instead of -1 / 0 / 1
as_numeric(hmp, to_r = TRUE, code_as = "012", verbose = FALSE)[, 1:8]
#>    snp allele chr pos cm L1 L2 L3
#> 1 snp1    C/T   1 100 NA  2  0  2
#> 2 snp2    A/G   1 500 NA  2  0  2
#> 3 snp3    G/T   2 200 NA  2  0  0

# Orient against a supplied reference allele rather than by frequency
as_numeric(hmp, to_r = TRUE, verbose = FALSE,
           method = "reference", ref_allele = c("C", "A", "G"))[, 1:8]
#>    snp allele chr pos cm L1 L2 L3
#> 1 snp1    C/T   1 100 NA  1 -1  1
#> 2 snp2    A/G   1 500 NA  1 -1  1
#> 3 snp3    G/T   2 200 NA  1 -1 -1
```

Reading from files (not run here, since it needs your data):

``` r
as_numeric("my_genotypes.hmp.txt", to_r = TRUE)
as_numeric("my_genotypes.vcf",     to_r = TRUE)
as_numeric("my_genotypes.bed",     to_r = TRUE)   # needs .bim and .fam alongside
```

## Filtering markers with `filter_geno()`

`filter_geno()` is one place to trim a genotype before simulating: it
drops monomorphic markers, filters by minor-allele frequency and by
heterozygosity, and prunes by linkage disequilibrium. A marker must pass
every requested filter.

``` r
data("SNP55K_maize282_maf04")
geno <- SNP55K_maize282_maf04

# MAF window + only markers that carry heterozygotes (e.g. before dominance)
g <- filter_geno(geno, maf_below = 0.48, hets = "include", verbose = FALSE)
nrow(g)
#> [1] 3884
```

Heterozygosity matters for dominance: a dominance deviation acts on
heterozygotes, so on a near-inbred panel `dominance()` errors unless the
loci have heterozygotes. Pre-filtering fixes that:

``` r
simulate_phenotype(g, seed = 60) |>
  additive(prop = 0.4, n_qtn = 5) |>
  dominance(prop = 0.1)
#> <phenotype_sim>  (realized · long format)
#>   Genotypes: g   Traits: 1   Architecture: independent   Seed: 60
#>   Variance partition (proportions of V_P):
#>     additive   0.40   (5 QTNs, geometric)
#>     dominance  0.10   (same QTNs as additive)
#>     residual   0.50
#>   Requested genetic share = 0.50   realized H² = 0.48
```

**LD pruning and haplotype blocks.** The three sliding-window pruners
take a `c(window, step, threshold)` specification; `blocks` groups
markers into Gabriel *et al.* (2002) haplotype blocks and keeps one tag
per block. (Shown on a single chromosome to keep the vignette quick.)

``` r
chr1 <- geno[geno$chr == geno$chr[1], ][1:400, ]

# pairwise composite-r2 pruning (r2 > 0.2 in 50-marker windows)
nrow(filter_geno(chr1, indep_pairwise = c(50, 5, 0.2), verbose = FALSE))
#> [1] 273

# phased (EM haplotype) r2 pruning
nrow(filter_geno(chr1, indep_pairphase = c(50, 5, 0.2), verbose = FALSE))
#> [1] 273

# variance-inflation-factor pruning (PLINK's --indep)
nrow(filter_geno(chr1, indep = c(50, 5, 2), verbose = FALSE))
#> [1] 264

# Gabriel haplotype blocks: keep one tag marker per block
nrow(filter_geno(chr1, blocks = TRUE, block_max_kb = 2000, verbose = FALSE))
#> [1] 329
```

Using the VIF pruner or the block definition prints a one-per-session
citation for the method it implements (PLINK’s VIF procedure, or Gabriel
*et al.* for the blocks); the pairwise pruners are generic r^2
operations and print nothing.

# 2. The genetic map

Meiosis needs distances in centiMorgans. `SNP55K_maize282_maf04` ships
with a **synthetic** map in `cm` — modelled from the physical positions,
*not* a published maize linkage map, and not to be quoted as measured
recombination distance.

``` r
chrom_length <- aggregate(cm ~ chr, data = geno, FUN = max)
chrom_length
#>    chr       cm
#> 1    1 219.3573
#> 2    2 172.2582
#> 3    3 168.4068
#> 4    4 175.4795
#> 5    5 158.9096
#> 6    6 123.2703
#> 7    7 128.5558
#> 8    8 128.1356
#> 9    9 114.2056
#> 10  10 108.5956

round(sum(chrom_length$cm))             # total genetic length, cM
#> [1] 1497
```

`synthetic_map()` builds such a map for your own data. Every argument:

``` r
# Default: pericentromeric suppression, length from cm_per_mb
cm_default <- synthetic_map(geno$chr, geno$pos)

# Uniform map: cM strictly proportional to bp
cm_uniform <- synthetic_map(geno$chr, geno$pos, suppression = 0)

# Fixed total length per chromosome
cm_fixed <- synthetic_map(geno$chr, geno$pos,
                          total_cm = c(200, 180, 175, 170, 160,
                                       130, 130, 125, 115, 110))

# Explicit centromere positions, and a wider/stronger suppressed region
chr1 <- geno$pos[geno$chr == 1]
mid <- aggregate(pos ~ chr, data = geno, FUN = median)
cm_custom <- synthetic_map(geno$chr, geno$pos,
                           centromere = mid$pos,
                           suppression = 0.95, width = 0.25)

round(c(default = max(cm_default[geno$chr == 1]),
        uniform = max(cm_uniform[geno$chr == 1]),
        fixed   = max(cm_fixed[geno$chr == 1]),
        custom  = max(cm_custom[geno$chr == 1])), 1)
#> default uniform   fixed  custom 
#>   219.4   219.4   200.0   219.4
```

The suppression model concentrates recombination on the chromosome arms:

``` r
along_chr1 <- data.frame(
  pos = chr1,
  cm  = cm_default[geno$chr == 1],
  bin = cut(chr1, breaks = 8, labels = FALSE)
)

span <- aggregate(cbind(cm, pos) ~ bin, data = along_chr1,
                  FUN = function(x) max(x) - min(x))
span$cM_per_Mb <- round(span$cm / (span$pos / 1e6), 2)

span[, c("bin", "cM_per_Mb")]     # low near the centromere, high on the arms
#>   bin cM_per_Mb
#> 1   1      1.06
#> 2   2      0.96
#> 3   3      0.65
#> 4   4      0.26
#> 5   5      0.26
#> 6   6      0.65
#> 7   7      0.96
#> 8   8      1.06
```

To use your own real map, just assign it:

``` r
geno$cm <- my_real_map_in_cM
```

# 3. The simulation grammar

A simulation starts with `simulate_phenotype()` and gains variance
components as layers. Each `prop` is a proportion of total phenotypic
variance, so heritability is the sum of the genetic proportions.

``` r
ph <- simulate_phenotype(geno, seed = 1) |>
  additive(prop = 0.5, n_qtn = 3)
ph
#> <phenotype_sim>  (realized · long format)
#>   Genotypes: geno   Traits: 1   Architecture: independent   Seed: 1
#>   Variance partition (proportions of V_P):
#>     additive   0.50   (3 QTNs, geometric)
#>     residual   0.50
#>   Requested genetic share = 0.50   realized H² = 0.48
```

## Every layer

``` r
full <- simulate_phenotype(geno, seed = 2) |>
  additive(prop = 0.3, n_qtn = 5) |>
  dominance(prop = 0.1) |>
  epistasis(prop = 0.1, n_pairs = 2, interaction = 2) |>
  vqtl(prop = 0.1)
#> In addition to citing:
#> Fernandes, S.B., Lipka, A.E. simplePHENOTYPES: SIMulation of pleiotropic, linked and epistatic phenotypes. BMC Bioinformatics 21, 491 (2020). https://doi.org/10.1186/s12859-020-03804-y
#> Please also cite:
#> Murphy et al. (2022), Heredity 129:93-102, doi:10.1038/s41437-022-00541-1, for plant vQTL simulation.
#> This message is displayed once per session.
full
#> <phenotype_sim>  (realized · long format)
#>   Genotypes: geno   Traits: 1   Architecture: independent   Seed: 2
#>   Variance partition (proportions of V_P):
#>     additive   0.30   (5 QTNs, geometric)
#>     dominance  0.10   (same QTNs as additive)
#>     epistasis  0.10   (2 pairs, 2-way)
#>     vqtl       0.10   (same QTNs as additive)
#>     residual   0.40
#>   Requested genetic share = 0.50   realized H² = 0.51
#>   (vqtl is a residual-heterogeneity component and is not counted in
#>    broad-sense heritability)
```

## Layer options

``` r
# Geometric effects with an explicit base (0.5, 0.25, 0.125, ...)
simulate_phenotype(geno, seed = 3) |>
  additive(prop = 0.5, n_qtn = 3, effect = 0.5)
#> <phenotype_sim>  (realized · long format)
#>   Genotypes: geno   Traits: 1   Architecture: independent   Seed: 3
#>   Variance partition (proportions of V_P):
#>     additive   0.50   (3 QTNs, geometric)
#>     residual   0.50
#>   Requested genetic share = 0.50   realized H² = 0.50

# Explicit effect series, one per QTN
simulate_phenotype(geno, seed = 3) |>
  additive(prop = 0.5, n_qtn = 3, effect = c(0.6, 0.3, 0.1))
#> <phenotype_sim>  (realized · long format)
#>   Genotypes: geno   Traits: 1   Architecture: independent   Seed: 3
#>   Variance partition (proportions of V_P):
#>     additive   0.50   (3 QTNs, geometric)
#>     residual   0.50
#>   Requested genetic share = 0.50   realized H² = 0.50

# Dominance and vQTL on their own loci rather than the additive ones
simulate_phenotype(geno, seed = 4) |>
  additive(prop = 0.3, n_qtn = 5) |>
  dominance(prop = 0.1, same_as_add = FALSE, n_qtn = 3) |>
  vqtl(prop = 0.1, same_as_add = FALSE, n_qtn = 2)
#> <phenotype_sim>  (realized · long format)
#>   Genotypes: geno   Traits: 1   Architecture: independent   Seed: 4
#>   Variance partition (proportions of V_P):
#>     additive   0.30   (5 QTNs, geometric)
#>     dominance  0.10   (3 QTNs)
#>     vqtl       0.10   (2 QTNs)
#>     residual   0.50
#>   Requested genetic share = 0.40   realized H² = 0.42
#>   (vqtl is a residual-heterogeneity component and is not counted in
#>    broad-sense heritability)

# Three-way epistasis
simulate_phenotype(geno, seed = 5) |>
  additive(prop = 0.3, n_qtn = 4) |>
  epistasis(prop = 0.2, n_pairs = 2, interaction = 3)
#> <phenotype_sim>  (realized · long format)
#>   Genotypes: geno   Traits: 1   Architecture: independent   Seed: 5
#>   Variance partition (proportions of V_P):
#>     additive   0.30   (4 QTNs, geometric)
#>     epistasis  0.20   (2 pairs, 3-way)
#>     residual   0.50
#>   Requested genetic share = 0.50   realized H² = 0.48
```

## Replications and per-trait proportions

``` r
reps <- simulate_phenotype(geno, n_traits = 2, n_reps = 3, seed = 6) |>
  additive(prop = c(0.6, 0.3), n_qtn = 4)     # different h2 per trait
table(phenotypes_long(reps)$rep)
#> 
#>   1   2   3 
#> 560 560 560
head(phenotypes_wide(reps))
#>      id rep    Trait_1    Trait_2
#> 1  4226   1 -0.3366227  1.4717462
#> 2  4722   1 -1.5119297  0.3993162
#> 3 33-16   1 -0.3033333  0.2544596
#> 4 38-11   1  0.9029617  0.9313731
#> 5  A188   1 -0.9249267  1.5941070
#> 6  A239   1  1.8810237 -1.9487912
```

## One call instead of a pipe

``` r
simulate_phenotype(geno, h2 = 0.5, n_qtn = 3, seed = 7)              # model "A"
#> <phenotype_sim>  (realized · long format)
#>   Genotypes: geno   Traits: 1   Architecture: independent   Seed: 7
#>   Variance partition (proportions of V_P):
#>     additive   0.50   (3 QTNs, geometric)
#>     residual   0.50
#>   Requested genetic share = 0.50   realized H² = 0.48
simulate_phenotype(geno, h2 = 0.6, n_qtn = 4, model = "AD", seed = 7)
#> <phenotype_sim>  (realized · long format)
#>   Genotypes: geno   Traits: 1   Architecture: independent   Seed: 7
#>   Variance partition (proportions of V_P):
#>     additive   0.30   (4 QTNs, geometric)
#>     dominance  0.30   (same QTNs as additive)
#>     residual   0.40
#>   Requested genetic share = 0.60   realized H² = 0.59
```

## Fine control: fixed loci, phase, means, subsets

Set the causal loci for one layer while the rest are drawn at random.
Give marker names or column indices:

``` r
simulate_phenotype(geno, h2 = 0.5, seed = 1) |>
  additive(qtn = c("ss196442916", "ss196439337", "ss196480535"))
#> <phenotype_sim>  (realized · long format)
#>   Genotypes: geno   Traits: 1   Architecture: independent   Seed: 1
#>   Variance partition (proportions of V_P):
#>     additive   0.50   (3 QTNs, geometric)
#>     residual   0.50
#>   Requested genetic share = 0.50   realized H² = 0.48
```

`phase = "repulsion"` alternates the effect signs so linked increasing
and decreasing alleles oppose one another (the realized variance is
still `prop`; only the sign structure changes):

``` r
ph_rep <- simulate_phenotype(geno, h2 = 0.5, seed = 1) |>
  additive(n_qtn = 6, phase = "repulsion")
qtn_table(ph_rep)$effect
#> [1]  0.500000 -0.250000  0.125000 -0.062500  0.031250 -0.015625
```

A per-trait intercept with `mean`, and a subset of individuals with
`individuals`:

``` r
ph_m <- simulate_phenotype(geno, n_traits = 2, h2 = 0.5, mean = c(10, 20),
                           seed = 1) |>
  additive(n_qtn = 5)
colMeans(phenotypes_wide(ph_m)[, c("Trait_1", "Trait_2")])
#> Trait_1 Trait_2 
#>      10      20
```

``` r
ph_s <- simulate_phenotype(geno, h2 = 0.5,
                           individuals = c("33-16", "38-11", "A188", "B73"),
                           seed = 1) |>
  additive(n_qtn = 3)
n_individuals_sim <- ph_s$n_ind
n_individuals_sim
#> [1] 4
```

`vary_qtn = TRUE` gives each replication its own QTNs, and
`distinct_chr = TRUE` puts each trait’s QTNs on separate chromosomes:

``` r
ph_v <- simulate_phenotype(geno, h2 = 0.5, n_reps = 3, vary_qtn = TRUE,
                           seed = 1) |>
  additive(n_qtn = 4)
# a different QTN set per replication
lengths(ph_v$layers[[1]]$qtn_reps)
#> [1] 1 1 1

ph_d <- simulate_phenotype(geno, n_traits = 2, h2 = 0.4, distinct_chr = TRUE,
                           seed = 1) |>
  additive(n_qtn = 3)
lapply(ph_d$layers[[1]]$qtn, function(i) unique(geno$chr[i]))
#> [[1]]
#> [1] 5 3
#> 
#> [[2]]
#> [1]  2  8 10
```

# 4. Genetic architectures

## Independent

Trait-specific QTNs, no enforced sharing.

``` r
ind <- simulate_phenotype(geno, n_traits = 2, h2 = 0.4, seed = 8) |>
  additive(n_qtn = 5)

qtn <- ind$layers[[1]]$qtn
head(qtn[[1]], 3)      # trait 1
#> [1] 9831 7675 8927
head(qtn[[2]], 3)      # trait 2: different loci
#> [1] 2171 7850  903
```

## Pleiotropy with a target genetic correlation

Shared QTNs whose effects are drawn from a multivariate normal, so the
genetic correlation is the target in expectation. Works for **any number
of traits**.

``` r
ph2 <- simulate_phenotype(geno, architecture = "pleiotropy", n_traits = 2,
                          cor = 0.6, seed = 9) |>
  additive(prop = 0.5, n_qtn = 200)
#> In addition to citing:
#> Fernandes, S.B., Lipka, A.E. simplePHENOTYPES: SIMulation of pleiotropic, linked and epistatic phenotypes. BMC Bioinformatics 21, 491 (2020). https://doi.org/10.1186/s12859-020-03804-y
#> Please also cite:
#> Prado et al. (in preparation), when controlling the correlation in the pleiotropic architecture.
#> This message is displayed once per session.

g <- genetic_values(ph2)
round(cor(g[, 1], g[, 2]), 3)
#> [1] 0.617
```

Any number of traits, with a full target matrix including negative
correlations:

``` r
R <- matrix(c( 1.0,  0.8, -0.4,
               0.8,  1.0, -0.2,
              -0.4, -0.2,  1.0), 3, 3)

# Average over 20 seeds, since any single simulation scatters around the target
total <- matrix(0, 3, 3)
for (s in 1:20) {
  sim <- simulate_phenotype(geno, architecture = "pleiotropy", n_traits = 3,
                            cor = R, h2 = 0.5, seed = s) |>
    additive(n_qtn = 300)
  total <- total + cor(genetic_values(sim))
}

round(total / 20, 2)      # compare with R above
#>         Trait_1 Trait_2 Trait_3
#> Trait_1    1.00    0.78   -0.39
#> Trait_2    0.78    1.00   -0.18
#> Trait_3   -0.39   -0.18    1.00
```

Partial pleiotropy — only some of each trait’s genetic variance is
shared — via `pi`:

``` r
part <- simulate_phenotype(geno, architecture = "pleiotropy", n_traits = 3,
                           cor = 0.5, pi = c(1, 0.8, 0.6), seed = 10) |>
  additive(prop = 0.5, n_qtn = 100)

# Every trait keeps n_qtn loci; a shared core carries the covariance
qtn_pi <- part$layers[[1]]$qtn
lengths(qtn_pi)
#> [1] 100 100 100

shared <- intersect(intersect(qtn_pi[[1]], qtn_pi[[2]]), qtn_pi[[3]])
length(shared)
#> [1] 80
```

Concentrating variance in major pleiotropic loci. This is available but
costs precision in the realized correlation, so both arguments default
to zero:

``` r
simulate_phenotype(geno, architecture = "pleiotropy", n_traits = 2,
                   cor = 0.6, n_pleio_major = 5, prop_var_major = 0.5,
                   seed = 11) |>
  additive(prop = 0.5, n_qtn = 100)
#> <phenotype_sim>  (realized · long format)
#>   Genotypes: geno   Traits: 2   Architecture: pleiotropy   Seed: 11
#>   Variance partition (proportions of V_P):
#>     additive   0.50   (100 QTNs, geometric)
#>     residual   [0.50, 0.50]
#>   Requested genetic share = [0.50, 0.50]   realized H² = [0.49, 0.57]
```

Impossible requests are an error rather than a silent approximation:

``` r
simulate_phenotype(geno, architecture = "pleiotropy", n_traits = 2,
                   cor = 0.9, pi = 0.5, seed = 12) |>
  additive(prop = 0.5, n_qtn = 20)
#> Error:
#> ! Biological constraint violated: cor^2 (0.81) cannot exceed pi_target * pi_secondary (0.25).
```

## Linkage disequilibrium

The causal variant sits near, but not on, the marker a study would test.

``` r
simulate_phenotype(geno, architecture = "ld", n_traits = 2,
                   ld_type = "indirect", seed = 13) |>
  additive(prop = 0.4, n_qtn = 3)
#> <phenotype_sim>  (realized · long format)
#>   Genotypes: geno   Traits: 2   Architecture: ld   Seed: 13
#>   Variance partition (proportions of V_P):
#>     additive   0.40   (3 QTNs, geometric)
#>     residual   [0.60, 0.60]
#>   Requested genetic share = [0.40, 0.40]   realized H² = [0.43, 0.43]

simulate_phenotype(geno, architecture = "ld", n_traits = 2,
                   ld_type = "direct", r2_max = 0.9, r2_min = 0.4, seed = 14) |>
  additive(prop = 0.4, n_qtn = 3)
#> <phenotype_sim>  (realized · long format)
#>   Genotypes: geno   Traits: 2   Architecture: ld   Seed: 14
#>   Variance partition (proportions of V_P):
#>     additive   0.40   (3 QTNs, geometric)
#>     residual   [0.60, 0.60]
#>   Requested genetic share = [0.40, 0.40]   realized H² = [0.39, 0.40]
```

## Combining architectures

``` r
pleio <- simulate_phenotype(geno, architecture = "pleiotropy", n_traits = 2,
                            cor = 0.5, seed = 15) |>
  additive(prop = 0.4, n_qtn = 5)
indep <- simulate_phenotype(geno, n_traits = 2, seed = 15) |>
  additive(prop = 0.3, n_qtn = 5)

complex_phenotypes(pleio, indep, h2 = 0.5)
#> <phenotype_sim>  (realized · long format)
#>   Genotypes: geno   Traits: 2   Architecture: complex   Seed: 15
#>   Variance partition (proportions of V_P):
#>     combined from: pleiotropy + independent
#>     genetic    [0.50, 0.50]
#>     residual   [0.50, 0.50]
#>   Requested genetic share = [0.50, 0.50]   realized H² = [0.53, 0.51]
```

# 5. Breeding populations

Meiosis is simulated directly, so a pedigree can be built from real
founders and phenotyped without leaving the package.

``` r
pop <- as_population(geno, individuals = c("33-16", "38-11", "A188", "B73"))
pop
#> <Population>
#>   Individuals: 4   Markers: 10650   Chromosomes: 10
#>   Genetic map: 1497 cM total (109-219 cM per chromosome)
#>   Origin: founder
#>   IDs: 33-16, 38-11, A188, B73
n_individuals(pop)
#> [1] 4
dim(dosages(pop))
#> [1] 10650     4
```

## Crossing, selfing, doubled haploids

``` r
f1  <- cross(pop[1], pop[2], n = 1, seed = 21)
#> In addition to citing:
#> Fernandes, S.B., Lipka, A.E. simplePHENOTYPES: SIMulation of pleiotropic, linked and epistatic phenotypes. BMC Bioinformatics 21, 491 (2020). https://doi.org/10.1186/s12859-020-03804-y
#> Please also cite:
#> Toledo, F.H., Perez-Rodriguez, P., Crossa, J. and Burgueno, J. (2019). isqg: A Binary Framework for in Silico Quantitative Genetics. G3 9(8):2425-2428, doi:10.1534/g3.119.400373, when using the crossing pipelines (cross, selfcross, double_haploid).
#> This message is displayed once per session.
f2  <- selfcross(f1, n = 200, seed = 22)
dh  <- double_haploid(f1, n = 200, seed = 23)
bc1 <- cross(f1, pop[1], n = 100, seed = 24)      # backcross to founder 1

c(F1 = n_individuals(f1), F2 = n_individuals(f2),
  DH = n_individuals(dh), BC1 = n_individuals(bc1))
#>  F1  F2  DH BC1 
#>   1 200 200 100
```

Checking the genetics came out right:

``` r
d <- dosages(pop)
inf <- which(d[, 1] != d[, 2] & d[, 1] != 0 & d[, 2] != 0)  # informative loci

# F1 is heterozygous everywhere the founders differ
all(dosages(f1)[inf, 1] == 0)
#> [1] TRUE

# F2 segregates 1:2:1
round(prop.table(table(dosages(f2)[inf, ])), 3)
#> 
#>    -1     0     1 
#> 0.254 0.491 0.255

# Doubled haploids are fully homozygous, ~1:1
any(dosages(dh) == 0)
#> [1] FALSE
round(prop.table(table(dosages(dh)[inf, ])), 3)
#> 
#>    -1     1 
#> 0.498 0.502

# Backcross to founder 1: half heterozygous, half like the recurrent parent
round(prop.table(table(dosages(bc1)[inf, ])), 3)
#> 
#>    -1     0     1 
#> 0.225 0.508 0.267
```

## Recurrent selfing

Single-seed descent: each generation one plant is selfed, and the next
generation descends from one of its progeny.

``` r
gen <- f1
het <- numeric(5)
for (i in 1:5) {
  gen <- selfcross(gen[1], n = 100, seed = 30 + i)
  het[i] <- mean(dosages(gen)[inf, ] == 0)
}
round(het, 3)     # heterozygosity halves each generation
#> [1] 0.502 0.221 0.085 0.068 0.029
```

## Reproducibility

``` r
identical(dosages(selfcross(f1, n = 10, seed = 99)),
          dosages(selfcross(f1, n = 10, seed = 99)))
#> [1] TRUE
```

## Phenotyping a pedigree

A `Population` is accepted anywhere `simulate_phenotype()` takes
genotypes.

``` r
simulate_phenotype(f2, seed = 40) |>
  additive(prop = 0.6, n_qtn = 5)
#> <phenotype_sim>  (realized · long format)
#>   Genotypes: <Population: selfcross(prog_1)>   Traits: 1   Architecture: independent   Seed: 40
#>   Variance partition (proportions of V_P):
#>     additive   0.60   (5 QTNs, geometric)
#>     residual   0.40
#>   Requested genetic share = 0.60   realized H² = 0.58
```

# 6. Inspecting a simulation

Nothing is written to disk unless you ask for it. Everything the
original `create_phenotypes()` wrote to `Log_Sim.txt` and its companion
files lives in the returned object instead.

`print()` shows the requested variance budget alongside what was
actually realized:

``` r
# dominance needs heterozygous loci; on this near-inbred panel filter to them
# first (otherwise dominance() errors -- see filter_geno() above)
gh <- filter_geno(geno, hets = "include", verbose = FALSE)
ph_i <- simulate_phenotype(gh, h2 = 0.5, seed = 60) |>
  additive(prop = 0.3, n_qtn = 5) |>
  dominance(prop = 0.2)
ph_i
#> <phenotype_sim>  (realized · long format)
#>   Genotypes: gh   Traits: 1   Architecture: independent   Seed: 60
#>   Variance partition (proportions of V_P):
#>     additive   0.30   (5 QTNs, geometric)
#>     dominance  0.20   (same QTNs as additive)
#>     residual   0.50
#>   Requested genetic share = 0.50   realized H² = 0.49
```

`genetic_values()` returns the genetic component, before the residual —
the equivalent of the old `Genetic_values.txt`:

``` r
gv <- genetic_values(ph_i)
head(gv)
#>           Trait_1
#> 4226  -0.77902269
#> 4722   0.04996093
#> 33-16  0.17749688
#> 38-11 -0.97032660
#> A188   0.36880079
#> A239  -0.20511095

# Realized heritability, computed by hand from the same pieces
round(var(gv[, 1]) / var(phenotypes_wide(ph_i)$Trait_1), 3)
#> [1] 0.491
```

The `var_explained` column of `qtn_table()` gives each
additive/dominance QTN’s marginal share of phenotypic variance (`NA` for
epistatic sets, where variance is not per-locus). `qtn_table()` names
the causal loci and their effects — the equivalent of
`Additive_QTNs.txt` and `QTN_effects_summary.txt`. Epistatic sets are
expanded one row per member, sharing a `set` number and effect:

``` r
qtn_table(ph_i)
#>      trait     layer set         snp chr       pos       maf  effect
#> 1  Trait_1  additive  NA ss196502704   4  26399385 0.4875000 0.50000
#> 2  Trait_1  additive  NA ss196505289   8 118005722 0.4589286 0.25000
#> 3  Trait_1  additive  NA ss196488121   8 163213942 0.4303571 0.12500
#> 4  Trait_1  additive  NA ss196469274   6   6514321 0.4482143 0.06250
#> 5  Trait_1  additive  NA ss196513900   6 147443964 0.4375000 0.03125
#> 6  Trait_1 dominance  NA ss196502704   4  26399385 0.4875000 0.50000
#> 7  Trait_1 dominance  NA ss196505289   8 118005722 0.4589286 0.25000
#> 8  Trait_1 dominance  NA ss196488121   8 163213942 0.4303571 0.12500
#> 9  Trait_1 dominance  NA ss196469274   6   6514321 0.4482143 0.06250
#> 10 Trait_1 dominance  NA ss196513900   6 147443964 0.4375000 0.03125
#>    var_explained QTN_t1 QTN_t2 ld_r2
#> 1   0.2600834274   <NA>   <NA>    NA
#> 2   0.0646210924   <NA>   <NA>    NA
#> 3   0.0159487281   <NA>   <NA>    NA
#> 4   0.0040225778   <NA>   <NA>    NA
#> 5   0.0010006474   <NA>   <NA>    NA
#> 6   0.1442923666   <NA>   <NA>    NA
#> 7   0.0360730916   <NA>   <NA>    NA
#> 8   0.0090182729   <NA>   <NA>    NA
#> 9   0.0022545682   <NA>   <NA>    NA
#> 10  0.0005636421   <NA>   <NA>    NA
```

Everything else the log recorded is a field on the object:

``` r
ph_i$seed
#> [1] 60
ph_i$n_traits
#> [1] 1
ph_i$var_budget
#>     trait component prop
#> 1 Trait_1  additive  0.3
#> 2 Trait_1 dominance  0.2
#> 3 Trait_1  residual  0.5
```

Note that a `vqtl()` layer is an explicit residual-heterogeneity
variance share, not part of the broad-sense `h2` budget. It varies the
conditional residual variance rather than adding genetic value;
`print()` says so when one is present:

``` r
simulate_phenotype(geno, h2 = 0.5, seed = 61) |>
  additive(prop = 0.3, n_qtn = 4) |>
  vqtl(prop = 0.2)
#> <phenotype_sim>  (realized · long format)
#>   Genotypes: geno   Traits: 1   Architecture: independent   Seed: 61
#>   Variance partition (proportions of V_P):
#>     additive   0.30   (4 QTNs, geometric)
#>     vqtl       0.20   (same QTNs as additive)
#>     residual   0.50
#>   Requested genetic share = 0.30   realized H² = 0.27
#>   (vqtl is a residual-heterogeneity component and is not counted in
#>    broad-sense heritability)
```

`plot()` gives a four-panel diagnostic summary – variance partition,
phenotype distribution, QTN effects, and (for two or more traits) the
genetic-value scatter:

``` r
ph_plot <- simulate_phenotype(geno, n_traits = 2, h2 = 0.5, seed = 62) |>
  additive(prop = 0.4, n_qtn = 8) |>
  dominance(prop = 0.1)
plot(ph_plot)
```

![](/Users/samuelbf/Library/CloudStorage/OneDrive-UniversityofArkansas/UARK/collaboration/Software/simplePHENOTYPES/docs/complete-reference_files/figure-gfm/inspect-plot-1.png)<!-- -->

# 7. Output

## In memory

``` r
ph3 <- simulate_phenotype(geno, n_traits = 2, seed = 50) |>
  additive(prop = 0.5, n_qtn = 3)

head(phenotypes_long(ph3))
#>      id   trait rep      value
#> 1  4226 Trait_1   1 -0.7711339
#> 2  4722 Trait_1   1  1.8858533
#> 3 33-16 Trait_1   1 -0.3975453
#> 4 38-11 Trait_1   1  0.7776585
#> 5  A188 Trait_1   1 -0.2648815
#> 6  A239 Trait_1   1 -1.1256463
head(phenotypes_wide(ph3))
#>      id rep    Trait_1    Trait_2
#> 1  4226   1 -0.7711339 -0.8819540
#> 2  4722   1  1.8858533  0.1406755
#> 3 33-16   1 -0.3975453  0.9748847
#> 4 38-11   1  0.7776585  0.2032506
#> 5  A188   1 -0.2648815  0.1948636
#> 6  A239   1 -1.1256463 -1.2485879
```

## Written to disk

Each call writes exactly the file you name — nothing else is created
alongside it.

``` r
out_dir <- file.path(tempdir(), "sp_demo")
dir.create(out_dir, showWarnings = FALSE)

write_phenotypes(ph3, file = file.path(out_dir, "phenotypes_long.txt"))
write_phenotypes(ph3, file = file.path(out_dir, "phenotypes_wide.csv"),
                 format = "wide", sep = ",")

list.files(out_dir)
#> [1] "phenotypes_long.txt" "phenotypes_wide.csv"
head(read.delim(file.path(out_dir, "phenotypes_long.txt")))
#>      id   trait rep      value
#> 1  4226 Trait_1   1 -0.7711339
#> 2  4722 Trait_1   1  1.8858533
#> 3 33-16 Trait_1   1 -0.3975453
#> 4 38-11 Trait_1   1  0.7776585
#> 5  A188 Trait_1   1 -0.2648815
#> 6  A239 Trait_1   1 -1.1256463
```

# 8. Original implementation of create_phenotypes

`create_phenotypes()` is the original interface, kept frozen so
published results stay reproducible. It receives bug fixes only; new
work should use the grammar above.

Its output always lands in a folder of its own, so a directory is never
littered with loose files:

``` r
legacy_home <- file.path(tempdir(), "sp_legacy")
dir.create(legacy_home, showWarnings = FALSE)

invisible(create_phenotypes(
  geno_obj = geno, add_QTN_num = 3, add_effect = 0.2,
  rep = 2, h2 = 0.5, model = "A",
  home_dir = legacy_home, verbose = FALSE
))
#> Simulation completed!
#> Results are saved at:/var/folders/k3/kdzmjr9j38q0gt1j2rjy7pl00000gp/T//Rtmpjkncsp/sp_legacy/simplePHENOTYPES_output

list.files(legacy_home)                                  # one folder
#> [1] "simplePHENOTYPES_output"
list.files(file.path(legacy_home, "simplePHENOTYPES_output"))
#> [1] "Additive_QTNs.txt"                   "Genetic_values.txt"                 
#> [3] "Log_Sim.txt"                         "Simulated_Data_2_Reps_Herit_0.5.txt"
```

Naming the folder yourself, and re-running without overwriting:

``` r
invisible(create_phenotypes(
  geno_obj = geno, add_QTN_num = 3, add_effect = 0.2,
  rep = 1, h2 = 0.5, model = "A",
  home_dir = legacy_home, output_dir = "my_run", verbose = FALSE
))
#> Simulation completed!
#> Results are saved at:/var/folders/k3/kdzmjr9j38q0gt1j2rjy7pl00000gp/T//Rtmpjkncsp/sp_legacy/my_run
invisible(create_phenotypes(
  geno_obj = geno, add_QTN_num = 3, add_effect = 0.2,
  rep = 1, h2 = 0.5, model = "A",
  home_dir = legacy_home, output_dir = "my_run", verbose = FALSE
))
#> Directory name provided by 'output_dir' alredy exists! 
#> Creating: /var/folders/k3/kdzmjr9j38q0gt1j2rjy7pl00000gp/T//Rtmpjkncsp/sp_legacy/my_run(1)
#> Simulation completed!
#> Results are saved at:/var/folders/k3/kdzmjr9j38q0gt1j2rjy7pl00000gp/T//Rtmpjkncsp/sp_legacy/my_run(1)

list.files(legacy_home)      # my_run and my_run(1): nothing is overwritten
#> [1] "my_run"                  "my_run(1)"              
#> [3] "simplePHENOTYPES_output"
```

# Session information

``` r
sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: aarch64-apple-darwin23
#> Running under: macOS Tahoe 26.6.2
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRblas.0.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
#> 
#> locale:
#> [1] C.UTF-8/C.UTF-8/C.UTF-8/C/C.UTF-8/C.UTF-8
#> 
#> time zone: America/Chicago
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] simplePHENOTYPES_1.4.0-9001 testthat_3.3.2             
#> 
#> loaded via a namespace (and not attached):
#>  [1] vctrs_0.7.3         knitr_1.51          cli_3.6.6          
#>  [4] xfun_0.60           rlang_1.3.0         otel_0.2.0         
#>  [7] purrr_1.2.2         pkgload_1.5.3       data.table_1.18.6.1
#> [10] glue_1.8.1          htmltools_0.5.9     rprojroot_2.1.1    
#> [13] SNPRelate_1.46.0    pkgbuild_1.4.8      brio_1.1.5         
#> [16] rmarkdown_2.32      evaluate_1.0.5      ellipsis_0.3.3     
#> [19] fastmap_1.2.0       yaml_2.3.12         mvtnorm_1.4-2      
#> [22] lifecycle_1.0.5     RhpcBLASctl_0.23-42 gdsfmt_1.48.2      
#> [25] memoise_2.0.1       compiler_4.6.1      fs_2.1.0           
#> [28] sessioninfo_1.2.4   rstudioapi_0.19.0   digest_0.6.39      
#> [31] R6_2.6.1            pillar_1.11.1       usethis_3.2.1      
#> [34] magrittr_2.0.5      tools_4.6.1         devtools_2.5.2     
#> [37] cachem_1.1.0        desc_1.4.3
```
