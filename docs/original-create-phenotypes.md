Original implementation of create_phenotypes
================
simplePHENOTYPES v 1.3 (2026-09-10)

This short tutorial presents some of the possible genetic settings one
could simulate, but it certainly does not explore all the possibilities.
For more information on specific input parameters, please check the help
documentation (?create_phenotypes).

# Load Sample Data set

Note that the data set used in all vignettes is already in numeric
format. In addition to the numeric format, simplePHENOTYPES’ parameter
`geno_obj` also takes an R object in HapMap format as input. Other input
options are VCF, GDS, and Plink bed/ped. These last formats should be
loaded from file with `geno_file` or `geno_path`.

``` r
library(simplePHENOTYPES)
data("SNP55K_maize282_maf04")
SNP55K_maize282_maf04[1:8, 1:10]
```

# Single Trait

The simplest option is the simulation of univariate traits. In the
example below, we are simulating ten single trait experiments with a
heritability of 0.7. In this setting, the simulated trait is controlled
by one large-effect QTN (`big_add_QTN_effect = 0.9`) and two small
effect QTNs. The additive effects of these last two QTNs follow a
geometric series starting with 0.2. Thus, the effect size of the first
of these two QTNs is 0.2, and the effect size of the second is
0.2<sup>2</sup>. Results are being saved at a temporary directory
(`home_dir = tempdir()`). Please see help files (?create_phenotypes) to
see which default values are being used.

``` r
create_phenotypes(
  geno_obj = SNP55K_maize282_maf04,
  add_QTN_num = 3,
  add_effect = 0.2,
  big_add_QTN_effect = 0.9,
  rep = 10,
  h2 = 0.7,
  model = "A",
  home_dir = tempdir())
```

**V2 equivalent.** The grammar has no `big_add_QTN_effect`; give the
large-effect QTN its value explicitly through `effect`. Heritability
emerges as the sum of the layer `prop`s, so a single additive layer
takes all of `h2 = 0.7`. (The realized phenotypes will not match v1
number-for-number — the v2 engine has its own clean seed scheme,
DECISION-009 — but the architecture is the same.)

``` r
simulate_phenotype(SNP55K_maize282_maf04, seed = 1, n_reps = 10) |>
  additive(prop = 0.7, n_qtn = 3, effect = c(0.9, 0.2, 0.2^2))
```

# Multiple Traits: Pleiotropy Architecture

simplePHENOTYPES provides three multi-trait simulation scenarios:
pleiotropy, partial pleiotropy, and spurious pleiotropy. In this
example, we are simulating three (`ntraits = 3`) pleiotropic
(`architecture = "pleiotropic"`) trait controlled by three additive and
four dominance QTNs. The effect size of the largest-effect additive QTN
is 0.3 for all traits (`big_add_QTN_effect = c(0.3, 0.3, 0.3)`), while
the additive and dominance effect sizes are 0.04, 0.2, and 0.1 for each
trait, respectively. Heritability for trait_1 is 0.2, while the
heritability of the two correlated traits is 0.4. Each replicate is
being recorded in a different file (`output_format = "multi-file"`) in a
folder named “Results_Pleiotropic”. In this setting, we do not specify
the correlation between traits; instead, the observed (realized)
correlation is an artifact of different allelic effects for each trait.
The same QTNs are used to generate phenotypes in all ten replications
(`vary_QTN = FALSE`)(default); alternatively, we could select different
QTNs in each replicate using `vary_QTN = TRUE`. As mentioned above, the
first QTN of each trait will get the effect provided by
big_add_QTN_effect; all other QTNs will have the effect size assigned by
`add_effect` and `dom_effect`. The vector `add_effect` contains one
allelic effect for each trait, and a geometric series (default) is being
used to generate allelic effects for each one of the two additive QTNs
(`add_QTN_num = 3`) and three dominance QTNs (`dom_QTN_num = 4`). All
results will be saved to file, and a data.frame with all phenotypes will
be assigned to an object called “test1” (to_r = TRUE).

``` r
 test1 <-  create_phenotypes(
    geno_obj = SNP55K_maize282_maf04,
    add_QTN_num = 3,
    dom_QTN_num = 4,
    big_add_QTN_effect = c(0.3, 0.3, 0.3),
    h2 = c(0.2, 0.4, 0.4),
    add_effect = c(0.04,0.2,0.1),
    dom_effect = c(0.04,0.2,0.1),
    ntraits = 3,
    rep = 10,
    vary_QTN = FALSE,
    output_format = "multi-file",
    architecture = "pleiotropic",
    output_dir = "Results_Pleiotropic",
    to_r = TRUE,
    seed = 10,
    model = "AD",
    sim_method = "geometric",
  home_dir = tempdir()
  )
```

Optionally, we may input a list of allelic effects
(`sim_method = "custom"`). In the example below, a geometric series
(custom_geometric) is being assigned and should generate the same
simulated data as the previous example (all.equal(test1, test2)). Notice
that since `big_add_QTN_effect` is non-NULL, we only need to provide
effects for two out of the three simulated additive QTNs. On the other
hand, all four dominance QTN must have an effect assigned on the
custom_geometric_d list. Importantly, the allelic effects are assigned
to each trait based on the order they appear in the list and not based
on the names, i.e., ‘trait_1’, ‘trait_2’, and ‘trait_3’.

``` r
 custom_geometric_a <- list(trait_1 = c(0.04, 0.0016),
                         trait_2 = c(0.2, 0.04),
                         trait_3 = c(0.1, 0.01))
 custom_geometric_d <- list(trait_1 = c(0.04, 0.0016, 6.4e-05, 2.56e-06),
                         trait_2 = c(0.2, 0.04, 0.008, 0.0016),
                         trait_3 = c(0.1, 0.01, 0.001, 1e-04))

 test2 <-  create_phenotypes(
   geno_obj = SNP55K_maize282_maf04,
   add_QTN_num = 3,
   dom_QTN_num = 4,
   big_add_QTN_effect = c(0.3, 0.3, 0.3),
   h2 = c(0.2,0.4, 0.4),
   add_effect = custom_geometric_a,
   dom_effect = custom_geometric_d,
   ntraits = 3,
   rep = 10,
   vary_QTN = FALSE,
   output_format = "multi-file",
   architecture = "pleiotropic",
   output_dir = "Results_Pleiotropic",
   to_r = T,
   sim_method = "custom",
   seed = 10,
   model = "AD",
  home_dir = tempdir()
 )
 
 all.equal(test1, test2)
```

**V2 equivalent.** `architecture = "pleiotropy"` shares the QTNs across
traits, and the big difference is that **`cor` now sets the genetic
correlation directly** — in v1 the correlation was an uncontrolled
by-product of the per-trait effect sizes. Per-layer `prop` replaces the
`h2` + `add_effect` / `dom_effect` bookkeeping; a custom effect series
still goes through `effect =`.

``` r
simulate_phenotype(SNP55K_maize282_maf04, architecture = "pleiotropy",
                   n_traits = 3, cor = 0.5, seed = 10) |>
  additive(prop = 0.3, n_qtn = 3) |>
  dominance(prop = 0.1, same_as_add = TRUE)
```

# Multiple Traits: Partial Pleiotropy Architecture

In this example, we simulate 20 replicates of three partially
pleiotropic traits (`architecture = "partially"`), which are
respectively controlled by seven, 13, and four QTNs. All QTNs will have
additive effects that follow a geometric series, where the effect size
of the i<sup>th</sup> QTN is add_effect^i. For instance, trait_2 is
controlled by three pleiotropic additive QTNs and ten trait-specific
additive QTNs; consequently, the first pleiotropic additive QTN will
have an additive effect of 0.33 and the 13<sup>th</sup> trait-specific
additive QTN will have an effect of 0.33<sup>13</sup>. Correlation among
traits is assigned to be equal to the cor_matrix object. All 20
replicates of these three simulated traits will be saved in one file,
specifically in a long format and with an additional column named “Rep”.
Results will be saved in a directory called “Results_Partially”. In this
example, the genotype file will also be saved in numeric format.

``` r
cor_matrix <- matrix(c(   1, 0.3, -0.9,
                        0.3,   1,  -0.5,
                       -0.9, -0.5,    1 ), 3)

sim_results <- create_phenotypes(
  geno_obj = SNP55K_maize282_maf04,
  ntraits = 3,
  pleio_a = 3,
  pleio_e = 2,
  same_add_dom_QTN = TRUE,
  degree_of_dom = 0.5,
  trait_spec_a_QTN_num = c(4, 10, 1),
  trait_spec_e_QTN_num = c(3, 2, 5),
  h2 = c(0.2, 0.4, 0.8),
  add_effect = c(0.5, 0.33, 0.2),
  epi_effect = c(0.3, 0.3, 0.3),
  epi_interaction = 2,
  cor = cor_matrix,
  rep = 20,
  output_dir = "Results_Partially",
  output_format = "long",
  architecture = "partially",
  out_geno = "numeric",
  to_r = TRUE,
  model = "AE",
  home_dir = tempdir()
)
```

**V2 equivalent.** There is no `"partially"` architecture; partial
pleiotropy is built by combining a **shared** (`"pleiotropy"`) model
with a **trait-specific** (`"independent"`) model through
`complex_phenotypes()`. A grammar layer uses one `n_qtn` per call, so
the per-trait `trait_spec_a_QTN_num = c(4, 10, 1)` collapses to a single
count here (build one independent model per trait and combine them if
you truly need different counts).

``` r
cor_matrix <- matrix(c(1, 0.3, -0.9, 0.3, 1, -0.5, -0.9, -0.5, 1), 3)
shared   <- simulate_phenotype(SNP55K_maize282_maf04, architecture = "pleiotropy",
                               n_traits = 3, cor = cor_matrix, seed = 10) |>
  additive(prop = 0.5, n_qtn = 3)
specific <- simulate_phenotype(SNP55K_maize282_maf04, architecture = "independent",
                               n_traits = 3, seed = 10) |>
  additive(prop = 0.5, n_qtn = 5)
complex_phenotypes(shared, specific, h2 = c(0.2, 0.4, 0.8))
```

# Multiple Traits: Spurious Pleiotropy Architecture

Another architecture implemented is Spurious Pleiotropy. In this case,
we have two options: direct or indirect LD (`type_of_ld = "indirect"`).
In the example below, we simulate a case of indirect LD with five
replicates of two traits controlled by three additive QTNs each. For
each QTN, a marker is first selected (intermediate marker), and then two
separate markers (one upstream and another downstream) are picked to be
QTNs for each of the two traits. This QTN selection is based on an
r<sup>2</sup> threshold of at most 0.8 (`ld_max=0.8`) with the
intermediate marker. The three QTNs will have additive effects that
follow a geometric series, where the effect size of the i<sup>th</sup>
QTN is 0.02<sup>i</sup> for one trait and 0.05<sup>i</sup> for the other
trait. Starting seed number is 200, and output phenotypes are saved in
one file, but in a “wide” format with each replicate of two traits being
added as additional columns. Plink fam, bim, and bed files are also
saved at Results_LD.

``` r
create_phenotypes(
  geno_obj = SNP55K_maize282_maf04,
  add_QTN_num = 3,
  h2 = c(0.2, 0.4),
  add_effect = c(0.02, 0.05),
  rep = 5,
  seed = 200,
  output_format = "wide",
  architecture = "LD",
  output_dir = "Results_LD",
  out_geno = "BED",
  remove_QTN = TRUE,
  ld_max =0.8,
  ld_min =0.2,
  model = "A",
  ld_method = "composite",
  type_of_ld = "indirect",
  home_dir = tempdir()
)
```

**V2 equivalent.** `architecture = "ld"` (two traits, each with its own
causal locus in linkage disequilibrium). `type_of_ld` → `ld_type` (now
defaulting to `"direct"`), and `ld_max` / `ld_min` → `r2_max` /
`r2_min`. `qtn_table()` reports the linked pair (`QTN_t1`, `QTN_t2`) and
their `ld_r2`.

``` r
simulate_phenotype(SNP55K_maize282_maf04, architecture = "ld", n_traits = 2,
                   ld_type = "indirect", r2_max = 0.8, r2_min = 0.2, seed = 200) |>
  additive(prop = c(0.2, 0.4), n_qtn = 3)
```

# Multiple Traits: Partial Pleiotropy Architecture with other useful parameters

The example below simulates five replicates of three traits. In each
replicate, different SNPs are selected to be the QTNs for each
experiment (`vary_QTN = TRUE`). These traits are controlled by three
pleiotropic (`pleio = 3`) additive and dominance QTNs
(`same_add_dom_QTN = TRUE` and `degree_of_dom = 1`); two pleiotropic
epistatic QTNs (`pleio_e = 2`); four, ten and one trait-specific
additive and dominance QTNs (`trait_spec_a_QTN_num = c(4, 10, 1)`); and
two, one and five epistatic trait-specific epistatic QTNs
(`trait_spec_e_QTN_num = c(2, 1, 5)`). In addition to the default
parameters, each genetic architecture may be simulated with many
auxiliary features. For instance, we may be interested in outputting the
amount of variance explained by each simulated QTN
(`QTN_variance = TRUE`) or setting a residual correlation between traits
(`cor_res = residual`) and thus, change the default option of
independent residuals. Notice that in this example, the heritability is
a 2x3 matrix (`h2 = heritability`). Each column of the matrix
“heritability” will be assigned to a different trait. In this case,
simplePHENOTYPES will loop over each row of `h2`, keeping all other
variables constant. Since rep = 5 and nrow(h2) = 2, ten experiments will
be simulated and saved in separate files. Simulated results will be
saved as “.fam” files used as GEMMA input. Simultaneously, one genotypic
file without the QTNs for the simulated traits will be saved for each
replication. Due to the option `vary_QTN = TRUE`, each experiment will
be simulated with different QTNs; thus, if we opt for
`remove_QTN = TRUE`, many potentially large files will be saved in the
output_dir folder. By default, simplePHENOTYPES will ask us if all these
files should be saved. To avoid this question, we may use
`warning_file_saver = FALSE`. In the present example, ten plink bed
files (which is also the input for GEMMA) are saved. Genotypic files for
rep one will be named `SNP55K_maize282_maf04_noQTN_rep_1.bed`,
`SNP55K_maize282_maf04_noQTN_rep_1.bim`, and
`SNP55K_maize282_maf04_noQTN_rep_1.fam`, whereas the phenotypic file
will be saved as `Simulated_Data__Rep1_Herit_0.2_0.8_0.7.fam`.
Importantly, the file `SNP55K_maize282_maf04_noQTN_rep_1.fam` does not
contain the phenotypic data and needs to be replaced by
`Simulated_Data__Rep1_Herit_0.2_0.8_0.7.fam` prior to its use by GEMMA
or other software that uses bed files. A parameter particularly useful,
especially when simulating dominance, is `constraints`. Here we only
“include” heterozygote SNPs to be used as QTNs (
`constraints = list(maf_above = 0.3, maf_below = 0.44, hets = "include")`).
Optionally, we may “remove” all the heterozygotes from consideration.
The other constrain options used here are to select only QTNs with minor
allele frequency between 0.3 and 0.44.

``` r
residual <- matrix(c(1, 0.1,-0.2,
                     0.1, 1,-0.1,-0.2,-0.1, 1), 3)
heritability <- matrix(c(0.2, 0.4, 0.8,
                         0.6, 0.7, 0.2), 2)
create_phenotypes(
  geno_obj = SNP55K_maize282_maf04,
  pleio_a = 3,
  pleio_e = 2,
  same_add_dom_QTN = TRUE,
  degree_of_dom = 1,
  trait_spec_a_QTN_num = c(4, 10, 1),
  trait_spec_e_QTN_num = c(2, 1, 5),
  epi_effect = c(0.01, 0.4, 0.2),
  add_effect = c(0.3, 0.2, 0.5),
  h2 = heritability,
  ntraits = 3,
  rep = 5,
  vary_QTN = TRUE,
  warning_file_saver = FALSE,
  output_dir = "Results_Partially_ADE",
  output_format = "gemma",
  architecture = "partially",
  model = "ADE",
  QTN_variance = TRUE,
  remove_QTN = TRUE,
  home_dir = tempdir(),
  constraints = list(
    maf_above = 0.3,
    maf_below = 0.44,
    hets = "include"
  ),
  cor_res = residual
)
```

**V2 equivalent.** Build partial pleiotropy as above (a `"pleiotropy"`
model plus an `"independent"` model, joined with
`complex_phenotypes()`), adding an `epistasis()` layer for the “E” of
ADE. `QTN_variance = TRUE` becomes `qtn_table()`, which always reports a
`var_explained` column; `vary_QTN` becomes `vary_qtn`; and constraining
QTNs by MAF or heterozygosity is now done up front with `filter_geno()`
(e.g. `filter_geno(geno, maf_above = 0.3, maf_below = 0.44, hets = "include")`).
Residual correlation (`cor_res`) has no grammar equivalent yet – it
lives only in the frozen `create_phenotypes()`.

# Using Selected Markers to be QTN

As of the version `1.2.15`, it is also possible to select what markers
will be used as QTNS. The arguments `QTN_list` takes a lists of
additive, dominance or epistatic QTNs. In the example below, the SNP
`ss196523212` is selected to be the additive QTN. Another argument that
is exemplified below is `epi_interaction`. It defines the number of
markers to be involved in an epistatic interaction. The parameter
`out_geno` is defining that the marker data used in the simulation will
be exported as plink BED files.

``` r
QTN_list <- list()
QTN_list$add[[1]] <- c("ss196523212")
QTN_list$dom[[1]] <- c("ss196510214", "ss196472187")
QTN_list$epi[[1]] <- c("ss196530605", "ss196475446")
create_phenotypes(
  geno_obj = SNP55K_maize282_maf04,
  add_QTN_num = 1,
  dom_QTN_num = 2,
  epi_QTN_num = 1,
  epi_interaction = 2,
  h2 = c(0.92, 0.4) ,
  add_effect = c(0.90, 0.2),
  dom_effect = c(0.01, 0.3),
  epi_effect = c(-0.3, 0.7),
  ntraits = 2,
  QTN_list = QTN_list,
  rep = 1,
  output_format = "gemma",
  out_geno = "BED",
  output_dir = "output_data",
  model = "ADE",
  home_dir = getwd()
)
```

**V2 equivalent.** `QTN_list` becomes each layer’s `qtn =` argument (the
other layers stay random); `epistasis()` takes an
`n_pairs x interaction` matrix.

``` r
simulate_phenotype(SNP55K_maize282_maf04, n_traits = 2, seed = 1) |>
  additive(prop = 0.5, qtn = "ss196523212") |>
  dominance(prop = 0.2, qtn = c("ss196510214", "ss196472187")) |>
  epistasis(prop = 0.2, qtn = matrix(c("ss196530605", "ss196475446"), nrow = 1))
```

# Using Multiple Marker Data Files

If files are saved by chromosome, they can be read directly into
create_phenotypes using options `geno_path` (recommendation: consider
having all marker data files in a separate folder). If multiple files
are saved in the same folder as the marker data, the parameter `prefix`
might be used to select only the marker data. For example, if your data
is saved as “WGS_chrm_1.hmp.txt”, …, “WGS_chrm_10.hmp.txt”, one would
use `prefix = "WGS_chrm_"` .

``` r
create_phenotypes(
  geno_path = "PATH/TO/FILE",
  prefix = "WGS_chrm_",
  add_QTN_num = 3,
  h2 = 0.2,
  add_effect = 0.02,
  rep = 5,
  seed = 200,
  output_format = "gemma",
  output_dir = "Results",
  model = "ADE",
  home_dir = tempdir()
)
```

**V2 equivalent.** `as_numeric()` reads one file at a time, so convert
the per-chromosome files and stack them, then pipe the numeric object
into `simulate_phenotype()`. (The multi-file `geno_path` / `prefix`
reader still lives in the frozen `create_phenotypes()` if you prefer
it.)

``` r
files <- list.files("PATH/TO/FILE", pattern = "^WGS_chrm_", full.names = TRUE)
num   <- do.call(rbind, lapply(files, as_numeric, to_r = TRUE))
simulate_phenotype(num, h2 = 0.2, n_qtn = 3, seed = 200) |>
  additive(prop = 0.2)
```

For a full, side-by-side migration reference see the **v1-to-v2**
vignette (`vignette("v1-to-v2")`).
