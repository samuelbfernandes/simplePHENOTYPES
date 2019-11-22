
<!-- README.md is generated from README.Rmd. Please edit that file -->

# simplePHENOTYPES

<!-- badges: start -->

<!-- badges: end -->

simplePHENOTYPES aims to make it simple to simulate traits under
pleiotropy, partial-pleiotropy and linkage disequilibrium.

This short tutorial presents some of the possible genetic architectures
that one could simulate, but it certainly does not explore all the
possibilities. For more information on specific input parameters, please
check the help documentation (?create\_phenotypes).

# Installation

In order to install simplePHENOTYPES, the following r packages will also
be installed:

  - From Bioconductor:
      - SNPRelate
      - gdsfmt
  - From CRAN:
      - mvtnorm  
      - lqmm  
      - data.table

<!-- end list -->

``` r
setRepositories(ind = 1:2)
devtools::install_github("samuelbfernandes/simplePHENOTYPES")
```

# Load sample dataset

Note that the example dataset is already numericalized. HapMap format
may also be used with options `geno_obj`, `genotypes_file` or
`geno_path`.

``` r
library(simplePHENOTYPES)
data("SNP55K_maize282_maf04")
SNP55K_maize282_maf04[1:8, 1:10]
#>            snp allele chr     pos cm 4226 4722 33-16 38-11 A188
#> 1: ss196422159    A/G   1  379844 NA    0    0     0     0    0
#> 2: ss196422171    G/A   1  613257 NA    0    0     0     0    0
#> 3: ss196422173    G/A   1  659354 NA    2    2     0     0    0
#> 4: ss196422186    G/C   1  992572 NA    2    0     0     0    0
#> 5: ss196500940    G/A   1 2044264 NA    2    2     2     2    2
#> 6: ss196422234    A/G   1 2044555 NA    0    0     0     0    0
#> 7: ss196422236    G/A   1 2045035 NA    2    2     2     2    2
#> 8: ss196422248    G/A   1 2428671 NA    0    2     2     0    2
```

# Single Trait

Simulating ten replicates of a single trait with a heritability of 0.7.
In this setting, the simulated trait is controlled by one large-effect
QTN (additive effect of 0.9) and two small effect QTNs with additive
effects following a geometric series starting with 0.2; thus the effect
size of the first of these two QTNs is 0.2, and the effect size of the
second is 0.2^2. Results are being saved at a temporary directory
(`home_dir = tempdir()`). Please see help files (?create\_phenotypes) to
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
  home_dir = getwd())
```

# Multiple Traits: Pleiotropic Architecture

Simulating three traits `ntraits = 3` controlled by the same three
additive QTN and four dominance QTN. The effect size of the
largest-effect additive QTN is 0.3 for all traits, while the additive
effect sizes are 0.04, 0.2 and 0.1 for each trait, respectively.
Heritability for trait\_1 is 0.2, while the heritability of the two
correlated traits is 0.4. Each replicate is being recorded in a
different file (`output_format = "multi-file"`) in a folder named
“Results\_Pleiotropic”. The correlation between traits is not
specified by the user; instead, the observed correlation is an artifact
of different allelic effects for each trait. The same QTNs are used to
generate phenotypes in all 10 replications (`vary_QTN = FALSE`);
alternatively, one could select different QTNs in each replicate using
`vary_QTN = TRUE` (default). The vector `add_effect` contains one
allelic effect for each trait and a geometric series (default) will be
used to generate allelic effects for each one of the three additive QTNs
(`add_QTN_num = 3`) and four dominance QTNs (`dom_QTN_num = 4`).

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
    home_dir = getwd()
  )
```

Optionally, users may input a list of allelic effects (`sim_method =
"custom"`). In the example below, a geometric series (custom\_geometric)
is being assigned and should generate the same simulated data as the
previous example (all.equal(test1, test2)). Notice that since
`big_add_QTN_effect` is provided, the user only needs to provide effects
for two out of the three additive QTNs being simulated. On the other
hand, all four dominance QTN must have an effect assigned on the
custom\_geometric\_d list.

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
   home_dir = getwd()
 )
 
 all.equal(test1, test2)
```

# Multiple Traits: Partially Pleiotropic Architecture

Simulating 20 replicates of three traits, which are respectively
controlled by seven, 13 and four QTNs. The first of the three
pleiotropic QTNs (`pleio_a = 3`) will have a large additive QTN effect
of 0.9. All other QTNs will have additive effects that follow a
geometric series, where the effect size of the i^th QTN has an effect
size of 0.3^i. Correlation among traits is assigned to be equal to the
cor\_matrix object. All 20 replicates of these three simulated traits
will be saved in one file, specifically in a long format and with an
additional column named “Rep”. Results will be saved in a directory
called “Results\_Partially”. In this example, the genotype file will be
saved in numeric format. Additionally, all phenotypes will be assigned
to an object called “sim\_results” (to\_r = TRUE) and loaded into R.

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
  cor = cor_matrix,
  rep = 20,
  output_dir = "Results_Partially",
  output_format = "long",
  architecture = "partially",
  out_geno = "numeric",
  to_r = TRUE,
  model = "AED",
  home_dir = getwd()
)
```

# Multiple Traits: Linkage Disequilibrium Architecture

Simulating five replicates of two linked traits controlled by three
additive QTNs each. For each QTN, a marker is first selected, and then
two separate markers that have an r^2 of at most 0.8 (`ld=0.8`) with the
selected marker will be the respective QTNs of the two traits. The three
QTNs will have additive effects that follow a geometric series, where
the effect size of the i^th QTN is 0.02^i for one trait and 0.05^i for
the other trait. Starting seed number is 200 and output phenotypes are
saved in one file, but in a “wide” format with each replicate of two
traits being added as additional columns. Plink fam, bim and bed files
are also saved at Results\_LD.

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
  out_geno = "plink",
  ld=0.8,
  model = "A",
  home_dir = getwd()
)
```

# Multiple Traits: Partially Pleiotropic Architecture with Epistatic effects

Simulating 20 replicates of three traits. In each replicate, different
SNPs are selected to be the QTNs for each trait. These traits are
controlled by three pleiotropic additive QTNs (`pleio = 3`); two
pleiotropic epistatic QTNs (`pleio_e = 2`); four, ten and one additive
trait-specific QTN (`specific_QTN_number = c(4,10, 1)`); and two, one
and five epistatic trait-specific QTNs (`trait_specific_e_QTN_number =
c(2,1, 5)`). Results will be saved as “.fam” files used as gemma input;
thus each replicate trait that is simulated will have its own unique
output file.

``` r
create_phenotypes(
  geno_obj = SNP55K_maize282_maf04,
  pleio_a = 3,
  pleio_e = 2,
  trait_spec_a_QTN_num = c(4,10, 1),
  trait_spec_e_QTN_num = c(2,1, 5),
  epi_effect = c(0.01, 0.4, 0.2),
  add_effect = c(0.3, 0.2, 0.5),
  h2 = c(0.2,0.4, 0.8),
  ntraits = 3,
  rep = 20,
  vary_QTN = TRUE,
  output_dir = "Results_Partially_E",
  output_format = "gemma",
  architecture = "partially",
  model = "AE",
  home_dir = getwd()
)
```

# Using your own data

If files are saved by chromosome, they can be read directly into
create\_phenotypes using options `geno_path` (recommendation: consider
having all marker data files in a separate folder). If multiple files
are saved in the same folder as the marker data, the parameter `prefix`
might be used to select only marker data. For example, if your data is
saved as “WGS\_chrm\_1.hmp.txt”, …, “WGS\_chrm\_10.hmp.txt”, one would
use `prefix = "WGS_chrm_"` .

``` r
create_phenotypes(
  geno_path = getwd(),
  prefix = "WGS_chrm_",
  add_QTN_num = 3,
  h2 = 0.2,
  add_effect = 0.02,
  rep = 5,
  seed = 200,
  output_format = "gemma",
  output_dir = "Results",
  model = "ADE",
  home_dir = getwd()
)
```

# Contact

Questions, suggestions, and bug reports are welcome and appreciated.

Author: Samuel B Fernandes and Alexander E Lipka

Contact: <samuelf@illinois.edu> or <fernandessb101@gmail.com>

Institution: University of Illinois at Urbana-Champaign
