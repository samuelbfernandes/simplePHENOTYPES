### create.simulated.data usage
Parameters:
genotypes   =   NULL   a dataframe with genotypes
file.G   =   NULL   genotypic dataset file. (Used if genotype is not provided)
file.Ext.G   =  NULL   genotypic dataset file extension. (Used if genotype is not provided)

default parameters for numericalization (mostly coming from GAPIT)
    file.from = 1
    file.to = 1
    SNP.effect = "Add"
    SNP.impute = "Middle"
    Create.indicator = FALSE
    Major.allele.zero = FALSE
    file.fragment = Inf   # reads entire file at once using data.table.
    maf_cutoff = NULL   # if a MAF cutoff is to be used.


Additive.QTN.number   =   NULL 
Epistatic.QTN.number   =   NULL if not provided, only Additive effects are simulated
additive.effect   =   NULL  if multiple traits are being simulated, it may have lenght = ntraits
epistatic.effect   =   NULL  if multiple traits are being simulated, it may have lenght = ntraits
big.additive.QTN.effect   =   NULL  if multiple traits are being simulated, it may have lenght = ntraits
model   =   c("pleiotropic", "partially")   Choose one of the two genetic models for multiple traits
overlap   =   NULL   number of additive QTNs in common. Used for the "partially" (pleiotropic) model
overlapE   =   NULL   number of epistatic QTNs in common. Used for the "partially" (pleiotropic) model
specific.QTN.number   =   NULL   number of additive trait specific QTNs. Used for the "partially" (pleiotropic) model
specific.E.QTN.number   =   NULL   number of epistatic trait specific QTNs. Used for the "partially" (pleiotropic) model
rep   =   NULL   number of simulations for each genetic settings
ntraits   =   1   number of traits 
h2   =   NULL   heritabilities (if a vector, simulation results will be generated for each h2 value)
h2_MT   =   NULL   heritability of correlated traits (length = ntraits - 1)
set.cor   =   TRUE   control for trait correlation?
correlation   =   NULL   correlation matrix of dimension ntraits x ntraits. (Used if set.cor = TRUE)
seed   =   NULL   Seed used for sampling and random variable generator.
home.dir   =   getwd()   home directory
output.dir   =   NULL   Name of folder where files will be saved
format   =   "multi-file"   output format. (one of: "long", "multi-file" or "wide")
out.geno   =   FALSE   if input is a hapmap file, should numericalized hapmap be saved?

If set.cor = FALSE and model = "pleiotropic", additive.effect and big.additive.QTN.effect should have multiple values

---
##  Examples
perform numericalization
notice that the file name is SNP55K_maize282_AGPv2_20100513_1.hmp.txt, the "1" will be extracted from the parameters
file.from = 1 and file.to = 1, for the case of one file for each chromosome.
a file named SNP55K_maize282_AGPv2_20100513_NUM.txt will be saved at the home directory

create.simulated.data(file.G = "SNP55K_maize282_AGPv2_20100513_", file.Ext.G = "hmp.txt", out.geno = TRUE)

---
## simulation of 10 replicates of a univariate trait:
source("auxiliar_functions.R")
library(data.table)
library(lqmm)
library(mvtnorm)
Additive.QTN.number <- 5
heritabilities.vector <- c(0.9)
replicates <- 10

Size of the additive effect of the largest QTL (must be (-1,1) but preferably (0,1))
big.additive.QTN.effect <- c(0.9)
additive.effect <- c(0.1)
geno <- fread("SNP55K_maize282_AGPv2_20100513_NUM.txt", data.table = F)

create.simulated.data(
  genotypes = geno,
  output.dir = "My_Simulated_Data",
  Additive.QTN.number = Additive.QTN.number,
  additive.effect = additive.effect,
  big.additive.QTN.effect = big.additive.QTN.effect,
  rep = replicates,
  h2 = heritabilities.vector,
  seed = 2019
)

---
## simulation of 50 replicates of a multi-variate trait:
source("auxiliar_functions.R")
library(data.table)
library(lqmm)
library(mvtnorm)
Additive.QTN.number <- 25
heritabilities.vector <- c(0.9, 0.2)
replicates <- 50

Size of the additive effect of the largest QTL (must be (-1,1) but preferably (0,1))
big.additive.QTN.effect <- c(0.9, 0.7, 0.5, 0.8)
additive.effect <- c(0.1, 0.2, 0.1, 0.5)

h2_MT = c(0.5, 0.9, 0.3)

overlap = 10
specific.QTN.number <- c(14, 3, 10, 5)

cor4 <- matrix(c(   1,  0.3,  0.9, -0.4,
                    0.3,    1, 0.6,  0.34, 
                    0.9, 0.6,    1, -0.2,
                    -0.4,  0.34, -0.2,    1), 4)

geno <- fread("SNP55K_maize282_AGPv2_20100513_NUM.txt", data.table = F)

create.simulated.data(
  genotypes = geno,
  output.dir = "My_Simulated_Data",
  additive.effect = additive.effect,
  big.additive.QTN.effect = big.additive.QTN.effect,
  set.cor = T,
  ntraits = 4,
  correlation = cor4,
  h2_MT = h2_MT,
  rep = replicates,
  h2 = heritabilities.vector,
  seed = 2019,
  format = "multi-file",
  model = "partially",
  overlap = overlap,
  specific.QTN.number = specific.QTN.number
)

