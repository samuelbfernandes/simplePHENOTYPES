setwd("~/Documents/UIUC2/multitrait-sim")
library(simplePHENOTYPES)
data(geno)
#----- Model 1 same QTNs ----
Additive.QTN.number <- 25

heritabilities.vector <- c(0.9, 0.2)
replicates <- 50

#Size of the additive effect of the largest QTL (must be (-1,1) but preferably (0,1))
big.additive.QTN.effect <- c(0.9, 0.7, 0.5, 0.8)
additive.effect <- c(0.1, 0.2, 0.1, 0.5)

#---- required for multi-traits ----
# ntraits (number of traits)
# h2_MT (heritability of each correlated trait)
# if set.cor = TRUE:
# correlation = correlation matrix
# if set.cor = FALSE set multiple values for the following objects:
# additive.effect
# big.additive.QTN.effect

h2_MT = c(0.5, 0.9, 0.3)
#---- model 2 only some QTNs in common -----
overlap = 10
specific.QTN.number <- c(14, 3, 10, 5)

cor4 <- matrix(c(   1,  0.3,  0.9, -0.4,
                  0.3,    1, 0.6,  0.34,
                  0.9, 0.6,    1, -0.2,
                 -0.4,  0.34, -0.2,    1), 4)


create.simulated.data(
  genotypes = geno,
  output.dir = "Test",
  Additive.QTN.number = Additive.QTN.number,
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

## include
old <- setwd(tempdir())
on.exit(setwd(old), add = TRUE)
