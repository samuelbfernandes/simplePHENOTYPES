# simplePHENOTYPES 1.3.0
## Major changes
Implemented the parameter "ld_max" (replacing "ld") and "ld_min".
## Minor changes
Fixed bug that changed the working directory after the simulation

# simplePHENOTYPES 1.2.16
## Major changes
Included the parameter 'mean', so traits can be simulated with the desired mean (intercept) value.
Included QTN_list option for the LD architecture. 
set default seed generator as RNGversion('3.5.1') to ensure reproducibility.
## Minor changes
Renamed some output QTN info files to make it standard across different architectures

# simplePHENOTYPES 1.2.15
## Major changes
Included the parameter QTN_list = list(add = NULL, dom = NULL, epi = NULL) to give the user the possibility to select the specific markers to be used as QTNs.
## Minor changes
Included 'Master Seed' in the log file to facilitate reproducibility. Now it only saves individual seed numbers when verbose = TRUE (default).

# simplePHENOTYPES 1.2.14
## Minor changes
check if 'out_geno' is either 'numeric', 'plink' or 'gds'
replaced the dependence lqmm and uses the function make_pd() to make cor matrix positive definite

# simplePHENOTYPES 1.2.13
## Minor changes
Fixed bug when reading multiple files using geno_path

# simplePHENOTYPES 1.2.12
## Minor changes
Fix bug that stopped simplePHENOTYPES when using geno_obj and architecture = "LD"

# simplePHENOTYPES 1.2.11
## Minor changes
Fix bug in the simulation of single trait using multiple h2 values 

# simplePHENOTYPES 1.2.10
## Minor changes
Fix bug that made the direct LD option stop running


# simplePHENOTYPES 1.2.9
## Minor changes
Fixed bug that also removed the cause of LD when remove_QTN = TRUE with architecture = "LD"
Fixed bug when reading multiple files using geno_path

# simplePHENOTYPES 1.2.8
## Minor changes
Set all additive parameters to NULL when model is dominance or epistasis.

# simplePHENOTYPES 1.2.7
## Minor changes
Fixed bug when more than 9 traits were simulated under the "partially" architecture
Fixed bug when saving file name with very large name due to a large number of traits

# simplePHENOTYPES 1.2.6
## Minor changes
Fixed bug in the QTN MAF calculation on the LD architecture
Fixed bug when importing VCF and exporting BED files (implemented by chr_prefix)


# simplePHENOTYPES 1.2.4
## Major changes
**Input**
1. Implemented options for input format as VCF, plink bed/ped files, GDS.
1. Changed dosage (numeric format) information from 0, 1, and 2 to -1 (aa), 0 (Aa) and 1 (AA).
1. Implemented a new type of spurious pleiotropy, direct LD (type\_of\_ld = "indirect").
1. Included the option for assigning a residual correlation among traits.
1. Implemented a constrain option to select only heterozygote or only homozygote QTNs.
1. Included the warning\_file\_saver option to skip asking if the user wants to save one genotype file for each rep when vary\_QTN = FALSE.

**Output**
1. Included a new output file with the summary linkage disequilibrium information on the selected spurious pleiotropy QTNs.
1. Included MAF in the outputted QTN information file.
1. Calculates the proportion of phenotypic variation explained by each QTN (QTN\_variance = TRUE).
1. Includes the option to remove QTNs from the genotype file (remove_QTN = TRUE).
1. Renamed <Taxa> by <Trait> in Tassel output format.


## Minor changes

Fixed bug that didn't recognize geno\_obj as HapMap.
Fixed bug when simulating dominance will all SNPs being homozygotes.
Fixed bug when reaching the end of the file while looking for SNPs in LD.
Fixed bug in importing geno\_file from other directories.
Fixed bug in selecting QTNs when marker data < 6 SNPs.
Included file removal when simulation does not complete.
Renamed file outputted as numeric.
Renamed constrain option.
Implemented check for biallelic markers.
Changed the QTN file name.
Implemented an interactive question before Check to remove QTNs with vary\_QTN = T.
Included check.names as FALSE in all data.frames.
Check if geno\_file and geno\_path are NULL.
Incorrect output name.

