# simplePHENOTYPES 1.2.4
## Major changes
**Input**
1. Implemented options for input format as VCF, plink bed/ped files, GDS.
1. Changed dosage (numeric format) information from 0, 1, and 2 to -1 (aa), 0 (Aa) and 1 (AA).
1. Implemented a now type of spurious pleiotropy, direct LD (type\_of\_ld = "indirect").
1. Included the option for assigning a residual correlation among traits.
1. Implemented a constrain option to select only heterozygote or only homozygote QTNs.
1. Included the warning\_file\_saver option to skip asking if the user wants to save one genotype file for each rep when vary\_QTN = FALSE.

**Output**
1. Included a new output file with the summary linkage disequilibrium information on the selected spurious pleiotropy QTNs.
1. Included MAF in the outputted QTN information file.
1. Calculates the proportion of phenotypic variation explained by each QTN (QTN\_variance = TRUE).
1. Includes the option to remove QTNs from genotype file (remove_QTN = TRUE).
1. Renamed <Taxa> by <Trait> in Tassel output format.


## Minor changes

Bug fix when geno\_obj is hapmap.
Bug fix when model includes dominance and there's no hets.
Bug fix when reaching the end of the file while looking for SNPs in LD.
Bug fix with residual correlation and h2=1.
Bug fix in importing geno_file from other folders.
Bug fix in numeric input format.
Bug fix in selecting QTNs when marker data < 6 SNPs.
Remove files when simulation does not complete.
Renamed file outputted as numeric.
Input missing data from.
Renamed constrain option.
Renamed gds files.
Check if QTNs are biallelic.
Changed QTN file name.
Check to remove QTNs with vary\_QTN = T.
Dont check names in data.frame.
Check if geno\_file and geno\_path are NULL.
Incorrect output name.

