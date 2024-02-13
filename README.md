This pipeline is for the procesing of the RNAseq data from "Translation-efficiency-driven-by-CNOT3-subunit-of-the-CCR4-NOT-complex-promotes-leukemogenesis".

Reads are processed by both Shell and R scripts. The Shell scripts are run first according to their number (which account for adaptor removal and alignment to the transcriptome) followed by the R scripts (DEseq2, feature analysis of mRNAs and gradient boosting)

In each directory, there is a common variables script that must be modified by the user

