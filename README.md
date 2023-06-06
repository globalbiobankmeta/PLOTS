# PLOTS

Scripts to make plots and generate lists of variants with different freq compared to gnomAD and having strand ambiguity

This repository is used for QC of individual GWAS summary statistis files in Global Biobank Meta-analysis Initiative. 
[WDL](https://github.com/openwdl/wdl) workflows and Google Compute Engine are used for computing. The workflows consist of cleaning/munging input files to the same format and running a meta-analysis.

1. Use [./readme](./readme) and [./Dockerfile](./Dockerfile) to generate the docker image

2. The WDL workflow in [https://github.com/globalbiobankmeta/META_ANALYSIS/blob/master/wdl/munge_sumstats_beforeQC_obtain_QClist.wdl](https://github.com/globalbiobankmeta/META_ANALYSIS/blob/master/wdl/munge_sumstats_beforeQC_obtain_QClist.wdl) and [wdl/munge_sumstats_beforeQC_obtain_QClist.json](munge_sumstats_beforeQC_obtain_QClist.json) is used to compare the allele frequency of genetic variants in individual GWAS summary statisic files. 
