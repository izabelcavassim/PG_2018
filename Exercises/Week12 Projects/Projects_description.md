# Analysis of GWAS summary statistics

## Summary statistics

To protect the privacy of the participants in GWA studies the raw genotype data is usually not made public. It has, however, become the norm that studies will publish the summary statistics (p-value, effect size, frequency etc.) for all the variants they have tested. For this reason there are many new methods designed to do further analyses based on summary statistics. In this project you will analyse summary statistcs from a large meta-analysis of BMI and:
- Do a per gene test and look for enriched gene sets.
- Estimate the heritability of BMI.

## Data
You can download summary statistics from a BMI study with ~700,000 individuals [here](http://portals.broadinstitute.org/collaboration/giant/images/6/63/Meta-analysis_Wood_et_al%2BUKBiobank_2018.txt.gz).

## Per gene test.
The gcta tool that you used in the last GWAS exercise can calculate a per gene test based on summary statistics (see [here](http://gcta.freeforums.net/thread/309/gcta-fastbat-based-association-analysis)). In order to see if the presence of multiple significant variants in the same gene is due to LD or multiple independent signals it is necessary to provide a data set in plink-format that the program can use to estimate the LD between variants. For this purpose you can use data from the 1000 genomes project that can be downloaded in plink format [here](https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz). The BMI study and LD data is based on hg19 so you should use the gene list file called "glist-hg19.txt". The genotype data consists of separate files for each chromosome so you can do a test per chromosome and then concatenate the results in the end.

### Gene set enrichment
Try to see if the genes affecting BMI are enriched in specific Gene Ontologies or Pathways or Tissue types. You can for example use the tool enrichR: http://amp.pharm.mssm.edu/Enrichr/

## Estimating heritability
The method called LD-score regression can be used to estimate heritability using summary statistics. The method is described in [this article](https://www.nature.com/articles/ng.3211). A description of how to use the software to estimate heritability can be found [here](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation). The software can be downloaded [here](https://github.com/bulik/ldsc).
