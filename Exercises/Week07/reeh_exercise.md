Background
----------

In this exercise section you will analyse positive selection in human populations. You will be looking at the X chromosome of individuals from the Simons genome diversity project. The advantage of X chromosomes in males is that they are haploid, and therefore, fully phased. We also have other reasons to belive that they often experience natural selection. Finally, they have typically not been investigated previously to the same extent as the autosomes.

You will perform many of the same analysis as in the [Sabeti et al. (2007)](https://www.nature.com/articles/nature06250) paper. They used these to find selection in the human HapMap data (SNP data) on the autosomes.

You will have data from the following populations:

| Population  | Individuals |
|-------------|:-----------:|
| Africa      |   26 males  |
| West Europe |   44 males  |
| South Asia  |   30 males  |
| East Asia   |   24 males  |

The data consists of **24198** SNPs from the region 73-81 Mb on the X chromosome and there is no missing data. The haplotype data for each population is found in separate files (**genotypes360\_400\_.**), whereas they use a common SNP identity file **snps360\_400\_filtered.snp**.

Package
-------

The package requires 2 inputs:

1.  The snp file (`snp_ID`, `Chromosome`, `position`, `reference allele` and `derived allele`).

2.  The haplotype file: haplotype of each individual encoded as reference allele, alternative allele and missing data N.

Analysis
--------

You will perform a genome wide scan and then focus on candidate SNPs. The package that you will be using on these analysis is `rehh`. And the vignette of the package can be found at blackboard or [here](https://cran.r-project.org/web/packages/rehh/vignettes/rehh.pdf). Please have a look at the document and get familiarized with the functions of the package.

### Reading in data in REHH format

``` r
# Install the package
#install.packages('rehh')
library(rehh)

# Define the directory of your working folder:
setwd("~/Dropbox/PG2018/exercises/rehh")

# Reading the data for each population:
hap360_400_AF <-data2haplohh(hap_file="genotypes360_400_AF",map_file="snps360_400_filtered",
                             recode.allele=TRUE, 
                             min_perc_geno.snp=100,
                             min_perc_geno.hap=80,
                             haplotype.in.columns=TRUE,
                             chr.name=1)
```

    ## Map file seems OK: 24198  SNPs declared for chromosome 1 
    ## Haplotype are in columns with no header
    ## Alleles are being recoded according to map file as:
    ##  0 (missing data), 1 (ancestral allele) or 2 (derived allele)
    ## Discard Haplotype with less than  80 % of genotyped SNPs
    ## No haplotype discarded
    ## Discard SNPs genotyped on less than  100 % of haplotypes
    ## No SNP discarded
    ## Data consists of 26 haplotypes and 24198 SNPs

    ## Map file seems OK: 24198  SNPs declared for chromosome 1 
    ## Haplotype are in columns with no header
    ## Alleles are being recoded according to map file as:
    ##  0 (missing data), 1 (ancestral allele) or 2 (derived allele)
    ## Discard Haplotype with less than  80 % of genotyped SNPs
    ## No haplotype discarded
    ## Discard SNPs genotyped on less than  100 % of haplotypes
    ## No SNP discarded
    ## Data consists of 30 haplotypes and 24198 SNPs

    ## Map file seems OK: 24198  SNPs declared for chromosome 1 
    ## Haplotype are in columns with no header
    ## Alleles are being recoded according to map file as:
    ##  0 (missing data), 1 (ancestral allele) or 2 (derived allele)
    ## Discard Haplotype with less than  80 % of genotyped SNPs
    ## No haplotype discarded
    ## Discard SNPs genotyped on less than  100 % of haplotypes
    ## No SNP discarded
    ## Data consists of 44 haplotypes and 24198 SNPs

    ## Map file seems OK: 24198  SNPs declared for chromosome 1 
    ## Haplotype are in columns with no header
    ## Alleles are being recoded according to map file as:
    ##  0 (missing data), 1 (ancestral allele) or 2 (derived allele)
    ## Discard Haplotype with less than  80 % of genotyped SNPs
    ## No haplotype discarded
    ## Discard SNPs genotyped on less than  100 % of haplotypes
    ## No SNP discarded
    ## Data consists of 24 haplotypes and 24198 SNPs

#### Q1. How many haplotypes and snps are found in each population?

Scan the region using iHS, Rsb and XP-EHH
-----------------------------------------

You should first perform a scan for each of the regions for extreme values of **i**ntegrated **h**aplotype **s**core (iHS), standardized to difference in allele frequency. This implies using the functions **scan\_hh**, followed by \*\*ihh2ihs and then plotting the results using either ihsplot or plotting by yourself

#### Q2. Try to also plot a histogram of the allelefrequncies of the SNPs in each population (part of the dataframe resulted from scan\_hh).

#### Do you find outliers with significant iHS? Then, record the SNP positions of the most significant SNPs for later analysis, using e.g. which.max() or which.min().

``` r
res.scanAF<-scan_hh(hap360_400_AF)
res.scanSA<-scan_hh(hap360_400_SA)
res.scanWE<-scan_hh(hap360_400_WE)
res.scanEA<-scan_hh(hap360_400_EA)

wg.ihsAF<-ihh2ihs(res.scanAF, freqbin = 0.05) 
ihsplot(wg.ihsAF, plot.pval = TRUE)

wg.ihsSA<-ihh2ihs(res.scanSA, freqbin = 0.05) 
ihsplot(wg.ihsSA, plot.pval = TRUE, ylim.scan = 2, main = "iHS South Asia")

wg.ihsWE<-ihh2ihs(res.scanWE, freqbin = 0.05) 
ihsplot(wg.ihsWE, plot.pval = TRUE)

wg.ihsEA<-ihh2ihs(res.scanEA, freqbin = 0.05) 
ihsplot(wg.ihsEA, plot.pval = TRUE)
```

Perform pairwise population tests
---------------------------------

There are two possible functions to be used, both require the dataframes of 'scan\_hh' class.

``` r
wg.rsbAFWE <- ies2rsb(res.scanAF,res.scanWE, popname1 = "Africa", popname2 = "W Europe", method = "bilateral")
rsbplot(wg.rsbAFWE, plot.pval = T)
wg.XPEHHAFWE <- ies2xpehh(res.scanAF,res.scanWE, popname1 = "Africa", popname2 = "W Europe", method = "bilateral")
xpehhplot(wg.XPEHHAFWE, plot.pval = T)
wg.XPEHHAFSA <- ies2xpehh(res.scanAF,res.scanSA, popname1 = "Africa", popname2 = "South Asia", method = "bilateral")

wg.XPEHHAFEA <- ies2xpehh(res.scanAF,res.scanEA, popname1 = "Africa", popname2 = "East Asia", method = "bilateral")
wg.XPEHHWESA <- ies2xpehh(res.scanWE,res.scanSA, popname1 = "WestEurope", popname2 = "South Asia", method = "bilateral")
```

Zooming in on interesting markers
---------------------------------

From the scan you can find SNPs that give extreme values of iEHS or of XPEHH for a set of populations. You can then analyse the haplotype structure around them. This is done by including the index position of the interested marker in the functions **calc\_ehhs** and **bifurcation.diagram**.

``` r
#For African Western Europe the most significant marker is 4901

calc_ehhs(hap360_400_WE, mrk=4901)
calc_ehh(hap360_400_WE,mrk = 4901)

calc_ehhs(hap360_400_AF, mrk=4901)
calc_ehh(hap360_400_AF,mrk = 4901)

layout(matrix(1:2,2,1))
bifurcation.diagram(hap360_400_WE,mrk_foc=4901,all_foc=1,nmrk_l=200,nmrk_r=200, refsize = 0.8,
main="Candidate X: Ancestral Allele")
bifurcation.diagram(hap360_400_AF,mrk_foc=4901,all_foc=1,nmrk_l=200,nmrk_r=200, refsize = 0.8,
main="Candidate X: Ancestral Allele")
```

### What is the biological function of the region around this snp?

Have a look at NCBI, remember that this dataset belongs to chromosome X and HG19 as reference.
