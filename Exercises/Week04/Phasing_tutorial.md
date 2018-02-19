# Phasing

The base calling exercise you identified the two bases at each position in each individual. However, you do not know which base goes on which of the two chromosomes. That means that you do not know if....

    -----A-----C------         -----T-----C------ 
    -----T-----G------    or   -----A-----G------


Jointly called SNPs in `/home/data/Allvariants_135_145_chr2.vcf`.


## Running Beagle

For additional information see the [Beagle 4.1 manual](https://faculty.washington.edu/browning/beagle/beagle_4.1_03Oct15.pdf)

    java -jar /home/shared/beagle.08Jun17.d8b.jar gt=/home/Data/Allvariants_135_145_chr2.vcf map=/home/shared/data/plink.chr2.GRCh37.map out=Allvariants_135_145_chr2_phased

## Browsing phased results

Download the phased VCF file:

    scp -P 8922 user_name@185.45.23.197:Allvariants_135_145_chr2_phased.vcf.gz dir_on_your_computer

Open the file in IGV (integrative genomics viewer): 
    
1. Choose Human hg19 as the reference genome.
2. Click `File > Load from File...` and select you phased VCF file.

Explore phases of haplotypes at two positions in the alignment:

Zoom all the way in and select a base in a European individual (ERR1025620 is English) close to position 135,000,000. Group alignments by genotype (right-click on the base in the alignment tracks to get a popup menu). Answer these questions:

1. What does the haplotypes look like?
2. Which haplotypes agree?
3. How wide is the region where they agree?

Now do the same for base 136608646 in a the same European individual and answer the same three questions.

What are the differences in the haplotype structure in alignment around the two different positions? What do you think produces these differences?

The SNP at position 136608646 is called 136608646. Try to search for 136608646 in the [UCSC genome browser](https://genome-euro.ucsc.edu/cgi-bin/hgGateway?redirect=manual&source=genome.ucsc.edu). Can you find anything that explains your observations?

# Estimating a recombination map

[LDhat  manual](https://github.com/auton1/LDhat/blob/master/manual.pdf)

## Format input data for LDhat

LDhat needs its input data in a particular format. We will use vcftools to produce these input files from the phased VCF file:

    vcftools --gzvcf Allvariants_135_145_chr2_phased.vcf.gz --chr 2 --ldhat --out recmap_data

Have a look at the two files produced using `less`. E.g.:

    less recmap_data.ldhat.sites 

> NB: press `q` to quit `less`

How do you think the information is encoded in these files?

## Running LDhat

To speed up computations you can make a lookup table first. That takes a while, so I did if for you. But it is done like this using the `complete` program that comes with LDhat:

    /usr/local/bin/complete -n 56 -rhomax 100 -n_pts 101 -theta 0.0001

- `-n 56`:the number of haplotypes (2 * 28).
- `-rhomax 100`: maximum rho ($4N_e r$) alowed.
- `-n_pts 101`: number of points in grid: 101 (recommended)
- `-theta 0.0001`: human theta ($4N_e \mu$).

That produces a file named `new_lk.txt`. The next step is to calculate the recombination map. It will take around 8 minutes for the entire dataset.

    rhomap -seq recmap_data.ldhat.sites -loc recmap_data.ldhat.locs -lk /home/shared/data/new_lk.txt -its 200000 -samp 2000 

- `-lk`: likelihood lookup table.
- `-its`: number of iterations of the MCMC chain.
- `-samp`: how often to sample from the MCMC chain.
- `-burn`: how many of the initial iterations to discard (should be at least 100,000 divided by the value of `-samp`).

Quoting the manual:

> Four files are generated as the Markov Chain progresses. `acceptance_rates.txt` details the acceptance rates. `hotspots.txt` contains details of the hotspots at each sample from the chain. 
When rhomap completes it writes three files:

- `acceptance_rates.txt`: acceptance rates of the MCMC. If they are lower than 1%. The program should be run with more iterations.
- `summary.txt`: summary of the recombination rates estimated.
- `rates.txt`: (quoting the manual) is the output from each sample detailing the recombination rate (expressed in $4N_e r$ per kb) between each SNP. summary.txt contains a summary of the samples from the chain detailing, for each SNP interval, the estimated genetic map position, the estimated recombination rate, and the hotspot density (the number of hotspots per kb per iteration). The rates.txt file can be summarized by use of the program stat.

## Analyze results

Open Rstudio from the directory where your output files are. Do that by navigating into the right directory and then type `rstudio` in the terminal.

> NB: To make the rstudio window pop up on your screen you need to use the `-Y` option when you log in to the server: `ssh -Y etc...`

Now paste this into Rstudio console:

```R
source("http://ldhat.sourceforge.net/R/coalescent.r")
```

That loads a lot of R functions written by the author of LDhat.

Now run this code:

```R
summary <- summarise.rhomap()
```

If you renamed any of the output files you can specify the new names like this:

```R
summary<-summarise.rhomap(rates.file = "<file_name>", burn.in = <float>, locs.file="<file_name>")
```

The summary produces two plots:

- A graph of the recombination rate across the sequence, along with confidence intervals.
- A plot showing how estimation of recombination rate has progressed with each MCMC sample. Notice that the initial run of MCMC samples are atypical. This is the "burn-in" of the MCMC. We want to remove that, so take notice of how many samples it corresponds to. If it is 50 they we can produce a new set of estimates that excludes this burn-in using the `stat` program that comes with LDhat:

		/usr/loca/bin/stat -input rates.txt -loc recmap_data.ldhat.locs -burn 50

This produces a file called `res.txt` that describes the confidence in the estimated recombination rate along the sequence.

Now try to plot the final results:

```R
library(ggplot2)
library(dplyr)
library(tidyr)
library(magrittr)

rates <- read.table('res.txt', header = T)
rates %>%
    # filter(Loci > 140000.000, Loci < 142000.000)
    ggplot(aes(x=Loci, y=Mean_rho, ymin=L95, ymax=U95)) +  
        geom_line(color='blue') +
        geom_ribbon(alpha=0.1) +
        theme_bw()
```

Look at the plots and ponder the following questions:

- Are there any recombination hotspots?
- Are there any regions where the estimated recombination rate is really low? If so why do you think that could be?


