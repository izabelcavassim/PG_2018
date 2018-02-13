# Phasing

The base calling exercise you identified the two bases at each position in each individual. However, you do not know which base goes on which of the two chromosomes. That means that you do not know if....

    -----A-----C------         -----T-----C------ 
    -----T-----G------    or   -----A-----G------

## Running Beagle

    java -jar /home/shared/beagle.08Jun17.d8b.jar gt=/home/data/Allvariants_135_145_chr2.vcf map=/home/data/plink.chr2.GRCh37.map out=Allvariants_135_145_chr2_phased

## Browsing phased results

    scp -P 8922 user_name@185.45.23.197:Allvariants_135_145_chr2_phased.vcf.gz dir_on_your_computer

Download the phased VCF file and open it in IGV. 

Explore phases of haplotypes at two positions in the alignment:

Select a base in a european individual close to position 130,000,000. Group alignments by base (right-click popup menu in alignment tracks). Answer these questions:

1. What does the haplotypes look like?
2. Which haplotypes agree?
3. How wide is the region where they agree?


Lactase persistence causal sap: rs4988235 at chr2:136608646-136608646

Do the same for base 136608646 in a European and answer the same three  questions.

Where are the differences in haplotype structure at the two different positions? Why do you think that is?


# Estimating a recombination map

[LDhat  manual](https://github.com/auton1/LDhat/blob/master/manual.pdf)

## Format input data for LDhat

vcftools --gzvcf Allvariants_135_145_chr2_phased.vcf.gz --chr 2 --ldhat --out recmap_data

## Running LDhat

To speed up computations you can make a lookup table first. That takes a while, so I did if for you. But it is done like this using the `complete` program that comes with LDhat:

    /usr/local/bin/complete -n 56 -rhomax 100 -n_pts 101 -theta 0.0001

- number of haplotypes: 56
- max rho: 100
- number of points in grid: 101 (recommended)
- human theta: 0.0001

That produces a file named `new_lk.txt`.

    rhomap -seq recmap_data.ldhat.sites -loc ldhat_input.ldhat.locs -lk new_lk.txt -its 60000000 -samp 40000 -burn 0


## Summarize results

    /usr/loca/bin/stat -input rates.txt -loc ldhat_input.ldhat.locs -burn 500

## Visualize results

Open Rstudio from the directory where your output files are. Do that by navigating into the right directory and then type `rstudio` in the terminal.

> NB: To make the rstudio window pop up on your screen you need to use the `-Y` option when you log in to the server: `ssh -Y etc...`

Now paste this into Rstudio console:

```R
source("http://ldhat.sourceforge.net/R/coalescent.r")
```

That loads a lot of R functions written by the author of LDhat.

Now run this code:

```R
summary<-rhomap(rates.file = "<file_name>", burn.in = <float>, locs.file="<file_name>")
```





