# Phasing using Beagle



## Running Beagle

    java -jar Phasing_tutorial_files/beagle.08Jun17.d8b.jar gt=Phasing_tutorial_files/Allvariants_135_145_chr2.vcf map=Phasing_tutorial/plink.chr2.GRCh37.map out=Allvariants_135_145_chr2_phased

## Browsing phased results

    scp -P 8922 <your_name>@185.45.23.197:<vcf_file> <dir_on_your_computer>

Download the phased VCF file and open it in IGV. 

Explore phases of haplotypes at two positions in the alignment:

Select a base in a european individual close to position 130,000,000. Group alignments by base (right-click popup menu in alignment tracks). Answer these questions:

1. What does the haplotypes look like?
2. Which haplotypes agree?
3. How wide is the region where they agree?

Do the same for base 136,576,577 in a European and answer the same three  questions.


