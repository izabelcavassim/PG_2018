# PSMC


    utils/fq2psmcfa -q20 diploid.fq.gz > diploid.psmcfa
    psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o diploid.psmc diploid.psmcfa
    utils/psmc2history.pl diploid.psmc | utils/history2ms.pl > ms-cmd.sh
    utils/psmc_plot.pl diploid diploid.psmc


../Data/sorted_ERR1019039_reads_135_145.fq


software/bcftools-1.6/misc/vcfutils.pl vcf2fq -d 5 -D 34 -Q 30 --indv ERR1025617_chr2_piece_dedup data/Allvariants_135_145_chr2.vcf





vcftools --vcf data/Allvariants_135_145_chr2.vcf --indv ERR1025617_chr2_piece_dedup --recode --recode-INFO-all --out single_indiv



Finally, we run our consensus calling pipeline, consisting of a linked set of samtools, bcftools, and vcfutils.pl commands:

samtools mpileup -Q 30 -q 30 -u -v \
-f /home/Data/Homo_sapiens.GRCh37.75.dna.chromosome.2.fa -r chr2 P964.bam |  
bcftools call -c |  
vcfutils.pl vcf2fq -d 5 -D 34 -Q 30 > P964.$CHR.fq

This takes as input an aligned bam file and a reference genome, generates an mpileup using samtools, calls the consensus sequence with bcftools, and then filters and converts the consensus to fastq format, writing the results for each chromosome to a separate fastq file. Some parameter explanations:

samtools:
-Q and -q in mpileup determine the cutoffs for baseQ and mapQ, respectively
-v tells mpileup to produce vcf output, and -u says that should be uncompressed
-f is the reference fasta used (needs to be indexed)
-r is the region to call the mpileup for (in this case, a particular chromosome based on the array task id)
P964.bam is the bam file to use
bcftools:
call -c calls a consensus sequence from the mpileup using the original calling method
vcfutils.pl:
-d 5 and -d 34 determine the minimum and maximum coverage to allow for vcf2fq, anything outside that range is filtered
-Q 30 sets the root mean squared mapping quality minimum to 30