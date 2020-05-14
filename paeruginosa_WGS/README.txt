Whole population genome sequencing analysis for P. aeruginosa populations

Fastq's can be found in NCBI at the following BioprojectID:

SubmissionID:	SUB6668180 
BioProject ID:	PRJNA595915

The code for trimming reads using trimmomatic and calling variants with breseq can be found in trimmomatic_shellscript.sh and breseq_shellscript.sh.

The resulting breseq outputs for all of the populations was combined into csv files and utf characters were removed (pa_tob_snps.csv and pa_tob_nje.csv). SNPs can be viewed in pa_tob_snps and new junction evidence can be viewed in pa_tob_nje. 

Code for removing ancestral sequences and applying frequency filters can be viewed in mutation_filtering.R. Reshaping data for generating frequency tables (figures 2 and S2), allele frequency plots (figure S4) and muller plots (figure 3) can be viewed in frequency_tables_mullerplots.R.

day12_ab_gene.csv and day12.csv represent outputs from the frequency_tables_mullerplots.R script that were subsequently used to create figures 2 and S2, respectively, using code in abaumannii/WGS.
