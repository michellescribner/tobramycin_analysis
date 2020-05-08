Whole population genome sequencing analysis for A. baumannii populations


The code for trimming reads using trimmomatic and calling variants is analogous to trimmomatic_shellscript.sh and breseq_shellscript.sh in /paeruginosa_WGS

The resulting breseq outputs for all of the populations was combined into csv files and utf characters were removed (pa_tob_snps.csv and pa_tob_nje.csv). SNPs can be viewed in ab_tob_snps and new junction evidence can be viewed in ab_tob_nje. 

Code for removing ancestral sequences and applying frequency filters can be viewed in mutation_filtering.R. Reshaping data for generating frequency tables (figures 2 and S2), allele frequency plots (figure S4) and muller plots (figure 3) can be viewed in frequency_tables_mullerplots.R.