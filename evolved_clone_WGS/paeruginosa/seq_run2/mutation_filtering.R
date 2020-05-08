## Filtering Clones from PA14 Tobi M9 experiment

# the headers were all wonky, used Nate's script
# /home/mrs186/scripts/BreseqCat_clones.py -d /home/mrs186/tobim9/breseqclones/clones
# scp -r mrs186@beagle.mmg.pitt.edu://home/mrs186/tobim9/breseqclones/clones/Breseq_Output.xlsx /Users/mrs/Documents/PA14_Tobi_M9/breseq2
# save as csv

system("iconv -c -f utf-8 -t ascii /Users/mrs/Documents/PA14_Tobi_M9/breseq2/Breseq_Output_snps.csv > /Users/mrs/Documents/PA14_Tobi_M9/breseq2/breseq_output_snps_noutf.csv")
snps <- read.csv("/Users/mrs/Documents/PA14_Tobi_M9/breseq2/breseq_output_snps_noutf.csv",header=TRUE)

ref_snps <- read.csv("/Users/mrs/Documents/PA14_Tobi_M63/Breseq_Output_noutf.csv",header=TRUE)
ancestor_snps <- subset(ref_snps,Sample == "78")
ancestor_snps <- ancestor_snps[,-5]

## Remove all common mutations based on the Position column
snps_noref <- snps[ !(snps$Position %in% ancestor_snps$Position), ]

## Remove all common mutations based on the Gene column
snps_noref_gene <- snps[ !(snps$Gene %in% ancestor_snps$Gene), ]

## Create the filtered csv ##
write.csv(snps_noref_gene,file=("snps_noref.csv"))
