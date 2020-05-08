#P. aeruginosa clones selected for sequencing from tob m9 media populations

#create breseq output document on beagle
# ssh mrs186@beagle.mmg.pitt.edu
# cd /home/mrs186/tobim9/clones_190516/breseq
# /home/mrs186/scripts/BreseqCatEdited.py -d /home/mrs186/tobim9/clones_190516/breseq

#copy to computer
#scp mrs186@beagle.mmg.pitt.edu://home/mrs186/tobim9/clones_190516/breseq/Breseq_Output.xlsx /Users/mrs/Documents/PA14_Tobi_M9/clones

#save SNP tab as csv

#Convert to noutf characters 
system("iconv -c -f utf-8 -t ascii /Users/mrs/Documents/PA14_Tobi_M9/clones/Breseq_Output.csv > /Users/mrs/Documents/PA14_Tobi_M9/clones/breseq_output_noutf.csv")

# remove SNPs already present in the ancestor 
snps <- read.csv("/Users/mrs/Documents/PA14_Tobi_M9/clones/breseq_output_noutf.csv",header=TRUE)
ref_snps <- read.csv("/Users/mrs/Documents/PA14_Tobi_M63/Breseq_Output_noutf.csv",header=TRUE)
ancestor_snps <- subset(ref_snps,Sample == "78")

## Remove all common mutations based on the Position column
snps_noref <- snps[ !(snps$Position %in% ancestor_snps$Position), ]
nrow(snps_noref)#192

## Remove all common mutations based on the Gene column
snps_noref_gene <- snps[ !(snps$Gene %in% ancestor_snps$Gene), ]
nrow(snps_noref_gene)#37

## Create the filtered csv ##
write.csv(snps_noref_gene,file=("snps_noref.csv"))
