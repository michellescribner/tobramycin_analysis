#Analysis of WGS reads for clones isolated from abaum populations in Scribner et al

#1 clone each was sequenced from biofilm populations 1 and 3 on days 4,7, and 12
#245 is pop 1 day 4, 246 is pop 1 day 7, 247 is pop1 day 12
#248 is pop 3 day 4, 249 is pop 3 day 7, 250 is pop3 day12

setwd("/Users/mrs/Documents/PA14_Tobi_M9/abaumtob")

#Trim and filter reads
#for i in 245 246 247 248 249 250 ; do trimmomatic PE -phred33 /home/mrs186/abaumtobseq/seq_200313/reads/CooperLabMRS/031320_$i/*R1_001.fastq.gz /home/mrs186/abaumtobseq/seq_200313/reads/CooperLabMRS/031320_$i/*R2_001.fastq.gz "$i"_forward_paired.fq.gz "$i"_forward_unpaired.fq.gz "$i"_reverse_paired.fq.gz "$i"_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70;done

#Breseq on populations for variant calling
#for i in 245 246 247 248 249 250 ; do time breseq -r /home/cwm47/ref_genomes/Abaumannii17978/NZ_CP012004_GCF_001077675.1_ASM107767v1_NC009084pAB1_009083pAB2_genomic.gbff /home/mrs186/abaumtobseq/seq_200313/trimmed/"$i"_forward_paired.fq.gz /home/mrs186/abaumtobseq/seq_200313/trimmed/"$i"_forward_unpaired.fq.gz /home/mrs186/abaumtobseq/seq_200313/trimmed/"$i"_reverse_paired.fq.gz -o /home/mrs186/abaumtobseq/seq_200313/breseq/"$i" -j 8; done 

# create breseq output excel file 
# This compiles each breseq output and forms them into an .xlsx file with tabs for SNPs, MC, and NJE ##
#/home/mrs186/scripts/BreseqCatEdited.py -d /home/mrs186/abaumtobseq/seq_200313/breseq

# Copy to computer =
#scp -r mrs186@beagle.mmg.pitt.edu://home/mrs186/abaumtobseq/seq_200313/breseq/Breseq_Output.xlsx /Users/mrs/sequences/seq_200313

# Save the snp tab as a csv format ##

## Convert to noutf characters ##
## Gets rid of spaces in all columns except description. convert arrows to '>'. get rid of all commas, remove Â, ‑,– ##
system("iconv -c -f utf-8 -t ascii /Users/mrs/sequences/seq_200313/Breseq_Output_snps.csv > /Users/mrs/sequences/seq_200313/Breseq_Output_snps_noutf.csv")

snps <- read.csv("/Users/mrs/sequences/seq_200313/Breseq_Output_snps_noutf.csv",header=TRUE)
snps$Position <- as.numeric(gsub(",", "", snps$Position))

# Remove SNPs already present in the ancestor
ref_snps <- read.csv("/Users/mrs/Documents/tob_paper_data/abaumannii_WGS/abaum_ancestor.csv",header=TRUE)
ref_snps$Position <- gsub(":.*", "", ref_snps$Position)
ref_snps$Position <- gsub(",", "", ref_snps$Position)

## Remove all common mutations based on the Position column
snps_noref <- snps[ !(snps$Position %in% ref_snps$Position), ]
## Or Remove all common mutations based on the Gene column
snps_gene_noref <- snps[ !(snps$Gene %in% ref_snps$Gene), ]
nrow(snps_noref)#9

## Create the filtered csv ##
write.csv(snps_noref_clean,file=("/Users/mrs/sequences/seq_200313/clones_snps_noref.csv"))

################### New Junction Evidence #########################

nje <- read.csv("/Users/mrs/sequences/seq_200313/Breseq_Output_nje.csv",header=TRUE)

evens <- seq(from = 2, to = nrow(nje), by = 2)
nje_even <- nje[evens, ]
nje_even <- nje_even[, c(3,4,9,10,11)]
colnames(nje_even) <- c("position2", "reads..cov2", "annotation2", "gene2", "product2")
nje_odd <- nje[-evens, ]
nje_done <- cbind(nje_odd, nje_even)

#remove ancestral nje
nje_done <- subset(nje_done, nje_done$Gene != "A1S_3469/A1S_3470")
nje_done <- subset(nje_done, nje_done$Gene != "ACX60_RS06475/ACX60_RS06480")
nje_done <- subset(nje_done, nje_done$Gene != "ACX60_RS13510/ACX60_RS13515")

write.csv(nje_done, "/Users/mrs/sequences/seq_200313/nje_done.csv")

#combine these manually with SNP mutations 

