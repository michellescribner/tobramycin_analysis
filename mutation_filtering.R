### Filtering for Evolution of PA14 in Tobramycin for Scribner et al  ###

library("reshape2")

setwd("/Users/mrs/Documents/tob_paper_data/paeruginosa_WGS")

# Mutations called in breseq as described in Scribner et al were converted to CSV files with no utf characters and informative sample names were substituted for those assigned for sequencing.
snps <- read.csv("/Users/mrs/Documents/tob_paper_data/paeruginosa_WGS/pa_tob_snps.csv",header=TRUE)

# Remove SNPs already present in the ancestor
ancestor_snps <- read.csv("/Users/mrs/Documents/tob_paper_data/paeruginosa_WGS/pa_ancestor.csv",header=TRUE)
snps_noref <- snps[ !(snps$Position %in% ancestor_snps$Position), ]
nrow(snps_noref)#1429

#add gene_annot column
snps_noref$gene_annot <- paste(snps_noref$Gene,snps_noref$Annotation,sep="::")

## remove manually filtered mutations
# The breseq documentation as of 4/30/2020 reads "Polymorphism prediction is still considered a somewhat experimental feature. It continues to be actively developed." 
# Therefore, we manually examined read pileups for mutations reported by breseq in polymorphism mode. Mutations that appear biologically implausible were removed from further analysis. 
# These mutations and the rationale motivating their removal are listed in manuallyfiltered_mutations and in the manuscript.
manuallyfiltered_mutations <- read.csv("/Users/mrs/Documents/tob_paper_data/paeruginosa_WGS/manuallyfiltered_mutations.csv", header=TRUE)
manuallyfiltered_mutations$gene_annot <- paste(manuallyfiltered_mutations$gene,manuallyfiltered_mutations$annotation,sep="::")
snps_noref <- snps_noref[ !(snps_noref$gene_annot %in% manuallyfiltered_mutations$gene_annot), ]

### Rename locus tages to gene names or old RefSeq locus tags
snps_noref_clean <- snps_noref
snps_noref_clean$Gene <- gsub( "fusA", "fusA1", as.character(snps_noref_clean$Gene))
snps_noref_clean$Gene <- gsub( "PA14_RS01795", "ptsP", as.character(snps_noref_clean$Gene))
snps_noref_clean$Gene <- gsub( "PA14_RS25825", "pmrB", as.character(snps_noref_clean$Gene))
snps_noref_clean$Gene <- gsub( "PA14_RS06610", "wspF", as.character(snps_noref_clean$Gene))
snps_noref_clean$Gene <- gsub( "PA14_RS06585", "wspA", as.character(snps_noref_clean$Gene))
snps_noref_clean$Gene <- gsub( "PA14_RS24870", "morA", as.character(snps_noref_clean$Gene))
snps_noref_clean$Gene <- gsub( "PA14_RS18655", "lasR", as.character(snps_noref_clean$Gene))
snps_noref_clean$Gene <- gsub( "PA14_RS26180", "PA14_64050", as.character(snps_noref_clean$Gene))
snps_noref_clean$Gene <- gsub( "PA14_RS22520", "lysR", as.character(snps_noref_clean$Gene))
snps_noref_clean$Gene <- gsub( "PA14_RS23785", "PA14_58320", as.character(snps_noref_clean$Gene))
snps_noref_clean$Gene <- gsub( "PA14_RS19965", "PA14_49160", as.character(snps_noref_clean$Gene))
snps_noref_clean$Gene <- gsub( "PA14_RS20335", "PA14_50060", as.character(snps_noref_clean$Gene))
snps_noref_clean$Gene <- gsub( "PA14_RS09370", "PA14_23130", as.character(snps_noref_clean$Gene))
snps_noref_clean$Gene <- gsub( "PA14_RS16730", "PA14_41270", as.character(snps_noref_clean$Gene))
snps_noref_clean$Gene <- gsub( "PA14_RS03595", "rplB", as.character(snps_noref_clean$Gene))
snps_noref_clean$Gene <- gsub( "PA14_RS09465", "orfK", as.character(snps_noref_clean$Gene))
snps_noref_clean$Gene <- gsub( "PA14_RS29915", "orfL", as.character(snps_noref_clean$Gene))
snps_noref_clean$Gene <- gsub( "PA14_RS09470", "orfH", as.character(snps_noref_clean$Gene))
snps_noref_clean$Gene <- gsub( "PA14_RS09495", "orfN", as.character(snps_noref_clean$Gene))

#remove % sign from frequency column
snps_noref_clean$Frequency <- gsub( "%", "", as.character(snps_noref_clean$Frequency))
#convert frequency values to numeric
snps_noref_clean$Frequency <- as.numeric(as.character(snps_noref_clean$Frequency))

#split sample column into components
sample <- colsplit(snps_noref_clean$Sample, ", ", c("day", "pop", "b_p", "xmic"))
sample$day <- gsub( "Day ", "", as.character(sample$day))
sample$pop <- gsub( "Pop ", "", as.character(sample$pop))
snps_noref_clean <- cbind(sample, snps_noref_clean)

snps <- snps_noref_clean

#combine gene and annotation in one column
snps$gene_annot_desc <-  paste(snps$Gene,snps$Annotation, snps$Description,sep="::")
snps$gene_annot_desc <- as.factor(snps$gene_annot_desc)
#remove % and convert frequency to numeric
snps$Frequency <- gsub( "%", "", as.character(snps$Frequency))
snps$Frequency <- as.numeric(as.character(snps$Frequency))
#melt 
melt_snp <- melt(snps, id=c("Sample", "day", "pop", "b_p", "xmic", "gene_annot_desc", "Gene", "Frequency","Description"),measure.vars = c("Position"))

# cast data frame - organizing with each mutation as the rows and the frequency of that mutation on a given day as the columns
snps_cast <- t(dcast(melt_snp,Sample~gene_annot_desc,mean, value.var = "Frequency",fill=0))
snps_cast <- as.data.frame(snps_cast,header=TRUE)
colnames(snps_cast) <- as.character(unlist(snps_cast[1,]))
#remove sample column
snps_cast <- snps_cast[-1,]
#convert frequency values to numeric class - start as "character"
snps_cast[,1:ncol(snps_cast)] <- apply(snps_cast[,1:ncol(snps_cast)], 2, function(x) as.numeric(as.character(x)))
nrow(snps_cast) #1009
snps_cast$Sums <- rowSums(snps_cast)

#filter only rows with greater than 25% total frequency
snps_cast_25 <- subset(snps_cast, snps_cast$Sums >= 25) #greater than 25%
nrow(snps_cast_25)#47

t_snps_cast_25 <- t(snps_cast_25)
write.csv(t_snps_cast_25, file="/Users/mrs/Documents/tob_paper_data/paeruginosa_WGS/t_snps_25.csv")

############### New Junction Evidence ####################

# read in NJE mutations
nje <- read.csv("/Users/mrs/Documents/tob_paper_data/paeruginosa_WGS/pa_tob_nje.csv",header=TRUE)

#remove nje mutations that were in the ancestor's breseq output
nje$gene_annot_desc <- paste(nje$Gene,nje$Annotation, nje$Product, nje$`2.Gene`, nje$`2.Annotation`, nje$`2.Product`,sep="::")
nje_noref <- subset(nje, nje$gene_annot_desc != "PA14_RS30510/PA14_RS14515::intergenic(146/+106)::hypothetical protein/hypothetical protein::PA14_RS14525/PA14_RS14530::intergenic(69/218)::hypothetical protein/Tn3 family transposase")
nje_noref <- subset(nje_noref, nje_noref$gene_annot_desc != "PA14_RS19845/PA14_RS19850::intergenic(45/+61)::tRNA 2thiocytidine(32) synthetase TtcA/integrase::PA14_RS31160/PA14_RS19915::intergenic(+317/59)::hypothetical protein/hypothetical protein")

#remove % sign from frequency column
nje_noref$Freq <- gsub( "%", "", as.character(nje_noref$Freq))
#convert frequency values to numeric
nje_noref$Freq <- as.numeric(as.character(nje_noref$Freq))

#cast
nje_cast <- t(dcast(nje_noref,Population~gene_annot_desc,mean, value.var = "Freq",fill=0))
nje_cast <- as.data.frame(nje_cast,header=TRUE)
colnames(nje_cast) <- as.character(unlist(nje_cast[1,]))
#remove sample column
nje_cast <- nje_cast[-1,]
#convert frequency values to numeric class - start as "character"
nje_cast[,1:ncol(nje_cast)] <-apply(nje_cast[,1:ncol(nje_cast)], 2, function(x) as.numeric(as.character(x)))
nrow(nje_cast) #44
nje_cast$Sums <- rowSums(nje_cast)

# Mutations in nje_cast were further examined manually by analyzing read pileups.

# We took a conservative approach to analyzing junction calls identified in polymorphism mode: 
# we included in further analysis only junction calls that occurred in loci (or loci in related pathways) that were confirmed to confer tobramycin resistance by clonal WGS and phenotyping
# or loci previously associated with aminoglycoside resistance (ptsP, HPr, and NADH quinone oxidoreductase for A. baumannii and 23S rRNA for P. aeruginosa).

# While several mutations in 23S rRNA were reported by breseq, only one was determined to be of high enough quality by examining read pileups to be included in subsequent analysis. 
# This mutation was estimated to occur at high frequencies at consecutive timepoints whereas other 23S rRNA other mutations were singletons or reached only low frequencies (ie, would not pass frequency filters applied to SNPs). 
# This mutation also occured at trajectories that complemented the allele frequencies of other contending mutations in that population.
# Thus, this mutation was included in subsequent analysis as its read pileups and dynamics both supported that it was an important player influencing the population dynamics within that population.

