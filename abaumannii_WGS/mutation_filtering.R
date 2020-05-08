#Example code for filtering WGS data for A. baumannii populations in Scribner et al

library(dplyr)
library(reshape2)

######## SNPs #############

# Mutations called in breseq as described in Scribner et al were converted to CSV files with no utf characters 

# read in csv file containing snp data for all populations
df_snp <- read.csv("/Users/mrs/Documents/tob_paper_data/abaumannii_WGS/ab_tob_snps.csv")

#Mutation data frames found in common reference = false positives, read in ancestor
df_1_10ref_snps <- read.csv("/Users/mrs/Documents/tob_paper_data/abaumannii_WGS/abaum_ancestor.csv", header=TRUE)
df_1_10ref_snps$chr_pos <- paste(df_1_10ref_snps$Seq.ID,df_1_10ref_snps$Position,sep="_")
df_1_10ref_snps$Position <- gsub(",", "", df_1_10ref_snps$Position)

#remove common mutations
df_snp_noref <- df_snp[ !(df_snp$Position %in% df_1_10ref_snps$Position), ]

#create a column for gene_desc
df_snp_noref$gene_description <- paste(df_snp_noref$Description, df_snp_noref$Gene, sep="::") 
#combine gene_description_annotation in one column
df_snp_noref$desc_gene_annot <-  paste(df_snp_noref$gene_description, df_snp_noref$Annotation, sep="::") 
df_snp_noref$desc_gene_annot <- as.factor(df_snp_noref$desc_gene_annot)
#remove '%' symbol
df_snp_noref$Frequency <- gsub( "%", "", as.character(df_snp_noref$Frequency))
df_snp_noref$Frequency <- as.numeric(as.character(df_snp_noref$Frequency))
#convert Position to numeric
df_snp_noref$Position <- gsub(":.*","",df_snp_noref$Position)
df_snp_noref$Position <- as.numeric(as.character(gsub(",","",df_snp_noref$Position)))
#melted data frame for future casting
m_df_snp_noref <- melt(df_snp_noref, id=c("sample_name","replicate","day","biofilm_plank","xmic","Annotation","Gene","desc_gene_annot"),measure.vars = c("Frequency"))
#cast data frame
c_df_snp_noref <- t(dcast(m_df_snp_noref,sample_name+xmic~desc_gene_annot,mean, value.var = "value",fill=0))
c_df_snp_noref <- as.data.frame(c_df_snp_noref,header=TRUE)
colnames(c_df_snp_noref) <- as.character(unlist(c_df_snp_noref[1,]))
c_df_snp_noref  <- c_df_snp_noref[-1,]

cast_df <- c_df_snp_noref
cast_df_rows <- rownames(cast_df)
cast_df <-sapply(cast_df[2:nrow(cast_df),], function(x) as.numeric(as.character(x)))

rownames(cast_df) <- cast_df_rows[2:length(cast_df_rows)]
cast_df <- as.data.frame(cast_df)
cast_df_sum <- transform(cast_df, count=rowSums(cast_df!=0.0), sum=rowSums(cast_df))
cast_df_sum$desc_gene_annot <- rownames(cast_df_sum)

cast_df_25 <- subset(cast_df_sum, cast_df_sum$sum > 25)

### Mutations in ptsP are prevalent in New Junction Evidence. 
# Given that we have isolated clones of ptsP mutants and know the mutation to be real, 
# I have added mutations in these genes and related (Hpr) into the analysis. 
# I have also added NADH mutations as these are related to ETC like cyoAB mutations.

######################## NJE #######################

df_nje <- as.data.frame(read.csv("/Users/mrs/Documents/tob_paper_data/abaumannii_WGS/ab_tob_nje.csv"))
df_nje$Sample <- as.character(df_nje$Sample)

# subset high quality mutations:
#ACX60_RS15180 phosphocarrier protein HPr, ACX60_RS16050 phosphoenolpyruvateprotein phosphotransferase, ACX60_RS14610 NADHquinone oxidoreductase subunit B
df_nje <- df_nje[df_nje$Gene == "ACX60_RS15180" | df_nje$Gene == "ACX60_RS16050" | df_nje$Gene == "ACX60_RS14610", ]

#create a column for desc_gene_annot
df_nje$desc_gene_annot <- paste(df_nje$Product, df_nje$Gene, df_nje$Annotation, sep="::") 
df_nje <- (merge(df_nje,sample_key,by="Sample"))
#remove '%' symbol
df_nje$Freq <- gsub( "%", "", as.character(df_nje$Freq))
df_nje$Freq <- as.numeric(as.character(df_nje$Freq))
#melted data frame for future casting
m_df_nje <- melt(df_nje, id=c("sample_name","replicate","day","biofilm_plank","xmic","Annotation","Gene","desc_gene_annot"),measure.vars = c("Freq"))

#merge with df_snps
m_df_nje <- rbind(m_df_nje, m_df_snp_noref)

#cast data frame
c_df_nje <- t(dcast(m_df_nje,sample_name~desc_gene_annot,mean, value.var = "value",fill=0))
c_df_nje <- as.data.frame(c_df_nje,header=TRUE)
colnames(c_df_nje) <- as.character(unlist(c_df_nje[1,]))
c_df_nje  <- c_df_nje[-1,]
cast_df_nje <- c_df_nje
#convert all values from factor to numeric
cast_df_nje_rows <- rownames(cast_df_nje)
cast_df_nje <-sapply(cast_df_nje[2:nrow(cast_df_nje),], function(x) as.numeric(as.character(x)))

rownames(cast_df_nje) <- cast_df_nje_rows[2:length(cast_df_nje_rows)]
cast_df_nje <- as.data.frame(cast_df_nje)
cast_df_nje_sum <- transform(cast_df_nje, count=rowSums(cast_df_nje!=0.0), sum=rowSums(cast_df_nje))
cast_df_nje_sum$desc_gene_annot <- rownames(cast_df_nje_sum)

cast_df_nje_25 <- subset(cast_df_nje_sum, cast_df_nje_sum$sum > 25)

### Biologically implausible mutations were removed manually below. 
# Mutations were removed according to the rationale described in Scribner et al to remove mutations suspected to be mapping errors on manual examinations of read pileups.
cast_df_nje_25 <- subset(cast_df_nje_25, cast_df_nje_25$desc_gene_annot != "hypothetical protein::ACX60_RS02305::G499G(GGCGGT)") #There are many snps in same gene indicating poor mapping, mutations occur inconsistently in populations.
cast_df_nje_25 <- subset(cast_df_nje_25, cast_df_nje_25$desc_gene_annot != "hypothetical protein::ACX60_RS02305::I437V(ATTGTT)") 
cast_df_nje_25 <- subset(cast_df_nje_25, cast_df_nje_25$desc_gene_annot != "hypothetical protein::ACX60_RS02305::I501I(ATTATC)") 
cast_df_nje_25 <- subset(cast_df_nje_25, cast_df_nje_25$desc_gene_annot != "hypothetical protein::ACX60_RS04030::G1324G(GGTGGC)") #Inconsistent populations
cast_df_nje_25 <- subset(cast_df_nje_25, cast_df_nje_25$desc_gene_annot != "hypothetical protein::ACX60_RS04030::G1326G(GGTGGC)") 
cast_df_nje_25 <- subset(cast_df_nje_25, cast_df_nje_25$desc_gene_annot != "hypothetical protein::ACX60_RS14465::A55A(GCCGCT)") #Inconsistent populations
cast_df_nje_25 <- subset(cast_df_nje_25, cast_df_nje_25$desc_gene_annot != "hypothetical protein::ACX60_RS14465::T52P(ACCCCC)") 
cast_df_nje_25 <- subset(cast_df_nje_25, cast_df_nje_25$desc_gene_annot != "hypothetical protein::ACX60_RS14465::T52T(ACCACA)") 
cast_df_nje_25 <- subset(cast_df_nje_25, cast_df_nje_25$desc_gene_annot != "hypothetical protein::ACX60_RS14655::G572G(GGTGGC)") 
cast_df_nje_25 <- subset(cast_df_nje_25, cast_df_nje_25$desc_gene_annot != "hypothetical protein::ACX60_RS14655::Y578Y(TATTAC)") 
cast_df_nje_25 <- subset(cast_df_nje_25, cast_df_nje_25$desc_gene_annot != "tRNAGly/tRNAGly::ACX60_RS17465/ACX60_RS17470::intergenic(+1/46)") #Only ends of reads
cast_df_nje_25 <- subset(cast_df_nje_25, cast_df_nje_25$desc_gene_annot != "peptidase, M23 family protein::ACX60_RS00605::coding(644751/2439nt)") #Junction, inconsistent pops, low coverage
cast_df_nje_25 <- subset(cast_df_nje_25, cast_df_nje_25$desc_gene_annot != "hypothetical protein::ACX60_RS00600::coding(8587/183nt)") # Low frequency junction in inconsistent populations near a region of repeats
cast_df_nje_25 <- subset(cast_df_nje_25, cast_df_nje_25$desc_gene_annot != "protein TolA::ACX60_RS04595::R272R(CGTCGA)")  #Lots of mutations within gene suggesting poor mapping, occurs twice in two different populations and never in consecutive timepoints
cast_df_nje_25 <- subset(cast_df_nje_25, cast_df_nje_25$desc_gene_annot != "regulatory protein LysR/::A1S_3470/::intergenic(245/)") #  Weird junction with mostly ends on one side, inconsistent populations

# These mutations are duplicate mutations to those reported in new junction evidence.
# They reflect the other side of the junction for the mutations reported above, therefore have been removed as these are redundant. 
cast_df_nje_25 <- subset(cast_df_nje_25, cast_df_nje_25$desc_gene_annot != "phosphocarrier protein HPr::ACX60_RS15180::coding(27/270nt)")
cast_df_nje_25 <- subset(cast_df_nje_25, cast_df_nje_25$desc_gene_annot != "NADHquinone oxidoreductase subunit B::ACX60_RS14610::coding(673/678nt)")
cast_df_nje_25 <- subset(cast_df_nje_25, cast_df_nje_25$desc_gene_annot != "phosphoenolpyruvateprotein phosphotransferase::ACX60_RS16050::coding(1251/2295nt)")
cast_df_nje_25 <- subset(cast_df_nje_25, cast_df_nje_25$desc_gene_annot != "phosphoenolpyruvateprotein phosphotransferase::ACX60_RS16050::coding(1252/2295nt)")
cast_df_nje_25 <- subset(cast_df_nje_25, cast_df_nje_25$desc_gene_annot != "phosphoenolpyruvateprotein phosphotransferase::ACX60_RS16050::coding(1742/2295nt)")
cast_df_nje_25 <- subset(cast_df_nje_25, cast_df_nje_25$desc_gene_annot != "phosphoenolpyruvateprotein phosphotransferase::ACX60_RS16050::coding(1788/2295nt)")
cast_df_nje_25 <- subset(cast_df_nje_25, cast_df_nje_25$desc_gene_annot != "phosphoenolpyruvateprotein phosphotransferase::ACX60_RS16050::coding(2232/2295nt)")
cast_df_nje_25 <- subset(cast_df_nje_25, cast_df_nje_25$desc_gene_annot != "phosphoenolpyruvateprotein phosphotransferase::ACX60_RS16050::coding(496/2295nt)")

write.csv(cast_df_nje_25, file = "/Users/mrs/Documents/tob_paper_data/abaumannii_WGS/nje_25_filtered.csv")

