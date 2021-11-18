#P. aeruginosa clones selected for sequencing from tob m9 media populations

#Trim and filter reads
for i in {63..84}; do trimmomatic PE -phred33 /home/dmux/190516/CooperLabMRS/051619_$i/*R1_001.fastq.gz /home/dmux/190516/CooperLabMRS/051619_$i/*R2_001.fastq.gz "$i"_forward_paired.fq.gz "$i"_forward_unpaired.fq.gz "$i"_reverse_paired.fq.gz "$i"_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70;done

#Breseq on populations for variant calling
for i in {63..84}; do time breseq -r /home/cwm47/ref_genomes/Paeruginosa/PA14_GCF_000014625.1_ASM1462v1_genomic.gbff  /home/mrs186/tobim9/clones/trim/"$i"_forward_paired.fq.gz /home/mrs186/tobim9/clones/trim/"$i"_forward_unpaired.fq.gz /home/mrs186/tobim9/clones/trim/"$i"_reverse_paired.fq.gz -o /home/mrs186/tobim9/clones/breseq/breseq_"$i" -j 8; done 

#END - moved analysis into R
