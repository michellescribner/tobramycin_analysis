#! bin/bash

## Trimmomatic for pa14 in tobramycin in M9 CLONES ##

## scp /Users/mrs/Documents/PA14_Tobi_M9/PA14_tobi_m9_clones_trim_shellscript.sh mrs186@cepacia.mmg.pitt.edu://home/mrs186/tobim9 ##

trimmomatic PE -phred33 -threads 4 -trimlog trim.log /home/dmux/181025/CooperLabMRS/102518_71/936_2_S71_R1_001.fastq.gz /home/dmux/181025/CooperLabMRS/102518_71/936_2_S71_R2_001.fastq.gz /home/mrs186/tobim9/trimmed/936_2_forward_paired.fq.gz /home/mrs186/tobim9/trimmed/936_2_forward_unpaired.fq.gz  /home/mrs186/tobim9/trimmed/936_2_reverse_paired.fq.gz /home/mrs186/tobim9/trimmed/936_2_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
trimmomatic PE -phred33 -threads 4 -trimlog trim.log /home/dmux/181025/CooperLabMRS/102518_72/936_3_S72_R1_001.fastq.gz /home/dmux/181025/CooperLabMRS/102518_72/936_3_S72_R2_001.fastq.gz /home/mrs186/tobim9/trimmed/936_3_forward_paired.fq.gz /home/mrs186/tobim9/trimmed/936_3_forward_unpaired.fq.gz  /home/mrs186/tobim9/trimmed/936_3_reverse_paired.fq.gz /home/mrs186/tobim9/trimmed/936_3_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
trimmomatic PE -phred33 -threads 4 -trimlog trim.log /home/dmux/181025/CooperLabMRS/102518_73/936_4_S73_R1_001.fastq.gz /home/dmux/181025/CooperLabMRS/102518_73/936_4_S73_R2_001.fastq.gz /home/mrs186/tobim9/trimmed/936_4_forward_paired.fq.gz /home/mrs186/tobim9/trimmed/936_4_forward_unpaired.fq.gz  /home/mrs186/tobim9/trimmed/936_4_reverse_paired.fq.gz /home/mrs186/tobim9/trimmed/936_4_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
trimmomatic PE -phred33 -threads 4 -trimlog trim.log /home/dmux/181025/CooperLabMRS/102518_74/997_1_S74_R1_001.fastq.gz /home/dmux/181025/CooperLabMRS/102518_74/997_1_S74_R2_001.fastq.gz /home/mrs186/tobim9/trimmed/997_1_forward_paired.fq.gz /home/mrs186/tobim9/trimmed/997_1_forward_unpaired.fq.gz  /home/mrs186/tobim9/trimmed/997_1_reverse_paired.fq.gz /home/mrs186/tobim9/trimmed/997_1_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
trimmomatic PE -phred33 -threads 4 -trimlog trim.log /home/dmux/181025/CooperLabMRS/102518_75/997_2_S75_R1_001.fastq.gz /home/dmux/181025/CooperLabMRS/102518_75/997_2_S75_R2_001.fastq.gz /home/mrs186/tobim9/trimmed/997_2_forward_paired.fq.gz /home/mrs186/tobim9/trimmed/997_2_forward_unpaired.fq.gz  /home/mrs186/tobim9/trimmed/997_2_reverse_paired.fq.gz /home/mrs186/tobim9/trimmed/997_2_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
trimmomatic PE -phred33 -threads 4 -trimlog trim.log /home/dmux/181025/CooperLabMRS/102518_76/997_3_S76_R1_001.fastq.gz /home/dmux/181025/CooperLabMRS/102518_76/997_3_S76_R2_001.fastq.gz /home/mrs186/tobim9/trimmed/997_3_forward_paired.fq.gz /home/mrs186/tobim9/trimmed/997_3_forward_unpaired.fq.gz  /home/mrs186/tobim9/trimmed/997_3_reverse_paired.fq.gz /home/mrs186/tobim9/trimmed/997_3_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
trimmomatic PE -phred33 -threads 4 -trimlog trim.log /home/dmux/181025/CooperLabMRS/102518_77/997_4_S77_R1_001.fastq.gz /home/dmux/181025/CooperLabMRS/102518_77/997_4_S77_R2_001.fastq.gz /home/mrs186/tobim9/trimmed/997_4_forward_paired.fq.gz /home/mrs186/tobim9/trimmed/997_4_forward_unpaired.fq.gz  /home/mrs186/tobim9/trimmed/997_4_reverse_paired.fq.gz /home/mrs186/tobim9/trimmed/997_4_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
trimmomatic PE -phred33 -threads 4 -trimlog trim.log /home/dmux/181025/CooperLabMRS/102518_78/997_5_S78_R1_001.fastq.gz /home/dmux/181025/CooperLabMRS/102518_78/997_5_S78_R2_001.fastq.gz /home/mrs186/tobim9/trimmed/997_5_forward_paired.fq.gz /home/mrs186/tobim9/trimmed/997_5_forward_unpaired.fq.gz  /home/mrs186/tobim9/trimmed/997_5_reverse_paired.fq.gz /home/mrs186/tobim9/trimmed/997_5_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
trimmomatic PE -phred33 -threads 4 -trimlog trim.log /home/dmux/181025/CooperLabMRS/102518_79/916_S79_R1_001.fastq.gz /home/dmux/181025/CooperLabMRS/102518_79/916_S79_R2_001.fastq.gz /home/mrs186/tobim9/trimmed/916_forward_paired.fq.gz /home/mrs186/tobim9/trimmed/916_forward_unpaired.fq.gz  /home/mrs186/tobim9/trimmed/916_reverse_paired.fq.gz /home/mrs186/tobim9/trimmed/916_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
trimmomatic PE -phred33 -threads 4 -trimlog trim.log /home/dmux/181025/CooperLabMRS/102518_80/938_S80_R1_001.fastq.gz /home/dmux/181025/CooperLabMRS/102518_80/938_S80_R2_001.fastq.gz /home/mrs186/tobim9/trimmed/938_forward_paired.fq.gz /home/mrs186/tobim9/trimmed/938_forward_unpaired.fq.gz  /home/mrs186/tobim9/trimmed/938_reverse_paired.fq.gz /home/mrs186/tobim9/trimmed/938_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
trimmomatic PE -phred33 -threads 4 -trimlog trim.log /home/dmux/181025/CooperLabMRS/102518_81/957_S81_R1_001.fastq.gz /home/dmux/181025/CooperLabMRS/102518_81/957_S81_R2_001.fastq.gz /home/mrs186/tobim9/trimmed/957_forward_paired.fq.gz /home/mrs186/tobim9/trimmed/957_forward_unpaired.fq.gz  /home/mrs186/tobim9/trimmed/957_reverse_paired.fq.gz /home/mrs186/tobim9/trimmed/957_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
trimmomatic PE -phred33 -threads 4 -trimlog trim.log /home/dmux/181025/CooperLabMRS/102518_82/978_S82_R1_001.fastq.gz /home/dmux/181025/CooperLabMRS/102518_82/978_S82_R2_001.fastq.gz /home/mrs186/tobim9/trimmed/978_forward_paired.fq.gz /home/mrs186/tobim9/trimmed/978_forward_unpaired.fq.gz  /home/mrs186/tobim9/trimmed/978_reverse_paired.fq.gz /home/mrs186/tobim9/trimmed/978_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
trimmomatic PE -phred33 -threads 4 -trimlog trim.log /home/dmux/181025/CooperLabMRS/102518_83/1027_4_S83_R1_001.fastq.gz /home/dmux/181025/CooperLabMRS/102518_83/1027_4_S83_R2_001.fastq.gz /home/mrs186/tobim9/trimmed/1027_4_forward_paired.fq.gz /home/mrs186/tobim9/trimmed/1027_4_forward_unpaired.fq.gz  /home/mrs186/tobim9/trimmed/1027_4_reverse_paired.fq.gz /home/mrs186/tobim9/trimmed/1027_4_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
trimmomatic PE -phred33 -threads 4 -trimlog trim.log /home/dmux/181025/CooperLabMRS/102518_84/1027_3_S84_R1_001.fastq.gz /home/dmux/181025/CooperLabMRS/102518_84/1027_3_S84_R2_001.fastq.gz /home/mrs186/tobim9/trimmed/1027_3_forward_paired.fq.gz /home/mrs186/tobim9/trimmed/1027_3_forward_unpaired.fq.gz  /home/mrs186/tobim9/trimmed/1027_3_reverse_paired.fq.gz /home/mrs186/tobim9/trimmed/1027_3_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
trimmomatic PE -phred33 -threads 4 -trimlog trim.log /home/dmux/181025/CooperLabMRS/102518_85/1027_2_S85_R1_001.fastq.gz /home/dmux/181025/CooperLabMRS/102518_85/1027_2_S85_R2_001.fastq.gz /home/mrs186/tobim9/trimmed/1027_2_forward_paired.fq.gz /home/mrs186/tobim9/trimmed/1027_2_forward_unpaired.fq.gz  /home/mrs186/tobim9/trimmed/1027_2_reverse_paired.fq.gz /home/mrs186/tobim9/trimmed/1027_2_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
trimmomatic PE -phred33 -threads 4 -trimlog trim.log /home/dmux/181025/CooperLabMRS/102518_86/1027_1_S86_R1_001.fastq.gz /home/dmux/181025/CooperLabMRS/102518_86/1027_1_S86_R2_001.fastq.gz /home/mrs186/tobim9/trimmed/1027_1_forward_paired.fq.gz /home/mrs186/tobim9/trimmed/1027_1_forward_unpaired.fq.gz  /home/mrs186/tobim9/trimmed/1027_1_reverse_paired.fq.gz /home/mrs186/tobim9/trimmed/1027_1_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
trimmomatic PE -phred33 -threads 4 -trimlog trim.log /home/dmux/181025/CooperLabMRS/102518_87/1017_4_S87_R1_001.fastq.gz /home/dmux/181025/CooperLabMRS/102518_87/1017_4_S87_R2_001.fastq.gz /home/mrs186/tobim9/trimmed/1017_4_forward_paired.fq.gz /home/mrs186/tobim9/trimmed/1017_4_forward_unpaired.fq.gz  /home/mrs186/tobim9/trimmed/1017_4_reverse_paired.fq.gz /home/mrs186/tobim9/trimmed/1017_4_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
trimmomatic PE -phred33 -threads 4 -trimlog trim.log /home/dmux/181025/CooperLabMRS/102518_88/1017_5_S88_R1_001.fastq.gz /home/dmux/181025/CooperLabMRS/102518_88/1017_5_S88_R2_001.fastq.gz /home/mrs186/tobim9/trimmed/1017_5_forward_paired.fq.gz /home/mrs186/tobim9/trimmed/1017_5_forward_unpaired.fq.gz  /home/mrs186/tobim9/trimmed/1017_5_reverse_paired.fq.gz /home/mrs186/tobim9/trimmed/1017_5_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
trimmomatic PE -phred33 -threads 4 -trimlog trim.log /home/dmux/181025/CooperLabMRS/102518_89/917_1_S89_R1_001.fastq.gz /home/dmux/181025/CooperLabMRS/102518_89/917_1_S89_R2_001.fastq.gz /home/mrs186/tobim9/trimmed/917_1_forward_paired.fq.gz /home/mrs186/tobim9/trimmed/917_1_forward_unpaired.fq.gz  /home/mrs186/tobim9/trimmed/917_1_reverse_paired.fq.gz /home/mrs186/tobim9/trimmed/917_1_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
trimmomatic PE -phred33 -threads 4 -trimlog trim.log /home/dmux/181025/CooperLabMRS/102518_90/917_2_S90_R1_001.fastq.gz /home/dmux/181025/CooperLabMRS/102518_90/917_2_S90_R2_001.fastq.gz /home/mrs186/tobim9/trimmed/917_2_forward_paired.fq.gz /home/mrs186/tobim9/trimmed/917_2_forward_unpaired.fq.gz  /home/mrs186/tobim9/trimmed/917_2_reverse_paired.fq.gz /home/mrs186/tobim9/trimmed/917_2_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
trimmomatic PE -phred33 -threads 4 -trimlog trim.log /home/dmux/181025/CooperLabMRS/102518_91/917_3_S91_R1_001.fastq.gz /home/dmux/181025/CooperLabMRS/102518_91/917_3_S91_R2_001.fastq.gz /home/mrs186/tobim9/trimmed/917_3_forward_paired.fq.gz /home/mrs186/tobim9/trimmed/917_3_forward_unpaired.fq.gz  /home/mrs186/tobim9/trimmed/917_3_reverse_paired.fq.gz /home/mrs186/tobim9/trimmed/917_3_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
trimmomatic PE -phred33 -threads 4 -trimlog trim.log /home/dmux/181025/CooperLabMRS/102518_92/917_4_S92_R1_001.fastq.gz /home/dmux/181025/CooperLabMRS/102518_92/917_4_S92_R2_001.fastq.gz /home/mrs186/tobim9/trimmed/917_4_forward_paired.fq.gz /home/mrs186/tobim9/trimmed/917_4_forward_unpaired.fq.gz  /home/mrs186/tobim9/trimmed/917_4_reverse_paired.fq.gz /home/mrs186/tobim9/trimmed/917_4_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
trimmomatic PE -phred33 -threads 4 -trimlog trim.log /home/dmux/181025/CooperLabMRS/102518_93/917_5_S93_R1_001.fastq.gz /home/dmux/181025/CooperLabMRS/102518_93/917_5_S93_R2_001.fastq.gz /home/mrs186/tobim9/trimmed/917_5_forward_paired.fq.gz /home/mrs186/tobim9/trimmed/917_5_forward_unpaired.fq.gz  /home/mrs186/tobim9/trimmed/917_5_reverse_paired.fq.gz /home/mrs186/tobim9/trimmed/917_5_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
trimmomatic PE -phred33 -threads 4 -trimlog trim.log /home/dmux/181025/CooperLabMRS/102518_94/947_1_S94_R1_001.fastq.gz /home/dmux/181025/CooperLabMRS/102518_94/947_1_S94_R2_001.fastq.gz /home/mrs186/tobim9/trimmed/947_1_forward_paired.fq.gz /home/mrs186/tobim9/trimmed/947_1_forward_unpaired.fq.gz  /home/mrs186/tobim9/trimmed/947_1_reverse_paired.fq.gz /home/mrs186/tobim9/trimmed/947_1_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
trimmomatic PE -phred33 -threads 4 -trimlog trim.log /home/dmux/181025/CooperLabMRS/102518_95/947_2_S95_R1_001.fastq.gz /home/dmux/181025/CooperLabMRS/102518_95/947_2_S95_R2_001.fastq.gz /home/mrs186/tobim9/trimmed/947_2_forward_paired.fq.gz /home/mrs186/tobim9/trimmed/947_2_forward_unpaired.fq.gz  /home/mrs186/tobim9/trimmed/947_2_reverse_paired.fq.gz /home/mrs186/tobim9/trimmed/947_2_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
trimmomatic PE -phred33 -threads 4 -trimlog trim.log /home/dmux/181025/CooperLabMRS/102518_96/947_3_S96_R1_001.fastq.gz /home/dmux/181025/CooperLabMRS/102518_96/947_3_S96_R2_001.fastq.gz /home/mrs186/tobim9/trimmed/947_3_forward_paired.fq.gz /home/mrs186/tobim9/trimmed/947_3_forward_unpaired.fq.gz  /home/mrs186/tobim9/trimmed/947_3_reverse_paired.fq.gz /home/mrs186/tobim9/trimmed/947_3_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
trimmomatic PE -phred33 -threads 4 -trimlog trim.log /home/dmux/181027/CooperLabMRS/102718_79/1017_1_S79_R1_001.fastq.gz /home/dmux/181027/CooperLabMRS/102718_79/1017_1_S79_R2_001.fastq.gz /home/mrs186/tobim9/trimmed/1017_1_forward_paired.fq.gz /home/mrs186/tobim9/trimmed/1017_1_forward_unpaired.fq.gz  /home/mrs186/tobim9/trimmed/1017_1_reverse_paired.fq.gz /home/mrs186/tobim9/trimmed/1017_1_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
trimmomatic PE -phred33 -threads 4 -trimlog trim.log /home/dmux/181027/CooperLabMRS/102718_80/1017_2_S80_R1_001.fastq.gz /home/dmux/181027/CooperLabMRS/102718_80/1017_2_S80_R2_001.fastq.gz /home/mrs186/tobim9/trimmed/1017_2_forward_paired.fq.gz /home/mrs186/tobim9/trimmed/1017_2_forward_unpaired.fq.gz  /home/mrs186/tobim9/trimmed/1017_2_reverse_paired.fq.gz /home/mrs186/tobim9/trimmed/1017_2_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
trimmomatic PE -phred33 -threads 4 -trimlog trim.log /home/dmux/181027/CooperLabMRS/102718_81/1017_3_S81_R1_001.fastq.gz /home/dmux/181027/CooperLabMRS/102718_81/1017_3_S81_R2_001.fastq.gz /home/mrs186/tobim9/trimmed/1017_3_forward_paired.fq.gz /home/mrs186/tobim9/trimmed/1017_3_forward_unpaired.fq.gz  /home/mrs186/tobim9/trimmed/1017_3_reverse_paired.fq.gz /home/mrs186/tobim9/trimmed/1017_3_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
trimmomatic PE -phred33 -threads 4 -trimlog trim.log /home/dmux/181027/CooperLabMRS/102718_82/10274_5_S82_R1_001.fastq.gz /home/dmux/181027/CooperLabMRS/102718_82/10274_5_S82_R2_001.fastq.gz /home/mrs186/tobim9/trimmed/10274_5_forward_paired.fq.gz /home/mrs186/tobim9/trimmed/10274_5_forward_unpaired.fq.gz  /home/mrs186/tobim9/trimmed/10274_5_reverse_paired.fq.gz /home/mrs186/tobim9/trimmed/10274_5_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
trimmomatic PE -phred33 -threads 4 -trimlog trim.log /home/dmux/181027/CooperLabMRS/102718_84/947_4_S84_R1_001.fastq.gz /home/dmux/181027/CooperLabMRS/102718_84/947_4_S84_R2_001.fastq.gz /home/mrs186/tobim9/trimmed/947_4_forward_paired.fq.gz /home/mrs186/tobim9/trimmed/947_4_forward_unpaired.fq.gz  /home/mrs186/tobim9/trimmed/947_4_reverse_paired.fq.gz /home/mrs186/tobim9/trimmed/947_4_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
trimmomatic PE -phred33 -threads 4 -trimlog trim.log /home/dmux/181027/CooperLabMRS/102718_85/947_5_S85_R1_001.fastq.gz /home/dmux/181027/CooperLabMRS/102718_85/947_5_S85_R2_001.fastq.gz /home/mrs186/tobim9/trimmed/947_5_forward_paired.fq.gz /home/mrs186/tobim9/trimmed/947_5_forward_unpaired.fq.gz  /home/mrs186/tobim9/trimmed/947_5_reverse_paired.fq.gz /home/mrs186/tobim9/trimmed/947_5_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
trimmomatic PE -phred33 -threads 4 -trimlog trim.log /home/dmux/181027/CooperLabMRS/102718_86/1021_1_S86_R1_001.fastq.gz /home/dmux/181027/CooperLabMRS/102718_86/1021_1_S86_R2_001.fastq.gz /home/mrs186/tobim9/trimmed/1021_1_forward_paired.fq.gz /home/mrs186/tobim9/trimmed/1021_1_forward_unpaired.fq.gz  /home/mrs186/tobim9/trimmed/1021_1_reverse_paired.fq.gz /home/mrs186/tobim9/trimmed/1021_1_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70
























