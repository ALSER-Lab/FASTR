#!/bin/bash

echo "=============================================="
echo "FASTR - minimap2 Benchmarking"
echo "=============================================="

#source ~/anaconda3/bin/activate root

#conda activate fastq_bench


#Original minimap2 - Illumina
sbatch -p mi2101x -N 1 -n 16  -t 12:00:00 -o "/work1/malser/mohammedalser/illumina-ERR15909551-results/FASTR-minimap2-NoQality-original-illumina-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/illumina-ERR15909551-results/FASTR-minimap2-NoQality-original-illumina-stats-time-memory.txt" /work1/malser/mohammedalser/minimap2-original/minimap2 -Q -a -x sr -t 16 /work1/malser/mohammedalser/GCF_000001405.40_GRCh38.p14_genomic.fna /work1/malser/mohammedalser/ERR15909551.fastq -o /work1/malser/mohammedalser/illumina-ERR15909551-results/FASTR-minimap2-NoQality-original-illumina.SAM"

#Original minimap2 - HiFi
sbatch -p mi2101x -N 1 -n 16  -t 12:00:00 -o "/work1/malser/mohammedalser/hifi-ERR13491966-results/FASTR-minimap2-NoQality-original-hifi-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/hifi-ERR13491966-results/FASTR-minimap2-NoQality-original-hifi-stats-time-memory.txt" /work1/malser/mohammedalser/minimap2-original/minimap2 -Q -a -x map-hifi -t 16 /work1/malser/mohammedalser/GCF_000001405.40_GRCh38.p14_genomic.fna /work1/malser/mohammedalser/ERR13491966.fastq -o /work1/malser/mohammedalser/hifi-ERR13491966-results/FASTR-minimap2-NoQality-original-hifi.SAM"

#Original minimap2 - ONT
sbatch -p mi2101x -N 1 -n 16  -t 12:00:00 -o "/work1/malser/mohammedalser/ont-SRR33464820-results/FASTR-minimap2-NoQality-original-ont-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/ont-SRR33464820-results/FASTR-minimap2-NoQality-original-ont-stats-time-memory.txt" /work1/malser/mohammedalser/minimap2-original/minimap2 -Q -a -x map-ont -t 16 /work1/malser/mohammedalser/GCF_000001405.40_GRCh38.p14_genomic.fna /work1/malser/mohammedalser/SRR33464820.fastq -o /work1/malser/mohammedalser/ont-SRR33464820-results/FASTR-minimap2-NoQality-original-ont.SAM"



#modified minimap2 - Illumina
sbatch -p mi2101x -N 1 -n 16  -t 12:00:00 -o "/work1/malser/mohammedalser/illumina-ERR15909551-results/FASTR-minimap2-modified-illumina-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/illumina-ERR15909551-results/FASTR-minimap2-modified-illumina-stats-time-memory.txt" /work1/malser/mohammedalser/minimap2/minimap2 -Q -a -x sr -t 16 /work1/malser/mohammedalser/GCF_000001405.40_GRCh38.p14_genomic.fna /work1/malser/mohammedalser/illumina-ERR15909551-results/ERR15909551.fastr_mode2.fastr -o /work1/malser/mohammedalser/illumina-ERR15909551-results/FASTR-minimap2-modified-illumina.SAM"

#modified minimap2 - HiFi
sbatch -p mi2101x -N 1 -n 16  -t 12:00:00 -o "/work1/malser/mohammedalser/hifi-ERR13491966-results/FASTR-minimap2-modified-hifi-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/hifi-ERR13491966-results/FASTR-minimap2-modified-hifi-stats-time-memory.txt" /work1/malser/mohammedalser/minimap2/minimap2 -Q -a -x map-hifi -t 16 /work1/malser/mohammedalser/GCF_000001405.40_GRCh38.p14_genomic.fna /work1/malser/mohammedalser/hifi-ERR13491966-results/ERR13491966.fastr_mode2.fastr -o /work1/malser/mohammedalser/hifi-ERR13491966-results/FASTR-minimap2-modified-hifi.SAM"

#modified minimap2 - ONT
sbatch -p mi2101x -N 1 -n 16  -t 12:00:00 -o "/work1/malser/mohammedalser/ont-SRR33464820-results/FASTR-minimap2-modified-ont-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/ont-SRR33464820-results/FASTR-minimap2-modified-ont-stats-time-memory.txt" /work1/malser/mohammedalser/minimap2/minimap2 -Q -a -x map-ont -t 16 /work1/malser/mohammedalser/GCF_000001405.40_GRCh38.p14_genomic.fna /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820.fastr_mode2.fastr -o /work1/malser/mohammedalser/ont-SRR33464820-results/FASTR-minimap2-modified-ont.SAM"



## Checking alignment results between FASTR and FASTQ

awk -F'\t' '!/^@/{print substr($1,13) , $4}' /work1/malser/mohammedalser/illumina-ERR15909551-results/FASTR-minimap2-NoQality-original-illumina.SAM | head -n 20


awk -F'\t' '!/^@/{print $1 , $4}' /work1/malser/mohammedalser/illumina-ERR15909551-results/FASTR-minimap2-modified-illumina.SAM | head -n 20


cmp original.txt modified.txt




#Original minimap2 - Illumina
sbatch -p mi2101x -N 1 -n 16  -t 12:00:00 -o "/work1/malser/mohammedalser/illumina-ERR15909551-results/FASTR-minimap2-original-illumina-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/illumina-ERR15909551-results/FASTR-minimap2-original-illumina-stats-time-memory.txt" /work1/malser/mohammedalser/minimap2-original/minimap2 -a -x sr -t 16 /work1/malser/mohammedalser/GCF_000001405.40_GRCh38.p14_genomic.fna /work1/malser/mohammedalser/ERR15909551.fastq -o /work1/malser/mohammedalser/illumina-ERR15909551-results/FASTR-minimap2-original-illumina.SAM"

#Original minimap2 - HiFi
sbatch -p mi2101x -N 1 -n 16  -t 12:00:00 -o "/work1/malser/mohammedalser/hifi-ERR13491966-results/FASTR-minimap2-original-hifi-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/hifi-ERR13491966-results/FASTR-minimap2-original-hifi-stats-time-memory.txt" /work1/malser/mohammedalser/minimap2-original/minimap2 -a -x map-hifi -t 16 /work1/malser/mohammedalser/GCF_000001405.40_GRCh38.p14_genomic.fna /work1/malser/mohammedalser/ERR13491966.fastq -o /work1/malser/mohammedalser/hifi-ERR13491966-results/FASTR-minimap2-original-hifi.SAM"

#Original minimap2 - ONT
sbatch -p mi2101x -N 1 -n 16  -t 12:00:00 -o "/work1/malser/mohammedalser/ont-SRR33464820-results/FASTR-minimap2-original-ont-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/ont-SRR33464820-results/FASTR-minimap2-original-ont-stats-time-memory.txt" /work1/malser/mohammedalser/minimap2-original/minimap2 -a -x map-ont -t 16 /work1/malser/mohammedalser/GCF_000001405.40_GRCh38.p14_genomic.fna /work1/malser/mohammedalser/SRR33464820.fastq -o /work1/malser/mohammedalser/ont-SRR33464820-results/FASTR-minimap2-original-ont.SAM"



#SAM decom minimap2 - Illumina
sbatch -p mi2101x -N 1 -n 16  -t 12:00:00 -o "/work1/malser/mohammedalser/illumina-ERR15909551-results/FASTR-minimap2-original-illumina-SAM_decom-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/illumina-ERR15909551-results/FASTR-minimap2-original-illumina-SAM_decom-stats-time-memory.txt" samtools fastq -@ 16 /work1/malser/mohammedalser/illumina-ERR15909551-results/FASTR-minimap2-original-illumina.SAM > /work1/malser/mohammedalser/illumina-ERR15909551-results/ERR15909551-SAM_decom.fastq"

#SAM decom minimap2 - HiFi
sbatch -p mi2101x -N 1 -n 16  -t 12:00:00 -o "/work1/malser/mohammedalser/hifi-ERR13491966-results/FASTR-minimap2-original-hifi-SAM_decom-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/hifi-ERR13491966-results/FASTR-minimap2-original-hifi-SAM_decom-stats-time-memory.txt" samtools fastq -@ 16 /work1/malser/mohammedalser/hifi-ERR13491966-results/FASTR-minimap2-original-hifi.SAM > /work1/malser/mohammedalser/hifi-ERR13491966-results/ERR13491966-SAM_decom.fastq"

#SAM decom minimap2 - ONT
sbatch -p mi2101x -N 1 -n 16  -t 12:00:00 -o "/work1/malser/mohammedalser/ont-SRR33464820-results/FASTR-minimap2-original-ont-SAM_decom-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/ont-SRR33464820-results/FASTR-minimap2-original-ont-SAM_decom-stats-time-memory.txt" samtools fastq -@ 16 /work1/malser/mohammedalser/ont-SRR33464820-results/FASTR-minimap2-original-ont.SAM > /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820-SAM_decom.fastq"








#BAM minimap2 - Illumina
sbatch -p mi2101x -N 1 -n 16  -t 12:00:00 -o "/work1/malser/mohammedalser/illumina-ERR15909551-results/FASTR-minimap2-original-illumina-BAM-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/illumina-ERR15909551-results/FASTR-minimap2-original-illumina-BAM-stats-time-memory.txt" samtools view -bS -@ 16 /work1/malser/mohammedalser/illumina-ERR15909551-results/FASTR-minimap2-original-illumina.SAM -o /work1/malser/mohammedalser/illumina-ERR15909551-results/FASTR-minimap2-original-illumina.BAM"

#BAM minimap2 - HiFi
sbatch -p mi2101x -N 1 -n 16  -t 12:00:00 -o "/work1/malser/mohammedalser/hifi-ERR13491966-results/FASTR-minimap2-original-hifi-BAM-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/hifi-ERR13491966-results/FASTR-minimap2-original-hifi-BAM-stats-time-memory.txt" samtools view -bS -@ 16 /work1/malser/mohammedalser/hifi-ERR13491966-results/FASTR-minimap2-original-hifi.SAM -o /work1/malser/mohammedalser/hifi-ERR13491966-results/FASTR-minimap2-original-hifi.BAM"

#BAM minimap2 - ONT
sbatch -p mi2101x -N 1 -n 16  -t 12:00:00 -o "/work1/malser/mohammedalser/ont-SRR33464820-results/FASTR-minimap2-original-ont-BAM-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/ont-SRR33464820-results/FASTR-minimap2-original-ont-BAM-stats-time-memory.txt" samtools view -bS -@ 16 /work1/malser/mohammedalser/ont-SRR33464820-results/FASTR-minimap2-original-ont.SAM -o /work1/malser/mohammedalser/ont-SRR33464820-results/FASTR-minimap2-original-ont.BAM"





#BAM decom minimap2 - Illumina
sbatch -p mi2101x -N 1 -n 16  -t 12:00:00 -o "/work1/malser/mohammedalser/illumina-ERR15909551-results/FASTR-minimap2-original-illumina-BAM_decom-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/illumina-ERR15909551-results/FASTR-minimap2-original-illumina-BAM_decom-stats-time-memory.txt" samtools fastq -@ 16 /work1/malser/mohammedalser/illumina-ERR15909551-results/FASTR-minimap2-original-illumina.BAM > /work1/malser/mohammedalser/illumina-ERR15909551-results/ERR15909551-BAM_decom.fastq"

#BAM decom minimap2 - HiFi
sbatch -p mi2101x -N 1 -n 16  -t 12:00:00 -o "/work1/malser/mohammedalser/hifi-ERR13491966-results/FASTR-minimap2-original-hifi-BAM_decom-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/hifi-ERR13491966-results/FASTR-minimap2-original-hifi-BAM_decom-stats-time-memory.txt" samtools fastq -@ 16 /work1/malser/mohammedalser/hifi-ERR13491966-results/FASTR-minimap2-original-hifi.BAM > /work1/malser/mohammedalser/hifi-ERR13491966-results/ERR13491966-BAM_decom.fastq"

#BAM decom minimap2 - ONT
sbatch -p mi2101x -N 1 -n 16  -t 12:00:00 -o "/work1/malser/mohammedalser/ont-SRR33464820-results/FASTR-minimap2-original-ont-BAM_decom-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/ont-SRR33464820-results/FASTR-minimap2-original-ont-BAM_decom-stats-time-memory.txt" samtools fastq -@ 16 /work1/malser/mohammedalser/ont-SRR33464820-results/FASTR-minimap2-original-ont.BAM > /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820-BAM_decom.fastq"





#CRAM minimap2 - Illumina
sbatch -p mi2101x -N 1 -n 16  -t 12:00:00 -o "/work1/malser/mohammedalser/illumina-ERR15909551-results/FASTR-minimap2-original-illumina-CRAM-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/illumina-ERR15909551-results/FASTR-minimap2-original-illumina-CRAM-stats-time-memory.txt" samtools view -@ 16 -C -T /work1/malser/mohammedalser/GCF_000001405.40_GRCh38.p14_genomic.fna -o /work1/malser/mohammedalser/illumina-ERR15909551-results/FASTR-minimap2-original-illumina.CRAM /work1/malser/mohammedalser/illumina-ERR15909551-results/FASTR-minimap2-original-illumina.BAM"

#CRAM minimap2 - HiFi
sbatch -p mi2101x -N 1 -n 16  -t 12:00:00 -o "/work1/malser/mohammedalser/hifi-ERR13491966-results/FASTR-minimap2-original-hifi-CRAM-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/hifi-ERR13491966-results/FASTR-minimap2-original-hifi-CRAM-stats-time-memory.txt" samtools view -@ 16 -C -T /work1/malser/mohammedalser/GCF_000001405.40_GRCh38.p14_genomic.fna -o /work1/malser/mohammedalser/hifi-ERR13491966-results/FASTR-minimap2-original-hifi.CRAM /work1/malser/mohammedalser/hifi-ERR13491966-results/FASTR-minimap2-original-hifi.BAM"

#CRAM minimap2 - ONT
sbatch -p mi2101x -N 1 -n 16  -t 12:00:00 -o "/work1/malser/mohammedalser/ont-SRR33464820-results/FASTR-minimap2-original-ont-CRAM-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/ont-SRR33464820-results/FASTR-minimap2-original-ont-CRAM-stats-time-memory.txt" samtools view -@ 16 -C -T /work1/malser/mohammedalser/GCF_000001405.40_GRCh38.p14_genomic.fna -o /work1/malser/mohammedalser/ont-SRR33464820-results/FASTR-minimap2-original-ont.CRAM /work1/malser/mohammedalser/ont-SRR33464820-results/FASTR-minimap2-original-ont.BAM"






#CRAM decom minimap2 - Illumina
sbatch -p mi2101x -N 1 -n 16  -t 12:00:00 -o "/work1/malser/mohammedalser/illumina-ERR15909551-results/FASTR-minimap2-original-illumina-CRAM_decom-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/illumina-ERR15909551-results/FASTR-minimap2-original-illumina-CRAM_decom-stats-time-memory.txt" samtools fastq -@ 16 -T /work1/malser/mohammedalser/GCF_000001405.40_GRCh38.p14_genomic.fna /work1/malser/mohammedalser/illumina-ERR15909551-results/FASTR-minimap2-original-illumina.CRAM > /work1/malser/mohammedalser/illumina-ERR15909551-results/FASTR-minimap2-original-illumina-CRAM_decom.fastq"

#CRAM decom minimap2 - HiFi
sbatch -p mi2101x -N 1 -n 16  -t 12:00:00 -o "/work1/malser/mohammedalser/hifi-ERR13491966-results/FASTR-minimap2-original-hifi-CRAM_decom-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/hifi-ERR13491966-results/FASTR-minimap2-original-hifi-CRAM_decom-stats-time-memory.txt" samtools fastq -@ 16 -T /work1/malser/mohammedalser/GCF_000001405.40_GRCh38.p14_genomic.fna /work1/malser/mohammedalser/hifi-ERR13491966-results/FASTR-minimap2-original-hifi.CRAM > /work1/malser/mohammedalser/hifi-ERR13491966-results/FASTR-minimap2-original-hifi-CRAM_decom.fastq"

#CRAM decom minimap2 - ONT
sbatch -p mi2101x -N 1 -n 16  -t 12:00:00 -o "/work1/malser/mohammedalser/ont-SRR33464820-results/FASTR-minimap2-original-ont-CRAM_decom-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/ont-SRR33464820-results/FASTR-minimap2-original-ont-CRAM_decom-stats-time-memory.txt" samtools fastq -@ 16 -T /work1/malser/mohammedalser/GCF_000001405.40_GRCh38.p14_genomic.fna /work1/malser/mohammedalser/ont-SRR33464820-results/FASTR-minimap2-original-ont.CRAM > /work1/malser/mohammedalser/ont-SRR33464820-results/FASTR-minimap2-original-ont-CRAM_decom.fastq"




