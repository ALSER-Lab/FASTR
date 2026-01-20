#!/bin/bash

echo "=============================================="
echo "FASTR Benchmarking"
echo "=============================================="


#source ~/anaconda3/bin/activate root

#conda activate fastq_bench


#XZ COMPRESSION
sbatch -p mi2101x -N 1 -n 16 -t 12:00:00 -o "/work1/malser/mohammedalser/illumina-ERR15909551-results/xz_fast-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/illumina-ERR15909551-results/xz_fast-stats-time-memory.txt" xz -1 -k -T 16 -c /work1/malser/mohammedalser/ERR15909551.fastq > /work1/malser/mohammedalser/illumina-ERR15909551-results/ERR15909551.xz_fast.xz"

sbatch -p mi2101x -N 1 -n 16 -t 12:00:00 -o "/work1/malser/mohammedalser/illumina-ERR15909551-results/xz_best-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/illumina-ERR15909551-results/xz_best-stats-time-memory.txt" xz -9 -k -T 16 -c /work1/malser/mohammedalser/ERR15909551.fastq > /work1/malser/mohammedalser/illumina-ERR15909551-results/ERR15909551.xz_best.xz"


#PIGZ COMPRESSION
sbatch -p mi2101x -N 1 -n 16 -t 12:00:00 -o "/work1/malser/mohammedalser/illumina-ERR15909551-results/pigz_fast-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/illumina-ERR15909551-results/pigz_fast-stats-time-memory.txt" pigz -1 -k -p 16 -c /work1/malser/mohammedalser/ERR15909551.fastq > /work1/malser/mohammedalser/illumina-ERR15909551-results/ERR15909551.pigz_fast.gz"

sbatch -p mi2101x -N 1 -n 16 -t 12:00:00 -o "/work1/malser/mohammedalser/illumina-ERR15909551-results/pigz_best-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/illumina-ERR15909551-results/pigz_best-stats-time-memory.txt" pigz -9 -k -p 16 -c /work1/malser/mohammedalser/ERR15909551.fastq > /work1/malser/mohammedalser/illumina-ERR15909551-results/ERR15909551.pigz_best.gz"

sbatch -p mi2101x -N 1 -n 16 -t 12:00:00 -o "/work1/malser/mohammedalser/illumina-ERR15909551-results/pigz_zopfli-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/illumina-ERR15909551-results/pigz_zopfli-stats-time-memory.txt" pigz -11 -k -p 16 -c /work1/malser/mohammedalser/ERR15909551.fastq > /work1/malser/mohammedalser/illumina-ERR15909551-results/ERR15909551.pigz_zopfli.gz"


# SPRING COMPRESSION
sbatch -p mi2101x -N 1 -n 16 -t 12:00:00 -o "/work1/malser/mohammedalser/illumina-ERR15909551-results/spring-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/illumina-ERR15909551-results/spring-stats-time-memory.txt" spring -c -t 16 -i /work1/malser/mohammedalser/ERR15909551.fastq -o /work1/malser/mohammedalser/illumina-ERR15909551-results/ERR15909551.spring.gz -w /work1/malser/mohammedalser/illumina-ERR15909551-results/tmp"


# RENANO COMPRESSION
sbatch -p mi2101x -N 1 -n 16 -t 12:00:00 -o "/work1/malser/mohammedalser/illumina-ERR15909551-results/renano_refless-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/illumina-ERR15909551-results/renano_refless-stats-time-memory.txt" renano -t 16 /work1/malser/mohammedalser/ERR15909551.fastq /work1/malser/mohammedalser/illumina-ERR15909551-results/ERR15909551.renano_refless.gz"


#FASTR Mode 0 COMPRESSION
sbatch -p mi2101x -N 1 -n 16 -t 12:00:00 -o "/work1/malser/mohammedalser/illumina-ERR15909551-results/fastr_mode0-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/illumina-ERR15909551-results/fastr_mode0-stats-time-memory.txt" python /work1/malser/mohammedalser/FASTR-main/src/to_fastr.py /work1/malser/mohammedalser/ERR15909551.fastq /work1/malser/mohammedalser/illumina-ERR15909551-results/ERR15909551.fastr_mode0.fastr --mode 0 --qual_scale log --seq_type illumina_sra --workers 16 --phred_alpha phred94"


#FASTR Mode 1 COMPRESSION
sbatch -p mi2101x -N 1 -n 16 -t 12:00:00 -o "/work1/malser/mohammedalser/illumina-ERR15909551-results/fastr_mode1-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/illumina-ERR15909551-results/fastr_mode1-stats-time-memory.txt" python /work1/malser/mohammedalser/FASTR-main/src/to_fastr.py /work1/malser/mohammedalser/ERR15909551.fastq /work1/malser/mohammedalser/illumina-ERR15909551-results/ERR15909551.fastr_mode1.fastr --mode 1 --qual_scale log --seq_type illumina_sra --workers 16 --phred_alpha phred94"

#FASTR Mode 2 COMPRESSION
sbatch -p mi2101x -N 1 -n 16 -t 12:00:00 -o "/work1/malser/mohammedalser/illumina-ERR15909551-results/fastr_mode2-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/illumina-ERR15909551-results/fastr_mode2-stats-time-memory.txt" python /work1/malser/mohammedalser/FASTR-main/src/to_fastr.py /work1/malser/mohammedalser/ERR15909551.fastq /work1/malser/mohammedalser/illumina-ERR15909551-results/ERR15909551.fastr_mode2.fastr --mode 2 --qual_scale log --seq_type illumina_sra --workers 16 --phred_alpha phred94"


#FASTR Mode 3 COMPRESSION
sbatch -p mi2101x -N 1 -n 16 -t 12:00:00 -o "/work1/malser/mohammedalser/illumina-ERR15909551-results/fastr_mode3-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/illumina-ERR15909551-results/fastr_mode3-stats-time-memory.txt" python /work1/malser/mohammedalser/FASTR-main/src/to_fastr.py /work1/malser/mohammedalser/ERR15909551.fastq /work1/malser/mohammedalser/illumina-ERR15909551-results/ERR15909551.fastr_mode3.fastr --mode 3 --qual_scale log --seq_type illumina_sra --workers 16 --phred_alpha phred94"


#XZ COMPRESSION on FASTR MODE2
sbatch -p mi2101x -N 1 -n 16 -t 12:00:00 -o "/work1/malser/mohammedalser/illumina-ERR15909551-results/xz_fast_fastr_mode2-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/illumina-ERR15909551-results/xz_fast_fastr_mode2-stats-time-memory.txt" xz -1 -k -T 16 -c /work1/malser/mohammedalser/illumina-ERR15909551-results/ERR15909551.fastr_mode2.fastr > /work1/malser/mohammedalser/illumina-ERR15909551-results/ERR15909551.xz_fast_fastr_mode2.fastr.xz"

sbatch -p mi2101x -N 1 -n 16 -t 12:00:00 -o "/work1/malser/mohammedalser/illumina-ERR15909551-results/xz_best_fastr_mode2-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/illumina-ERR15909551-results/xz_best_fastr_mode2-stats-time-memory.txt" xz -9 -k -T 16 -c /work1/malser/mohammedalser/illumina-ERR15909551-results/ERR15909551.fastr_mode2.fastr > /work1/malser/mohammedalser/illumina-ERR15909551-results/ERR15909551.xz_best_fastr_mode2.fastr.xz"


#PIGZ COMPRESSION on FASTR MODE2
sbatch -p mi2101x -N 1 -n 16 -t 12:00:00 -o "/work1/malser/mohammedalser/illumina-ERR15909551-results/pigz_fast_fastr_mode2-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/illumina-ERR15909551-results/pigz_fast_fastr_mode2-stats-time-memory.txt" pigz -1 -k -p 16 -c /work1/malser/mohammedalser/illumina-ERR15909551-results/ERR15909551.fastr_mode2.fastr > /work1/malser/mohammedalser/illumina-ERR15909551-results/ERR15909551.pigz_fast_fastr_mode2.gz"

sbatch -p mi2101x -N 1 -n 16 -t 12:00:00 -o "/work1/malser/mohammedalser/illumina-ERR15909551-results/pigz_best_fastr_mode2-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/illumina-ERR15909551-results/pigz_best_fastr_mode2-stats-time-memory.txt" pigz -9 -k -p 16 -c /work1/malser/mohammedalser/illumina-ERR15909551-results/ERR15909551.fastr_mode2.fastr > /work1/malser/mohammedalser/illumina-ERR15909551-results/ERR15909551.pigz_best_fastr_mode2.gz"

sbatch -p mi2101x -N 1 -n 16 -t 12:00:00 -o "/work1/malser/mohammedalser/illumina-ERR15909551-results/pigz_zopfli_fastr_mode2-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/illumina-ERR15909551-results/pigz_zopfli_fastr_mode2-stats-time-memory.txt" pigz -11 -k -p 16 -c /work1/malser/mohammedalser/illumina-ERR15909551-results/ERR15909551.fastr_mode2.fastr > /work1/malser/mohammedalser/illumina-ERR15909551-results/ERR15909551.pigz_zopfli_fastr_mode2.gz"


echo "=============================================="
echo "Total time and peak memory results"
echo "=============================================="

cd /work1/malser/mohammedalser/illumina-ERR15909551-results/
awk -F': ' 'BEGIN{print "Filename,Total_CPU_Sec,Peak_Mem_KB"} /User time/{u=$2} /System time/{s=$2} /Maximum resident/{print FILENAME "," u+s "," $2}' *.txt > Illumina-ERR15909551-performance_summary.csv

echo "Filename,SizeBytes" > Illumina-ERR15909551-size_summary.csv && find . -maxdepth 1 -type f ! -name "*.txt" -printf "%f,%s\n" >> Illumina-ERR15909551-size_summary.csv

python /work1/malser/mohammedalser/FASTR-main/merger.py Illumina-ERR15909551-size_summary.csv Illumina-ERR15909551-performance_summary.csv Illumina-ERR15909551-results.csv ERR15909551 illumina