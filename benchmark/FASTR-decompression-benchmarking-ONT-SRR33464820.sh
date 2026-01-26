#!/bin/bash

echo "=============================================="
echo "FASTR Decompression Benchmarking"
echo "=============================================="


#XZ DECOMPRESSION
sbatch -p mi2101x -N 1 -n 16 -t 12:00:00 -o "/work1/malser/mohammedalser/ont-SRR33464820-results/xz_fast_decom-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/ont-SRR33464820-results/xz_fast_decom-stats-time-memory.txt" xz -d -k -T 16 -c /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820.xz_fast.xz > /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820.xz_fast_decom.fastq"

sbatch -p mi2101x -N 1 -n 16 -t 12:00:00 -o "/work1/malser/mohammedalser/ont-SRR33464820-results/xz_best_decom-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/ont-SRR33464820-results/xz_best_decom-stats-time-memory.txt" xz -d -k -T 16 -c /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820.xz_best.xz > /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820.xz_best_decom.fastq"


#PIGZ DECOMPRESSION
sbatch -p mi2101x -N 1 -n 16 -t 12:00:00 -o "/work1/malser/mohammedalser/ont-SRR33464820-results/pigz_fast_decom-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/ont-SRR33464820-results/pigz_fast_decom-stats-time-memory.txt" pigz -d -k -p 16 -c /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820.pigz_fast.gz > /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820.pigz_fast_decom.fastq"

sbatch -p mi2101x -N 1 -n 16 -t 12:00:00 -o "/work1/malser/mohammedalser/ont-SRR33464820-results/pigz_best_decom-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/ont-SRR33464820-results/pigz_best_decom-stats-time-memory.txt" pigz -d -k -p 16 -c /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820.pigz_best.gz > /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820.pigz_best_decom.fastq"

sbatch -p mi2101x -N 1 -n 16 -t 12:00:00 -o "/work1/malser/mohammedalser/ont-SRR33464820-results/pigz_zopfli_decom-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/ont-SRR33464820-results/pigz_zopfli_decom-stats-time-memory.txt" pigz -d -k -p 16 -c  /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820.pigz_zopfli.gz > /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820.pigz_zopfli_decom.fastq"


# SPRING DECOMPRESSION
sbatch -p mi2101x -N 1 -n 16 -t 12:00:00 -o "/work1/malser/mohammedalser/ont-SRR33464820-results/spring_decom-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/ont-SRR33464820-results/spring_decom-stats-time-memory.txt" spring -d -l -t 16 -i /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820.spring.gz -o /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820.spring_decom.fastq -w /work1/malser/mohammedalser/ont-SRR33464820-results/tmp"


# RENANO DECOMPRESSION
sbatch -p mi2101x -N 1 -n 16 -t 12:00:00 -o "/work1/malser/mohammedalser/ont-SRR33464820-results/renano_decom-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/ont-SRR33464820-results/renano_decom-stats-time-memory.txt" renano -d -t 16 /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820.renano_refless.gz /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820.renano_decom.fastq"


#FASTR Mode 0 DECOMPRESSION
sbatch -p mi2101x -N 1 -n 16 -t 12:00:00 -o "/work1/malser/mohammedalser/ont-SRR33464820-results/fastr_mode0_decom-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/ont-SRR33464820-results/fastr_mode0_decom-stats-time-memory.txt" python /work1/malser/mohammedalser/FASTR-main/src/toFASTQ.py /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820.fastr_mode0.fastr /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820.fastr_mode0_decom.fastq --mode 0 --num_workers 16 --phred_alpha phred94"


#FASTR Mode 1 DECOMPRESSION
sbatch -p mi2101x -N 1 -n 16 -t 12:00:00 -o "/work1/malser/mohammedalser/ont-SRR33464820-results/fastr_mode1_decom-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/ont-SRR33464820-results/fastr_mode1_decom-stats-time-memory.txt" python /work1/malser/mohammedalser/FASTR-main/src/toFASTQ.py /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820.fastr_mode1.fastr /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820.fastr_mode1_decom.fastq --mode 1 --num_workers 16 --phred_alpha phred94"


#FASTR Mode 2 DECOMPRESSION
sbatch -p mi2101x -N 1 -n 16 -t 12:00:00 -o "/work1/malser/mohammedalser/ont-SRR33464820-results/fastr_mode2_decom-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/ont-SRR33464820-results/fastr_mode2_decom-stats-time-memory.txt" python /work1/malser/mohammedalser/FASTR-main/src/toFASTQ.py /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820.fastr_mode2.fastr /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820.fastr_mode2_decom.fastq --mode 2 --num_workers 16 --phred_alpha phred94"


#FASTR Mode 3 DECOMPRESSION
sbatch -p mi2101x -N 1 -n 16 -t 12:00:00 -o "/work1/malser/mohammedalser/ont-SRR33464820-results/fastr_mode3_decom-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/ont-SRR33464820-results/fastr_mode3_decom-stats-time-memory.txt" python /work1/malser/mohammedalser/FASTR-main/src/toFASTQ.py /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820.fastr_mode3.fastr /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820.fastr_mode3_decom.fastq --mode 3 --num_workers 16 --phred_alpha phred94 --headers_file /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820.fastr_mode3_headers.txt"


#XZ DECOMPRESSION on FASTR MODE2
sbatch -p mi2101x -N 1 -n 16 -t 12:00:00 -o "/work1/malser/mohammedalser/ont-SRR33464820-results/xz_fast_fastr_mode2_decom-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/ont-SRR33464820-results/xz_fast_fastr_mode2_decom-stats-time-memory.txt" xz -d -k -T 16 -c /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820.xz_fast_fastr_mode2.fastr.xz > /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820.xz_fast_fastr_mode2_decom.fastr"

sbatch -p mi2101x -N 1 -n 16 -t 12:00:00 -o "/work1/malser/mohammedalser/ont-SRR33464820-results/xz_best_fastr_mode2_decom-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/ont-SRR33464820-results/xz_best_fastr_mode2_decom-stats-time-memory.txt" xz -d -k -T 16 -c /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820.xz_best_fastr_mode2.fastr.xz > /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820.xz_best_fastr_mode2_decom.fastr"


#PIGZ DECOMPRESSION on FASTR MODE2
sbatch -p mi2101x -N 1 -n 16 -t 12:00:00 -o "/work1/malser/mohammedalser/ont-SRR33464820-results/pigz_fast_fastr_mode2_decom-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/ont-SRR33464820-results/pigz_fast_fastr_mode2_decom-stats-time-memory.txt" pigz -d -k -p 16 -c /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820.pigz_fast_fastr_mode2.gz > /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820.pigz_fast_fastr_mode2_decom.fastr"

sbatch -p mi2101x -N 1 -n 16 -t 12:00:00 -o "/work1/malser/mohammedalser/ont-SRR33464820-results/pigz_best_fastr_mode2_decom-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/ont-SRR33464820-results/pigz_best_fastr_mode2_decom-stats-time-memory.txt" pigz -d -k -p 16 -c /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820.pigz_best_fastr_mode2.gz > /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820.pigz_best_fastr_mode2_decom.fastr"

sbatch -p mi2101x -N 1 -n 16 -t 12:00:00 -o "/work1/malser/mohammedalser/ont-SRR33464820-results/pigz_zopfli_fastr_mode2_decom-stats-command-out.txt" --wrap="\time -v --output "/work1/malser/mohammedalser/ont-SRR33464820-results/pigz_zopfli_fastr_mode2_decom-stats-time-memory.txt" pigz -d -k -p 16 -c  /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820.pigz_zopfli_fastr_mode2.gz > /work1/malser/mohammedalser/ont-SRR33464820-results/SRR33464820.pigz_zopfli_fastr_mode2_decom.fastr"
