# FASTQ Compression Benchmark
Reproducible scripts for benchmarking compression and decompression of FASTQ files across general-purpose, fastq-specific and fastq-alternative tools.


---

## Create Conda environment with all tools that we installed:

```
conda env create -f environment.yml
```

## Activate the environment

```
source ~/anaconda3/bin/activate root

conda activate fastq_bench
```
## Download the data
```
fasterq-dump SRR33464820

fasterq-dump ERR15909551

fasterq-dump ERR13491966
```

## Download Human genome
https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/

## Run compression tools
Use the sbatch jobs inside the 3 compression files, be careful with dependencies (e.g., FASTR mode 2 must be generated before running pigz on FASTR mode 2)

## Run decompression tools
Use the sbatch jobs inside the 3 decompression files, be careful with dependencies (e.g., pigz on FASTR mode 2 must be generated before running pigz decompression)

## Run minimap2, SAM, BAM, CRAM tools
Use the minimap2 shell script and be careful with dependencies

## Collect all results
To build a nice summary table of all results:
```

cd /work1/malser/mohammedalser/ont-SRR33464820-results/
awk -F': ' 'BEGIN{print "Filename,Total_CPU_Sec,Peak_Mem_KB"} /User time/{u=$2} /System time/{s=$2} /Maximum resident/{print FILENAME "," u+s "," $2}' *.txt > ONT-SRR33464820-performance_summary.csv

echo "Filename,SizeBytes" > ONT-SRR33464820-size_summary.csv && find . -maxdepth 1 -type f ! -name "*.txt" -printf "%f,%s\n" >> ONT-SRR33464820-size_summary.csv

python /work1/malser/mohammedalser/FASTR-main/merger.py ONT-SRR33464820-size_summary.csv ONT-SRR33464820-performance_summary.csv ONT-SRR33464820-results.csv SRR33464820 ont
```

Enjoy it!


