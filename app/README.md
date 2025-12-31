# Using FASTR in minimap2:
We integrated a FASTR interpreter in the latest minimap2 (version: 2.30-r1290-dirty, cloned on Dec 12, 2025) as a proof-of-concept of how easy to adapt our new file format.
```
make
./minimap2 -x sr GCF_000001405.40_GRCh38.p14_genomic.fna.gz fastr-mode2-illumina-utdir_1.fastr -o minimap2-mapping-FASTR-illumina.SAM
```
