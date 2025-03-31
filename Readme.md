# Personal Genomes
This package creates personal genomes for alignment of samples with known genetic variation, converts coordiantes of files aligned to the personal genome back to GRCh38, and counts alleles of the known variants.

## Dependencies
Personal Genomes Requires
```
- python (>=3.9.12)
- numpy (>=1.20.3)
- pysam (>=0.15.3)
```

## Usage

1) Use `write_genome_and_all_vcf.py` to create a personal genome.
2) Align fasta files to the personal genome.
3) Use `convert_to_hg38.py` to convert SAM file from personal genome coordinates to GRCh38 coordinates.
4) Sort reads in the GRCh SAM file. 
5) Use `count_alleles.py` to count alleles of known variants in the sorted SAM file.
