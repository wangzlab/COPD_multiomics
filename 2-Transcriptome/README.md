## 1. Software and database required

Software:

- Cutadapt (v2.5): https://cutadapt.readthedocs.io/en/stable/
- Hisat2 (v2.1.0): https://daehwankimlab.github.io/hisat2/
- RSEM (v1.3.1): https://deweylab.github.io/RSEM/
- DESeq2 (v1.20): https://bioconductor.org/packages/release/bioc/html/DESeq2.html

Database:

- Human genome database (hg38): https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/

## 2. Data preparation

Create directories and copy raw sequencing files

```shell
mkdir 00_rawdata 01_cutadapt 03_hisat 04_rsem 05_normalize
mv *.fastq 00_rawdata
cd 00_rawdata/
ls *_1.fastq.gz > ../filelist
cd ..
```

## 3. Quality filtering
