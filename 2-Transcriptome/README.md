

## 1. Software and database required

**Software:**

- Cutadapt (v2.5): https://cutadapt.readthedocs.io/en/stable/
- Hisat2 (v2.1.0): https://daehwankimlab.github.io/hisat2/
- RSEM (v1.3.3): https://deweylab.github.io/RSEM/
- DESeq2 (v1.20): https://bioconductor.org/packages/release/bioc/html/DESeq2.html

**Database:**

- RSEM reference database（GRch38）https://www.gencodegenes.org/human/

## 2. Data preparation

**Create directories and copy raw sequencing files：**

```shell
mkdir 00_rawdata 01_fastp 02_cutadapt 03_rsem
mv *.fastq 00_rawdata
cd 00_rawdata/
ls *_1.fastq.gz > ../filelist
cd ..
```

## 3. Quality filtering

**Adapter sequence detection：**

```shell
for i in `cat filelist`
do
	j=${i/_1./_2.}
	fastp \
	--detect_adapter_for_pe \
	-i 00_rawdata/$i \
	-I 00_rawdata/$j \
	-o 01_fastp/$i \
	-O 01_fastp/$j \
	-z 4 -q 20 -u 40 -n 6
done
```

**Remove adapter：**

```shell
for i in `cat filelist`
do
	j=${i/_1./_2.}
	cutadapt \
	--discard-trimmed \
	-O 18 \
	-j 20 \
	-b AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
	-B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
	-o 02_cutadapt/$i \
	-p 02_cutadapt/$j \
	00_rawdata/$i \
	00_rawdata/$j 
done
```

## 4. Sequence alignment and gene expression calculation

**Install RSEM:**

```shell
wget https://github.com/deweylab/RSEM/archive/v1.3.3.tar.gz
tar -zxvf v1.3.3.tar.gz
cd rsem
make install prefix=/bigdata/gaojy/biosoft/rsem/ 
echo 'PATH=$PATH:/bigdata/gaojy/biosoft/rsem/bin/' >> ~/.bashrc
source ~/.bashrc
```

**Build index：**

```shell
rsem-prepare-reference \
-p 20 \
--gtf \
./gencode.v39.annotation.gtf \
--hisat2-hca \
./GRCh38.primary_assembly.genome.fa \
./grch38
```

**Align RNA-seq against human genome and calculate gene expression:**

```shell
for i in `cat filelist`
do
	j=${i/_1./_2.}
	k=${i%%_*}
	rsem-calculate-expression \
	--paired-end \
	-p 20 \
	--hisat2-hca \
	02_cutadapt/$i \
	02_cutadapt/$j \
	/bigdata/gaojy/biosoft/rsem/database/grch38/grch38 \
	03_rsem/$k
done
```

The count and FPKM data of all genes in each sample should be present in [SampleID].genes.results

Concatenate the results of each sample to generate the file count.txt

## 5. Variance stabilizing transformation of count data

```R
data<-read.table("count.txt",row.names=1,sep="\t",header=T)
group<-read.table("group.txt",sep="\t",header=T) ## optional metadata
dds <- DESeqDataSetFromMatrix(countData = data, colData = group, design= ~ Group)
vst<-assay(varianceStabilizingTransformation(dds))
write.table(vst,"transcriptome.txt",append=FALSE,sep="\t",quote=FALSE)
```

