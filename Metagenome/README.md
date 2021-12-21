## 1. Data preparation

Install sunbeam (https://sunbeam.readthedocs.io/en/latest/), create directories, and copy raw sequencing files

```shell
conda activate sunbeam
mkdir 00_rawdata 01_cutadapt 02_trimmomatic 03_komplexity 04_decontam 05_taxonomy 06_assembly 07_annotation 08_binning 09_dimred 10_diffabund
mv *.fastq 00_rawdata
cd 00_rawdata/
ls *_1.fastq.gz > ../filelist
cd ..
```

## 2. Adaptor trimming

```shell
for i in `cat filelist`
	do
	j = ${i/_1./_2.}
	cutadapt --discard-trimmed -O 17 -b GTTTCCCAGTCACGATC -b GTTTCCCAGTCACGATCNNNNNNNNNGTTTCCCAGTCACGATC -B GTTTCCCAGTCACGATC -B GTTTCCCAGTCACGATCNNNNNNNNNGTTTCCCAGTCACGATC -o 01_cutadapt/$i -p 01_cutadapt/$j  00_rawdata/$i 00_rawdata/$j --quiet -j 5  
	done
```

## 3. Quality trimming

```shell
mkdir 02_trimmomatic/unpaired

for i in `cat filelist`
	do 
	j = ${i/_1./_2.}
	trimmomatic PE -threads 5 -phred33 01_cutadapt/$i 01_cutadapt/$j 02_trimmomatic/$i 02_trimmomatic/unpaired/$i 02_trimmomatic/$j 02_trimmomatic/unpaired/$j ILLUMINACLIP:/bigdata/yangjh/miniconda3/envs/sunbeam/share/trimmomatic-0.36-6/adapters/NexteraPE-PE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	done
```

## 4. Low complexity reads removal

```shell
for i in `cat filelist`
	do 
	j = ${i/_1./_2.}
	k = ${i%%_*}
	for rp in 02_trimmomatic/$i 02_trimmomatic/$j; do gzip -dc $rp |kz |awk '{if ($4<0.55) print $1}'>> 03_komplexity/filtered_ids/${k}_0.55; done 
	done
	
for i in `cat filelist`
	do 
	j=${i/_1./_2.}
	k=${i%%_*}
	gzip -dc 02_trimmomatic/$i  |rbt fastq-filter 03_komplexity/filtered_ids/${k}_0.55 |gzip > 03_komplexity/$i 
	gzip -dc 02_trimmomatic/$j  |rbt fastq-filter 03_komplexity/filtered_ids/${k}_0.55 |gzip > 03_komplexity/$j 
	done
```

## 5. Decontamination

Download human reference genome hg38 (https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/)

Download the script decontam.py from this github folder to your current working directory

```shell
mkdir 04_decontam/bam

for i in `cat filelist`
	do 
	j = ${i/_1./_2.}
	k = ${i%%_*}
	bwa mem -M -t 5 /software/databases/hg38/hg38.fasta 03_komplexity/$i 03_komplexity/$j |samtools view -bSF4 - > 04_decontam/bam/${k}.bam
	done

cd 04_decontam/bam
for i in *.bam
	do
	python ../../decontam.py $i > ${i}.hostlist
	done

cd ../..

for i in `cat filelist`
    do
	j = ${i/_1./_2.}
	k = ${i%%_*}
	gzip -dc 03_komplexity/$i |rbt fastq-filter 04_decontam/bam/${k}.bam.hostlist |gzip > 04_decontam/$i 
	gzip -dc 03_komplexity/$j |rbt fastq-filter 04_decontam/bam/${k}.bam.hostlist |gzip > 04_decontam/$j
	done
```

## 6. Kraken2 taxonomic annotation

Download and build standard kraken2 database (https://github.com/DerrickWood/kraken2/wiki/Manual)

```shell
for i in `cat filelist`
	do 
	j = ${i/_1./_2.}
	k = ${i%%_*}
	kraken2 --db /software/databases/kraken/kraken2_standard --report 05_taxonomy/${k}.taxa.tsv --use-mpa-style --use-name --thread 1 --paired 04_decontam/$i 04_decontam/$j  > /dev/null 
	done
```

## 7. Megahit assembly

Co-assembly for samples whose clean sequencing data are < 500M

```shell
gzip 04_decontam/*.fastq.gz

for i in `cat filelist |sed -e 's/.gz//'`
	do
	j = ${i/_1./_2.}
	k = ${i%_*}
	sum = `du -sm ${k}*|awk '{sum += $1}END{print sum}'`
    	if [$sum -lt 500]
    	then
    		echo $i >> assemblelist
    	else
    		cat $i >> 04_decontam/remain_1.fastq
    		cat $j >> 04_decontam/remain_2.fastq
    	fi
    done
echo "remain_1.fastq" >> assemblelist
cd ..

for i in `cat assemblelist`
	do 
	j = ${i/_1./_2.}
	k = ${i%%_*}
	megahit -l 04_decontam/$i 04_decontam/$j -o 06_assembly/${k} -t 20
	done
```

Use a customized script to rename, combine and select contigs greater than 300bp

For example, save below perl script as rename_contigs.pl, which use bioperl to filter sequences below 300bp

```perl
#!bin/perl

use strict;
use lib "/home/wangzhang/perl5/lib/perl5/";
use Bio::SeqIO;

open (IN, "assemblelist");
while (<IN>) {
	my $input = ($_ =~ /(\S+)_1.fastq/);
	my $seq_in = Bio::SeqIO -> new( -format => 'fasta', -file => "06_assembly/$input/final.contigs.fa"); 
	my $seq_out = Bio::SeqIO -> new (-format => 'fasta', -file => ">06_assembly/$input/final.contigs.rename.300bp.fa");
	my %length = (); 
	my %sequence = ();
	while (my $seq = $seq_in -> next_seq()) { #read each sequence into $seq and push them into array
		my $seqid = $seq -> id;
		$seqid = $input."-".$seqid;
		my $seqlength = $seq -> length;
		my $sequence = $seq -> seq;
		$length{$seqid} = $seqlength;
		$sequence{$seqid} = $sequence;
	}
	for my $key (sort {$length{$b} <=> $length{$a}} keys %length) {
		next if ($length{$key} <= 300);
		my $seq2 = new Bio::Seq ('-id' => $key, '-seq' => $sequence{$key});
		$seq_out -> write_seq ($seq2);
	}
}
```

```shell
perl rename_contigs.pl 
cat 06_assembly/*/final.contigs.rename.300bp.fa > 06_assembly/all_contigs.300bp.fa
```

## 8. Functional annotation

CD-HIT for de-redundancy, prodigal for ORF identification, diamond for alignment to KEGG database

