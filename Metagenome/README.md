## 1. Software and database required

Software:

- Sunbeam (v2.1): https://sunbeam.readthedocs.io/en/latest/
- Kraken 2 (v2.0.8): https://ccb.jhu.edu/software/kraken2/
- Megahit (v1.1.3): https://github.com/voutcn/megahit
- Prodigal (v2.6.3): https://github.com/hyattpd/Prodigal
- NCBI-Blast (v2.9.0): https://ftp.ncbi.nlm.nih.gov/blast/executables/
- CD-HIT (v4.8.1): http://weizhong-lab.ucsd.edu/cd-hit/
- BBMap (v38.44): https://sourceforge.net/projects/bbmap/
- MetaWRAP (v1.2.1): https://github.com/bxlab/metaWRAP
- Bioperl module (v5.26.2): https://bioperl.org/



Database: 

- Human genome database (hg38): https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/
- Kraken 2 standard database: https://benlangmead.github.io/aws-indexes/k2
- KEGG protein sequence database: https://www.genome.jp/dbget/
- 6,530 representative bacterial genome database (manually curated, available upon request)

## 2. Data preparation

Create directories, and copy raw sequencing files

```shell
conda activate sunbeam
mkdir 00_rawdata 01_cutadapt 02_trimmomatic 03_komplexity 04_decontam 05_taxonomy 06_assembly 07_annotation 08_abundcalc 09_binning
mv *.fastq 00_rawdata
cd 00_rawdata/
ls *_1.fastq.gz > ../filelist
cd ..
```

## 3. Adaptor trimming

```shell
for i in `cat filelist`
	do
	j = ${i/_1./_2.}
	cutadapt --discard-trimmed -O 17 -b GTTTCCCAGTCACGATC -b GTTTCCCAGTCACGATCNNNNNNNNNGTTTCCCAGTCACGATC -B GTTTCCCAGTCACGATC -B GTTTCCCAGTCACGATCNNNNNNNNNGTTTCCCAGTCACGATC -o 01_cutadapt/$i -p 01_cutadapt/$j  00_rawdata/$i 00_rawdata/$j --quiet -j 5  
	done
```

## 4. Quality trimming

```shell
mkdir 02_trimmomatic/unpaired

for i in `cat filelist`
	do 
	j = ${i/_1./_2.}
	trimmomatic PE -threads 5 -phred33 01_cutadapt/$i 01_cutadapt/$j 02_trimmomatic/$i 02_trimmomatic/unpaired/$i 02_trimmomatic/$j 02_trimmomatic/unpaired/$j ILLUMINACLIP:/bigdata/yangjh/miniconda3/envs/sunbeam/share/trimmomatic-0.36-6/adapters/NexteraPE-PE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	done
```

## 5. Low complexity reads removal

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

## 6. Decontamination

Download human reference genome hg38

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

## 7. Kraken2 taxonomic annotation

Download and build standard kraken2 database

Kraken2 taxonomic annotation (read-based)

```shell
for i in `cat filelist`
	do 
	j = ${i/_1./_2.}
	k = ${i%%_*}
	kraken2 --db /software/databases/kraken/kraken2_standard --report 05_taxonomy/${k}.taxa.tsv --use-mpa-style --use-name --thread 1 --paired 04_decontam/$i 04_decontam/$j  > /dev/null 
	done
```

Use below perl script to merge Kraken 2 output to a matrix file

```perl
#!bin/perl

use strict;
open (IN, "filelist");
while (<IN>) {
	chop;
	my $id = ($_ =~ /(\S+)_1.fastq/);
	open (IN1, "05_taxonomy/$id-taxa.tsv");
	my $dump=<IN1>;
	while (<IN1>) {
		chop;
		next if (/d__Eukaryota/);
		next unless (/s__/);
		@a=split("\t",$_);
		$taxa{$id}{$a[0]}=$a[1];
		$trans_taxa{$a[0]}{$id}=$a[1];
	}
}
open (OUT, ">05_taxonomy/taxonomy.txt");
print OUT "#SampleID";
for my $key (sort keys %taxa) {
	print OUT "\t$key";
}
print OUT "\n";
for my $key (sort keys %trans_taxa) {
	print OUT $key;
	for my $key2 (sort keys %taxa) {
		if (exists $trans_taxa{$key}{$key2}) {
			print OUT "\t$trans_taxa{$key}{$key2}";
		}
		else {
			print OUT "\t0.0";
		}
	}
	print OUT "\n";
}
```

The file generated in this step was uploaded as taxonomy.txt in this github folder

## 8. Megahit assembly

Uni-assembly for samples whose are data >= 500M, co-assembly for samples whose data are < 500M

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

for i in `cat assemblelist`
	do 
	j = ${i/_1./_2.}
	k = ${i%%_*}
	megahit -l 04_decontam/$i 04_decontam/$j -o 06_assembly/${k} -t 20
	done
```

## 9. Functional annotation

prodigal for ORF identification

```shell
for i in `cat assemblelist`
	do 
	j = ${i/_1./_2.}
	k = ${i%%_*}
	prodigal -i ${k}/final.contigs.fa -a ${k}/proteins.faa -d ${k}/nucleotides.fna -o ${k}/coord.gff -f gff -p meta -m
	done
```

Use a customized script to rename, combine and select ORF nucleotide sequences greater than 300bp

For example, save below perl script as rename_contigs.pl, which use bioperl to filter sequences below 300bp

```perl
#!bin/perl

use strict;
use lib "/home/wangzhang/perl5/lib/perl5/";
use Bio::SeqIO;

open (IN, "assemblelist");
while (<IN>) {
	my $input = ($_ =~ /(\S+)_1.fastq/);
	my $seq_in = Bio::SeqIO -> new( -format => 'fasta', -file => "06_assembly/$input/nucleotides.fna"); 
	my $seq_out = Bio::SeqIO -> new (-format => 'fasta', -file => ">06_assembly/$input/nucleotides.rename.300bp.fa");
	my %length = (); 
	my %sequence = ();
	while (my $seq = $seq_in -> next_seq()) {
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
cat 06_assembly/*/nucleotides.rename.300bp.fa > 07_annotation/all.orf.300bp.fa
```

CD-HIT to generate non-redundant gene catalogue

```shell
cd-hit-est -i 07_annotation/all.orf.300bp.fa -o 07_annotation/all.orf.300bp.fa.95.90 -c 0.95 -n 8 -M 10200 -aS 0.9 -T 20 -d 0 -G 0 -g 1
```

functional annotation by alignment to KEGG database

```shell
diamond blastx -p 10 -d /software/databases/KEGG_database/KEGG_database_all_protein.dmnd -q 07_annotation/all.orf.300bp.fa.95.90 -e 1e-4 -k 1 --sensitive -o 07_annotation/all.orf.300bp.fa.95.90.kegg.annotation.out -f 6 qseqid sseqid stitle pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp
```

The functional annotation file is saved as: all.orf.300bp.fa.95.90.kegg.annotation.out

taxonomic annotation (gene-based) by BLASTn to a manually curated bacterial genome database (6,530 representative genomes, available upon request)

```shell
blastn -query 07_annotation/all.orf.300bp.fa.95.90 -db /software/databases/bacterial_genomes/6530set.all.fa -outfmt 6 -out 07_annotation/all.orf.300bp.6530set_cov0.8_id0.65.btab -max_target_seqs 99999 -perc_identity 0.65 -qcov_hsp_perc 0.8 -evalue 0.01
```

Download 6530set_taxonomy.txt and save it and below perl script as gene_level_taxa_assign.pl in 07_annotation folder

```perl
#!bin/perl
use strict;

my %count=();
my %species=();
my %phylum=();
my %count=();
my %subcount=();
my %hit=();

open (IN, "6530set_taxonomy.txt");
while (<IN>) {
	chop;
	my @a=split("\t",$_);
	$species{$a[0]}=$a[2];
	$genus{$a[0]}=$a[3];
	$phylum{$a[0]}=$a[4];
}

open (IN, "all.orf.300bp.6530set_cov0.8_id0.65.btab");
while (<IN>) {
	chop;
	my @a=split("\t",$_);
	$count{$a[0]} ++;
	($genome)=($a[1] =~ /(GCA_\d+\.\d)/);
	if (! exists $hit{$a[0]}{$genome} or $a[11]>$hit{$a[0]}{$genome}) {
		$hit{$a[0]}{$genome}=$a[11];
	}
}
for my $key (keys %count) {
	my $consensus=sprintf('%.0f',$count{$key}*0.5);
	if ($consensus < 1) {
		$subcount{$a[0]}=1;
	}
	else {
		$subcount{$a[0]}=$consensus;
	}

}
open (OUT, ">all.orf.300bp.taxa.annotation.out");
for my $key (keys %hit) {
	my $line=1;
	my %allspecies=();
	my %allgenus=();
	my %allphylum=();
	for my $key2 (sort {$hit{$key}{$b}<=>$hit{$key}{$a}} keys %{$hit{$key}}) {
		last if ($line>$subcount{$key});
		$allspecies{$species{$key2}}++;
		$allgenus{$genus{$key2}}++;
		$allphylum{$phylum{$key2}}++;
		$line++;
	}
	@species=sort {$allspecies{$b}<=>$allspecies{$a}} keys %allspecies;
	@genus=sort {$allgenus{$b}<=>$allgenus{$a}} keys %allgenus;
	@phylum=sort {$allphylum{$b}<=>$allphylum{$a}} keys %allphylum;
	if ($allspecies{$species[0]}/$subcount{$key}>=0.95) {
		print OUT $key."\t"."Species"."\t".$species[0]."\n";
	}
	elsif ($allgenus{$genus[0]}/$subcount{$key}>=0.85) {
		print OUT $key."\t"."Genus"."\t".$genus[0]."\n";
	}
	elsif ($allphylum{$phylum[0]}/$subcount{$key}>=0.65) {
		print OUT $key."\t"."Phylum"."\t".$phylum[0]."\n";
	}
	else {
		print OUT $key."\tUnclassified\n";
	}
}
```

```shell
perl gene_level_taxa_assign.pl
```

The gene-level taxonomic annotation was saved as all.orf.300bp.taxa.annotation.out

## 10. Gene abundance calculation

Mapping clean reads of each sample to the non-redundant gene catalogue

Downsize all bam to at most 3M reads to adjust for sequencing depth

```shell
for i in `cat filelist`
	do 
	j = ${i/_1./_2.}
	k = ${i%%_*}
	/software/bbmap/bbmap.sh ref=07_annotation/all.orf.300bp.fa.95.90 in=04_decontam/$i in2=04_decontam/$j out=08_abundcalc/${k}.bam k=13 minid=0.90 thread=20 nodisk=true rpkm=08_abundcalc/${k}.rpkm
	samtools view -F 4 08_abundcalc/${k}.bam -o 08_abundcalc/${k}.sam
	grep '^\@' 08_abundcalc/${k}.sam 08_abundcalc/${k}_3M.sam
	grep -v '^\@' 08_abundcalc/${k}.sam 08_abundcalc/${k}_noheader.sam
	/software/sample-master/sample -k 3000000 -o -p 08_abundcalc/${k}_nohead.sam >> 08_abundcalc/${k}_3M.sam
	samtools view -bS 08_abundcalc/${k}_3M.sam -o 08_abundcalc/${k}_3M.bam
	rm 08_abundcalc/*sam
	samtools sort 08_abundcalc/${k}_3M.bam -o 08_abundcalc/${k}.sorted.bam -@ 10
	/software/jgi_summarize_bam_contig_depths --outputDepth 08_abundcalc/${k}.depth.txt 08_abundcalc/${k}.sorted.bam
	done
```

Aggregate gene depth file to KEGG orthologues using below perl script

```perl
#!bin/perl
use strict;
my %ko=();
my %depth=();
my %sampleid=();

open (IN, "07_annotation/all.orf.300bp.fa.95.90.kegg.annotation.out");
while (<IN>) {
	chop;
	my @a=split("\t",$_);
	if ($a[2] =~ /(K\d+)/) {
		($ko{$a[0]}) = ($a[2] =~ /(K\d+)/);
		$ko{$ko{$a[0]}}=1;
	}
}
close IN;

open (IN, "filelist");
while (<IN>) {
	chop;
	my $id=$_;
	$sampleid{$id}=1;
	open (IN1, "08_abundcalc/$id.depth.txt");
	my $dump=<IN1>;
	while (<IN1>) {
		chop;
		s/^(\S+)[^\t]{0,}\t/$1\t/g;
		my @a=split("\t",$_);
		if (exists $ko{$a[0]}) {
			$depth{$ko{$a[0]}}{$id}+=$a[2];
		}
	}
}

open (OUT, ">metagenome.txt");
print OUT "#SampleID";
for my $key (sort keys %sampleid) {
	print OUT "\t$key";
}
print OUT "\n";
for my $key (sort keys %ko) {
	print OUT $key;
	for my $key2 (sort keys %sampleid) {
		if (exists $depth{$key}{$key2}) {
			print OUT "\t$depth{$key}{$key2}";
		}
		else {
			print OUT "\t0.0";
		}
	}
	print OUT "\n";
}
```

The file generated in this step was uploaded as metagenome.txt in this github folder.

## 11. Binning

Binning using metawrap pipeline

```shell
mkdir 04_decontam/remain 
mv 04_decontam/remain*.fastq 04_decontam/remain/ ## move the redundant remain_1.fastq remain_2.fastq to a separate folder
cat 06_assembly/*/final_contigs.fa > 06_assembly/all_contigs.fa
metawrap binning -o 09_binning/binning -t 10 -a 06_assembly/all_contigs.fa --metabat2 --maxbin2 --concoct 04_decontam/*.fastq
metawrap bin_refinement -o 09_binning/bin_refine -t 10 -A 09_binning/binning/metabat2_bins/ -B 09_binning/binning/maxbin2_bins/ -C 09_binning/binning/concoct_bins/ -c 0 -x 100 ## remove contamination/completeness cutoff to retain as many refined bins as possible
mkdir 09_binning/bin_refine/metawrap_50_10_bins/
cat 09_binning/bin_refine/metawrap_0_100_bins.stats |gawk '$2>=50 && $3<=10 {print "cp 09_binning/bin_refine/metawrap_0_100_bins/"$1".fa 09_binning/bin_refine/metawrap_50_10_bins/."}' ## store completeness >= 50% & contamination <= 10% bins to a separate folder for further analyses
metawrap blobology -a 06_assembly/all_contigs.fa -t 30 -o 09_binning/blobology --bins 09_binning/bin_refine/metawrap_50_10_bins 04_decontam/*fastq
metawrap quant_bins -b 09_binning/bin_refine/metawrap_50_10_bins -o 09_binning/quant_bins -a 06_assembly/all_contigs.fa 04_decontam/*fastq
metawrap reassemble_bins -o 09_binning/reassembly -1 04_decontam/*1.fastq -2 04_decontam/*2.fastq -t 30 -m 800 -c 50 -x 10 -b 09_binning/bin_refine/metawrap_50_10_bins ##check to see whether bin stats improve after re-assembly
metawrap classify_bins -b 09_binning/bin_refine/metawrap_50_10_bins -o 09_binning/orig_taxa -t 30 ## original bins
metawrap classify_bins -b 09_binning/reassembly/reassembled_bins -o 09_binning/reassemble_taxa -t 30 ## reassembled bins
```
