#!bin/perl
## This script was used to extract metabolite-CIDm and CIDm-host target mapping information from STITCH database ##

open (IN, "chemical.sources.v5.0.tsv"); #### stitch id mapping file ####
while (<IN>) {
	chop;
	next if (/^#/);
	s/CID\S//g;
	@a=split("\t",$_);
	if (/PC\t(\d+)/) {
		$pc{$1}="CIDm".$a[0];
	}
	if (/PS\t(\d+)/) {
		$ps{$1}="CIDm".$a[0];
	}
	$merge{$a[1]}="CIDm".$a[0];
	if (/CHEBI:(\d+)/) {
		$chebi{$1}="CIDm".$a[0];
	}
	if (/KEGG/) {
		$kegg{$a[3]}="CIDm".$a[0];
	}
}

open (IN1, "metabolite_IDs.txt"); ### metabolite IDs list ###
open (OUT1, ">metabo2CIDm.txt"); ## output ##
while (<IN>) {
	chop;
	@a=split("\t",$_);
	if ($chebi{$a[5]} ne $merge{$a[4]}) {
		print STDERR "chebi and pubchem IDs for $a[0] do not match\n"; ### raise a flag if IDs do not match based on different mapping rules
	}
	else {
		print OUT1 $a[0]."\t".$ps{$a[3]}."\t".$merge{$a[3]}."\t".$pc{$a[3]}."\t".$chebi{$a[4]}."\t".$kegg{$a[0]}."\t".$kegg{$a[5]}."\n";
	}	
}

open (IN, "9606.protein_chemical.links.detailed.v5.0.tsv");
while (<IN>) {
	chop;
	s/9606\.//g;
	@a=split("\t",$_);
	next if ($a[6]<700);
	#$hash{$a[0]}{$a[1]}=$a[5];
	$score{$a[0]}{$a[1]}=$a[2]."\t".$a[3]."\t".$a[4]."\t".$a[5]."\t".$a[6];
}
open (IN, "9606.actions.v5.0.tsv");
while (<IN>) {
	chop;
	s/9606\.//g;
	@a=split("\t",$_);
	next if ($a[5]<700);
	$int{$a[0]}{$a[1]}=$a[2];
	$int{$a[1]}{$a[0]}=$a[2];
}

open (IN, "human_gene_ids.txt");
while (<IN>) {
	chop;
	@a=split("\t",$_);
	$geneid{$a[2]}=$a[9]."\t".$a[8];
}

open (IN2, "all_CIDm.txt"); #### stitch CIDm compound ID list #####
open (OUT2, ">all_CIDm_targets.txt"); ## output ##
while (<IN>) {
	chop;
#	@a=split("\t",$_);
	$id=$_;
	if (exists $int{$id}) {
		for my $key (keys %{$int{$id}}) {
			print OUT2 $id."\t".$key."\t".$geneid{$key}."\t".$int{$id}{$key}."\t".$score{$id}{$key}."\n";
		}
	}
}

