#!bin/perl
## This script was used to extract ID mapping information between metabolites in metabolomic data and MetaCyc compounds ##


open (IN, "compounds.tab");
open (OUT1, ">metacyc_compounds_IDs_SMILEs.tab");
while (<IN>) {
	chop;
	$id=();
	$pubchem=();
	$chebi=();
	$smile=();
	$kegg=();
	$hmdb=();
	#@a=split("\t",$_);
	if ($_=~ /^UNIQUE-ID - (\S+)/) {
		$id=$1;
	}
	if ($_=~ /DBLINKS - \(PUBCHEM "(\d+)"/) {
		$pubchem=$1;
	}
	if ($_=~ /DBLINKS - \(CHEBI "(\d+)"/) {
		$chebi = $1;
	}
        if ($_=~ /DBLINKS - \(LIGAND-CPD "(C\d+)"/) {
                $kegg = $1;
        }
        if ($_=~ /DBLINKS - \(HMDB "(HMDB\d+)"/) {
                $hmdb = $1;
        }
	if ($_=~ /SMILES - (\S+)/) {
		$smile= $1;
	}
	print OUT1 $id."\t".$pubchem."\t".$chebi."\t".$kegg."\t".$hmdb."\t".$smile."\n";
}

open (IN, "metacyc_compounds_IDs_SMILEs.tab");
while (<IN>) {
	chop;
	@a=split("\t",$_);
	if ($a[3] ne '') {
		$kegg{$a[3]}=$a[0];
	}
	if ($a[4] ne '') {
		($tmp)=($a[4] =~ /HMDB(\d+)/);
		$hmdb = "HMDB00".$tmp;
		$hmdb{$hmdb}=$a[0];
	}
	if ($a[1] ne '') {
		$pubchem{$a[1]}=$a[0];
	}
	if ($a[2] ne '') {
		$chebi{$a[2]}=$a[0];
	}
}

open (IN, "cmpd_description.txt");
open (OUT2, ">cmpd2metabo_IDmatch.txt");
$dump=<IN>;
while (<IN>) {
	chop;
	my %match=();
	@a=split("\t",$_);
	if ($a[0] =~ /^C/) {
		$match{$kegg{$a[0]}}=1;
		$match{$hmdb{$a[3]}}=1;
		$match{$pubchem{$a[4]}}=1;
		$match{$chebi{$a[5]}}=1;
		$match{$kegg{$a[6]}}=1;
		#print $a[0]."\t".$kegg{$a[0]}."\t".$hmdb{$a[3]}."\t".$pubchem{$a[4]}."\t".$chebi{$a[5]}."\t".$kegg{$a[6]}."\n";
	}
	else {
		$match{$hmdb{$a[0]}}=1;
		$match{$hmdb{$a[3]}}=1;
		$match{$pubchem{$a[4]}}=1;
		$match{$chebi{$a[5]}}=1;
		$match{$kegg{$a[6]}}=1;
		#print $a[0]."\t".$hmdb{$a[0]}."\t".$hmdb{$a[3]}."\t".$pubchem{$a[4]}."\t".$chebi{$a[5]}."\t".$kegg{$a[6]}."\n";
	}
	for my $key (sort keys %match) {
		if ($key ne '') {
			print OUT2 $a[0]."\t".$key."\n";
		}
	}
}
