#!/usr/bin/perl -w
use strict;

# script to make peak.bed for homer peak calling for peaks associated with list of genes in each module from WGCNA
# txt file output from Excel will not work
# cat -v testmodcol.txt | sed s/'\^M//g'  this does not work
# save as Unix(LF) in TextWrangler 


# perl ~/bin/select_module_faire_cutoffv1.pl gene_moduleColors.txt only_T144high_vs_C24.txt outbed2


# 24.04.2014  janet.higgins@tgac.ac.uk
my %module = ();
my @genelist = ();
#my %transID3 = ();

# Check arguments
unless (@ARGV == 3) {
  print "Usage: select_module_faire.pl <module_list.txt> < homer_annotation.txt > <output peak.bed> \n";
  exit;
}

# Open FHs
open (MOD, "<$ARGV[0]") or die "Couldn't open input module list file($!)\n";
open (ANN, "<$ARGV[1]") or die "Couldn't open annotation file($!)\n";
open (OUT, ">$ARGV[2]") or die "Couldn't open input bed file($!)\n";

while (<MOD>) {
 	push(@genelist,$_);  
}


while (<ANN>) {
    my $line = $_;
 	chomp $line;
 	next if ($line =~ /^PeakID/);
 	my @ann = split /\t/,$line; 	
  	my $chr = $ann[1];
  	my $id = $ann[0];
  	my $start = $ann[2];
  	my $end = $ann[3];
  	my $sig   =   $ann[5];
    my $strand    =   $ann[4];
  	my $gene = $ann[15];
  	my $TSS = $ann[9];
  	
  	if ($TSS < 100000 && $TSS > -100000) {
  	if ($sig > 0.9) {
  	
  	print "@genelist\n";
  	my %search = map { $_ => 1 } @genelist;

		if(exists($search{$gene})) {
		#print "$gene\n";
		#print "$gene\t$TSS\t$sig\n";
 		print OUT "$chr\t$start\t$end\t$id\t$sig\t$strand\n";	
 }
}
}
}

exit;