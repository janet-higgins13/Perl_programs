#!/usr/bin/perl -w
use strict;

# script to find  out how many genes in a list are in each WGCNA module 
# perl ~/bin/hub_distributionv1.pl tophubs_signed.txt Sarah_significanGenes.txt hubs.txt


# 01.05.2014  janet.higgins@tgac.ac.uk
my %module = ();
my %newlist = ();
my %search = ();
my $gene = ();
my $Sgene = ();
my $col = ();

# Check arguments
unless (@ARGV == 3) {
  print "Usage: module.pl <module_colour.txt> < diff_exp_gene.txt> <out.txt> \n";
  exit;
}

# Open FHs
open (MOD, "<$ARGV[0]") or die "Couldn't open input module list file($!)\n";
open (SIGGENES, "<$ARGV[1]") or die "Couldn't open input gene file($!)\n";
open (OUT, ">$ARGV[2]") or die "Couldn't open input bed file($!)\n";

while (<MOD>) {
    my $l = $_;
 	chomp $l;
 	my @ann = split /\t/,$l; 	
 	my $gene = $ann[0];
 	my $col = $ann[1];
 	$gene =~ s/"//g;
 	$col =~ s/"//g;
 	#print "$col\t$gene\n";
 	push(@{$module{$col}},$gene); 
 	push(@{$newlist{$col}},$gene); 
}

my $size = keys%module;
print "The number of modules is $size\n";

# This prints out list and number of genes in each module
foreach my $col(keys %module) {	
		#print "$col\t";
		my $geneno = scalar(@{$module{$col}});
		#print out number of genes in each module
		#print "$geneno \n";
		#print out list of genes in each module
		#print "@{$module{$col}}\n";
}	

foreach my $col(keys %newlist) {	
   @{$newlist{$col}} = undef;

}


   
while (<SIGGENES>) {
	my $l = $_;
	chomp $l;
	my @gann = split /\t/,$l; 	
 	my $expgene = $gann[1];   
 	#for Sarah list
 	#my $expgene = $gann[1];  
   		foreach my $col(keys %module) {	
		
	#foreach my $element(@genelist) {	
		my @genelist = @{$module{$col}};
  		my %search = map { $_ => 1 } @genelist;
		if(exists($search{$expgene})) {
		#print "$col\t$expgene\n";
 		 push(@{$newlist{$col}},$expgene);  
 		} 
 	}
}


foreach my $col(keys %newlist) {	
		print OUT "$col\t";
		my $geneno = scalar(@{$newlist{$col}});
		#print out number of genes in each module
		my $num = $geneno -1;
		print OUT "$num \t";
		#print out list of genes in each module
		print OUT "@{$newlist{$col}}\n";

}




exit;