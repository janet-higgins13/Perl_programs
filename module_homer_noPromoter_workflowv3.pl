#!/usr/bin/perl -w
use strict;

# script to run homer on genes in each module
# use peaks which are not found in promoter region
# use control peaks as background
#bsub -o outtest -q Test128 -J perl "perl ~/bin/module_homer_noPromoter_workflowv3.pl"

#source homer-4.6
#source blat-35
#source samtools-0.1.19
#source weblogo-2.8.2
#source ghostscript-9.10


# 07.07.2014  janet.higgins@tgac.ac.uk
my %module = ();
#my %newlist = ();
my %search = ();
my $gene = ();
#my $Sgene = ();
my $col = ();
my $ann = ();
my @list = ();

# for genes down reg at 144 and modules containing downreg genes
#@list = `ls *Cgene_list.txt*`;
#@list = `ls midnightblue_module_Cgene_list.txt*`;
#$ann = "control24.Zpeaks.txt";

# for genes up reg at 144 and modules containing upreg genes
@list = `ls *Tgene_list.txt*`;
$ann = "treated144.Zpeaks.txt";

foreach my $filename(@list) {	
my $name = $filename;
chomp $filename;
$filename =~ s/_list.txt//;
open (GENE, "<$name") or die "Couldn't open gene list file($!)\n";
    print "$filename\n";
    while (<GENE>) {
   		my $gene = $_;
   		chomp $gene;
   		#print "$gene\n";
 		push(@{$module{$filename}},$gene); 
 		}	
}

my $size = keys%module;
print "The number of modules is $size\n";

# this produces a peak.bed file for each module using annotated FAIRE peaks
foreach my $col(keys %module) {	
	print "$col\n";
	open (ANN, "<$ann") or die "Couldn't open annotation file($!)\n";
	open (OUTPEAK, ">$col\_peak.bed") or die "Couldn't open output bed file($!)\n";
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
  		if (defined $gene) {
  			my $TSS = $ann[9];
  			if ($TSS < 100000 && $TSS > -100000) {
  		#if ($sig > 0.9) {
  				my @genelist = @{$module{$col}};
  				my %search = map { $_ => 1 } @genelist;
					if(exists($search{$gene})) {
 					print OUTPEAK "$chr\t$start\t$end\t$id\t$sig\t$strand\n";	
 						
 				}
			}
		}
	}
}


# this produces a peak.bed file for each module using annotated FAIRE peaks using only peaks with a zinba value of >0.9
foreach my $col(keys %module) {	
	print "$col\n";
	open (ANN, "<$ann") or die "Couldn't open annotation file($!)\n";
	open (OUTPEAK, ">$col\_9peak.bed") or die "Couldn't open output bed file($!)\n";
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
  		if (defined $gene) {
  			my $TSS = $ann[9];
  			if ($TSS < 100000 && $TSS > -100000) {
  				if ($sig > 0.9) {
  					my @genelist = @{$module{$col}};
  					my %search = map { $_ => 1 } @genelist;
						if(exists($search{$gene})) {
 						print OUTPEAK "$chr\t$start\t$end\t$id\t$sig\t$strand\n";							
 					}
				}
			}
		}
	}
}

# this produces a peak.bed file for each module using annotated FAIRE peaks removing peaks in promoter regions
foreach my $col(keys %module) {	
	print "$col\n";
	open (ANN, "<$ann") or die "Couldn't open annotation file($!)\n";
	open (OUTPEAK, ">$col\_noPro_peak.bed") or die "Couldn't open output bed file($!)\n";
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
    	my $pro    =   $ann[7];
    	next if ($pro =~ /^promoter-TSS/);
  		my $gene = $ann[15];
  		if (defined $gene) {
  			my $TSS = $ann[9];
  			if ($TSS < 100000 && $TSS > -100000) {
  					my @genelist = @{$module{$col}};
  					my %search = map { $_ => 1 } @genelist;
						if(exists($search{$gene})) {
 						print OUTPEAK "$chr\t$start\t$end\t$id\t$sig\t$strand\n";							
				}
			}
		}
	}
}


# this produces a peak.bed file for each module using annotated FAIRE peaks removing peaks in promoter regions
foreach my $col(keys %module) {	
	print "$col\n";
	open (ANN, "<$ann") or die "Couldn't open annotation file($!)\n";
	open (OUTPEAK, ">$col\_noPro_50_peak.bed") or die "Couldn't open output bed file($!)\n";
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
    	my $pro    =   $ann[7];
    	next if ($pro =~ /^promoter-TSS/);
  		my $gene = $ann[15];
  		if (defined $gene) {
  			my $TSS = $ann[9];
  			if ($TSS < 50000 && $TSS > -50000) {
  					my @genelist = @{$module{$col}};
  					my %search = map { $_ => 1 } @genelist;
						if(exists($search{$gene})) {
 						print OUTPEAK "$chr\t$start\t$end\t$id\t$sig\t$strand\n";							
				}
			}
		}
	}
}

# this produces a peak.bed file for each module using annotated FAIRE peaks removing peaks in promoter regions
foreach my $col(keys %module) {	
	print "$col\n";
	open (ANN, "<$ann") or die "Couldn't open annotation file($!)\n";
	open (OUTPEAK, ">$col\_noPro_25_peak.bed") or die "Couldn't open output bed file($!)\n";
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
    	my $pro    =   $ann[7];
    	next if ($pro =~ /^promoter-TSS/);
  		my $gene = $ann[15];
  		if (defined $gene) {
  			my $TSS = $ann[9];
  			if ($TSS < 25000 && $TSS > -25000) {
  					my @genelist = @{$module{$col}};
  					my %search = map { $_ => 1 } @genelist;
						if(exists($search{$gene})) {
 						print OUTPEAK "$chr\t$start\t$end\t$id\t$sig\t$strand\n";							
				}
			}
		}
	}
}

exit;