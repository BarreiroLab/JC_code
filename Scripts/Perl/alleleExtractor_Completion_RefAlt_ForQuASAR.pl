#!/usr/bin/perl -w

#use diagnostics;
use strict;
use Getopt::Long qw (:config gnu_getopt);

my $infile;
my $pos;


my $parms = GetOptions (
        "infile=s" => \$infile,
        "pos=s"=>\$pos,
) or die();

my %errors = (
        usg => "Usage: # $0 --infile=<file.alleleCount.txt> --pos=<chr_pos_ref.file.txt>\n\n"
 #      nofile => "Please specify a good entry file\n",
);



if(!$infile && !$pos ){
        DIE ($errors{"usg"});
}

my $infile_prefix = $infile;$infile_prefix=~s/(.*)\.\w{1,}/$1/;


open( POS, "$pos" ) or die("can't open $pos!\n");

my @positions;
my %withAlt;
while(my $line = <POS>){
        chomp $line;
	my @line_content = split(/[\t\s]/,$line);
	my $chr = $line_content[0];
	my $pos = $line_content[1];
	my $ref = $line_content[2];
	my $alt = $line_content[3];
	my $id = $line_content[4];
	my $maf = $line_content[5];
	push(@positions,"$chr;$pos;$ref");
	$withAlt{"$chr;$pos;$ref"}="$alt;$id;$maf";
}
close POS;


open( OUTFILE, ">$infile_prefix.alleleCountCorrectedRefAlt.txt" ) or dieWithUnexpectedError("can't write to $infile_prefix.alleleCountCorrectedRefAlt.txt!");

open( INPUT, "$infile") or die("can't open $infile!\n");
my $index = 0;
while(my $line = <INPUT>){
	chomp $line;
        my @line_content = split(/\t/,$line);
        my $chr = $line_content[0];
        my $pos = $line_content[1];
	my $ref = $line_content[2];
	my $bases = $line_content[3];
        my @infos = split(/;/,$withAlt{"$chr;$pos;$ref"});               
	my $alt = $infos[0];
        my $id = $infos[1];
        my $maf = $infos[2];


	if($positions[$index] ne "$chr;$pos;$ref"){
		while($positions[$index] ne "$chr;$pos;$ref"){
			my @posRef = split(/;/,$positions[$index]);
		 	my @infos = split(/;/,$withAlt{"$chr;$pos;$ref"});
        	        my $alt = $infos[0];
	                my $id = $infos[1];
                	my $maf = $infos[2];
			my $postmp = $posRef[1];
			#do not print for QuASAR
			#print OUTFILE $posRef[0] . "\t" . ($postmp-1) . "\t$postmp\t$ref\t$alt\t$id\t$maf\t0\t0\t0\n";
			$index +=1;
		}
	}
	if($positions[$index] eq "$chr;$pos;$ref"){
		#GET COUNT TOTAL :
		my @b = split(/:/,$bases);
		my @counts = split(/,/,$b[1]);
		my $ref_count;
		my $alt_count;
		my $total_count;
		my $other_count;
		my @alleles= ("A","C","G","T","N");
		my %counts_Allele = (
			"A" => $counts[0],
			"C" => $counts[1],
			"G" => $counts[2],
			"T" => $counts[3],
			"N" => $counts[4],
		);
		
		$ref_count = $counts_Allele{$ref};
		$alt_count = $counts_Allele{$alt};
		
		$total_count = $counts_Allele{"A"}+$counts_Allele{"C"}+$counts_Allele{"G"}+$counts_Allele{"T"}+$counts_Allele{"N"};
		$other_count = $total_count - ($ref_count+$alt_count);
		
		if(($ref_count+$alt_count>0) && $maf > 0 ){
			print OUTFILE $chr . "\t" . ($pos-1) . "\t" . $pos . "\t$ref\t$alt\t$id\t$maf\t$ref_count\t$alt_count\t$other_count\n";
		}
	}
	$index +=1;
}
if($index < scalar(@positions)){
	for(my $j = $index ; $j < scalar(@positions) ; $j+=1 ){
		my @posRef = split(/;/,$positions[$j]);
		my $chr = $posRef[0];
		my $pos = $posRef[1];
		my $ref = $posRef[2];
	        my @infos = split(/;/,$withAlt{"$chr;$pos;$ref"});
        	my $alt = $infos[0];
	        my $id = $infos[1];
        	my $maf = $infos[2];
		#Do not print for QuASAR
		#print OUTFILE $posRef[0] . "\t" . ($posRef[1]-1) . "\t" . $posRef[1] . "\t$ref\t$alt\t$id\t$maf\t0\t0\t0\n";
	}
}

close INPUT;
close OUTFILE;
sub DIE {
        my ($msg) = @_;

        TRACE ($msg, 0);
        exit (1);
}
sub TRACE {
        my ($msg, $level) = @_;
        print $msg;
}
