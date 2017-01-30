#!/usr/bin/perl -w

#use diagnostics;
use strict;
use Data::Dumper;

my $infile;
my $outfile;

if(@ARGV != 1 ){
        die "Usage: $0 <InputFile>";
}
if( @ARGV == 1 ){
        $infile = shift;
}

my $infile_prefix = $infile;$infile_prefix=~s/(.*)\..*/$1/; 
$outfile = $infile_prefix . ".Geno_out";
my @header_array = ();

open( INPUT, "$infile" ) or die("can't open $infile!");
open( OUTPUT, ">$outfile" ) or die("can't write to $outfile!");

#Hashes for rs ID
my %rs_hash;

while(my $line = <INPUT>){
	chomp $line;
	if($. == 1 ){
		@header_array = split( "\t", $line );

		print OUTPUT $header_array[0] . "\t" . $header_array[1] . "\t" . $header_array[2] . "\t" . $header_array[3] . "\t0_genotype\t1_genotype\t2_genotype";
		for(my $i=4 ; $i < scalar(@header_array) ; $i++){
			print OUTPUT "\t" . $header_array[$i];
		}
		next;
	}
	my @line_content = split( "\t", $line );

	my $SNP_id = $line_content[0];
	my $chromo = $line_content[1];
	my $pos_on_chromo = $line_content[2];
	my $direction = $line_content[3];


	my %genotype_hash;
	for(my $i=4 ; $i < scalar(@line_content); $i++){
		foreach my $j(split("/",$line_content[$i])){
			if($j ne "I" && $j ne "D" && $j ne "-" && $j ne "N"){
				$genotype_hash{"$j"} +=1;
			}
		}
	}

	#Define the most frequent allele
	my $zero;
	my $one; 
	my $two;
	my $first_key = "";
	if( scalar(keys %genotype_hash) > 1 ){
		foreach my $key (sort keys %genotype_hash){
			if($key ne "N" && $key ne "-" && $key ne "I" && $key ne "D"){
				if($first_key eq ""){
					$first_key = $key;
				}
				else{
					if($genotype_hash{$key} > $genotype_hash{$first_key}){
						$zero = "$key/$key";
						$one = "$first_key/$key";
						$two = "$first_key/$first_key";
					}
					elsif($genotype_hash{$key} < $genotype_hash{$first_key}){
						$zero = "$first_key/$first_key";
	                                	$one = "$first_key/$key";
	        	                        $two = "$key/$key";
					}
					#If they are equals take the first one
					else{
						$zero = "$first_key/$first_key";
                	                	$one = "$first_key/$key";
	                        	        $two = "$key/$key";
					}	
				}
			}
		}
	}
	elsif(scalar(keys %genotype_hash) == 1 ){
		my @tmp = keys %genotype_hash;
		if($tmp[0]=~/[ACGT]/){
			$zero = $one = $two = $tmp[0]."/".$tmp[0];
		}
		else{
			$zero = $one = $two = "NA";
		}
	}
	else{
		$zero = $one = $two = "NA";
	}
	#Make the output
	print OUTPUT "\n";
	print OUTPUT "$SNP_id\t$chromo\t$pos_on_chromo\t$direction\t$zero\t$one\t$two";

	for(my $i=4 ; $i < scalar(@line_content); $i++){
		my @alleles = split("/",$line_content[$i]);
                my $geno1 = $alleles[0] . "/" . $alleles[1];
		my $geno2 = $alleles[1] . "/" . $alleles[0];
		

                print OUTPUT "\t0" if( ($geno1 eq $zero || $geno2 eq $zero) && ($geno1 ne $one && $geno2 ne $one) && ($geno1 ne $two && $geno2 ne $two) );
		print OUTPUT "\t0" if ($geno1 eq $zero && $geno1 eq $one && $geno1 eq $two);
		print OUTPUT "\t1" if( ($geno1 eq $one || $geno2 eq $one) && ($geno1 ne $zero && $geno2 ne $zero) && ($geno1 ne $two && $geno2 ne $two) );
		print OUTPUT "\t2" if( ($geno1 eq $two || $geno2 eq $two) && ($geno1 ne $one && $geno2 ne $one) && ($geno1 ne $zero && $geno2 ne $zero) );
		print OUTPUT "\tNA" if($line_content[$i] eq "N/N" || $line_content[$i] eq "-/-");
		print OUTPUT "\tNA" if($line_content[$i] eq "I/I" || $line_content[$i] eq "D/D" || $line_content[$i] eq "I/D" || $line_content[$i] eq "D/I");
        }
}

close INPUT;
close OUTPUT;
