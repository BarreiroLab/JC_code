#!/usr/bin/perl -w

#use diagnostics;
use strict;

my $infile;
my $indivFile;
my $outfile;

if(@ARGV != 2 ){
        die "Usage: $0 <InputFile.haps> <InputFile.sample>";
}
if( @ARGV == 2 ){
        $infile = shift;
	$indivFile = shift;
}

my $infile_prefix = $infile;$infile_prefix=~s/(.*)\..*/$1/; 

my @indivTable;
print "Begin Reading Individual File\n";
open( INDIV, "$indivFile") or die("can't open $indivFile!");
while(my $line = <INDIV>){
	chomp $line;
	if($. == 1){next;}
	if($. == 2){next;}
	my @line_content = split(/[\s\t]/,$line);
	push(@indivTable,$line_content[1]);
}
close INDIV;
print "End Reading Individual File\n";

open( INPUT, "$infile" ) or die("can't open $infile!");
open( OUTPUT, ">$infile.shapeIT.inp_hapguess_switch.out") or die("can't write to $infile.shapeIT.inp_hapguess_switch.out!"); 

print OUTPUT 
"********************************************
*                                          *
*      Output from shapeIT2fastphase.pl    *
*      Code by J-C Grenier                 *
*      jean.christophe.grenier\@gmail.com  *
*                                          *
********************************************
BEGIN GENOTYPES\n";

my %rsHash;
my @rsOrder;

print "Begin Reading Haplotype File\n";

while(my $line = <INPUT>){
	chomp $line;
	my @line_content = split(/[\s\t]/,$line);
	push(@rsOrder,$line_content[0]);
	my $allelicCount = 0;
	#Determine major and minor alleles :
	for(my $i = 5 ; $i < scalar(@line_content) ; $i += 2){
                my $A = $line_content[$i];
                my $B = $line_content[$i+1];
		$allelicCount += $A + $B;
        }
	my $majorZero = 0;
	#IF ALLELE SUM >= NB_alleles/2
	## THe MAJOR ALLELE IS 1
	if($allelicCount >= scalar(@indivTable)){
		$majorZero = 0;
	}
	## THE MAJOR ALLELE IS 0
	else{
		$majorZero = 1;
	}

	for(my $i = 5 ; $i < scalar(@line_content) ; $i += 2){
		my $A;
		my $B;
		if($majorZero == 1){
			$A = $line_content[$i] eq 0 ? '1' : '0';
			$B= $line_content[$i+1] eq 0 ? '1' : '0';
		}
		else{
			$A = $line_content[$i] eq 0 ? '0' : '1';
                        $B= $line_content[$i+1] eq 0 ? '0' : '1';
		}
		if(!defined $rsHash{$line_content[0]}){
			my @tmp = ();
			push(@tmp,"$A:$B");
			$rsHash{$line_content[0]} = \@tmp;
		}
		else{
			push(@{$rsHash{$line_content[0]}},"$A:$B");
		}
	}

}
close INPUT;



print "End Reading Haplotype File\n";
print "Begin Writing OUTPUT file\n";
my $indiv_count = 0;
foreach my $indiv( @indivTable ){
	print OUTPUT "# ID $indiv\n";
	my $counter = 0;
	foreach my $rs (@rsOrder){
		my @alleles = split(/:/,${$rsHash{$rs}}[$indiv_count]);
		if($counter == 0){
			print OUTPUT $alleles[0];
		}
		else{
			print OUTPUT " " . $alleles[0];
		}
		$counter += 1;
	}
	print OUTPUT "\n";
	$counter =0;
	foreach my $rs (@rsOrder){
                my @alleles = split(/:/,${$rsHash{$rs}}[$indiv_count]);
                if($counter == 0){
                        print OUTPUT $alleles[1];
                }
                else{
                        print OUTPUT " " . $alleles[1];
                }
		$counter += 1;
        }
	print OUTPUT "\n";
	$indiv_count += 1;
}
print OUTPUT "END GENOTYPES\n";
close OUTPUT;
print "END\n";
