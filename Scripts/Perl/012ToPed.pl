#!/usr/bin/perl -w

#use diagnostics;
use strict;

my $infile;

if(@ARGV != 1 ){
        die "Usage: $0 <InputFile>";
}
if( @ARGV == 1 ){
        $infile = shift;
}

my $infile_prefix = $infile;$infile_prefix=~s/.*\/(.*)\.\w*/$1/;


open( INPUT, "$infile" ) or dieWithUnexpectedError("can't open $infile!");

open( OUTFILE, ">$infile_prefix.ped" ) or dieWithUnexpectedError("can't write to $infile_prefix.ped!");
open( OUTFILEMAP, ">$infile_prefix.map" ) or dieWithUnexpectedError("can't write to $infile_prefix.map!");

my @samples;
my %genotypes;

while(my $line = <INPUT>){
        chomp $line;
	my @line_content = split(/[\s\t]/,$line);
	if($.==1){
		for( my $i=6; $i < scalar(@line_content) ; $i+=1){
			push(@samples,$line_content[$i]);
			my @tmp;
			$genotypes{$line_content[$i]} = \@tmp;
		}
	}
	else{
		my $id = $line_content[0];
		my $chr = $line_content[1];
		my $pos = $line_content[2];
		my $g1 = $line_content[3].$line_content[3]; 
		my $g2 = $line_content[4].$line_content[4];
		my $g3 = $line_content[5];
		print OUTFILEMAP "$chr\t$id\t0\t$pos\n";
		my $j = 0;
		for (my $i=6; $i < scalar(@line_content) ; $i+=1){
			my $geno;
			if($line_content[$i] eq '0'){
				$geno = $g1;
			}
			elsif($line_content[$i] eq '1'){
                                $geno = $g2;
                        }
			elsif($line_content[$i] eq '2'){
                                $geno = $g3;
                        }
                        push(@{$genotypes{$samples[$j]}},$geno);
			$j += 1;
                }
	}
	

}
close INPUT;

foreach my $s(@samples){
	print OUTFILE "$s\t$s\t0\t0\t0\t-9";
	foreach my $g(@{$genotypes{$s}}){
		my @a = split(//,$g);
		print OUTFILE "\t" . $a[0] ."\t". $a[1];
	}
	print OUTFILE "\n";
}

close OUTFILE;
close OUTFILEMAP;
