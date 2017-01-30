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

open( INPUT, "$infile" ) or die("can't open $infile!");

my @samples;
my %genes;
my (@no_feature,@ambiguous,@too_low_aQual,@not_aligned,@alignment_not_unique);

while(my $line = <INPUT>){
        chomp $line;

	my @line_content = split(/[\s\t]+/,$line);
	if($.==1){
		foreach my $s(1 .. $#line_content){
			push(@samples,$line_content[$s]);
			$genes{$line_content[$s]}=0;
		}
	}
	elsif($line_content[0]eq 'no_feature'){
		foreach my $s(1 .. $#line_content){
			push(@no_feature,$line_content[$s]);
		}
	} 
	elsif($line_content[0]eq 'ambiguous'){
                foreach my $s(1 .. $#line_content){
                        push(@ambiguous,$line_content[$s]);
                }

	}
	elsif($line_content[0]eq 'too_low_aQual'){
                foreach my $s(1 .. $#line_content){
                        push(@too_low_aQual,$line_content[$s]);
                }

	}
	elsif($line_content[0]eq 'not_aligned'){
                foreach my $s(1 .. $#line_content){
                        push(@not_aligned,$line_content[$s]);
                }

	}
	elsif($line_content[0]eq 'alignment_not_unique'){
                foreach my $s(1 .. $#line_content){
                        push(@alignment_not_unique,$line_content[$s]);
                }

	}
	else{
		#push(@genes,$line_content[0]);
		
		#Sum all genes
		foreach my $s(1 .. $#line_content){
               	        $genes{$samples[$s-1]}+=$line_content[$s];
        	}
	}
	
}
close INPUT;

foreach my $s(0 .. $#samples){
	system("mkdir $samples[$s]");
	my $output=$samples[$s]."/counts.txt.summary";
	open (OUTFILE,">$output");
	
	print OUTFILE "Status\t".$samples[$s].".bam\n";
	print OUTFILE "Assigned\t".$genes{$samples[$s]}."\n";
	print OUTFILE "Unassigned_Ambiguity\t".$ambiguous[$s]."\n";

	print OUTFILE "Unassigned_MultiMapping\t".$alignment_not_unique[$s]."\n";
	print OUTFILE "Unassigned_NoFeatures\t".$no_feature[$s]."\n";
	print OUTFILE "Unassigned_Unmapped\t".$not_aligned[$s]."\n";
	print OUTFILE "Unassigned_MappingQuality\t".$too_low_aQual[$s]."\n";
print OUTFILE "Unassigned_FragmentLength\t0
Unassigned_Chimera\t0
Unassigned_Secondary\t0
Unassigned_Nonjunction\t0
Unassigned_Duplicate\t0\n";


	close OUTFILE;	
}
