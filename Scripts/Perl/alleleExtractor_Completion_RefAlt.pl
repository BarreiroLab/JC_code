#!/usr/bin/perl -w

#use diagnostics;
use strict;
use Getopt::Long qw (:config gnu_getopt);

my $infile;
my $pos;
my $posFilter;


my $parms = GetOptions (
        "infile=s" => \$infile,
        "pos=s"=>\$pos,
	"posFilter=s"=>\$posFilter,
) or die();

my %errors = (
        usg => "Usage: # $0 --infile=<file.alleleCount.txt> --pos=<chr_pos_ref.file.txt> --posFilter=<pos_from_VCF.txt>\n\n"
 #      nofile => "Please specify a good entry file\n",
);



if(!$infile && !$pos ){
        DIE ($errors{"usg"});
}

my $infile_prefix = $infile;$infile_prefix=~s/(.*)\.\w{1,}/$1/;


open( POS, "$pos" ) or die("can't open $pos!\n");

my @positions;
while(my $line = <POS>){
        chomp $line;
	my @line_content = split(/[\s\t]/,$line);
	my $chr = $line_content[0];
	my $pos = $line_content[1];
	my $ref = $line_content[2];

	push(@positions,"$chr\_$pos\_$ref");
}
close POS;


my %posF = ();

if($posFilter){
	open( POS, "$posFilter" ) or die("can't open $posFilter!\n");


	while(my $line = <POS>){
        	chomp $line;
	        my @line_content = split(/[\s\t]/,$line);
	        my $chr = $line_content[0];
        	my $pos = $line_content[1];
		$posF{"$chr\_$pos"} = 1;
	}
	close POS;
}


open( OUTFILE, ">$infile_prefix.alleleCountCorrectedRefAlt.txt" ) or dieWithUnexpectedError("can't write to $infile_prefix.alleleCountCorrectedRefAlt.txt!");

open( INPUT, "$infile") or die("can't open $infile!\n");
my $index = 0;
while(my $line = <INPUT>){
	chomp $line;
        my @line_content = split(/[\s\t]/,$line);
        my $chr = $line_content[0];
        my $pos = $line_content[1];
	my $ref = $line_content[2];
	my $bases = $line_content[3];

	if($positions[$index] ne "$chr\_$pos\_$ref"){
		while($positions[$index] ne "$chr\_$pos\_$ref"){
			my @posRef = split(/_/,$positions[$index]);

			if(defined $posFilter){		
				if(defined $posF{"$chr\_$pos"}){
					print OUTFILE $posRef[0] . "\t" . $posRef[1] . "\t" . $posRef[2] . "\t" . "Ref,Alt,Filter:0,0,PASS\n";
				}
				else{
					print OUTFILE $posRef[0] . "\t" . $posRef[1] . "\t" . $posRef[2] . "\t" . "Ref,Alt,Filter:0,0,FILTERED\n";
				}
			}
			else{
				print OUTFILE $posRef[0] . "\t" . $posRef[1] . "\t" . $posRef[2] . "\t" . "Ref,Alt:0,0\n";
			}
			$index +=1;
		}
	}
	if($positions[$index] eq "$chr\_$pos\_$ref"){
		#GET COUNT TOTAL :
		my @b = split(/:/,$bases);
		my @counts = split(/,/,$b[1]);
		my $total = 0;
		my $ref_count;
		my $alt_count=0;
		my $indexa= 0;
		my $indexRef = -1;
		if($ref eq 'A'){
			$ref_count = $counts[0];
			$indexRef = 0;
		}
		elsif($ref eq 'C'){
			$ref_count = $counts[1];
			$indexRef = 1;
		}
		elsif($ref eq 'G'){
			$ref_count = $counts[2];
			$indexRef = 2;
                }
		elsif($ref eq 'T'){
			$ref_count = $counts[3];
			$indexRef =3;
                }
		foreach my $i(@counts){
			if($indexa == 4){
                                last;
                        }
			if($indexa == $indexRef){}
			else{
				
				if($i > $alt_count){
					$alt_count = $i;
				}
			}
#			$total += $i;
			$indexa +=1;
		}
		if(defined $posFilter ){
			print OUTFILE "$chr\t$pos\t$ref\tRef,Alt,Filter:".$ref_count.",".$alt_count;
		}
		else{
			print OUTFILE "$chr\t$pos\t$ref\tRef,Alt:".$ref_count.",".$alt_count;
		}
		
	}

	if(defined $posF{"$chr\_$pos"} && $posFilter ne ""){
		print OUTFILE ",PASS";
	}
	elsif(!defined $posF{"$chr\_$pos"} && !$posFilter){
	}
	else{
		print OUTFILE ",FILTERED";
	}
	print OUTFILE "\n";
	$index +=1;
}
if($index < scalar(@positions)){
	for(my $j = $index ; $j < scalar(@positions) ; $j+=1 ){
		my @posRef = split(/_/,$positions[$j]);

		if(defined $posFilter){
			print OUTFILE $posRef[0] . "\t" . $posRef[1] . "\t" . $posRef[2] . "\t" . "Ref,Alt,Filter:0,0";
		}
		else{
			print OUTFILE $posRef[0] . "\t" . $posRef[1] . "\t" . $posRef[2] . "\t" . "Ref,Alt:0,0";
		}
		my $chr = $posRef[0];
		my $pos = $posRef[1];
		if(defined $posF{"$chr\_$pos"} && $posFilter ne ""){
	                print OUTFILE ",PASS";
	        }
        	elsif(!defined $posF{"$chr\_$pos"} && !$posFilter){
	        }
        	else{
	                print OUTFILE ",FILTERED";
        	}
		print OUTFILE "\n";
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
