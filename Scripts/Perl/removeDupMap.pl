#!/usr/bin/perl -w

#use diagnostics;
use strict;

my $ofh = select LOG;
$| = 1;
select $ofh;


my $infile;

if(@ARGV != 1 ){
        die "Usage: $0 <InputFile>";
}
if( @ARGV == 1 ){
        $infile = shift;
}

my $infile_prefix = $infile;$infile_prefix=~s/(.*)\.\w*/$1/;


open( INPUT, "$infile" ) or die("can't open $infile!");

open( OUTFILE, ">$infile_prefix.noDupToExtract.txt" ) or die("can't write to $infile_prefix.noDupToExtract.txt!");

my %snp;
my $presPos = "NA";
my @array;
while(my $line = <INPUT>){
        chomp $line;
	my @line_content = split(/[\s\t]/,$line);
	my $chr = $line_content[0];
	my $id = $line_content[1];
	my $pos = $line_content[0]."_".$line_content[3];
	if($pos ne $presPos && !defined ($snp{$pos})){
		push(@array,$pos);
		$presPos = $pos;
	}
	else{}
	if(!defined $snp{$pos}){
		$snp{$pos} = $line;
		next;
	}
	else{
		if($id=~/^rs/){
			$snp{$pos} = $line;
		}
		elsif((split(/[\s\t]/,$snp{$pos}))[1] =~ /^rs/){
			next;
		}
		elsif($id=~/^kgp/){
                        $snp{$pos} = $line;
                }
                elsif((split(/[\s\t]/,$snp{$pos}))[1] =~ /^kgp/){
                        next;
                }
		elsif($id=~/^GA/){
                        $snp{$pos} = $line;
                }
                elsif((split(/[\s\t]/,$snp{$pos}))[1] =~ /^GA/){
                        next;
                }
		elsif($id=~/^exm/){
                        $snp{$pos} = $line;
                }
                elsif((split(/[\s\t]/,$snp{$pos}))[1] =~ /^exm/){
                        next;
                }
		else{
			$snp{$pos}=$line;
		}
	}
}
close INPUT;

foreach my $s(@array){
	my $ofh = select LOG;
	$| = 1;
	select $ofh;

	print OUTFILE $snp{$s} . "\n";
}

close OUTFILE;
