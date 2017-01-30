#!/usr/bin/perl
# script by Thibault de Malliard
#
# Usage :
# format_plink_lgen_for_fst.pl <individual list> <lgen file>
#    -individual list: an individual by line with it city with a num code like : Cd_001,1


if (@ARGV != 2) {
	print STDERR "not enough arguments\n";
	print STDERR "  $0 <individual list> <lgen file>\n";
	exit;
}

open INDIV, $ARGV[0];

while (<INDIV>) {
	chomp;
	split(',');
	$data{$_[0]} = $_[1];
}

close INDIV;

open SNPS, $ARGV[1];

while (<SNPS>) {
	chomp;
	split(/\s{1,}/);
	if($_[3] == 0){
		my $indiv = $_[0].":".$_[1];
		print "$_[2],$_[0]:$_[1],$data{$indiv},-,-\n";
	}
	else{
		$_[3] = $_[3]==1 ? 'A' : 'B'; 
		$_[4] = $_[4]==1 ? 'A' : 'B';
		my $indiv = $_[0].":".$_[1];
		print "$_[2],$_[0]:$_[1],$data{$indiv},$_[3],$_[4]\n";
	}
}

close SNPS;
