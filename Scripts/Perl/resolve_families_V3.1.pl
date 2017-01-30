#!/usr/bin/perl -w 

use strict;

my $infile;

my $numK;

if(@ARGV != 1 ){
        die "Usage: $0 <InputFile> ";
}
if( @ARGV == 1 ){
        $infile = shift;
}

my $infile_prefix = $infile;$infile_prefix=~s/(.*)\..*/$1/;

my %indiv_rel;

open(FILE,"$infile") || die "Could not open $infile!\n";

while(my $line = <FILE>){
	my @line_content = split(/[\s\t]{1,}/,$line);
	#The first family
	if(exists $indiv_rel{$line_content[1].":".$line_content[2]}){
		${$indiv_rel{$line_content[1].":".$line_content[2]}}{$line_content[3].":".$line_content[4]}=$line_content[10];
	}
	else{
		my %tmp;
		$tmp{$line_content[3].":".$line_content[4]}=$line_content[10];
		$indiv_rel{$line_content[1].":".$line_content[2]} = \%tmp;
	}
	if(exists $indiv_rel{$line_content[3].":".$line_content[4]}){
                ${$indiv_rel{$line_content[3].":".$line_content[4]}}{$line_content[1].":".$line_content[2]}=$line_content[10];
        }
	else{
		my %tmp ;
                $tmp{$line_content[1].":".$line_content[2]}=$line_content[10];
                $indiv_rel{$line_content[3].":".$line_content[4]} = \%tmp;
	}
}
close FILE;

open(OUTREM,">$infile_prefix.indivToRem.txt") || die "Cannot write to $infile_prefix.indivToRem.txt!\n";
open(OUTFAM,">$infile_prefix.groupByIndiv.txt") || die "Cannot write to $infile_prefix.groupByIndiv.txt!\n";

foreach my $key(sort keys %indiv_rel){
	print OUTFAM "Related to : $key\n";
	foreach my $key2(sort keys %{$indiv_rel{$key}}){
		print OUTFAM "$key2\n";
	}
	print OUTFAM "\n";
}
close OUTFAM;

#PRINT A MATRIX OF RELATEDNESS :
open(OUTFAMMAT,">$infile_prefix.relIndivMatrix.txt") || die "Cannot write to $infile_prefix.relIndivMatrix.txt!\n";

my $i = 0;
foreach my $key(sort keys %indiv_rel){
	if($i==0){
		print OUTFAMMAT "\t$key";
		$i =1;
	}
	else{
		print OUTFAMMAT "\t$key";
	}
}

foreach my $key(sort keys %indiv_rel){
	print OUTFAMMAT "\n$key";
	foreach my $key2(sort keys %indiv_rel){
		if(defined ${$indiv_rel{$key}}{$key2}){
			print OUTFAMMAT "\t".${$indiv_rel{$key}}{$key2};
		}
		elsif($key2 eq $key){
			print OUTFAMMAT "\t1";
		}
		else{
			print OUTFAMMAT "\t0";
		}
	}
}

close OUTFAMMAT;

#PRINT THE REJECTED INDIVIDUALS
#ESTABLISH A LIST
my %toReject;

foreach my $key(sort { scalar(keys %{$indiv_rel{$b}}) <=> scalar(keys %{$indiv_rel{$a}}) } keys %indiv_rel){
	#IF THE RELATED INDIV TO THE PRESENT ONE ONLY HAS 1 RELATED (so the key we're in)
	#Keep only $key Individual..
	my $found =0;
	if( scalar(keys %{$indiv_rel{$key}}) > 0){
		my @list_of_related = sort { scalar(keys %{$indiv_rel{$a}}) <=> scalar(keys %{$indiv_rel{$b}}) } keys %{$indiv_rel{$key}};
		if(scalar(@list_of_related) == 1 && scalar( keys %{$indiv_rel{$list_of_related[0]}}) ){
			$toReject{$list_of_related[0]} = 1;
			delete $indiv_rel{$list_of_related[0]};
		}
		else{
			#IF THERE ARE MULTIPLE INDIV RELATED TO THAT ONE:
			#DELETE IF KEY IS PRESENT IN MORE THAN 1 GROUP
			foreach my $key2(sort keys %indiv_rel){
				if($key2 ne $key && scalar(keys %{$indiv_rel{$key2}} > 0) ){
					if(defined ${$indiv_rel{$key2}}{$key}){
						$found +=1;
					}
				}
			}
			if($found > 1){
				$toReject{$key}=1;
				delete $indiv_rel{$key};
			}
		}
	}
}


foreach my $key(sort keys %toReject){
	my @key2 = split(/:/,$key);
	print OUTREM $key2[0] . "\t" . $key2[1] . "\n";
}
close OUTREM;
