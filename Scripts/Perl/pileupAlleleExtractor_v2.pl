#!/usr/bin/perl -w

#use diagnostics;
use strict;

my $infile;

my $minQ = 17;

if(@ARGV != 1 ){
        die "Usage: $0 <InputFile>";
}
if( @ARGV == 1 ){
        $infile = shift;
}

my $infile_prefix = $infile;$infile_prefix=~s/.*\/(.*)\.\w*/$1/;


open( INPUT, "$infile" ) or dieWithUnexpectedError("can't open $infile!");

open( OUTFILE, ">$infile_prefix.alleleCount.txt" ) or dieWithUnexpectedError("can't write to $infile_prefix.alleleCount.txt!");

while(my $line = <INPUT>){
        chomp $line;
	my @line_content = split(/[\s\t]/,$line);
	my $chr = $line_content[0];
	my $pos = $line_content[1];
	my $ref = uc($line_content[2]);
	my $dp = $line_content[3];
	my $bases = uc($line_content[4]);
	my $score = $line_content[5];

	#Remove positional characters
	$bases =~ s/[\^\]\$\[]//g;	
	my @bases = split(//,$bases);
	
	my %allele_count;

	#Treat the indels
	my $indel_status = 0;
	my $indelSize = 0;
	foreach my $base (@bases){
		if($base eq '+' || $base eq '-'){
			$indel_status = 1;
			$base = "!";
		}
		elsif($indel_status && $base =~ m/\d/){
			if($indelSize == 0){
				$indelSize = $base;
			}
			else{
				$indelSize .=$base;
			}
			$base = "!";
		}
		elsif($indel_status && $indelSize >0 && $base !~ m/\d/){
			$base = "!"; #Transform in string again after and remove those
			$indelSize = $indelSize -1;
		}
	}
	
	my $bases2 = join('',@bases);
	$bases2 =~ s/!//g;

	#Phase2: REMOVE THE BASES WITH QUALITY LESS THAN $minQ
	@bases = split(//,$bases2);
	my @score = split(//,$score);	

	my @nuc = ('A','C','G','T','*');
	foreach my $t(@nuc){
		$allele_count{$t}=0;
	}

	for (my $i = 0 ; $i < scalar(@score) ; $i++){
		if( ord($score[$i])-33 < $minQ ){}
		else{
			if( $bases[$i] =~ m/[\.,]/){
				if(!defined $allele_count{$ref}){
					$allele_count{$ref} = 0;
				}
				$allele_count{$ref} +=1; 
			}
			elsif($bases[$i] =~ m/[\<\>]/){}
			else{
				if(!defined $allele_count{$bases[$i]}){
					$allele_count{$bases[$i]} = 0;
				}
				$allele_count{$bases[$i]} +=1;
			}
		}
	}
	
	#my @altKeysSorted = sort {$allele_count_alt{$b} <=> $allele_count_alt{$a} } keys %allele_count_alt;

=comment
	if(!defined $allele_count_ref{$ref}){
		$allele_count_ref{$ref} = 0;
	}

	if(scalar(keys %allele_count_alt) == 0){
		print OUTFILE "$chr\t$pos\tref,.:dpRef,.\t$ref,.:" . $allele_count_ref{$ref} . ",.\n";
	}
	elsif(scalar(keys %allele_count_alt) == 1){
                print OUTFILE "$chr\t$pos\tref,alt:dpRef,dpAlt\t$ref";
		foreach my $keys (@altKeysSorted){
                        print OUTFILE ",$keys";
                }
                print OUTFILE ":" . $allele_count_ref{$ref};
                foreach my $keys (@altKeysSorted){
                        print OUTFILE ",".$allele_count_alt{$keys};
                }
                print OUTFILE "\n";
        }
	elsif(scalar(keys %allele_count_alt) == 2){
                print OUTFILE "$chr\t$pos\tref,alt1,alt2:dpRef,dpAlt1,dpAlt2\t$ref";
		foreach my $keys (@altKeysSorted){
			print OUTFILE ",$keys";
		}
		print OUTFILE ":" . $allele_count_ref{$ref};
		foreach my $keys (@altKeysSorted){
                        print OUTFILE ",".$allele_count_alt{$keys};
                }
		print OUTFILE "\n";
		
        }
	elsif(scalar(keys %allele_count_alt) == 3){
                print OUTFILE "$chr\t$pos\tref,alt1,alt2,alt3:dpRef,dpAlt1,dpAlt2,dpAlt3\t$ref";
		foreach my $keys (@altKeysSorted){
                        print OUTFILE ",$keys";
                }
                print OUTFILE ":" . $allele_count_ref{$ref};
                foreach my $keys (@altKeysSorted){
                        print OUTFILE ",".$allele_count_alt{$keys};
                }
                print OUTFILE "\n";
        }
	else{
		print OUTFILE "$chr\t$pos\tref,alt:dpRef,dpAlt\t$ref";
		foreach my $keys (@altKeysSorted){
                        print OUTFILE ",$keys";
                }
                print OUTFILE ":" . $allele_count_ref{$ref};
                foreach my $keys (@altKeysSorted){
                        print OUTFILE ",".$allele_count_alt{$keys};
                }
                print OUTFILE "\n";
	}
=cut
	print OUTFILE "$chr\t$pos\t$ref\t";
	my $index_tmp = 0;
	foreach my $keys (@nuc){
		if($index_tmp==0){
			print OUTFILE $keys;
		}
		else{
	        	print OUTFILE ",".$keys;
		}
		$index_tmp+=1;
        }
	print OUTFILE ":";
	$index_tmp = 0;
	foreach my $keys (@nuc){
		if($index_tmp==0){
	                print OUTFILE $allele_count{$keys};
		}
		else{
			print OUTFILE ",".$allele_count{$keys};
		}
		$index_tmp+=1;
        }
	print OUTFILE "\n";
}
close INPUT;
close OUTFILE;
