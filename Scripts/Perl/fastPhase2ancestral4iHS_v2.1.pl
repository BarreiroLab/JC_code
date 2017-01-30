#!/usr/bin/perl -w

######################
#Script wrote by Jean-Christophe Grenier
#jean.christophe.grenier@gmail.com
######################

##### LOG #####
# Modifications
# June 5th 2012 : Don't consider a snp if no ancestral information
# March 15th 2013 : consider the ancestral information coming from Alan's file (alignment consensus)
#

#use diagnostics;
use strict;
use IO::Handle;

$| = 1;

my $infile;
my $map;
my $ancestral;
my $recombination;
my $freqFile;
my $outfile;
my $printBad;

if(@ARGV != 6 ){
        die "Usage: $0 <InputFile.fastPhase> <mapFile> <FreqFile> <ancestralInfoFile> <recombinationMap> <1 or 0 (1=print unknown Ancestral as Major)>";
}
if( @ARGV == 6 ){
        $infile = shift;
	$map = shift;
	$freqFile = shift;
	$ancestral = shift;
	$recombination = shift;
	$printBad = shift;
}

my $infile_prefix = $infile;$infile_prefix=~s/(.*)\..*/$1/; 



my %snpPos;
my %rs_ID;
print "Begin Read Map file\n";
readMapFile($map);
print "End Read Map file\n";

#Will be made in readFreqFile
my %majorAllele;
my %minorAllele;
my @snpOK;
my @snpOrder; 
print "Begin Read Freq file\n";
readFreqFile($freqFile);
print "End Read Freq file\n";

my %ancestralA;

open( OUTPUTFLAG,">$infile.flags.txt") or die ("can't write to $infile.flags.txt!");
my %ancestralAllele;
my %derivedAllele;
#Will be made in readAncestr function
print "Begin Read Ancestral file\n";
readAncestr($ancestral);
close( OUTPUTFLAG);
print "End Read Ancestral file\n";


my %genetic;
print "Begin Read Genetic file\n";
getGeneticPos($recombination);
print "End Read Genetic file\n";


open( INPUT, "$infile" ) or die("can't open $infile!");
open( OUTPUTHAP, ">$infile.ihshap") or die("can't write to $infile.ihshap!"); 
open( OUTPUTMAP, ">$infile.ihsmap") or die ("can't write to $infile.ihsmap!");

my $READ_Begin = 0;
while(my $line = <INPUT>){
	chomp $line;
	#Phased data from fastPHASE (*hapguess_switch.out)
	#IN THIS FORMAT : 1 = Major Allele, 0 = Minor Allele
	if($line =~ /BEGIN GENOTYPES/){
		$READ_Begin = 1;
	}
	elsif($line =~ /END GENOTYPES/){
		last;
	}
	elsif($line =~ /^#.*/ && $READ_Begin eq 1){
	}
	elsif($READ_Begin == 1){
#		$count ++;
		my @line_content = split(/[\s\t]{1,}/,$line);
		my $lineToPrint = "";
		my $count =0 ;
		foreach my $snp(@snpOrder){
			if($snpOK[$count] eq '0'){$count++; next;}
			elsif($snpOK[$count] eq '1' && !$ancestralAllele{$snp}){
				$snpOK[$count] = '0';
			}
			#Major
			elsif($line_content[$count] == "1"){
				#Look for the ancestral
				if($ancestralAllele{$snp} eq $majorAllele{$snp}){
	                 	       $lineToPrint .= "1\t";
				}
				else{
					$lineToPrint .= "0\t";
				}
                        }
			#Minor
                        else{
				#Look for the ancestral
				if($ancestralAllele{$snp} eq $minorAllele{$snp}){
	                	        $lineToPrint .= "1\t";
				}
				else{
					$lineToPrint .= "0\t";
				}
                        }
			$count ++;
		}
		print OUTPUTHAP $lineToPrint . "\n";
	}
}
close INPUT;
close OUTPUTHAP;

my $count = 0;
foreach my $snpPres (@snpOrder){
	#MOD BY JCG 20/03/13
	
	my $pos = $snpPos{$snpPres};
	if($snpOK[$count] eq '0'){$count++; next;}
	elsif( $snpPos{$snpPres} && $ancestralAllele{$snpPres} ){
		#print "$snpPres " . $pos . " " . $genetic{$snpPres} . " " . $derivedAllele{$snpPres} . " " . $ancestralAllele{$snpPres} . "\n";
		print OUTPUTMAP "$snpPres " . $pos . " " . $genetic{$snpPres} . " " . $derivedAllele{$snpPres} . " " . $ancestralAllele{$snpPres} . "\n";
	}
	$count ++;
}
close OUTPUTMAP;

sub getGeneticPos{
	my $file = shift;
	open( INPUT, "$file" ) or die("can't open $file!");
	my $indexSNP = 0;
	my @sortPos = sort {$a <=> $b} keys %rs_ID;
	my $lastPos = 0;
	my $lastGen = 0;
        while( <INPUT> ){
		my $line = $_;
		if($.==1){next;}
		if($indexSNP >= scalar(@sortPos)){last;}
		chomp $line;
                my @line_content = split(/[\t\s]{1,}/,$line);
		if($sortPos[$indexSNP] < $line_content[1] && $line_content[3] == 0 ){
			#IF NO INFO AT THE BEGINNING OF A CHROMOSOME, PRINT 0
			if($printBad == 0){
                                $snpOK[$indexSNP] = 0;
                        }
			$genetic{ $rs_ID{$sortPos[$indexSNP]} } = 0;
			$lastPos = $line_content[1];
			$lastGen = $line_content[3];
			$indexSNP +=1;
			redo;
		}
		elsif($sortPos[$indexSNP] == $line_content[1]){
			#IF WE HAVE A PERFECT MATCH
			$genetic{ $rs_ID{$sortPos[$indexSNP]} } = $line_content[3];
			$lastPos = $line_content[1];
			$lastGen = $line_content[3];
			$indexSNP +=1;
		}
		elsif($sortPos[$indexSNP] > $line_content[1]){
			#IF WE'RE NOT IN A MIDDLE OF 2 POSITIONS WITH INFO GO TO NEXT
			$lastPos = $line_content[1];
			$lastGen = $line_content[3];
		}
		elsif($sortPos[$indexSNP] < $line_content[1] && $lastPos != 0){
			#PRINT THE CURRENT POSITION BASED ON KNOWN INFO
			my $diffPos = $line_content[1] - $lastPos;
			my $diffGen = $line_content[3] - $lastGen;
			my $perPos = $diffGen/$diffPos;
			$genetic{ $rs_ID{$sortPos[$indexSNP]} } = $lastGen + (($sortPos[$indexSNP]-$lastPos)*$perPos);
			
		#	$lastPos = $line_content[0];
		#	$lastGen = $line_content[2];
			$indexSNP +=1;
			redo;
		}
	}
	while($indexSNP < scalar(@sortPos)){
		if($printBad == 0){
                	$snpOK[$indexSNP] = 0;
                }
		$genetic{ $rs_ID{$sortPos[$indexSNP]} } = $lastGen;
		$indexSNP +=1;
	}
	
}

sub readMapFile{
	my $file = shift;
        open( INPUT, "$file" ) or die("can't open $file!");
        while(my $line = <INPUT>){
                chomp $line;
                my @line_content = split(/[\t\s]{1,}/,$line);
                my $rs_id = $line_content[1];
		if($rs_ID{$line_content[3]}){
			push(@snpOK,0);
			push(@snpOrder,$rs_id);
			$snpPos{$rs_id}=$line_content[3];
	#		next;
		}
		else{
			push(@snpOK,1);
			push(@snpOrder,$rs_id);
		        $rs_ID{$line_content[3]} = $rs_id;
			$snpPos{$rs_id}=$line_content[3];
		}
        }

        close INPUT;
}

sub readFreqFile{
	my $file = shift;
	open( INPUT, "$file" ) or die("can't open $file!");
	while(my $line = <INPUT>){
		if($.==1){next;}
	        chomp $line;
		my @line_content = split(/[\t\s]{1,}/,$line);
		my $rs_id = $line_content[1];
		if($snpPos{$rs_id}){
			#push(@snpOK,1);			#RECONSIDER
			my $minor = $line_content[2];
			my $major = $line_content[3];
			$majorAllele{$rs_id}=$major;
			$minorAllele{$rs_id} = $minor;
		}
		else{
	#		push(@snpOK,0);
		
		}
	}
	close INPUT;
}

sub readAncestr{
        my $file = shift;
	my $index = 0;
        open( INPUT, "$file" ) or die("can't open $file!");
        while(my $line = <INPUT>){
#		if($.==1){next;}
                chomp $line;

                my @line_content = split(/[\t\s]{1,}/,$line);
		my $chr = $line_content[0];
		my $pos = $line_content[1];
		my $ancestral = $line_content[2];
		if($ancestral =~ m/\w/){
			$ancestral = uc($ancestral);
		}
		my $rs = $rs_ID{$pos};
		
		my $derived;
		
		#MEANS THAT THE FIRST POS IN THE MAP DOESN'T HAVE INFO ABOUT THE ANCESTRAL
		if($snpPos{$snpOrder[$index]} < $pos){
			while($snpPos{$snpOrder[$index]} < $pos && $snpOrder[$index] ){
				print OUTPUTFLAG $line . "\tNotInAnc\n";
                	        $ancestral = $majorAllele{$rs};
                        	$derived = $minorAllele{$rs};
				if($printBad == 0){
					$snpOK[$index] = 0;
				}
				$ancestralAllele{$rs} = $ancestral;
		                $derivedAllele{$rs} = $derived;
        		        $index +=1;
			}
		}
		if($ancestral eq 'N' || $ancestral eq '.' || $ancestral eq '-' || $snpPos{$snpOrder[$index]} > $pos){
			print OUTPUTFLAG $line . "\tNotEnoughInfo\n";
			$ancestral = $majorAllele{$rs};
			$derived = $minorAllele{$rs};
			if($printBad == 0){
                                $snpOK[$index] = 0;
                        }
		}
		else{
			if($ancestral ne $majorAllele{$rs} && $ancestral ne $minorAllele{$rs}){
				my $ancestral_old = $ancestral;
				$ancestral = flip($ancestral);
				
				if($ancestral ne $majorAllele{$rs} && $ancestral ne $minorAllele{$rs}){
					$ancestral = $majorAllele{$rs};
                                        $derived = $minorAllele{$rs};
                                        print OUTPUTFLAG $line . "\tMismatchWithAnc\n";
					if($printBad == 0){
	                	                $snpOK[$index] = 0;
        		                }	
				}
			}
			
			if( ($majorAllele{$rs} eq $ancestral && $majorAllele{$rs} eq 'A' && $minorAllele{$rs} eq 'T') || 
				($majorAllele{$rs} eq $ancestral && $majorAllele{$rs} eq 'T' && $minorAllele{$rs} eq 'A') ||
				($minorAllele{$rs} eq $ancestral && $majorAllele{$rs} eq 'A' && $minorAllele{$rs} eq 'T') ||
				($minorAllele{$rs} eq $ancestral && $majorAllele{$rs} eq 'T' && $minorAllele{$rs} eq 'A') ){
				
				if($ancestral eq $minorAllele{$rs}){
					$ancestral = $majorAllele{$rs};
					$derived = $minorAllele{$rs};
					print OUTPUTFLAG $line . "\tATtype_majorChoosen\n";
					if($printBad == 0){
                         		       $snpOK[$index] = 0;
		                        }
				}
				else{
					$ancestral = $majorAllele{$rs};
					$derived = $minorAllele{$rs};
				}
			}
			elsif( ($majorAllele{$rs} eq $ancestral && $majorAllele{$rs} eq 'C' && $minorAllele{$rs} eq 'G') || 
                                ($majorAllele{$rs} eq $ancestral && $majorAllele{$rs} eq 'G' && $minorAllele{$rs} eq 'C') ||
                                ($minorAllele{$rs} eq $ancestral && $majorAllele{$rs} eq 'C' && $minorAllele{$rs} eq 'G') ||
                                ($minorAllele{$rs} eq $ancestral && $majorAllele{$rs} eq 'G' && $minorAllele{$rs} eq 'C') ){

				if($ancestral eq $minorAllele{$rs}){
                                        $ancestral = $majorAllele{$rs};
					$derived = $minorAllele{$rs};
					print OUTPUTFLAG $line . "\tGCtype_majorChoosen\n";
					if($printBad == 0){
                	         	       $snpOK[$index] = 0;
		                        }
                                }
				else{
                                        $ancestral = $majorAllele{$rs};
                                        $derived = $minorAllele{$rs};
                                }
                        }
			else{
				#Find the derived allele
				if($ancestral eq $majorAllele{$rs}){
					 $derived = $minorAllele{$rs};
				}
				elsif($ancestral eq $minorAllele{$rs}){
					 $derived = $majorAllele{$rs};
				}
			}
			
		}
		$ancestralAllele{$rs} = $ancestral;
		$derivedAllele{$rs} = $derived;
		$index +=1;
        }

        close INPUT;
}

sub flip{
	my $a = shift;
	if($a eq 'A'){return 'T';}
	if($a eq 'C'){return 'G';}
	if($a eq 'G'){return 'C';}
	if($a eq 'T'){return 'A';}
}
