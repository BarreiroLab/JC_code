#!/usr/bin/perl -w

#use diagnostics;
use strict;

my $infile;


if(@ARGV != 1 ){
        die "Usage: $0 <InputFile.gtf>";
}
if( @ARGV == 1 ){
        $infile = shift;
}

my $infile_prefix = $infile;$infile_prefix=~s/(.*)\..*/$1/;


open( INPUT, "$infile" ) or dieWithUnexpectedError("can't open $infile!");

open( OUTFILE, ">$infile_prefix.1transcriptPerGene.gtf" ) or dieWithUnexpectedError("can't write to $infile_prefix.1transcriptPerGene.gtf!");

my %genes;
my %geneOrder;
my %strand;
my %chromosome;

my $counter = 1;
while(my $line = <INPUT>){
        chomp $line;
	if($line !~ "exon"){next;}
	my @line_content = split(/[\t\s]/,$line);
	my $chr = $line_content[0];
	my $strd = $line_content[6];
	my $start = $line_content[3];
	my $stop = $line_content[4];
	my $gene = $line_content[9];
	my @geneArray = split(/\"/,$gene);
	my $geneName = $geneArray[1];
	if(!defined $geneOrder{"$geneName:$chr"}){
		$geneOrder{"$geneName:$chr"} = $counter;
		$counter +=1;
	}
	else{}
	$chromosome{"$geneName:$chr"} = $chr;
	$strand{"$geneName:$chr"}=$strd;	

	if(!defined $genes{"$geneName:$chr"}){
		my %tmp = ();
		for( my $i=$start ; $i <= $stop ; $i++){
			$tmp{$i} = 1;
		}
		$genes{"$geneName:$chr"}= \%tmp;
	}
	else{
		for( my $i=$start ; $i <= $stop ; $i++){
                        ${$genes{"$geneName:$chr"}}{$i} = 1;
                }
	}
}
close INPUT;

#Print the file
foreach my $gene_tmp (sort { $geneOrder{$a}<=>$geneOrder{$b} } keys %geneOrder){
	my $prevPos = 0;
	my $start;
	my $end;
	my $gene_Ok = (split(/:/,$gene_tmp))[0];
	foreach my $pos(sort {$a<=>$b} keys %{$genes{$gene_tmp}}){
		if($prevPos == 0){
			$start = $pos;
			$prevPos = $pos;
		}
		else{
			#Print the exon
			if($pos != ($prevPos+1)){
				$end = $prevPos;
				print OUTFILE $chromosome{$gene_tmp} . "\tunknown\texon\t$start\t$end\t.\t" .  $strand{$gene_tmp} . "\t.\tgene_id \"$gene_Ok\"; gene_name \"$gene_Ok\"; transcript_id \"$gene_Ok\";\n";

				$start = $pos;
				$prevPos = $pos;
				
			}
			else{
				$prevPos = $pos;
			}
		}
	}
	$end = $prevPos;
	#Print last exon
	print OUTFILE $chromosome{$gene_tmp} . "\tunknown\texon\t$start\t$end\t.\t" .  $strand{$gene_tmp} . "\t.\tgene_id \"$gene_Ok\"; gene_name \"$gene_Ok\"; transcript_id \"$gene_Ok\";\n";;
}

close OUTFILE;
