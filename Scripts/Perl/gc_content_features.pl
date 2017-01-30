#!/usr/bin/perl
# script by tibo
#
# get gc content from feature list and fasta file
#
# USage: script.pl <features (gtf)> <fasta file>

print STDERR "Get GC content by features\n  $ARGV[0]\n  $ARGV[1]\n";

# GTF
open GTF, $ARGV[0];

while (<GTF>) {
	if (/^#/) { next } #headers
	chomp;
	@cols = split;
	$cols[9] =~ s/[";]//g;
	$gtf{$cols[0]}{$cols[9]}{$cols[3]} = $cols[4];
	
}


print STDERR '' . (keys %gtf) . " features\n";

# read fasta file chr by chr, extract GC data; 
open FASTA, $ARGV[1];
$first = 1;

while (<FASTA>) {
	chomp;
	if (/^>(\w+)/) { #new chr
		if (!$first) { # not 1st chromosome : data to compute
		
			# for each feat, get GC content
			gc_content();
		}
		$chr = $1;
		$first = 0;
		$seq = '';
		print STDERR "$chr\n";
	}
	$seq .= $_;
}

# last chr: for each feat, get GC content
gc_content();

print STDERR "DONE\n";

sub gc_content {
	foreach $feat (keys %{ $gtf{$chr} }) {
		$feat_seq = '';
		
		# get all feature sequence
		foreach $start (keys %{ $gtf{$chr}{$feat} }) {
			$feat_seq .= substr($seq, $start, ($gtf{$chr}{$feat}{$start} - $start));
		}
		
		# compute GC content
		$gc = ($feat_seq =~ tr/GC/GC/) / length($feat_seq);
		
		$result{$feat} = length($feat_seq);
		$result{$feat} = $gc;
	}
}

$str1 = '';
$str2 = '';
foreach $feat (sort keys %result) {
	$str1 .= "$feat\t";
	$str2 .= "$result{$feat}\t";
}
chop $str1;
chop $str2;

print "$str1\n$str2\n";