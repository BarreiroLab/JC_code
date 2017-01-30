#!/usr/bin/perl -w

=for comment
#####HISTORIC#####
#Made by Jean-Christophe Grenier
#Script to make bash files to launch aligning pipeline on Exome Sequencing (paired-end) data from Illumina
=========
8/03/2013 : First draft of the script (JCG)

=========
=cut



#use diagnostics;
use strict;
use Getopt::Long qw (:config gnu_getopt);

my $usr = $ENV{'LOGNAME'};
my $currentDir = $ENV{'PWD'};

my $picard = "/home/apps/Logiciels/picard_tools/picard-tools-1.56";

my $center = "CHUM";

my $pair1;
my $pair2;
my $single;
my $phred = 33;
my $rg_id = "DefaultID";
my $rg_sample = "DefaultSample";
my $rg_library = "DefaultLib";
my $rg_platform = "Illumina";
my $rg_center = "Innovation_Center";
my $p = 12;
my $ref = "hg19";
my $ref_path_fasta = "";
my $ref_path_bwa = "";
my $gtf_path = "";

my $adapt1 = "CGGTTCAGCAGGAATGCCGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGAC";
my $adapt2 = "CGGCATTCCTGCTGAACCGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAG";

my $ref_hg18_bwa = "";
my $ref_hg18_fasta = "";
my $gtf_hg18 = "";

my $ref_hg19_bwa = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/hg19/bwa_0.7.7/hs_ref_GRCh37.fasta";
my $ref_hg19_fasta = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/hg19/hs_ref_GRCh37.fasta";
#my $ref_hg19CEU_bwa = "/exec5/GROUP/awadalla/grenier2/references/hg19/CEU/bwa/CEUref_hg19.fasta";
#my $ref_hg19CEU_fasta = "/exec5/GROUP/awadalla/grenier2/references/hg19/CEU/CEUref_hg19.fasta";
my $gtf_hg19 = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/hg19/GTF/Homo_sapiens.GRCh37.75.withCHR.ENSEMBL.stdCHR.gtf";

my $trim_method = "trim_galore";
my $trimQ = 20;
my $id_prefix = "";

my $parms = GetOptions (
	"single=s" => \$single,
	"pair1=s" => \$pair1,
	"pair2=s" => \$pair2,
	"phred=s" => \$phred,
	"a1=s" => \$adapt1,
	"a2=s" => \$adapt2,
	"rg_sample=s" => \$rg_sample,
	"rg_id=s" => \$rg_id,
	"rg_library=s" => \$rg_library,
	"rg_platform=s" => \$rg_platform,
	"rg_center=s" => \$rg_center,
	"ref=s"=>\$ref,
	"multithreading=s" => \$p,
	"trim=s" => \$trim_method,
	"trimq=s" => \$trimQ,
	"prefix=s"=> \$id_prefix,
) or die();

my %errors = (
        usg => "Usage: # $0 (--single=<input.fastq.gz> or (--pair1=<inputfile1.fastq.gz> --pair2=inputfile2.fastq.gz)) --a1=adapter1 --a2=adapter2 [--phred=33|64 (default=33)] [--ref=hg19 (default=hg19)] [--trim=trim_galore|fastx] [--trimq=VALUE (default=0)] [--rg_sample=STRING (default=DefaultSample)] [--rg_id=STRING (default=DefaultID)] [--rg_library=STRING (default=DefaultLib)] [--rg_platform=STRING (default=Illumina)] [--rg_center=STRING (default=Innovation_Center)] [--multithreading=INT (default=12 (max=12))] [--prefix=id_prefix]\n\n
For rg_id we suggest you to put : run#_lane#_sampleID\n"
 #      nofile => "Please specify a good entry file\n",
);
if ( (!$pair1 || !$pair2) && !$single ){
        DIE ($errors{"usg"});
}


if($ref eq "hg18"){
	print STDERR "HG18 is not yet available : TODO\n";
	exit;
}
elsif($ref eq "hg19"){
	$ref_path_bwa = $ref_hg19_bwa;
	$ref_path_fasta = $ref_hg19_fasta;
	$gtf_path = $gtf_hg19;
}
#elsif($ref eq "hg19CEU"){
#	$ref_path_bwa = $ref_hg19CEU_bwa;
#	$ref_path_fasta = $ref_hg19CEU_fasta;
##	$gtf_path = $gtf_hg19;
#}


print STDOUT "Parameters in input : \n";
if($single){
	print STDOUT "infile = $single\n";
}
else{
	print STDOUT "infiles = $pair1 $pair2\n";
}
print STDOUT "ref = $ref_path_bwa\n";
print STDOUT "rg_sample = $rg_sample\n";
print STDOUT "rg_id = $rg_id\n";
print STDOUT "rg_library = $rg_library\n";
print STDOUT "rg_platform = $rg_platform\n";
print STDOUT "rg_center = $rg_center\n";

my ($pair1_prefix,$pair2_prefix,$single_prefix);

if($single){
	$single =~ /([^\/]+)\.fastq\.gz$/;
        $single_prefix = $1;

        $id_prefix = $single_prefix;
}
else{
	$pair1 =~ /([^\/]+)\.fastq\.gz$/;
	$pair1_prefix = $1;

	$pair2 =~ /([^\/]+)\.fastq\.gz$/;
	$pair2_prefix = $1;

	if($id_prefix eq ""){
		$pair1 =~ /([^\/]+)_1\.fastq\.gz$/;
		$id_prefix = $1;
	}
	else{
	#	$id_prefix = $1;
	}
}

#PRINT ALL THIS IN A BASH FILE
my $outfile = "$id_prefix.trim$trimQ.$ref.bwa";
open(OUT,">$outfile.sh") or die("Cannot write to $outfile.sh");
print OUT "#/bin/sh\n";
print OUT "cd $currentDir\n";
print OUT "mkdir $outfile\n";

if($phred == 64){
        if($single){
                print OUT "gunzip -c $single | /opt/EMBOSS/bin/seqret fastq-illumina::stdin fastq::stdout | gzip -c > $single_prefix.phred33.fastq.gz\n";

                $single_prefix = $single_prefix . ".phred33";
                $single = "$single_prefix.fastq.gz";
        }
        else{   
                print OUT "gunzip -c $pair1 | /opt/EMBOSS/bin/seqret fastq-illumina::stdin fastq::stdout | gzip -c > $pair1_prefix.phred33.fastq.gz &\n";
                print OUT "gunzip -c $pair2 | /opt/EMBOSS/bin/seqret fastq-illumina::stdin fastq::stdout | gzip -c > $pair2_prefix.phred33.fastq.gz &\n";
                print OUT "wait\n";

                $pair1_prefix = $pair1_prefix . ".phred33";
                $pair2_prefix = $pair2_prefix . ".phred33";

                $pair1="$pair1_prefix.fastq.gz";
                $pair2="$pair2_prefix.fastq.gz";
        }
}

###FIRST STEP : REMOVE ADAPTORS
print OUT "#ADAPTORS REMOVAL\n";
if($trim_method eq "trim_galore"){
	if($single){
		print OUT "time trim_galore -q $trimQ -o $outfile --phred33 -a $adapt1 $single\n";
	}
	else{
		print OUT "time trim_galore -q $trimQ -o $outfile --phred33 -a $adapt1 -a2 $adapt2 --paired $pair1 $pair2\n";
	}
	
}
else{
	if($single){
		print OUT "gzip -c $single | fastx_clipper -v -M11 -C -a GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG | fastx_clipper -v -M8 -C -a GATCGGAAGAGCACACGTCT | fastx_clipper -v -M11 -C -a ACACTCTTTCCCTACACGACGCTCTTCCGATCT | fastq_quality_trimmer -v -t 30 -l 30 | gzip -c > $outfile/$single_prefix\_trimmed.fq.gz &\n";
	}
	else{
		print OUT "gzip -c $pair1 | fastx_clipper -v -M11 -C -a GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG | fastx_clipper -v -M8 -C -a GATCGGAAGAGCACACGTCT | fastx_clipper -v -M11 -C -a ACACTCTTTCCCTACACGACGCTCTTCCGATCT | fastq_quality_trimmer -v -t 30 -l 30 | gzip -c > $outfile/$pair1_prefix\_trimmed.fq.gz &\n";
		print OUT "gzip -c $pair2 | fastx_clipper -v -M11 -C -a GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG | fastx_clipper -v -M8 -C -a GATCGGAAGAGCACACGTCT | fastx_clipper -v -M11 -C -a ACACTCTTTCCCTACACGACGCTCTTCCGATCT | fastq_quality_trimmer -v -t 30 -l 30 | gzip -c > $outfile/$pair2_prefix.trimmed.fq.gz &\n";
		print OUT "wait\n";
	}
}

print OUT "cd $outfile\n";

if($single){
	print OUT "mv $single_prefix\_val_1.fq.gz $single_prefix\_trimmed.fq.gz\n";
#	print OUT "gzip $single_prefix\_trimmed.fq\n";
}
else{
	print OUT "mv $pair1_prefix\_val_1.fq.gz $pair1_prefix\_trimmed.fq.gz\n";
	print OUT "mv $pair2_prefix\_val_2.fq.gz $pair2_prefix\_trimmed.fq.gz\n";

#	print OUT "gzip $pair1_prefix\_trimmed.fq &\n";
#	print OUT "gzip $pair2_prefix\_trimmed.fq &\n";
#	print OUT "wait\n";
}

###2nd STEP : ALIGN WITH BWA
if($single){
	print OUT "time bwa aln -t 6 $ref_path_bwa $single_prefix\_trimmed.fq.gz > $single_prefix\_trimmed.sai\n";
        print OUT "time bwa samse $ref_path_bwa $single_prefix\_trimmed.sai $single_prefix\_trimmed.fq.gz | samtools view -bS - > $outfile.bam\n";
}
else{
	print OUT "time bwa aln -t 6 $ref_path_bwa $pair1_prefix\_trimmed.fq.gz > $pair1_prefix\_trimmed.sai &\n";
	print OUT "time bwa aln -t 6 $ref_path_bwa $pair2_prefix\_trimmed.fq.gz > $pair2_prefix\_trimmed.sai &\n";
	print OUT "wait\n";
	print OUT "time bwa sampe $ref_path_bwa $pair1_prefix\_trimmed.sai $pair2_prefix\_trimmed.sai $pair1_prefix\_trimmed.fq.gz $pair2_prefix\_trimmed.fq.gz | samtools view -bS - > $outfile.bam\n";
}

###3rd STEP : REORDERING + INDEX
print OUT "#REORDERING + INDEX\n";
my $bam_prefix = $outfile;
print OUT "java -Xmx10g -jar $picard/SortSam.jar TMP_DIR=\$LSCRATCH VALIDATION_STRINGENCY=LENIENT INPUT=$bam_prefix.bam OUTPUT=$bam_prefix.PicardSort.bam SORT_ORDER=coordinate\n";
$bam_prefix = "$bam_prefix.PicardSort";
print OUT "samtools index $bam_prefix.bam\n";

print OUT "java -Xmx10g -jar $picard/ReorderSam.jar TMP_DIR=\$LSCRATCH VALIDATION_STRINGENCY=LENIENT I=$bam_prefix.bam O=$bam_prefix.PicardReorder.bam REFERENCE=$ref_path_fasta\n";
print OUT "rm $bam_prefix.bam\n";

$bam_prefix = "$bam_prefix.PicardReorder";
print OUT "samtools index $bam_prefix.bam\n";

print OUT "echo \"samtools view -c -F 0x4 $bam_prefix.bam > $bam_prefix.stats.txt\" > $bam_prefix.stats.sh\n";
print OUT "qsub -d `pwd` -l walltime=1:00:00,nodes=1:ppn=1 $bam_prefix.stats.sh\n";


if(!$single){
	###4th STEP : REMOVE PCR DUPLICATES
	print OUT "#REMOVE PCR DUPLICATES\n";

	print OUT "time java -jar $picard/MarkDuplicates.jar TMP_DIR=\$LSCRATCH REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT INPUT=$bam_prefix.bam OUTPUT=$bam_prefix.woPCRdup.bam METRICS_FILE=$outfile.picard.metrics.txt\n";
	$bam_prefix = "$bam_prefix.woPCRdup";
	print OUT "samtools index $bam_prefix.bam\n";
	
	print OUT "echo \"samtools view -c -F 0x4 $bam_prefix.bam > $bam_prefix.stats.txt\" > $bam_prefix.stats.sh\n";
	print OUT "qsub -d `pwd` -l walltime=1:00:00,nodes=1:ppn=1 $bam_prefix.stats.sh\n";

	###5th STEP : KEEP ONLY PROPERLY PAIRED + INDEX
	print OUT "#KEEP ONLY PROPERLY PAIRED + INDEX\n";
	print OUT "samtools view -f 0x0002 -b -o $bam_prefix.PP.bam $bam_prefix.bam\n";
	$bam_prefix = "$bam_prefix.PP";
	print OUT "samtools index $bam_prefix.bam\n";

	print OUT "echo \"samtools view -c -F 0x4 $bam_prefix.bam > $bam_prefix.stats.txt\" > $bam_prefix.stats.sh\n";
	print OUT "qsub -d `pwd` -l walltime=1:00:00,nodes=1:ppn=1 $bam_prefix.stats.sh\n";

}
###6th STEP : 
print OUT "#Replace read group information in order to make GATK works well\n";
print OUT "java -Xmx10g -jar $picard/AddOrReplaceReadGroups.jar INPUT=$bam_prefix.bam OUTPUT=$bam_prefix.1RG.bam RGID=$rg_id RGLB=$rg_library RGPL=$rg_platform RGPU=barcode RGSM=$id_prefix RGCN=$rg_center VALIDATION_STRINGENCY=LENIENT\n";
$bam_prefix = "$bam_prefix.1RG";
print OUT "samtools index $bam_prefix.bam\n";

close OUT;

system("chmod u+x $outfile.sh");

sub DIE {
        my ($msg) = @_;

        TRACE ($msg, 0);
        exit (1);
}
sub TRACE {
        my ($msg, $level) = @_;
        print $msg;
}
