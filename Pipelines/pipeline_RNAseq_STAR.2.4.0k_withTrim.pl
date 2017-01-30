#!/usr/bin/perl -w

=for comment
#####HISTORIC#####

=========
2/10/2011 : First draft of the script (JCG)
5/1/2012 : GTF CHANGED : human_refGene.CurrentRelease.hg19.gtf
30/3/2012 : Single-end reads modification (AP)
6/1/2015: Changed STAR version to 2.3.1.z15 (was .z14)
=========
=cut



#use diagnostics;
use strict;
use Getopt::Long qw (:config gnu_getopt);

my $usr = $ENV{'LOGNAME'};
my $currentDir = $ENV{'PWD'};

my $picard = "/home/apps/Logiciels/picard_tools/picard-tools-1.56";

my $center = "CHUSJ";

my $single;
my $pair1;
my $pair2;
my $phred = 33;
my $rg_id = "DefaultID";
my $rg_sample = "DefaultSample";
my $rg_library = "DefaultLib";
my $rg_platform = "Illumina";
my $rg_center = "Innovation_Center";
my $p = 12;
my $max_multihit = 10;
my $ref = "hg19"; ##TODO change default ref
my $ref_path_fasta = "";
my $ref_path_star = "";
my $gtf_path = "";
my $mm = 4;
my $star = "/exec5/GROUP/barreiro/barreiro/barreiro_group/Programs/STAR_2.4.0k/STAR/bin/Linux_x86_64_static/STAR";
#my $tophat2 = "tophat2.0.11";

#CLASSIC TRUSEQ ADAPTERS
my $adapt1 = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
my $adapt2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";

my $ref_hg18_star = "";
my $ref_hg18_fasta = "";
my $gtf_hg18 = "";


########CHANGE ALL THE LINKS
#print "CHANGE THE LINKS\n";exit;
my $ref_hg19_star = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/hg19/STAR/";
my $ref_hg19_fasta = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/hg19/hs_ref_GRCh37.fasta";
my $ref_hg19CEU_star = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/hg19/CEU/STAR_ENSEMBL75";
my $ref_hg19CEU_fasta = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/hg19/CEU/CEUref_hg19.fasta";
my $gtf_hg19_old = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/hg19/GTF/human_refGene.CurrentRelease.hg19.gtf";
my $gtf_hg19 = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/hg19/GTF/human_refGene.UCSC2013-03-06-11-23-03.hg19.gtf";
my $gtf_hg19_ENSEMBL75 = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/hg19/GTF/Homo_sapiens.GRCh37.75.withCHR.ENSEMBL.stdCHR.gtf";
my $gtf = "ENSEMBL75";

my $ref_panTro4_star = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/Pan_troglodytes/Ensembl/CHIMP2.1.4/Sequence/Bowtie2/chimp2.1.4";
my $ref_panTro4_fasta = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/Pan_troglodytes/Ensembl/CHIMP2.1.4/Sequence/WholeGenomeFasta/chimp2.1.4.fa";
my $gtf_panTro4 = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/Pan_troglodytes/Ensembl/CHIMP2.1.4/Annotation/Archives/archive-2013-03-06-10-19-20/Genes/genes.gtf";

my $ref_mmul_1_star = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/Macaca_mulatta/Ensembl_iGenomes/Macaca_mulatta/Ensembl/Mmul_1/Sequence/STAR";
my $ref_mmul_1_fasta = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/Macaca_mulatta/Ensembl_iGenomes/Macaca_mulatta/Ensembl/Mmul_1/Sequence/WholeGenomeFasta/Mmul_1.fa";
my $gtf_mmul_1 = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/Macaca_mulatta/Ensembl_iGenomes/Macaca_mulatta/Ensembl/Mmul_1/Annotation/Archives/archive-2013-03-06-19-53-22/Genes/genes.gtf";

my $ref_MacaM_star ="/exec5/GROUP/barreiro/barreiro/barreiro_group/references/Macaca_mulatta/rheMac4_MacaM/STAR";
my $ref_MacaM_fasta = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/Macaca_mulatta/rheMac4_MacaM/Ref/MacaM_Rhesus_Genome_v7.fasta";
my $gtf_MacaM = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/Macaca_mulatta/rheMac4_MacaM/GTF/MacaM_Rhesus_Genome_Annotation_v7.6.8.gtf";

my $ref_mm10_star = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/mouse_mm10/Mus_musculus/UCSC/mm10/Sequence/STAR/genome";
my $ref_mm10_fasta = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/mouse_mm10/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa";
my $gtf_mm10 = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/mouse_mm10/Mus_musculus/UCSC/mm10/Annotation/Archives/archive-2013-03-06-15-06-02/Genes/genes.gtf";


my $ref_micMur1_star = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/Microcebus_murinus/ENSEMBL/STAR";
my $ref_micMur1_fasta = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/Microcebus_murinus/ENSEMBL/Microcebus_murinus.micMur1.dna_sm.toplevel.fa";
my $gtf_micMur1 = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/Microcebus_murinus/ENSEMBL/GTF/Microcebus_murinus.micMur1.78.gtf";


my $trim_method = "";
my $trimQ = 20;

my $parms = GetOptions (
	"single=s" => \$single,
	"pair1=s" => \$pair1,
	"pair2=s" => \$pair2,
	"rg_sample=s" => \$rg_sample,
	"rg_id=s" => \$rg_id,
	"rg_library=s" => \$rg_library,
	"rg_platform=s" => \$rg_platform,
	"rg_center=s" => \$rg_center,
	"max_multihit=s" => \$max_multihit,
	"ref=s"=>\$ref,
	"multithreading=s" => \$p,
        "trim=s" => \$trim_method,
        "trimq=s" => \$trimQ,
	"a1=s" => \$adapt1,
	"a2=s" => \$adapt2,
	"phred=s" => \$phred,
	"gtf=s" => \$gtf,
	"mismatches=s"=>\$mm,
	
) or die();

my %errors = (
        usg => "Usage: # $0 --single=<inputfile.fastq.gz> --pair1=<inputfile1.fastq.gz> --pair2=<inputfile2.fastq.gz> --a1=adapter1 --a2=adapter2 [--phred=33|64 (default=33)] [--trim=trim_galore|fastx] [--trimq=VALUE (default=0)] [--ref=hg18|hg19|hg19CEU (default=hg19)] [--gtf=UCSC|ENSEMBL (default=ENSEMBL for human)] [--rg_sample=STRING (default=DefaultSample)] [--rg_id=STRING (default=DefaultID)] [--rg_library=STRING (default=DefaultLib)] [--rg_platform=STRING (default=Illumina)] [--rg_center=STRING (default=Innovation_Center)] [--max_multihit=INT (default=10)] [--multithreading=INT (default=12 (max=12))] [--mismatches=INT (default=4)]\n\n
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
	$ref_path_star = $ref_hg19_star;
	$ref_path_fasta = $ref_hg19_fasta;
	($gtf eq "UCSC2013") ? $gtf_path = $gtf_hg19 : $gtf_path = $gtf_hg19_ENSEMBL75;
}
elsif($ref eq "hg19CEU"){
	$ref_path_star = $ref_hg19CEU_star;
	$ref_path_fasta = $ref_hg19CEU_fasta;
	($gtf eq "UCSC2013") ? $gtf_path = $gtf_hg19 : $gtf_path = $gtf_hg19_ENSEMBL75;
}
elsif($ref eq "micMur1"){
	$ref_path_star = $ref_micMur1_star;
	$ref_path_fasta = $ref_micMur1_fasta;
	$gtf_path = $gtf_micMur1;
}
elsif($ref eq "chimp"){
        $ref_path_star = $ref_panTro4_star;
        $ref_path_fasta = $ref_panTro4_fasta;
        $gtf_path = $gtf_panTro4;
}
elsif($ref eq "macaque"){
        $ref_path_star = $ref_mmul_1_star;
        $ref_path_fasta = $ref_mmul_1_fasta;
        $gtf_path = $gtf_mmul_1;
	$gtf = "ENSEMBL70";
}
elsif($ref eq "MacaM"){
        $ref_path_star = $ref_MacaM_star;
        $ref_path_fasta = $ref_MacaM_fasta;
        $gtf_path = $gtf_MacaM;
	$gtf = "MacaM_v7.6.8";
}
elsif($ref eq "mouse"){
        $ref_path_star = $ref_mm10_star;
        $ref_path_fasta = $ref_mm10_fasta;
        $gtf_path = $gtf_mm10;
}


print STDOUT "Parameters in input : \n";
if(!$single){print STDOUT "infiles = $pair1 $pair2\n";}
else{
	print STDOUT "infiles = $single\n";
}
print STDOUT "ref = $ref_path_star\n";
print STDOUT "rg_sample = $rg_sample\n";
print STDOUT "rg_id = $rg_id\n";
print STDOUT "rg_library = $rg_library\n";
print STDOUT "rg_platform = $rg_platform\n";
print STDOUT "rg_center = $rg_center\n";
print STDOUT "max_multihit = $max_multihit\n";

my ($pair1_prefix,$pair2_prefix,$id_prefix,$single_prefix);

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

        $pair1 =~ /([^\/]+)_R1.*\.fastq\.gz$/;
        $id_prefix = $1;
}


my $bam_prefix = "";

#NEED TO EXTRACT COLUMNS FROM SAMPLE INFO
#PRINT ALL THIS IN A BASH FILE
##TODO MODIFY FOR SINGLE END
my $outfile;

if(!$single){
	if($trim_method ne ""){
		$outfile = "$rg_id.g$max_multihit.$ref.trimQ$trimQ";
	}
	else{
		$outfile = "$rg_id.g$max_multihit.$ref";
	}
	open(OUT,">$outfile.sh") or die("Cannot write to $outfile.sh");
	print OUT "#/bin/sh\n";
	print OUT "cd $currentDir\n";
}
else{
	if($trim_method ne ""){
		$outfile = "$rg_id.$ref.trimQ$trimQ";
	}
	else{
		$outfile = "$rg_id.$ref";
	}
	open(OUT,">$outfile.sh") or die("Cannot write to $outfile.sh");
	print OUT "#/bin/sh\n";
	print OUT "cd $currentDir\n";
}

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

print OUT "#ADAPTORS REMOVAL\n";
if($trim_method ne ""){
	if($trim_method eq "trim_galore"){
        	if($single){
			print OUT "mkdir $outfile\n";
        	        print OUT "time trim_galore -q $trimQ -o $outfile --phred33 -a $adapt1 $single\n";
	        }
        	else{
			print OUT "mkdir $outfile\n";
        	        print OUT "time trim_galore -q $trimQ -o $outfile --phred33 -a $adapt1 -a2 $adapt2 --paired $pair1 $pair2\n";
	        }
	
	}
	else{
        	if($single){
	                print OUT "gzip -c $single | fastx_clipper -v -M11 -C -a GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG | fastx_clipper -v -M8 -C -a GATCGGAAGAGCACACGTCT | fastx_clipper -v -M11 -C -a ACACTCTTTCCCTACACGACGCTCTTCCGATCT | fastq_quality_trimmer -v -t 30 -l 30 | gzip -c > $outfile/$single_prefix\_trimmed.fastq.gz &\n";
        	}
	        else{
        	        print OUT "gzip -c $pair1 | fastx_clipper -v -M11 -C -a GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG | fastx_clipper -v -M8 -C -a GATCGGAAGAGCACACGTCT | fastx_clipper -v -M11 -C -a ACACTCTTTCCCTACACGACGCTCTTCCGATCT | fastq_quality_trimmer -v -t 30 -l 30 | gzip -c > $outfile/$pair1_prefix\_trimmed.fq.gz &\n";
	                print OUT "gzip -c $pair2 | fastx_clipper -v -M11 -C -a GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG | fastx_clipper -v -M8 -C -a GATCGGAAGAGCACACGTCT | fastx_clipper -v -M11 -C -a ACACTCTTTCCCTACACGACGCTCTTCCGATCT | fastq_quality_trimmer -v -t 30 -l 30 | gzip -c > $outfile/$pair2_prefix\_trimmed.fq.gz &\n";
        	        print OUT "wait\n";
	        }
	}


	if($single){
        	$single_prefix = $single_prefix . "_trimmed";
	        $single = "$single_prefix.fq.gz";
	}	
	else{
        	$pair1_prefix = $pair1_prefix . "_val_1";
	        $pair2_prefix = $pair2_prefix . "_val_2";

        	$pair1="$pair1_prefix.fq.gz";
	        $pair2="$pair2_prefix.fq.gz";
	}
	print OUT "cd $outfile\n";
}



###FIRST STEP : STAR 1pass
print OUT "#STAR\n";
#print OUT "module load STAR/2.3.1y\n";
#FIRST PASS
if(!$single){
	print OUT "time $star --genomeDir $ref_path_star --readFilesIn $pair1 $pair2 --runThreadN $p --readFilesCommand zcat --outSAMstrandField intronMotif --outFileNamePrefix $outfile.$gtf.$mm"."MM. --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFilterMismatchNmax $mm \n";

# -p $p --rg-id $rg_id --rg-sample $rg_sample --rg-library $rg_library --rg-platform-unit $rg_platform --rg-center $rg_center -G $gtf_path $ref_path_star $pair1 $pair2\n";
}
else{
	print OUT "time $star --genomeDir $ref_path_star --readFilesIn $single --runThreadN $p --readFilesCommand zcat --outSAMstrandField intronMotif --outFileNamePrefix $outfile.$gtf.$mm"."MM.  --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFilterMismatchNmax $mm \n";
}

#REF MAKING
#Treat the SJDBfile
print OUT "awk 'BEGIN {OFS=\"\\t\"; strChar[0]=\".\"; strChar[1]=\"+\"; strChar[2]=\"-\";} {if(\$5>0){print \$1,\$2,\$3,strChar[\$4]}}' $outfile.$gtf.$mm"."MM.SJ.out.tab > $outfile.$gtf.$mm"."MM.SJ.strand.out.tab\n";

if($ref eq "micMur1"){
	print OUT "mkdir genomeStar_2pass_$mm"."MM ; $star --runMode genomeGenerate --genomeDir genomeStar_2pass_$mm"."MM --genomeFastaFiles $ref_path_fasta --sjdbFileChrStartEnd $outfile.$gtf.$mm"."MM.SJ.strand.out.tab --sjdbGTFfile $gtf_path --sjdbOverhang 99 --runThreadN 12 --genomeChrBinNbits 12\n";
}
else{
	print OUT "mkdir genomeStar_2pass_$mm"."MM ; $star --runMode genomeGenerate --genomeDir genomeStar_2pass_$mm"."MM --genomeFastaFiles $ref_path_fasta --sjdbFileChrStartEnd $outfile.$gtf.$mm"."MM.SJ.strand.out.tab --sjdbGTFfile $gtf_path --sjdbOverhang 99 --runThreadN 12\n";
}
print OUT "rm $outfile.$gtf.$mm"."MM.Aligned.out.sam\n";

#Align 2 pass
if(!$single){
	print OUT "time $star --genomeDir genomeStar_2pass_$mm"."MM --readFilesIn $pair1 $pair2 --runThreadN $p --readFilesCommand zcat --outSAMstrandField intronMotif --outFileNamePrefix $outfile.2pass.$gtf.$mm"."MM. --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFilterMismatchNmax $mm --outSAMattributes NH nM NM MD --outSAMattrRGline  ID:$rg_id PU:$rg_platform PL:$rg_platform LB:$rg_library SM:$rg_sample CN:$rg_center --outSAMtype BAM SortedByCoordinate\n";
}
else{
	print OUT "time $star --genomeDir genomeStar_2pass_$mm"."MM --readFilesIn $single --runThreadN $p --readFilesCommand zcat --outSAMstrandField intronMotif --outFileNamePrefix $outfile.2pass.$gtf.$mm"."MM.  --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFilterMismatchNmax $mm --outSAMattributes NH nM MD --outSAMattrRGline  ID:$rg_id PU:$rg_platform PL:$rg_platform LB:$rg_library SM:$rg_sample CN:$rg_center --outSAMtype BAM SortedByCoordinate\n";
}

print OUT "rm -rf genomeStar_2pass_$mm"."MM\n";

$bam_prefix= "$outfile.2pass.$gtf.$mm"."MM";
###SECOND STEP : REMOVE PCR DUPLICATES
if(!$single){
	print OUT "#REMOVE PCR DUPLICATES\n";
	print OUT "time java -jar $picard/MarkDuplicates.jar REMOVE_DUPLICATES=true INPUT=$bam_prefix.bam OUTPUT=$outfile.woPCRdup.bam METRICS_FILE=$outfile.picard.metrics.txt\n";
	$bam_prefix = "$outfile.woPCRdup";
	print OUT "samtools index $bam_prefix.bam\n";
}



###THIRD STEP : REORDERING + INDEX
#if(!$single){
#	print OUT "#REORDERING + INDEX\n";
#	print OUT "java -Xmx10g -jar $picard/ReorderSam.jar I=$bam_prefix.bam O=$bam_prefix.PicardReorder.bam REFERENCE=$ref_path_fasta\n";
#	$bam_prefix = "$bam_prefix.PicardReorder";
#	print OUT "samtools index $bam_prefix.bam\n";
#}
#else{
#	$bam_prefix = $outfile;
#        print OUT "#REORDERING + INDEX\n";
#        print OUT "java -Xmx10g -jar $picard/ReorderSam.jar I=accepted_hits.bam O=$bam_prefix.PicardReorder.bam REFERENCE=$ref_path_fasta\n";
#        $bam_prefix = "$bam_prefix.PicardReorder";
#        print OUT "samtools index $bam_prefix.bam\n"
#}


###FOURTH STEP : KEEP ONLY PROPERLY PAIRED + INDEX
if(!$single){
	print OUT "#KEEP ONLY PROPERLY PAIRED + INDEX\n";
	print OUT "samtools view -f 0x0002 -b -o $bam_prefix.PP.bam $bam_prefix.bam\n";
	$bam_prefix = "$bam_prefix.PP";
	print OUT "samtools index $bam_prefix.bam\n";
}
else{}

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
