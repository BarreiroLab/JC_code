#!/usr/bin/perl
##Variant filtration pipeline
#script by thibault de Malliard and Jean-Christophe Grenier
#
#input: following arguments:
#   -root_dir <dir>              #directory where the output directory snp_calling, indel_calling or snp_indel_calling is created. Usualy it is dir where reads are present
#   -bam <bam file>              #path to the bam file
#   -calling <SNP>|<indel>|<both> #perform SNP or indel calling or both
#   -outdup                     #NOT REQUIRED: If present, output duplicates in side file
#   -woPCRdup                    #NOT REQUIRED: will remove PCR duplicates if this option is presents
#   -markPCRdup                  #NOT REQUIRED: will add PCR duplicates flag (0x400) to bam file if this option is presents
#   -q<number>                   #ONLY WHEN CALLING SNP: gets SNP with quality equal and above <number> (example: -q31   gets SNPs with quality >= 31)
#   -mpileup                    #NOT REQUIRED: will use mpileup for snp calling. Otherwise use pileup.
#   -norecal                   #NOT REQUIRED: skip the GATK recalibration
#   -reference <hg19>|<hg19CEU>|<hg18> #NOT REQUIRED : will use hg18 by default.
#   -ref <reference genome>    #NOT REQUIRED: use this fasta reference genome instead of human default one. Only use if the ref you want to use is not in the list.
#   -rnaseq			#NOT REQUIRED: use this if you are working on RNAseq data (will perform mapping qual reassignment, remove abundant sequences and only keep the uniquely mapped reads.
#   -fastq1			#ONLY WHEN rnaseq is used : give the path to the R1 fastq file for your sample
#   -fastq2			#ONLY WHEN rnaseq is used : give the path to the R2 fastq file for your sample
#   -indelRealignment 1       #NOT REQUIRED: preform GATK indel realigment before SNP or INDEL calling
#   -o <string>                #NOT REQUIRED: add a subdir to calling directory. Usefull when calling with different bam files but same root dir
#
#output:
#    * SNPs and indels Files will be all located in snp_indel_calling directory, created in the root_dir.
#    * an intermediate dir is created in this dir with all intermediate files 
#
# changelog:
# 14 dec 2010: 
#  - change intermediate dir
#  - qsub not working when launching job on node.
# 04 jan 2011:
#  - add markPCRdup option
#  - qsub working on node
# 05 jan 2011:
#  - add outdup option
# 07 jan 2011
#  - consider snp calling with no quality set (q0)
# 23 jan 2011
#  - add exome_cover option. Give exome coverage with sureselect files with snp_calling
# 28 jan 2011
#  - use /shared_scratch for some intermediate files
# 04 march 2011
#  - add -o option: add a subdir to snp_calling directory. Usefull when calling with different bam files but samr root dir
# 08 april 2011
#  - added super argument tibo to do whatever tibo wants
# 15 april 2011
#  - delete some intermediary files: SNP: pileup q0 + bam woPCRdup (keep woPCRdup_recalibrate); Indels: woPCRdup + woPCRdup_recalibrate (keep woPCRdup_recalibrate_baq)
# 20 mai 2011
#  - bioscope pipeline provide a unique bam for both snp and indel calling ==> use only 1 folder for intermediate files up to the calling
#  - added "-calling both" to call snp and indels at the same time.
# 30 mai 2011
#  - calculate the exome cover with -exome_cover argument
#  - upgrade GATK: from 1.0.4705 to 1.0.5777
#  - change dbsnp file (ROD file not supported anymore) for vcf file (v132 hg18)
# 1er juin 2011 - by JCG
#  - added GATK RealignerTargetCreator and IndelRealigner
#  - using special GATK variable ($gatk_real)
#  - added : removing unmapped reads using samtools 
# 03 june 2011
#  - add -mpileup command. Currently only perform SNP calling.
#  - remove .withPCRdup suffix when keeping PCR duplicates.
# 06 june 2011
#  - correct output path for mpileup
#  - add norecal argument (no value) to skip the reclibration
#  - add ref argument. Followed by a reference genome, will use this ref instead of default human.
# 14 june 2011
#  - update documentation
#  - remove exit command when mpileup is finished: calling both works!
# 13 october 2011
#  - update the code with the updated gatk version : use dbsnp vcf instead of rod.
#  - use gatk fixed by tibo (cigar and splice junction). This is a temp solution
#    waiting for a gatk update.
# 16 november 2011
#  - update the code for the usage of already made intervals file for the IndelRealigner
#  - update some part of the code : reference parameter now possible : specify either hg18, hg19 or hg19CEU
# 19 january 2012
#  - update for the RQCHP briarée cluster
#  - update of the reference path files and programs path files
# 24 january 2012
#  - add option to only Keep uniquely mapped reads for RNAseq
# 01 February 2012
#  - for RNAseq : put all the Mapping quality to 60 for the recalibration to work
# 10 February 2012
#  - Add part to remove the abundant sequences
# 01 march 2012
#  - Modified the options for mpileup and bcftools (added -g option for obtaining genotype likelihood and the option -p 1 to do a genotype call to everything.)
#30 march 2012
#  - Modified the options for bcftools to get the DP (DP4) in the column for every individual.
#May 1st 2012
#  - Indel calling part modified in order to launch it on the RQCHP (windows partitioning)
#Sept 26 2012
#  - Add option to compute Mito SNPs individually
#March 13th 2013
#  - Add dbsnp 137 option of hg19
#  - Add Capture kit bed files option
#  - Add 1000Genomes fasta reference, dbsnp and capture kit
#October 4th 2013
#  - Changing of GATK version (2.7.1)
#  - Add of trimming option
#MAI 9th 2014
#  - Changing of GATK version (3.1.0)
#  - Adapting script to new variant calling pipeline.
#  - Add dbsnp 138 option for hg19 and hg19CEU
#  - Add new 1000G indel ref
#  - Put the indelRealignment Reference back on
#  - Using now HC with GVCF option
#MAY 28th 2014
#  - Adding BED for Exome with Nimblegen kit
#APRIL 14th 2015
#  - Adding an option to work with WGS
#JUNE 2015
#  - Added an option for ENSEMBL GRCH37 reference

use strict;
use warnings;

use Data::Dumper;
use File::Find;
use File::Basename;
use File::Path qw(mkpath rmtree);
#$| = 1; # flush, to have stderr and stdout in same file

print "###########Parse arguments:\n";
my %args = get_args(@ARGV);
my $cpu=1;
my $mem=2;
if (exists($args{'cpu'})) {
        $cpu = $args{'cpu'};
	$mem = int($cpu*1.8);
	if($cpu==12){
		$mem=22;
	}
}

my $usr = $ENV{'LOGNAME'};
my $currentDir = $ENV{'PWD'};

my $noexec = 0;

my $cluster = 1;
my $tibo = 0;
my $noAbund = 0;
my $indelRealign = 0;
my $rnaseq = 0;
my $wgs = 0;
my $mito = 0;
my $kit = "50";
my $um = 0;
my $casava = 0;
my $trim = 0;
my $dbsnp_version = "138";
my $star = 0;
my $tophat2 = 0;
my $tophat = 0;

my $dbsnp;

########TODO########
##CHANGE ALL THE LINKS

#print "CHANGE THE LINKS\n";exit;

my $dbsnp_hg18 = '';
my $dbsnp_hg19 = "/RQusagers/$usr/barreiro_group/references/hg19/dbSNP/dbsnp_132.hg19.sort.VCF";
my $dbsnp_hg19_138 = "/RQusagers/$usr/barreiro_group/references/hg19/dbSNP/dbsnp_138.b37.excluding_sites_after_129.vcf";
my $dbsnp_hg19_138_1000G = "/RQusagers/grenier1/barreiro_group/references/hg19/dbSNP/dbsnp_138.b37.excluding_sites_after_129.NOCHR.1000G.vcf";
my $indel_1000G_138_1000G = "/RQusagers/$usr/barreiro_group/references/hg19/dbSNP/1000G_phase1.indels.b37.NOCHR.vcf";

my $dbsnp_hg19CEU = "/RQusagers/$usr/barreiro_group/references/hg19/dbSNP/dbsnp_132.hg19CEU.sort.VCF";
my $dbsnp_hg19CEU_138 = "/RQusagers/$usr/barreiro_group/references/hg19/dbSNP/dbsnp_138.b37.excluding_sites_after_129.CEU_converted.forJC.vcf";
my $indel_1000G_138_CEU = "/exec5/GROUP/awadalla/grenier2/references/hg19/CEU/Intervals/1000G_phase1.indels.b37.vcf";

my $indel_1000G_138 = "/RQusagers/$usr/barreiro_group/references/hg19/dbSNP/1000G_phase1.indels.b37.vcf";
my $indel_1000G = $indel_1000G_138_CEU;

my $kit_bed;
my $hg19_bedFile = "/RQusagers/$usr/barreiro_group/references/hg19/GTF/1transcriptPerGene/human_refGene.UCSC2011-08-30.hg19.1transcriptPerGene.sort.bed";
my $hg19_bed_50mb_1000G = "/RQusagers/$usr/barreiro_group/references/hg19/Exome/20110225.called_exome_targets.consensus.annotation.bed";
my $hg19_bed_50mb = "/RQusagers/$usr/barreiro_group/references/hg19/Exome/SureSelect_All_Exon_50mb_without_annotation.v37.bed";
my $hg19_bed_38mb = "/RQusagers/$usr/barreiro_group/references/hg19/Exome/SureSelect_All_Exon_38Mb_for_pipeline_Hg19.bed";
my $hg19_bed_nimblegen = "/RQusagers/$usr/barreiro_group/references/hg19/Exome/120430_HG19_ExomeV3_UTR_EZ_HX1_sorted_merged.bed";
my $hg19_bedFile_mito = "/RQusagers/$usr/barreiro_group/references/hg19/GTF/1transcriptPerGene/human_refGene.UCSC2012-09-18.hg19.mito.1transcriptPerGene.sort.bed";

my $hg18_bed_50mb = "/RQusagers/$usr/barreiro_group/references/hg18/Exome/SureSelect_All_Exon_TargetedRegions_50Mb_for_pipeline.Sort.bed";
my $hg18_bed_38mb = "/RQusagers/$usr/barreiro_group/references/hg18/Exome/SureSelect_All_Exon_38Mb_for_pipeline.bed";
my $abundSeqRef = "/RQusagers/$usr/barreiro_group/references/hg19/AbundantSeq/AbundantSeq_bowtie2";
my $hiMem = 0;

#my $intervals_hg18 = '/data/reference/hg18/GATK_dbSNP_ref/dbsnp_132.hg18.tibo.vcf.intervals';
my $intervals_hg19 = "/RQusagers/$usr/barreiro_group/references/hg19/Intervals/hg19_1000Genomes_indels.intervals";
my $intervals_hg19CEU = "/RQusagers/$usr/barreiro_group/references/hg19/CEU/Intervals/hg19CEU_1000Genomes_indels.intervals";


my $intervals = $intervals_hg19CEU;
#my $dbsnp = '/data/reference/hg18/GATK_dbSNP_ref/dbSNP_130_hg18.tibo.rod';
my $scratch_dir = "/RQusagers/$usr/barreiro_group/tmp"; #temp files on scratch

my $ref;
my $ref_hg18 = "/RQusagers/$usr/barreiro_group/references/hg18/Ref_samtool_pileup/hg18_c.tibo.fasta";
my $ref_hg19 = "/RQusagers/$usr/barreiro_group/references/hg19/hs_ref_GRCh37.fasta";
my $ref_hg19CEU = "/RQusagers/$usr/barreiro_group/references/hg19/CEU/CEUref_hg19.fasta";
my $ref_hg191000G = "/RQusagers/$usr/barreiro_group/references/hg19/1000Genomes/human_g1k_v37.fasta";
my $ref_hg19ENSEMBL = "/exec5/GROUP/barreiro/barreiro/barreiro_group/Yohann/Projects/MLS/RNA-seq/Data/Homo_sapiens.GRCh37.75.dna.primary_assembly.reordered.fa";

my $gatk3_1_1 = "java -Djava.io.tmpdir=\$LSCRATCH -Xmx$mem"."g -jar /home/apps/Logiciels/GATK/3.1-1/GenomeAnalysisTK.jar -l INFO";
my $gatk3_1_1_real = "java -Djava.io.tmpdir=\$LSCRATCH -Xmx$mem"."g -jar /home/apps/Logiciels/GATK/3.1-1/GenomeAnalysisTK.jar";

my $gatk3 = "java -Djava.io.tmpdir=\$LSCRATCH -Xmx$mem"."g -jar /home/apps/Logiciels/GATK/3.3-0/GenomeAnalysisTK.jar -l INFO";
my $gatk3_real = "java -Djava.io.tmpdir=\$LSCRATCH -Xmx$mem"."g -jar /home/apps/Logiciels/GATK/3.3-0/GenomeAnalysisTK.jar";

my $gatk2 = "java -Djava.io.tmpdir=\$LSCRATCH -Xmx$mem"."g -jar /home/apps/Logiciels/GATK/2.7-1/GenomeAnalysisTK.jar -l INFO --default_platform ILLUMINA";
my $gatk2_real = "java -Djava.io.tmpdir=\$LSCRATCH -Xmx$mem"."g -jar /home/apps/Logiciels/GATK/2.7-1/GenomeAnalysisTK.jar";

my $gatk1 = "java -Xmx$mem"."g -jar /home/apps/Logiciels/GATK/GenomeAnalysisTK-1.3-14-g348f2db/GenomeAnalysisTK.jar -l INFO --default_platform ILLUMINA --default_read_group DefaultReadGroup";
my $gatk1_real = "java -Xmx$mem"."g -jar /home/apps/Logiciels/GATK/GenomeAnalysisTK-1.3-14-g348f2db/GenomeAnalysisTK.jar";

my $gatk = $gatk3;
my $gatk_real = $gatk3_real;

my $picard = "java -Xmx$mem"."g -jar /home/apps/Logiciels/picard_tools/picard-tools-1.56 TMP_DIR=\$LSCRATCH CREATE_INDEX=true";
my $picard_dup = "java -Xmx$mem"."g -jar /home/apps/Logiciels/picard_tools/MarkDuplicates_dup.jar TMP_DIR=\$LSCRATCH CREATE_INDEX=true";
my $dindel = 'dindel';
my $bowtie2 = "/exec5/GROUP/barreiro/barreiro/barreiro_group/bin/bowtie2";
my $bamUtil = "/RQusagers/$usr/barreiro_group/Programs/bamUtil_1.0.2/bamUtil/bin/bam";
my $output_dir;
my $intermediate_dir;

my $cmd = '';
#################


# print Dumper \%args;
print "DONE\n";

####################
print "###########verify arguments\n";

if (exists($args{'help'})) {
	print_help();
}
if(exists($args{'himem'})){
	$hiMem = 1;
}
if (exists($args{'tibo'})) {
	print "**** Super argument <tibo> is set ****\n";
	$tibo = 1;
}
if( exists($args{'noabund'})){
	if(!exists($args{'fastq1'}) || !exists($args{'fastq2'})){
		print STDERR "Bad Fastq files!\n";
	        exit 0;
	}
	else{
		my @fastq1 = split(/,/,$args{'fastq1'});
		my @fastq2 = split(/,/,$args{'fastq2'});
		my @fastq12 = (@fastq1,@fastq2);
		
		foreach my $fastq( @fastq12 ){
			if( !-e $fastq ){
				print STDERR "Bad Fastq files : $fastq!\n";
		                exit 0;
			}
		}
		print "**** Will remove the abundant sequences first ****\n";
	        $noAbund = 1;
	}
}
if( !exists($args{'reference'}) && !exists($args{'ref'}) ){
	$ref = $ref_hg19;
	$dbsnp = $dbsnp_hg19;
}
if (!exists($args{'calling'}) && exists($args{'nocalling'})) {
	print STDOUT "No calling set, preparation of the files will be done instead.\n";
	$output_dir = 'variant_calling_preparation';
        if(exists($args{'mito'})){
                $output_dir = 'variant_calling_preparation_mito';
        }
}
if( $args{'reference'} eq 'hg18'){
        $ref = $ref_hg18;
        $dbsnp = $dbsnp_hg18;
	if( $args{'kit'} eq '38'){
		$kit_bed = $hg18_bed_38mb;
	}
	else{
		$kit_bed = $hg18_bed_50mb;
	}
#	$intervals = $intervals_hg18;
}
if($args{'reference'} eq '1000G'){
	$ref = $ref_hg191000G;
	$dbsnp = $dbsnp_hg19_138_1000G;
	$kit_bed = $hg19_bed_50mb_1000G;
	$indel_1000G = $indel_1000G_138_1000G;
}
if( $args{'reference'} eq 'hg19ENSEMBL'){
        $ref = $ref_hg19ENSEMBL;
        $dbsnp = $dbsnp_hg19_138_1000G;
	$kit_bed = $hg19_bed_50mb_1000G;
        $indel_1000G = $indel_1000G_138_1000G;
}
if( $args{'reference'} eq 'hg19'){
        $ref = $ref_hg19;
	if($args{'dbsnp_version'} eq "138"){
		$dbsnp = $dbsnp_hg19_138;
	}
	else{
	        $dbsnp = $dbsnp_hg19;
	}
	if( lc($args{'kit'}) eq '38'){
                $kit_bed = $hg19_bed_38mb;
        }
        elsif(lc($args{'kit'}) eq '50'){
                $kit_bed = $hg19_bed_50mb;
        }
	elsif(lc($args{'kit'}) eq 'nimblegen'){
                $kit_bed = $hg19_bed_nimblegen;
        }
	$intervals = $intervals_hg19;
	$indel_1000G = $indel_1000G_138;
}
if( $args{'reference'} eq 'hg19CEU'){
        $ref = $ref_hg19CEU;
	if($args{'dbsnp_version'} eq "138"){
                $dbsnp = $dbsnp_hg19CEU_138;
        }
	else{
	        $dbsnp = $dbsnp_hg19CEU;
	}
	if( lc($args{'kit'}) eq '38'){
                $kit_bed = $hg19_bed_38mb;
        }
        elsif(lc($args{'kit'}) eq '50'){
                $kit_bed = $hg19_bed_50mb;
        }
	elsif(lc($args{'kit'}) eq 'nimblegen'){
                $kit_bed = $hg19_bed_nimblegen;
        }
	$intervals = $intervals_hg19CEU;
	$indel_1000G = $indel_1000G_138_CEU;
}
if (exists($args{'calling'}) && $args{'calling'} eq 'snp') {
	$output_dir = 'snp_calling';
	if(exists($args{'mito'})){
		$output_dir = 'snp_calling_mito';
	}
}
#if (exists($args{'calling'}) && $args{'calling'} eq 'indel') {
#	$output_dir = 'indel_calling';
#	if(exists($args{'mito'})){
#                $output_dir = 'indel_calling_mito';
#        }
#}
if (exists($args{'calling'}) && $args{'calling'} eq 'both') {
	$output_dir = 'snp_indel_calling';
	if(exists($args{'mito'})){
                $output_dir = 'snp_indel_calling_mito';
        }
}
if (exists($args{'indelRealignment'})) {
        print "**** IndelLocalRealignment with GATK will be done! ****\n";
        $indelRealign = 1;
}
if (exists($args{'wgs'})){
	print "*** Bam will be treated as a WGS ***\n";
	$wgs = 1;
}
if (exists($args{'rnaseq'})){
	print "*** Bam will be treated as RNAseq ***\n";
	$rnaseq = 1;
	if (exists($args{'star'})){
		$star =1;
	}
	if (exists($args{'tophat2'})){
		$tophat2 =1;
	}
	elsif(exists($args{'tophat'})){
		$tophat =1;
	}
}
if (exists($args{'mito'})){
        print "*** Computation of mitochondrial SNPs only ***\n";
        $mito = 1;
	$hg19_bedFile = $hg19_bedFile_mito;
}
if (exists($args{'o'})) {
	$output_dir .= "/$args{'o'}";
}
if (exists($args{'casava'})){
        print "**** Input was aligned and will be treated like CASAVA bam file for UM filter****\n";
        $casava = 1;
}
if (exists($args{'um'})){
        $um = 1;
}
if ( exists($args{'trim'})){
        $trim = 1;
}

if (exists($args{'ref'}) && !exists($args{'reference'}) ) {
	$ref = $args{'ref'};
}

if (!-e $args{'bam'}) {
	print STDERR "Bam file not found\n";
	exit 0;
}
if (!-r $args{'bam'}) {
	print STDERR "Could not read bam file\n";
	exit 0;
}
$args{'bam'} =~ /([^\/]+)\.bam$/;
my $bam_name = $1;
	
if (!-d $args{'root_dir'}) {
	print STDERR "Root directory not found\n";
	exit 0;
}
if (!-w $args{'root_dir'}) {
	print STDERR "Could not write in root directory\n";
	exit 0;
}


# if (!exists($args{'q'}) and $args{'calling'} eq 'SNP') {
# 	print STDERR "no quality value set\n";
# 	exit 0;
# }

print "DONE\n";
#######################
print "###########Verify needed ressources\n";

if (!-e $ref) {
	print STDERR "Reference file not found at $ref\n";
	exit 0;
}
if (!-r $ref) {
	print STDERR "Could not open $ref\n";
	exit 0;
}

if (!-e $dbsnp) {
	print STDERR "dbSNP file not found at $dbsnp\n";
	exit 0;
}
if (!-r $dbsnp) {
	print STDERR "Could not open $dbsnp\n";
	exit 0;
}

print "DONE\n";
###############################
print "###########create result and intermediary directories\n";

$output_dir = "$args{'root_dir'}/$output_dir";
# $intermediate_dir = "$intermediate_dir/$bam_name\_$^T\_$$"; #use custom and unique directory name
$intermediate_dir = "$output_dir/intermediate"; #use custom and unique directory name

if (!-d $output_dir) {
	mkpath($output_dir) or die "could not create $output_dir\n";
}
if (!-d $intermediate_dir) {
	mkpath($intermediate_dir) or die "could not create $intermediate_dir\n";
}

print "results saved in $output_dir\n";
print "Intermediate filed will be in $intermediate_dir\n";
print "DONE\n";

#### log to screen AND output dir
#open (STDOUT, "| tee -ai $output_dir/log.txt");
#open (STDERR, "| tee -ai $output_dir/log.txt");
print "###########################\n";
print `date +'%Y/%m/%d %H:%M:%S'`;
print "@ARGV\n";
print "###########################\n";
if( $args{'bam'} =~ m/woPCRdup\.PicardReorder\.PP\.noAbundant\.UM\.MQ60\.recalibrate\.bam/g){

}
else{
if($noAbund == 1 && $mito == 0){
	chdir "$args{'root_dir'}";
	print "####Remove abundant sequences!\n";
	$cmd = "$bowtie2 -p $cpu $abundSeqRef -1 $args{'fastq1'} -2 $args{'fastq2'} -S $bam_name.abundRef.sam";
	print "\t$cmd\n";
	print `$cmd` unless $noexec ==1;

	$cmd = "samtools view -bS -F4 -o $args{'root_dir'}/$bam_name.abundRef.bam $args{'root_dir'}/$bam_name.abundRef.sam";
	print "\t$cmd\n";
        print `$cmd` unless $noexec ==1;

	$cmd = "rm $args{'root_dir'}/$bam_name.abundRef.sam";
	print "\t$cmd\n";
        print `$cmd` unless $noexec ==1;

	$cmd = "samtools sort $args{'root_dir'}/$bam_name.abundRef.bam $args{'root_dir'}/$bam_name.abundRef.sort";
        print "\t$cmd\n";
        print `$cmd` unless $noexec ==1;

	$cmd = "samtools index $args{'root_dir'}/$bam_name.abundRef.sort.bam";
	print "\t$cmd\n";
        print `$cmd` unless $noexec ==1;
	
	$cmd = "samtools view $args{'root_dir'}/$bam_name.abundRef.sort.bam | awk '{print \$1}' > $args{'root_dir'}/ids_to_rem.txt";
	print "\t$cmd\n";
        print `$cmd` unless $noexec ==1;

	if($hiMem == 0){
		$cmd = "samtools view -h $args{'bam'} | grep -Fvf $args{'root_dir'}/ids_to_rem.txt | samtools view -bS -> $args{'root_dir'}/$bam_name.noAbundant.bam ; samtools index $args{'root_dir'}/$bam_name.noAbundant.bam ; rm $args{'root_dir'}/ids_to_rem.txt";
	}
	else{
		my $nb_line_split = (split(/ /,`wc -l $args{'root_dir'}/ids_to_rem.txt`))[0];
		$nb_line_split = int($nb_line_split/2)+1;
		system("split -d -l $nb_line_split $args{'root_dir'}/ids_to_rem.txt $args{'root_dir'}/ids_to_rem.");
		$cmd = "samtools view -h $args{'bam'} | grep -Fvf $args{'root_dir'}/ids_to_rem.00 | samtools view -bS -> $args{'root_dir'}/$bam_name.noAbundant00.bam ; samtools index $args{'root_dir'}/$bam_name.noAbundant00.bam ; rm $args{'root_dir'}/ids_to_rem.{00,txt} ; samtools view -h $args{'root_dir'}/$bam_name.noAbundant00.bam | grep -Fvf $args{'root_dir'}/ids_to_rem.01 | samtools view -bS -> $args{'root_dir'}/$bam_name.noAbundant.bam ; samtools index $args{'root_dir'}/$bam_name.noAbundant.bam ; rm $args{'root_dir'}/ids_to_rem.01 ; rm $args{'root_dir'}/$bam_name.noAbundant00.bam $args{'root_dir'}/$bam_name.noAbundant00.bam.bai";
	}
	print "\t$cmd\n";
        print `$cmd` unless $noexec ==1;

	$args{'bam'} = "$args{'root_dir'}/$bam_name.noAbundant.bam";
	$bam_name = "$bam_name.noAbundant";
}
if($mito == 1){
	print "####Extract mito uniquely mapped reads\n";
        $cmd = "samtools view -hb $args{'bam'} chrM > $intermediate_dir/$bam_name.mito.bam";
        print "\t$cmd\n";
        print `$cmd` unless $noexec == 1;

	$args{'bam'} = "$intermediate_dir/$bam_name.mito.bam";
        $bam_name = "$bam_name.mito";

}
if($rnaseq == 1){
	if($star == 1){
	#SPLIT N TRIM PART FOR RNASEQ
		print "####SPLIT N TRIM and REASSIGN MAPPING QUALITY TO 60 for reads (STAR)\n";
		$cmd = "$gatk -T SplitNCigarReads -R $ref -I  $args{'bam'} -o $intermediate_dir/$bam_name.SNT.MQ60.bam -rf ReassignMappingQuality -DMQ 60 -U ALLOW_N_CIGAR_READS";
		print "\t$cmd\n";
		print `$cmd` unless $noexec == 1;

		$args{'bam'} = "$intermediate_dir/$bam_name.SNT.MQ60.bam";
		$bam_name = "$bam_name.SNT.MQ60";
	}
	if($tophat == 1){
		print "#####Put all the MQ to 60 for GATK to be working properly (tophat 1.4)\n";
		$cmd = "$gatk -T PrintReads -I $args{'bam'} -rf ReassignMappingQuality -DMQ 60 -o $intermediate_dir/$bam_name.MQ60.bam -R $ref -U ALLOW_N_CIGAR_READS";
		
	#	print "\t$cmd\n";
	#	print `$cmd` unless $noexec == 1;
	#	$args{'bam'} = "$intermediate_dir/$bam_name.MQ60.bam";
		$bam_name = "$bam_name.MQ60";
	}
	if($um == 1){
	        print "####Remove non-uniquely mapped reads\n";  ###MODIFY FOR STAR OR TOPHAT2
        	$cmd = "samtools view -h $args{'bam'} | grep -P \"NH:i:1\\t|^@\" | samtools view -bS - > $intermediate_dir/$bam_name.UM.bam";
	        print "\t$cmd\n";
        	print `$cmd` unless $noexec == 1;
        
	        $args{'bam'} = "$intermediate_dir/$bam_name.UM.bam";
        	$bam_name = "$bam_name.UM";

	        print "###########Make the index of the new bam file\n";
        	$cmd = "samtools index $args{'bam'}";
	        print "\t$cmd\n";
        	print `$cmd` unless $noexec == 1;
	}
}
else{
	if($um == 1){
		print "####Remove non-uniquely mapped reads\n";
		if($casava == 1){
                        $cmd = "samtools view -b -F256 $args{'bam'} > $intermediate_dir/$bam_name.UM.bam";
                }
                else{
                        $cmd = "samtools view -h $args{'bam'} | grep -P \"XT:A:U\$|^@\" | samtools view -bS - > $intermediate_dir/$bam_name.UM.bam";
                }
	        print "\t$cmd\n";
	        print `$cmd` unless $noexec == 1;

        	$args{'bam'} = "$intermediate_dir/$bam_name.UM.bam";
	        $bam_name = "$bam_name.UM";

	        print "###########Make the index of the new bam file\n";
        	$cmd = "samtools index $args{'bam'}";
	        print "\t$cmd\n";
        	print `$cmd` unless $noexec == 1;
	}


}
###############################
if (!exists($args{'wopcrdup'}) && !exists($args{'removepcrdup'}) && $rnaseq!=1) {
        print "###########Keep PCR duplicates\n";
#       $bam_name = "$bam_name.withPCRdup";
} 
else{
        if (exists($args{'wopcrdup'}) && $rnaseq!=1) {
################################
                print "###########Remove PCR duplicate\n";

                if (exists($args{'outdup'})) {
                        $picard = "$picard_dup OUTDUP=$output_dir/$bam_name.PCRdup.bam";
                }

                $cmd = "$picard REMOVE_DUPLICATES=true INPUT=$args{'bam'} OUTPUT=$intermediate_dir/$bam_name.woPCRdup.bam METRICS_FILE=$intermediate_dir/$bam_name.picard.metrics.txt";

                print "\t$cmd\n";
                print `$cmd` unless $noexec == 1;

                $args{'bam'} = "$intermediate_dir/$bam_name.woPCRdup.bam";
                $bam_name = "$bam_name.woPCRdup";

                print "DONE\n";
        } 
	else{
		if($rnaseq !=1){
###############################
                print "###########mark PCR duplicate\n";

                $cmd = "$picard REMOVE_DUPLICATES=false INPUT=$args{'bam'} OUTPUT=$intermediate_dir/$bam_name.markPCRdup.bam METRICS_FILE=$intermediate_dir/$bam_name.picard.metrics.txt";

                print "\t$cmd\n";
                print `$cmd` unless $noexec == 1;

                $args{'bam'} = "$intermediate_dir/$bam_name.markPCRdup.bam";
                $bam_name = "$bam_name.markPCRdup";

                print "DONE\n";
		}
        }
}

if($indelRealign == 1){
	if($rnaseq!=1){
	        print "###########Remove the unmapped reads\n";
        	$cmd = "samtools view -b -F4 $args{'bam'} > $intermediate_dir/$bam_name.onlyMapped.bam";
	        print "\t$cmd\n";
        	print `$cmd` unless $noexec == 1;

	        $args{'bam'} = "$intermediate_dir/$bam_name.onlyMapped.bam";
        	$bam_name = "$bam_name.onlyMapped";

	        print "###########Make the index of the new bam file\n";
        	$cmd = "samtools index $args{'bam'}";
	        print "\t$cmd\n";
        	print `$cmd` unless $noexec == 1;
	}
	print "###########Make the Target Intervals first\n";
		
	if($intervals eq 'NA'){
	       	$cmd = "$gatk_real -R $ref -T RealignerTargetCreator -o $intermediate_dir/$bam_name.intervals -I $args{'bam'} --validation_strictness LENIENT";
		$intervals = "$intermediate_dir/$bam_name.intervals";
	      	$cmd .= " -known $dbsnp" if $ref =~ /hg18_c/;
		$cmd .= " -known $indel_1000G_138" if $ref =~ /GRCh37/;
		$cmd .= " -known $indel_1000G_138" if $ref =~ /CEU/;
		$cmd .= " -U ALLOW_N_CIGAR_READS" if $rnaseq ==1;
	        #$cmd .= " -known $dbsnp_pfalci" if $ref !~ /hg18_c/;
	        #$cmd .= " -L \"$targetIntervals\"" if($targetIntervals ne "NA"); 
	        print "\t$cmd\n";
       		print `$cmd` unless $noexec == 1;
	}

        print "###########Make the IndelLocalRealignment!\n";

	if($rnaseq == 1){
		$cmd = "$gatk_real -I $args{'bam'} -R $ref -T IndelRealigner -o $intermediate_dir/$bam_name.realignedGATK.bam -targetIntervals $intervals -known $indel_1000G --consensusDeterminationModel KNOWNS_ONLY -LOD 0.4 --validation_strictness LENIENT -U ALLOW_N_CIGAR_READS";
	}
	else{
	        $cmd = "$gatk_real -I $args{'bam'} -R $ref -T IndelRealigner -o $intermediate_dir/$bam_name.realignedGATK.bam -targetIntervals $intervals --consensusDeterminationModel KNOWNS_ONLY -LOD 0.4 --validation_strictness LENIENT";
	}
	$cmd .= " -known $dbsnp" if $ref =~ /hg18_c/;
	$cmd .= " -known $indel_1000G" if $ref =~ /hg19/;
        #$cmd .= " -known $dbsnp_pfalci" if $ref !~ /hg18_c/;
        print "\t$cmd\n";
        print `$cmd` unless $noexec == 1;

        $args{'bam'} = "$intermediate_dir/$bam_name.realignedGATK.bam";
        $bam_name = "$bam_name.realignedGATK";
}
###############################
if (!exists($args{'norecal'})) {
	#####################
	print "###########Recalibrate alignements\n";
		
	$cmd = "$gatk -T BaseRecalibrator -I $args{'bam'} -R $ref -knownSites $dbsnp -o $intermediate_dir/$bam_name.recal_data.table -nct $cpu";
	$cmd .= " -U ALLOW_N_CIGAR_READS" if $rnaseq ==1;
	print "\t$cmd\n";
	print `$cmd` unless $noexec == 1;
	$cmd = "$gatk_real -T PrintReads -R $ref -I $args{'bam'} -BQSR $intermediate_dir/$bam_name.recal_data.table -o $intermediate_dir/$bam_name.recalibrate.bam -nct $cpu";
	$cmd .= " -U ALLOW_N_CIGAR_READS" if $rnaseq ==1;
	print "\t$cmd\n";
	print `$cmd` unless $noexec == 1;
	if ( ($tibo == 1) && ($noexec != 1) ) {
		unlink($args{'bam'});
		unlink("$args{'bam'}.bai");
	}
		
	$args{'bam'} = "$intermediate_dir/$bam_name.recalibrate.bam";
	$bam_name = "$bam_name.recalibrate";


	print "DONE\n";
}
if($trim == 1){
        $cmd = "$bamUtil trimBam $args{'bam'} $intermediate_dir/$bam_name.trim2baseEnd.bam 2 ; samtools index $intermediate_dir/$bam_name.trim2baseEnd.bam";
        print "\t$cmd\n";
        print `$cmd` unless $noexec == 1;
        $args{'bam'} = "$intermediate_dir/$bam_name.trim2baseEnd.bam";
        $bam_name = "$bam_name.trim2baseEnd";
}
}	### END ELSE FOR THE MATCH "woPCRdup.PicardReorder.PP.noAbundant.UM.MQ60.recalibrate.bam" 	############BE AWARE############
################################
if($args{'nocalling'}){
        exit;
}


if ( ($args{'calling'} eq 'snp') or ($args{'calling'} eq 'both') ) {

	print "############GATK HC in GVCF MODE indel + snps\n";
	if($rnaseq!=1 && $wgs!=1){	
		$cmd = "$gatk -pairHMM VECTOR_LOGLESS_CACHING -T HaplotypeCaller -I $args{'bam'} -R $ref --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o $output_dir/$bam_name.g.vcf --dbsnp $dbsnp -stand_call_conf 30.0 -stand_emit_conf 0.0 -nct $cpu -L $kit_bed -dontUseSoftClippedBases -mmq 1";
	}
	elsif($wgs==1 && $rnaseq!=1){
		$cmd = "$gatk -pairHMM VECTOR_LOGLESS_CACHING -T HaplotypeCaller -I $args{'bam'} -R $ref --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o $output_dir/$bam_name.g.vcf --dbsnp $dbsnp -stand_call_conf 30.0 -stand_emit_conf 0.0 -nct $cpu  -mmq 1 -U ALLOW_N_CIGAR_READS";
	}
	else{
		$cmd = "$gatk -pairHMM VECTOR_LOGLESS_CACHING -T HaplotypeCaller -I $args{'bam'} -R $ref --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o $output_dir/$bam_name.g.vcf --dbsnp $dbsnp -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -nct $cpu  -mmq 1 -U ALLOW_N_CIGAR_READS";
	}
	print "$cmd\n";
               print `$cmd` unless $noexec == 1;

	print "DONE\n";
}
#########################################
close (STDOUT);
close (STDERR);
################################

sub get_args{
	my @args = @_;
	my $current_arg;
	my %result_args;
	
	foreach my $arg (@args) {
		if ($arg =~ /^-(\w+)$/ ) {
			$current_arg = lc $1;
			if ($1 =~ /^q(\d+)$/) {
				push @{ $result_args{'q'} }, $1; 
				next;
			}
			if(uc($current_arg) eq 'STAR'){
				$result_args{'star'}=1;
			}
			if(uc($current_arg) eq 'TOPHAT'){
                                $result_args{'tophat'}=1;
                        }
			if(uc($current_arg) eq 'TOPHAT2'){
                                $result_args{'tophat2'}=1;
                        }
			if (uc($current_arg) eq 'TRIM'){
                                $result_args{'trim'} = 1;
                        }
                        if(uc($current_arg) eq 'NOCALLING'){
                                $result_args{'nocalling'} =1;
                        }
			if (uc($current_arg) eq 'UM'){
				$result_args{'um'} = 1;
			}
			if (uc($current_arg) eq 'CASAVA'){
                                $result_args{'casava'} = 1;
                        }
			if (uc($current_arg) eq 'NOABUND'){
                                $result_args{'noabund'} = 1;
                        }
			if (uc($current_arg) eq 'MITO'){
                                $result_args{'mito'} = 1;
                        }
			if (uc($current_arg) eq 'HIMEM'){
                                $result_args{'himem'} = 1;
                        }
			if (uc($current_arg) eq 'RNASEQ'){
				$result_args{'rnaseq'} = 1;
			}
			if (uc($current_arg) eq 'WGS'){
                                $result_args{'wgs'} = 1;
                        }
			if (uc($current_arg) eq 'REFERENCE'){
				$result_args{'reference'} = 1;
			}
			if (uc($current_arg) eq 'INDELREALIGNMENT') {
				$result_args{'indelRealignment'} = 1;
				next;
			}
			if (uc($current_arg) eq 'WOPCRDUP') {
				$result_args{'wopcrdup'} = 1;
				next;
			} 
			if (uc($current_arg) eq 'MARKPCRDUP') {
				$result_args{'markpcrdup'} = 1;
				next;
			} 
			if (uc($current_arg) eq 'OUTDUP') {
				$result_args{'outdup'} = 1;
				next;
			} 
			if (uc($current_arg) eq 'EXOME_COVER') {
				$result_args{'exome_cover'} = 1;
				next;
			} 
			if (uc($current_arg) eq 'MPILEUP') {
				$result_args{'mpileup'} = 1;
				next;
			} 
			if (uc($current_arg) eq 'NORECAL') {
				$result_args{'norecal'} = 1;
				next;
			} 
			if ( (uc($current_arg) eq 'HELP') or (uc($current_arg) eq 'H') ) {
				$result_args{'help'} = 1;
				next;
			} 
		} else {
			$arg =~ /^'?([^']+)'?$/;
			$result_args{$current_arg} = $1;
		}
	}
	
	return %result_args;
}
