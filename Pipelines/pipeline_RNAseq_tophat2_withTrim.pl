#!/usr/bin/perl -w

=for comment
#####HISTORIC#####

=========
2/10/2011 : First draft of the script (JCG)
5/1/2012 : GTF CHANGED : human_refGene.CurrentRelease.hg19.gtf
30/3/2012 : Single-end reads modification (AP)

=========
=cut



#use diagnostics;
use strict;
use Getopt::Long qw (:config gnu_getopt);

my $usr = $ENV{'LOGNAME'};
my $currentDir = $ENV{'PWD'};

my $picard = "/home/apps/Logiciels/picard_tools/picard-tools-1.56";

my $center = "CHUM";

my $single;
my $pair1;
my $pair2;
my $phred = 33;
my $rg_id = "DefaultID";
my $rg_sample = "DefaultSample";
my $rg_library = "DefaultLib";
my $rg_platform = "Illumina";
my $rg_platform_unit = "Illumina";
my $rg_center = "Innovation_Center";
my $p = 12;
my $max_multihit = 20;
my $mate_inner_dist = 200;
my $ref = "hg19"; ##TODO change default ref
my $ref_path_fasta = "";
my $ref_path_bowtie = "";
my $gtf_path = "";
my $cuff = 0;

#CLASSIC TRUSEQ ADAPTERS
my $adapt1 = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
my $adapt2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";

my $ref_hg18_bowtie = "";
my $ref_hg18_fasta = "";
my $gtf_hg18 = "";

my $ref_hg19_bowtie = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/hg19/bowtie2/hs_ref_GRCh37";
my $ref_hg19_fasta = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/hg19/hs_ref_GRCh37.fasta";
my $ref_hg19CEU_bowtie = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/hg19/CEU/bowtie2/CEUref_hg19";
my $ref_hg19CEU_fasta = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/hg19/CEU/CEUref_hg19.fasta";
my $gtf_hg19_old = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/hg19/GTF/human_refGene.CurrentRelease.hg19.gtf";
my $gtf_hg19_UCSC = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/hg19/GTF/human_refGene.UCSC2013-03-06-11-23-03.hg19.gtf";
my $gtf_hg19 = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/hg19/GTF/Homo_sapiens.GRCh37.75.withCHR.ENSEMBL.stdCHR.gtf";

my $ref_panTro4_bowtie = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/Pan_troglodytes/Ensembl/CHIMP2.1.4/Sequence/Bowtie2/chimp2.1.4";
my $ref_panTro4_fasta = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/Pan_troglodytes/Ensembl/CHIMP2.1.4/Sequence/WholeGenomeFasta/chimp2.1.4.fa";
my $gtf_panTro4 = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/Pan_troglodytes/Ensembl/CHIMP2.1.4/Annotation/Archives/archive-2013-03-06-10-19-20/Genes/genes.gtf";

my $ref_mmul_1_bowtie = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/Macaca_mulatta/Ensembl_iGenomes/Macaca_mulatta/Ensembl/Mmul_1/Sequence/Bowtie2/Mmul_1";
my $ref_mmul_1_fasta = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/Macaca_mulatta/Ensembl_iGenomes/Macaca_mulatta/Ensembl/Mmul_1/Sequence/WholeGenomeFasta/Mmul_1.fa";
my $gtf_mmul_1 = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/Macaca_mulatta/Ensembl_iGenomes/Macaca_mulatta/Ensembl/Mmul_1/Annotation/Archives/archive-2013-03-06-19-53-22/Genes/genes.gtf";

my $ref_mm10_bowtie = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/mouse_mm10/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome";
my $ref_mm10_fasta = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/mouse_mm10/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa";
my $gtf_mm10 = "/exec5/GROUP/barreiro/barreiro/barreiro_group/references/mouse_mm10/Mus_musculus/UCSC/mm10/Annotation/Archives/archive-2013-03-06-15-06-02/Genes/genes.gtf";

my $ref_sheep_bowtie="/exec5/GROUP/barreiro/barreiro/barreiro_group/references/sheep_oviAri3/NCBI/bowtie2/oar_ref_Oar_v3.1.masked.goodNames";
my $ref_sheep_fasta="/exec5/GROUP/barreiro/barreiro/barreiro_group/references/sheep_oviAri3/NCBI/oar_ref_Oar_v3.1.masked.goodNames.fa";
my $gtf_sheep="/exec5/GROUP/barreiro/barreiro/barreiro_group/references/sheep_oviAri3/ENSEMBL/GTF/Ovis_aries.Oar_v3.1.75.gtf";

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
	"rg_platform-unit"=>\$rg_platform_unit,
	"rg_center=s" => \$rg_center,
	"max_multihit=s" => \$max_multihit,
	"mate_inner_dist=s" => \$mate_inner_dist,
	"ref=s"=>\$ref,
	"multithreading=s" => \$p,
	"cuff=s"=>\$cuff,
        "trim=s" => \$trim_method,
        "trimq=s" => \$trimQ,
	"a1=s" => \$adapt1,
	"a2=s" => \$adapt2,
	"phred=s" => \$phred,
	
) or die();

##TODO MODIFY FOR SINGLE END
my %errors = (
        usg => "Usage: # $0 --single=<inputfile.fastq.gz> --pair1=<inputfile1.fastq.gz> --pair2=<inputfile2.fastq.gz> --a1=adapter1 --a2=adapter2 [--phred=33|64 (default=33)] [--trim=trim_galore|fastx] [--trimq=VALUE (default=0)] [--ref=hg18|hg19|hg19CEU (default=hg19)] [--cuff=0|1 (default=0)] [--rg_sample=STRING (default=DefaultSample)] [--rg_id=STRING (default=DefaultID)] [--rg_library=STRING (default=DefaultLib)] [--rg_platform=STRING (default=Illumina)] [--rg_center=STRING (default=Innovation_Center)] [--max_multihit=INT (default=20)] [--mate_inner_dist=INT (default=200)] [--multithreading=INT (default=12 (max=12))]\n\n
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
	$ref_path_bowtie = $ref_hg19_bowtie;
	$ref_path_fasta = $ref_hg19_fasta;
	$gtf_path = $gtf_hg19;
}
elsif($ref eq "hg19CEU"){
	$ref_path_bowtie = $ref_hg19CEU_bowtie;
	$ref_path_fasta = $ref_hg19CEU_fasta;
	$gtf_path = $gtf_hg19;
}
elsif($ref eq "chimp"){
        $ref_path_bowtie = $ref_panTro4_bowtie;
        $ref_path_fasta = $ref_panTro4_fasta;
        $gtf_path = $gtf_panTro4;
}
elsif($ref eq "macaque"){
        $ref_path_bowtie = $ref_mmul_1_bowtie;
        $ref_path_fasta = $ref_mmul_1_fasta;
        $gtf_path = $gtf_mmul_1;
}
elsif($ref eq "mouse"){
        $ref_path_bowtie = $ref_mm10_bowtie;
        $ref_path_fasta = $ref_mm10_fasta;
        $gtf_path = $gtf_mm10;
}
elsif($ref eq "sheep"){
        $ref_path_bowtie = $ref_sheep_bowtie;
        $ref_path_fasta = $ref_sheep_fasta;
        $gtf_path = $gtf_sheep;
}



print STDOUT "Parameters in input : \n";
if(!$single){print STDOUT "infiles = $pair1 $pair2\n";}
else{
	print STDOUT "infiles = $single\n";
}
print STDOUT "ref = $ref_path_bowtie\n";
print STDOUT "rg_sample = $rg_sample\n";
print STDOUT "rg_id = $rg_id\n";
print STDOUT "rg_library = $rg_library\n";
print STDOUT "rg_platform = $rg_platform\n";
print STDOUT "rg_center = $rg_center\n";
print STDOUT "max_multihit = $max_multihit\n";
print STDOUT "mate_inner_dist = $mate_inner_dist\n";

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

        $pair1 =~ /([^\/]+)_R1\.fastq\.gz$/;
        $id_prefix = $1;
}


my $bam_prefix = "accepted_hits";

#NEED TO EXTRACT COLUMNS FROM SAMPLE INFO
#PRINT ALL THIS IN A BASH FILE
##TODO MODIFY FOR SINGLE END
my $outfile;

if(!$single){
	if($trim_method ne ""){
		$outfile = "$rg_id.r$mate_inner_dist.g$max_multihit.$ref.trimQ$trimQ";
	}
	else{
		$outfile = "$rg_id.r$mate_inner_dist.g$max_multihit.$ref";
	}
	open(OUT,">$outfile.sh") or die("Cannot write to $outfile.sh");
	print OUT "#/bin/sh\n";
	print OUT "cd $currentDir\n";
}
else{
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
        	$pair1_prefix = $pair1_prefix . "_trimmed";
	        $pair2_prefix = $pair2_prefix . "_trimmed";

        	$pair1="$pair1_prefix.fq.gz";
	        $pair2="$pair2_prefix.fq.gz";
	}
	print OUT "cd $outfile\n";
}



###FIRST STEP : TOPHAT2
print OUT "#TOPHAT2\n";
if(!$single){
	print OUT "time tophat2 -o $outfile -r $mate_inner_dist -p $p -g $max_multihit --rg-id $rg_id --rg-sample $rg_sample --rg-library $rg_library --rg-platform-unit $rg_platform_unit --rg-platform $rg_platform --rg-center $rg_center -G $gtf_path $ref_path_bowtie $pair1 $pair2\n";
}
else{
#TA LIGNE
	print OUT "time tophat2 -o $outfile -p $p -g $max_multihit --rg-id $rg_id --rg-sample $rg_sample --rg-library $rg_library --rg-platform-unit $rg_platform_unit --rg-platform $rg_platform --rg-center $rg_center -G $gtf_path $ref_path_bowtie $single\n";
}


print OUT "cd $outfile\n";
###SECOND STEP : REMOVE PCR DUPLICATES
if(!$single){
	print OUT "#REMOVE PCR DUPLICATES\n";
	print OUT "time java -jar $picard/MarkDuplicates.jar REMOVE_DUPLICATES=true INPUT=$bam_prefix.bam OUTPUT=$outfile.woPCRdup.bam METRICS_FILE=$outfile.picard.metrics.txt\n";
	$bam_prefix = "$outfile.woPCRdup";
	print OUT "samtools index $bam_prefix.bam\n";
}



###THIRD STEP : REORDERING + INDEX
if(!$single){
	print OUT "#REORDERING + INDEX\n";
	print OUT "java -Xmx10g -jar $picard/ReorderSam.jar I=$bam_prefix.bam O=$bam_prefix.PicardReorder.bam REFERENCE=$ref_path_fasta\n";
	$bam_prefix = "$bam_prefix.PicardReorder";
	print OUT "samtools index $bam_prefix.bam\n";
}
else{
	$bam_prefix = $outfile;
        print OUT "#REORDERING + INDEX\n";
        print OUT "java -Xmx10g -jar $picard/ReorderSam.jar I=accepted_hits.bam O=$bam_prefix.PicardReorder.bam REFERENCE=$ref_path_fasta\n";
        $bam_prefix = "$bam_prefix.PicardReorder";
        print OUT "samtools index $bam_prefix.bam\n"
}


###FOURTH STEP : KEEP ONLY PROPERLY PAIRED + INDEX
if(!$single){
	print OUT "#KEEP ONLY PROPERLY PAIRED + INDEX\n";
	print OUT "samtools view -f 0x0002 -b -o $bam_prefix.PP.bam $bam_prefix.bam\n";
	$bam_prefix = "$bam_prefix.PP";
	print OUT "samtools index $bam_prefix.bam\n";
}
else{}

###FIFTH STEP : CUFFLINKS
if($cuff == 1){
	print OUT "#CUFFLINKS\n";
	print OUT "cufflinks -p $p -o $outfile"."_cufflinks -N -G $gtf_path -u $gtf_path $bam_prefix.bam\n";
}

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
