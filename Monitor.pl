#!/usr/bin/perl -w
use strict;
use Cwd qw/abs_path/;
use File::Basename qw/basename dirname/;
use FindBin qw/$RealScript $RealBin/;
use Getopt::Long;
use threads;
use threads::shared;

my ($CONFIG,$OUT_DIR,$HELP);

GetOptions(
	"-c:s"    => \$CONFIG,
	"-o:s"    => \$OUT_DIR,
	"-h|help" => \$HELP
);

if(!($HELP)) {
	if(!($CONFIG && abs_path($CONFIG))) {
		print "Invalid configure file!\n";
	}

	if(!($OUT_DIR)) {
		print "Invalid output path!\n";
	}
}

die "#" x 80 . "\n" . "Usage:
\tperl $RealScript [Options]
Version:
\tV1.0 at 2013-11-11
Options:
\t-c   [s]  Config File of RRBS_Kit pipeline. <required>
\t-o   [s]  Directory which will store all results. <required>
\t-h        Display this help info.
Author:
\tShengjie Gao (gaoshengjie\@genomics.org.cn)
\t    Quan Zhou (zhouquan\@genomics.org.cn)
\t Wenlong Jia (jiawenlong\@genomics.org.cn)
\t     Dan Zou (zoudan\@genomics.cn)	
" . "#" x 80 . "\n" if($HELP || !($CONFIG && &abs_path($CONFIG)) || !($OUT_DIR));

if($OUT_DIR){
	`mkdir -p $OUT_DIR` unless(-d $OUT_DIR);
}

my $lib=abs_path($CONFIG);
my $outdir=abs_path($OUT_DIR);

my ($Bin,$rawref,$refc2t,$refg2a,$adaptor,$rscript,$emycut,$blank,$mpileupregion,$dbsnp);
my ($alignsoft,$javasoft,$samtools,$bcftools,$bissnpsoft,$vcfutils,$picardmergesam,$picardsort,$picardaddread,$insert1,$insert2);
my ($sgesel,$bptsel,$cpunum,$Project,$queue,$dirsample);
my @element;
my %sample;
my %samplese;

unless(-d "$outdir"){
	`mkdir $outdir`;
}
	
open LIB,$lib or die "Unable to open configure file!\n";
while(<LIB>){
	chomp;
	next if(/\#/);
	
######################### Variables
	if(/Region/){
		my @region=split /\=/;
        $region[1]=~s/^\s+//;
        my @ins=split /\s+/,$region[1];
        $insert1=$ins[0];
        $insert2=$ins[1];
    }
	
	if(/SGEsel/){
		my @tmp=split /\=/;
        $tmp[1]=~s/^\s+//;
        $sgesel=$tmp[1];
    }
	
	if(/BPTsel/){
		my @tmp=split /\=/;
        $tmp[1]=~s/^\s+//;
        $bptsel=$tmp[1];
    }
	
	if(/Project/){
		my @tmp=split /\=/;
		$tmp[1]=~s/^\s+//g;
		if($tmp[1] eq "") {
			$Project="";
		}
		else {
			$Project=$tmp[1];
		}
	}

	if(/Queue/){
		my @tmp=split /\=/;
		$tmp[1]=~s/^\s+//g;
		if($tmp[1] eq "") {
			$queue="";
		} 
		else {
			$queue=$tmp[1];	
		}
	}
	
######################### Softwares
	if(/Alignsoft/){
		my @tmp=split /\=/;
		$alignsoft=$tmp[1];
		$alignsoft=~s/^\s+//g;
	}

	if(/Javasoft/){
		my @tmp=split /\=/;
		$javasoft=$tmp[1];
		$javasoft=~s/\s+//g;
	}
	
	if(/Rscript/){
		my @rsp=split /\=/;
		$rscript=$rsp[1];
		$rscript=~s/\s+//g;
	}

######################### Sample
	if(/SplDir/){
		my @tmp=split /\=/;
		$dirsample=$tmp[1];
		$dirsample=~s/\s+//g;
		$dirsample=abs_path($dirsample);
	}
	
	if(/Adapter/){
		my @tmp=split /\=/;
		$adaptor=$tmp[1];
		$adaptor=~s/\s+//g;
		$adaptor="$dirsample/$adaptor";
		$adaptor=abs_path($adaptor);
		die "Error! Invalid adaptor file!\n" unless($adaptor && -e $adaptor);
	}
	
	if(/^Sample/){
		my @sp=split /\=/;
		$sp[1]=~s/^\s+//;
		my @saminfo=split /\s+/,$sp[1];
		$saminfo[5] = "$dirsample/$saminfo[5]";
		$saminfo[5] = abs_path($saminfo[5]);
		
		$saminfo[6] = "$dirsample/$saminfo[6]";
		$saminfo[6] = abs_path($saminfo[6]);
		if($saminfo[2] eq "PE"){
            die "Error! Invalid sample path of $saminfo[5]!\n" unless(-e $saminfo[5]);
            die "Error! Invalid sample path of $saminfo[6]!\n" unless(-e $saminfo[6]);
			$sample{$saminfo[0]}{$saminfo[1]}{$saminfo[3]}{$saminfo[4]}="$saminfo[2]\t$saminfo[5]\t$saminfo[6]\t$saminfo[7]\t$saminfo[8]\t$saminfo[9]";
		}elsif($saminfo[2] eq "SE"){
            die "Error! Invalid sample path of $saminfo[5]!\n" unless(-e $saminfo[5]);
			$sample{$saminfo[0]}{$saminfo[1]}{$saminfo[3]}{$saminfo[4]}="$saminfo[2]\t$saminfo[5]\t$saminfo[6]\t$saminfo[7]\t$saminfo[8]";
		}	
	}
}
close LIB;

$dbsnp ||=" ";
if($sgesel == 0) {
	$cpunum = 2;
}
else {
	$cpunum = 4;
}
###########################################################
# Default settings and setting check
###########################################################
# {{
my @chkTmp;
my $valTmp;
###############################
# Binpath
###############################
$Bin="./bin";
$Bin=abs_path($Bin);
die "Error! Invalid bin path!\n" unless($Bin && -d $Bin);
###############################
# Region
###############################
die "Error! Invalid region configuration!\n" unless($insert1 && $insert2);
my $region_enzyme="$insert1\,$insert2";
###############################
# Softwares
###############################
# Java 
@chkTmp=split(/\s+/, $javasoft);
die "Error! Invalid javasoft configuration!\n" unless($javasoft && -e $javasoft);
# R
@chkTmp=split(/\s+/, $rscript);
die "Error! Invalid rscript configuration!\n" unless($rscript && -e $rscript);
# Align software
@chkTmp = split(/\s+/, $alignsoft);
$alignsoft = shift @chkTmp;
if($alignsoft=~/bowtie/) {
	$alignsoft="$Bin/bowtie2";
}
elsif($alignsoft=~/bsmap/) {
	$alignsoft="$Bin/bsmap";
}elsif($alignsoft=~/bismark/i){
	$alignsoft="$Bin/bismark";
}
my $alignval = join(" ", @chkTmp);
# picard MergeSamFiles
$picardmergesam="$Bin/picard-tools-1.114/MergeSamFiles.jar";
die "Error! Invalid picardmergesam configuration!\n" unless($picardmergesam && -e $picardmergesam);
$picardmergesam="$javasoft -Xmx4g -jar " . "$picardmergesam";
# picard MergeSamFiles
$picardsort="$Bin/picard-tools-1.114/SortSam.jar";
die "Error! Invalid picardsort configuration!\n" unless($picardsort && -e $picardsort);
$picardsort="$javasoft -Xmx4g -jar " . "$picardsort";
# picard SortSam
$picardsort="$Bin/picard-tools-1.114/SortSam.jar";
die "Error! Invalid picardsort configuration!\n" unless($picardsort && -e $picardsort);
$picardsort="$javasoft -Xmx4g -jar " . "$picardsort";
# picard AddOrReplaceReadGroups
$picardaddread="$Bin/picard-tools-1.114/AddOrReplaceReadGroups.jar";
die "Error! Invalid picardaddread configuration!\n" unless($picardaddread && -e $picardaddread);
$picardaddread="$javasoft -Xmx4g -jar " . "$picardaddread";
# bissnpsoft
$bissnpsoft="$Bin/BisSNP-0.82.2.jar";
die "Error! Invalid bissnpsoft configuration!\n" unless($bissnpsoft && -e $bissnpsoft);
$bissnpsoft="$javasoft -Xmx10g -jar " . "$bissnpsoft";
# others
$samtools=$Bin."/samtools";
$bcftools=$Bin."/bcftools";
$vcfutils=$Bin."/vcfutils.pl";
###############################
# Common Data
###############################
# directory
my $dircommon="./data/common";
$dircommon=abs_path($dircommon);
# raw ref
$rawref="$dircommon/hg19.fa";
die "Error! Invalid rawref file!\n" unless($rawref && -e $rawref);
# ref c2t
$refc2t="$dircommon/C2T/c2t.trans.fa";
die "Error! Invalid refc2t file!\n" unless($refc2t && -e $refc2t);
# ref g2a
$refg2a="$dircommon/G2A/g2a.trans.fa";
die "Error! Invalid refg2a file!\n" unless($refg2a && -e $refg2a);
# blank
$blank="$dircommon/chr_len.bed";
die "Error! Invalid blank file!\n" unless($blank && -e $blank);
# emycut
$emycut="$dircommon/hg19.fragment.bed.frag.40-220";
die "Error! Invalid emycut file!\n" unless($emycut && -e $emycut);
# mpileupregion
$mpileupregion="$dircommon/hg19.fragment.bed.mpileup";
die "Error! Invalid mpileupregion file!\n" unless($mpileupregion && -e $mpileupregion);
###############################
# Element Data
###############################
# directory
my $direlement="./data/element";
$direlement=abs_path($direlement);
# elements
@element=glob("$direlement/*");
# }}

my @chr;
for my $i(1..22,'X','Y','M'){
	push @chr,"chr$i";
}
push @chr,'Control';

my %hashchr;
open CHR_E,$blank or die "no blankfile\n";
my @chrname;
while(<CHR_E>){
	chomp;
    next unless(/\w/);
	my @a=split;
    	push @chrname,$a[0];
	$hashchr{$a[0]}=1;	
}
close CHR_E;
my %count; 
 @chrname = grep { ++$count{ $_ } < 2; } @chrname;

#die "@chrname";

open PRE,">$outdir/RRBS_prepare.pl" or die "can't create RRBS_prepare\n";
open DMR,">$outdir/RRBS_dmr.pl" or die "can't create RRBS_dmr.pl\n";
open ASM,">$outdir/RRBS_asm.pl" or die "can't create RRBS_asm.pl\n";
open REPT,">$outdir/RRBS_report.sh" or die "can't create RRBS_report.sh\n";

print PRE "#!/usr/bin/perl\nuse strict;\nuse lib \'$Bin/lib/\'\;\nuse threads;\nuse threads::shared;\nuse Parallel::ForkManager;\nmy \$pm = new Parallel::ForkManager(30);\nmy \@cmd;\n";
print ASM "#!/usr/bin/perl\nuse strict;\nuse lib \'$Bin/lib/\'\;\nuse threads;\nuse threads::shared;\nuse Parallel::ForkManager;\nmy \$pm = new Parallel::ForkManager(30);\nmy \@cmd;\n"; 
print DMR "#!/usr/bin/perl\nuse strict;\nuse lib \'$Bin/lib/\'\;\nuse threads;\nuse threads::shared;\nuse Parallel::ForkManager;\nmy \$pm = new Parallel::ForkManager(30);\nmy \@cmd_1;\nmy \@cmd_2;\n";
print REPT "perl $Bin/report_v3.pl $Bin $outdir report.table\n";

foreach my $case(keys %sample){
	system ("mkdir -p $outdir/Partshell/$case") unless (-d "$outdir/Partshell/$case");
	system ("mkdir -p $outdir/$case/DMR") unless (-d "$outdir/$case/DMR");
	my @type;
	my $normal;
	
	foreach my $type(keys %{$sample{$case}}){
		if($type =~ /Norm/) {
			$normal=$type;
        } 
		else {
            push @type,$type;
        }
	}

	foreach my $type(keys %{$sample{$case}}){
		my $coutregion="$outdir/$case/$type/cout\_$insert1\-$insert2";
		system ("mkdir -p $outdir/$case/$type/qmap");
		system ("mkdir -p $outdir/$case/$type/cout");
		system ("mkdir -p $outdir/$case/$type/cout\_$insert1\-$insert2");
		system ("mkdir -p $outdir/$case/$type/cout\_$insert1\-$insert2/cout_cov");
		print PRE "push \@cmd, \"$outdir/Partshell/$case/$type/rrbs_prepare.sh 2>$outdir/Partshell/$case/$type/rrbs_prepare.log && echo -n \'$case prepare analysis has been done at \'\; date +\'\%Y-\%m-\%d \%H:%M:%S\'\"\;\n";	
		print ASM "push \@cmd, \"$outdir/Partshell/$case/$type/rrbs_asm.sh 2>$outdir/Partshell/$case/$type/rrbs_asm.log && echo -n \'$case asm analysis has been done at \'\; date +\'\%Y-\%m-\%d \%H:%M:%S\'\"\;\n";
		
		
		print DMR "push \@cmd_1, \"$outdir/Partshell/$case/$type/rrbs_dmr_1.sh 2>$outdir/Partshell/$case/$type/rrbs_dmr.log && echo -n \'$case dmr analysis part 1 has been done at \'\; date +\'\%Y-\%m-\%d \%H:%M:%S\'\"\;\n";
		print DMR "push \@cmd_2, \"$outdir/Partshell/$case/$type/rrbs_dmr_2.sh 2>>$outdir/Partshell/$case/$type/rrbs_dmr.log && echo -n \'$case dmr analysis part 2 has been done at \'\; date +\'\%Y-\%m-\%d \%H:%M:%S\'\"\;\n";

		
		system ("mkdir -p $outdir/$case/$type/ASM/");
		system ("mkdir -p $outdir/Partshell/$case/$type");
		open RPRE, ">$outdir/Partshell/$case/$type/rrbs_prepare.sh" or die $!;
		open RASM, ">$outdir/Partshell/$case/$type/rrbs_asm.sh" or die $!;
		open RDMR_1, ">$outdir/Partshell/$case/$type/rrbs_dmr_1.sh" or die $!;
		open RDMR_2, ">$outdir/Partshell/$case/$type/rrbs_dmr_2.sh" or die $!;

		if($alignsoft=~/bowtie/){
			open PRE_AD,">$outdir/Partshell/$case/$type/RRBS-0-dealadptor.sh" or die "can't create $outdir/Partshell/RRBS-0-dealadapter.sh\n";
			open PRE_DFQ,">$outdir/Partshell/$case/$type/RRBS-1-dealfq.sh" or die "can't create $outdir/Partshell/RRBS-1-dealfq.sh\n";
			open PRE_Align, ">$outdir/Partshell/$case/$type/RRBS-2-align.sh" or die "can't create $outdir/Partshell/bowtie2.sh\n";
			open PRE_Sort, ">$outdir/Partshell/$case/$type/RRBS-3-sort.sh" or die "can't create RRBS-3-sort.sh\n";
			open PRE_RTN, ">$outdir/Partshell/$case/$type/RRBS-4-sam_return.sh" or die $!;
			open PRE_MG, ">$outdir/Partshell/$case/$type/RRBS-5-mergerepective.sh" or die $!;
			open PRE_UNQ, ">$outdir/Partshell/$case/$type/RRBS-6-sam_uniq.sh" or die $!;
			open PRE_AFDV, ">$outdir/Partshell/$case/$type/RRBS-7-sam2align.sh" or die $!;
		
			open DMR_Uniq, ">$outdir/Partshell/$case/$type/DMR-8-Uniq.sh" or die $!;
			open DMR_sort, ">$outdir/Partshell/$case/$type/DMR-9-sort.sh" or die $!;
			open DMR_Uniqtotal, ">$outdir/Partshell/$case/$type/DMR-10-Uniqtotal.sh" or die $!;
			open DMR_Stat, ">$outdir/Partshell/$case/$type/DMR-11-Stat.sh" or die $!;
			open DMR_Pos, ">$outdir/Partshell/$case/$type/DMR-12-posord.sh" or die $!;
			open DMR_Cat, ">$outdir/Partshell/$case/$type/DMR-13-cat.sh" or die $!;
			open DMR_sortpos, ">$outdir/Partshell/$case/$type/DMR-14-sortpos.sh" or die $!;
			open DMR_Qmap, ">$outdir/Partshell/$case/$type/DMR-15-qmap.sh" or die $!;
			open DMR_Cout, ">$outdir/Partshell/$case/$type/DMR-16-cout.sh" or die $!;
			open DMR_Getdmr, ">$outdir/Partshell/$case/$type/DMR-17-getdmr.sh" or die $!;
			open DMR_Draw, ">$outdir/Partshell/$case/$type/DMR-18-draw.sh" or die $!;
			open DMR_ANNO,">$outdir/Partshell/$case/$type/DMR-19-anno.sh" or die $!;
			open ASM_MGT,">$outdir/Partshell/$case/$type/ASM-8-mergetotal.sh" or die $!;	
			open ASM_FILM, ">$outdir/Partshell/$case/$type/ASM-9-filtermapread.sh" or die $!;
			open ASM_SORT, ">$outdir/Partshell/$case/$type/ASM-10-sort.sh" or die $!;
			open ASM_MERGE, ">$outdir/Partshell/$case/$type/ASM-11-mergesam.sh" or die $!;
			open ASM_MPIP, ">$outdir/Partshell/$case/$type/ASM-12-mpileup.sh" or die $!;
			open ASM_BCF, ">$outdir/Partshell/$case/$type/ASM-13-bcf.sh" or die $!;
			open ASM_CAT,">$outdir/Partshell/$case/$type/ASM-14-cat.sh" or die $!;
			open ASM_SNP, ">$outdir/Partshell/$case/$type/ASM-15-get_snp.sh" or die $!;
			open ASM_FILSNP, ">$outdir/Partshell/$case/$type/ASM-16-Filtersnp.sh" or die $!;
			open ASM_ASM, ">$outdir/Partshell/$case/$type/ASM-17-ASM.sh" or die $!;
			open ASM_Fil,">$outdir/Partshell/$case/$type/ASM-18-ASM_fil.sh" or die $!;
		
			my @samuniq; 
			my @ctsam; 
			my @gasam;

			foreach my $lib(keys %{$sample{$case}{$type}}){
				foreach my $flow(keys %{$sample{$case}{$type}{$lib}}){
					my @fq=split /\t/,$sample{$case}{$type}{$lib}{$flow};
					system ("mkdir -p $outdir/$case/$type/$lib/$flow");

					if($fq[0]=~/PE/){
						my $adp=$adaptor;
						print PRE_AD "perl $Bin/Dealadaptor.pl $adp $fq[1] $fq[2] $outdir/$case/$type/$lib/$flow/t_clean1.fq.gz $outdir/$case/$type/$lib/$flow/a_clean2.fq.gz  >$outdir/$case/$type/$lib/$flow/adaptor.info \n";
						print PRE_DFQ "perl $Bin/Dealfq.pl $outdir/$case/$type/$lib/$flow/t_clean1.fq.gz $outdir/$case/$type/$lib/$flow/a_clean2.fq.gz -o $outdir/$case/$type/$lib/$flow\n";				
						print PRE_Align "$alignsoft $alignval -I $fq[-2] -X $fq[-1] -x $refc2t -1 $outdir/$case/$type/$lib/$flow/ori.t.gz  -2 $outdir/$case/$type/$lib/$flow/ori.a.gz -S $outdir/$case/$type/$lib/$flow/C2T.sam  \n$alignsoft $alignval -I $fq[-2] -X $fq[-1]  -x $refg2a -1 $outdir/$case/$type/$lib/$flow/ori.t.gz  -2 $outdir/$case/$type/$lib/$flow/ori.a.gz -S $outdir/$case/$type/$lib/$flow/G2A.sam\n";
						print PRE_Sort "sort -s -k 1,1 -T $outdir/$case/$type/$lib/$flow/  $outdir/$case/$type/$lib/$flow/v1.fq -o $outdir/$case/$type/$lib/$flow/v1.fq \n";
						print PRE_Sort "perl $Bin/Sortsam.pl  $outdir/$case/$type/$lib/$flow/\n";
						print PRE_RTN "perl $Bin/Sam_return.pl --mode W $outdir/$case/$type/$lib/$flow/C2T.sam $outdir/$case/$type/$lib/$flow/v1.fq $outdir/$case/$type/$lib/$flow/C2T.pe.sam $outdir/$case/$type/$lib/$flow/C2T.se.sam\nperl $Bin/Sam_return.pl --mode C $outdir/$case/$type/$lib/$flow/G2A.sam $outdir/$case/$type/$lib/$flow/v1.fq $outdir/$case/$type/$lib/$flow/G2A.pe.sam $outdir/$case/$type/$lib/$flow/G2A.se.sam\n";
						print PRE_MG  "$picardmergesam I=$outdir/$case/$type/$lib/$flow/C2T.pe.sam I=$outdir/$case/$type/$lib/$flow/G2A.pe.sam O=$outdir/$case/$type/$lib/$flow/Merge.pe.sam VALIDATION_STRINGENCY=STRICT SO=queryname TMP_DIR=$outdir/$case/$type/$lib/$flow\n$picardmergesam I=$outdir/$case/$type/$lib/$flow/C2T.se.sam I=$outdir/$case/$type/$lib/$flow/G2A.se.sam O=$outdir/$case/$type/$lib/$flow/Merge.se.sam VALIDATION_STRINGENCY=STRICT SO=queryname TMP_DIR=$outdir/$case/$type/$lib/$flow\n";
						print PRE_UNQ "perl $Bin/Sam_uniq.pl -mode se $outdir/$case/$type/$lib/$flow/Merge.se.sam $outdir/$case/$type/$lib/$flow/Uniq.se.sam 2>$outdir/$case/$type/$lib/$flow/someerroinuniqse.sam\nperl $Bin/Sam_uniq.pl -mode pe $outdir/$case/$type/$lib/$flow/Merge.pe.sam $outdir/$case/$type/$lib/$flow/Uniq.pe.sam 2>$outdir/$case/$type/$lib/$flow/someerroinuniqpe.sam\n";
						print PRE_AFDV "perl $Bin/Get_methy.pl --mode pe  $rawref $outdir/$case/$type/$lib/$flow/Uniq.pe.sam $outdir/$case/$type/$lib/$flow/pe.align.gz $outdir/$case/$type/$lib/$flow/Uniq.pe.fil.sam\nperl $Bin/Get_methy.pl --mode se $rawref $outdir/$case/$type/$lib/$flow/Uniq.se.sam $outdir/$case/$type/$lib/$flow/se.align.gz $outdir/$case/$type/$lib/$flow/Uniq.se.fil.sam\n";
				
						print DMR_Uniq "perl $Bin/UniLabel_combine.pl --mode pe $outdir/$case/$type/$lib/$flow/pe.align.gz $outdir/$case/$type/$lib/$flow/pe.uniID 2>$outdir/$case/$type/$lib/$flow/pe.uniID.error\nperl $Bin/UniLabel_combine.pl --mode se $outdir/$case/$type/$lib/$flow/se.align.gz $outdir/$case/$type/$lib/$flow/se.uniID 2>$outdir/$case/$type/$lib/$flow/se.uniID.error\n";
						print DMR_sort "cat $outdir/$case/$type/$lib/$flow/se.uniID  $outdir/$case/$type/$lib/$flow/pe.uniID >$outdir/$case/$type/$lib/$flow/all.uniID \&\& sort -s -k 1,1 -T $outdir/$case/$type/$lib/$flow/ $outdir/$case/$type/$lib/$flow/all.uniID -o $outdir/$case/$type/$lib/$flow/all.uniID\n";
						print DMR_Uniqtotal "perl $Bin/UniLabel_combine.pl  --mode total $outdir/$case/$type/$lib/$flow/all.uniID  $outdir/$case/$type/$lib/$flow/withDupli.labelOrder.gz\n";
						print DMR_Stat "perl $Bin/Statistics.pl $emycut $outdir/$case/$type/$lib/$flow/withDupli.labelOrder.gz  $outdir/$case/$type/$lib/$flow/withDupli.labelOrder.enzyme.gz --out2 $outdir/$case/$type/$lib/$flow/withDupli.labelOrder.notinenzyme.gz --insert $insert1\:$insert2\n";
						print DMR_Pos "perl $Bin/PosOrd_combine.pl $outdir/$case/$type/$lib/$flow/withDupli.labelOrder.enzyme.gz  $outdir/$case/$type/$lib/$flow/posOrder.enzyme\n"; 			
				
						print ASM_MGT "$picardmergesam I=$outdir/$case/$type/$lib/$flow/Uniq.se.fil.sam I=$outdir/$case/$type/$lib/$flow/Uniq.pe.fil.sam O=$outdir/$case/$type/$lib/$flow/Tmp.total.bam VALIDATION_STRINGENCY=STRICT SO=queryname TMP_DIR=$outdir/$case/$type/$lib/$flow\n";
						print ASM_FILM "perl $Bin/Filtermapread.pl $outdir/$case/$type/$lib/$flow/Tmp.total.bam >$outdir/$case/$type/$lib/$flow/Uniq.total.sam\n";
						push @samuniq,"I=$outdir/$case/$type/$lib/$flow/Uniq.total.sam";
						push @ctsam, "I=$outdir/$case/$type/$lib/$flow/C2T.sam";
						push @gasam, "I=$outdir/$case/$type/$lib/$flow/G2A.sam";
					}
					
					if($fq[0]=~/SE/){
						my $adp=$adaptor;
						print PRE_AD "perl $Bin/Dealadaptor.pl $adp $fq[1]  $outdir/$case/$type/$lib/$flow/t_clean1.fq.gz >$outdir/$case/$type/$lib/$flow/adaptor.info -mode se \n";
						print PRE_DFQ "perl $Bin/Dealfq.pl $outdir/$case/$type/$lib/$flow/t_clean1.fq.gz  -o $outdir/$case/$type/$lib/$flow -mode se\n";
						print PRE_Align "$alignsoft $alignval -x $refc2t  $outdir/$case/$type/$lib/$flow/ori.t   -S $outdir/$case/$type/$lib/$flow/C2T.sam  \n$alignsoft $alignval -x $refg2a  $outdir/$case/$type/$lib/$flow/ori.t  -S $outdir/$case/$type/$lib/$flow/G2A.sam\n";
						print PRE_Sort "sort -s -k 1,1 -T $outdir/$case/$type/$lib/$flow/  $outdir/$case/$type/$lib/$flow/v1.fq -o $outdir/$case/$type/$lib/$flow/v1.fq \n";
						print PRE_Sort "perl $Bin/Sortsam.pl  $outdir/$case/$type/$lib/$flow/\n";
						print PRE_RTN "perl $Bin/Sam_return.pl --mode W --type se $outdir/$case/$type/$lib/$flow/C2T.sam $outdir/$case/$type/$lib/$flow/v1.fq $outdir/$case/$type/$lib/$flow/C2T.pe.sam $outdir/$case/$type/$lib/$flow/C2T.se.sam\nperl $Bin/Sam_return.pl --mode C --type se  $outdir/$case/$type/$lib/$flow/G2A.sam $outdir/$case/$type/$lib/$flow/v1.fq $outdir/$case/$type/$lib/$flow/G2A.pe.sam $outdir/$case/$type/$lib/$flow/G2A.se.sam\n";
						print PRE_MG  "$picardmergesam I=$outdir/$case/$type/$lib/$flow/C2T.se.sam I=$outdir/$case/$type/$lib/$flow/G2A.se.sam O=$outdir/$case/$type/$lib/$flow/Merge.se.sam VALIDATION_STRINGENCY=STRICT SO=queryname TMP_DIR=$outdir/$case/$type/$lib/$flow\n";
						print PRE_UNQ "perl $Bin/Sam_uniq.pl -mode se $outdir/$case/$type/$lib/$flow/Merge.se.sam $outdir/$case/$type/$lib/$flow/Uniq.se.sam 2>$outdir/$case/$type/$lib/$flow/someerroinuniqse.sam\n";
						print PRE_AFDV "perl $Bin/Get_methy.pl --mode se $rawref $outdir/$case/$type/$lib/$flow/Uniq.se.sam $outdir/$case/$type/$lib/$flow/se.align.gz $outdir/$case/$type/$lib/$flow/Uniq.se.fil.sam\n";
				
						print DMR_Uniq "perl $Bin/UniLabel_combine.pl --mode se $outdir/$case/$type/$lib/$flow/se.align.gz $outdir/$case/$type/$lib/$flow/se.uniID 2>$outdir/$case/$type/$lib/$flow/se.uniID.error\n";
						print DMR_sort "sort -s -k 1,1 -T $outdir/$case/$type/$lib/$flow/ $outdir/$case/$type/$lib/$flow/se.uniID -o $outdir/$case/$type/$lib/$flow/all.uniID\n";
						print DMR_Uniqtotal "perl $Bin/UniLabel_combine.pl  --mode total $outdir/$case/$type/$lib/$flow/all.uniID  $outdir/$case/$type/$lib/$flow/withDupli.labelOrder.gz\n";
						print DMR_Stat "perl $Bin/Statistics.pl $emycut $outdir/$case/$type/$lib/$flow/withDupli.labelOrder.gz  $outdir/$case/$type/$lib/$flow/withDupli.labelOrder.enzyme.gz --out2 $outdir/$case/$type/$lib/$flow/withDupli.labelOrder.notinenzyme.gz --insert $insert1\:$insert2\n";
						print DMR_Pos "perl $Bin/PosOrd_combine.pl $outdir/$case/$type/$lib/$flow/withDupli.labelOrder.enzyme.gz  $outdir/$case/$type/$lib/$flow/posOrder.enzyme\n";
						print ASM_MGT "$picardsort I=$outdir/$case/$type/$lib/$flow/Uniq.se.fil.sam  O=$outdir/$case/$type/$lib/$flow/Tmp.total.bam VALIDATION_STRINGENCY=STRICT SO=queryname TMP_DIR=$outdir/$case/$type/$lib/$flow\n";
						print ASM_FILM "perl $Bin/Filtermapread.pl $outdir/$case/$type/$lib/$flow/Tmp.total.bam >$outdir/$case/$type/$lib/$flow/Uniq.total.sam\n";
						push @samuniq,"I=$outdir/$case/$type/$lib/$flow/Uniq.total.sam";
						push @ctsam, "I=$outdir/$case/$type/$lib/$flow/C2T.sam";
						push @gasam, "I=$outdir/$case/$type/$lib/$flow/G2A.sam";
					}	
				}
			}

			print  DMR_Cat  "cat $outdir/$case/$type/*/*/posOrder.enzyme >$outdir/$case/$type/qmap/posOrder.enzyme\n";
			print  DMR_sortpos "sort -s -k 1,1 -k 2,2n -k 3,3n -T $outdir/$case/$type/qmap/ $outdir/$case/$type/qmap/posOrder.enzyme -o $outdir/$case/$type/qmap/posOrder.enzyme \&\& perl $Bin/filter_enzyme.pl $blank  $emycut  $outdir/$case/$type/qmap/posOrder.enzyme --insert $insert1\:$insert2\n";
            foreach my $blkchr(@chrname){
                print DMR_Qmap "perl $Bin/Maplist_combine.pl $outdir/$case/$type/qmap/$blkchr\.gz $blank\n";
                
            }
            #print DMR_Qmap "$Bin/Maplist_combine $outdir/$case/$type/qmap/posOrder.enzyme.filter $outdir/$case/$type/qmap/ $blank\n";
		
			my $cout_cov="$coutregion/cout_cov";
			foreach my $name(@chr){
				if(exists($hashchr{$name})){
					print DMR_Cout "perl $Bin/cout_9/Met_combine.pl $outdir/$case/$type/qmap/statistic\.$name\.gz $rawref 0.006 $outdir/$case/$type/cout/$name\.cout\.gz \&\&  perl $Bin/cout_9/cut_cout_cov.pl $emycut $outdir/$case/$type/cout/$name\.cout\.gz $coutregion/$name\.cout.cov 4 $insert1 $insert2 $coutregion/$name\.cout.gz $cout_cov/$name\.cout.gz\n";
					
					if($type=~/Norm/){
						print DMR_Getdmr "echo \"I am normal.\"\n";
						print DMR_ANNO "echo \"I am normal.\"\n";
					} else{
						print DMR_Getdmr "perl $Bin/Advance/dmr_core.pl -model core -i1 $outdir/$case/$type/cout\_$insert1\-$insert2/$name.cout.gz -i2  $outdir/$case/$normal/cout\_$insert1\-$insert2/$name\.cout.gz  -o $outdir/$case/DMR/$type\_normal_$name\.dmr \n";
						foreach my $elementfile(@element){	
							print DMR_ANNO "perl $Bin/anno_dmr.pl $elementfile  $outdir/$case/DMR/$type\_normal_$name\.dmr \n";
						}	
					}
				}
			}
		
			foreach my $elementone(@element){
				print DMR_Draw "perl $Bin/intersection_v2.pl $outdir/$case/$type/cout\_$insert1\-$insert2/ $elementone $emycut $rscript $region_enzyme\n";
			}

			print ASM_SORT "$picardmergesam ".join ("\t",@samuniq)." O=$outdir/$case/$type/ASM/Uniq.total.sort.bam  SO=coordinate TMP_DIR=$outdir/$case/$type/ASM \&\& $samtools index $outdir/$case/$type/ASM/Uniq.total.sort.bam\n";
			print ASM_MERGE "$picardmergesam ".join ("\t",@ctsam)." O=$outdir/$case/$type/ASM/C2T.sort.bam  SO=coordinate TMP_DIR=$outdir/$case/$type/ASM \&\& $samtools index $outdir/$case/$type/ASM/C2T.sort.bam\n";
			print ASM_MERGE "$picardmergesam ".join ("\t",@gasam)." O=$outdir/$case/$type/ASM/G2A.sort.bam  SO=coordinate TMP_DIR=$outdir/$case/$type/ASM \&\& $samtools index $outdir/$case/$type/ASM/G2A.sort.bam\n";
			print ASM_MPIP "$samtools mpileup -I  -gBDS -Q 20 -q 20 -l $mpileupregion -f $refc2t $outdir/$case/$type/ASM/C2T.sort.bam >$outdir/$case/$type/ASM/C2T.emycut.mpileup\n";
			print ASM_MPIP "$samtools mpileup -I  -gBDS -Q 20 -q 20 -l $mpileupregion -f $refg2a $outdir/$case/$type/ASM/G2A.sort.bam >$outdir/$case/$type/ASM/G2A.emycut.mpileup\n";	
			print ASM_BCF "$bcftools view -vcgA $outdir/$case/$type/ASM/C2T.emycut.mpileup \| perl $vcfutils varFilter -Q 15 >$outdir/$case/$type/ASM/C2T.emycut.mpileup.vcf\n";
			print ASM_BCF "$bcftools view -vcgA $outdir/$case/$type/ASM/G2A.emycut.mpileup \| perl $vcfutils varFilter -Q 15 >$outdir/$case/$type/ASM/G2A.emycut.mpileup.vcf\n";
			print ASM_CAT "perl $Bin/ASM_cat.pl $rawref $outdir/$case/$type/ASM/C2T.emycut.mpileup.vcf $outdir/$case/$type/ASM/G2A.emycut.mpileup.vcf $outdir/$case/$type/ASM/Uniq.total.emycut.mpileup.vcf\n";
			print ASM_SNP "perl $Bin/ASM_get_snp.pl $outdir/$case/$type/ASM/Uniq.total.emycut.mpileup.vcf >$outdir/$case/$type/ASM/Uniq.total.emycut.mpileup.vcf.snp\n";
			print ASM_FILSNP "perl $Bin/Filtersnp.pl $outdir/$case/$type/ASM/Uniq.total.emycut.mpileup.vcf.snp >$outdir/$case/$type/ASM/Uniq.total.emycut.mpileup.vcf.snp.fil\n";
			print ASM_ASM "perl  $Bin/ASM_get_CpG.pl $rawref $outdir/$case/$type/ASM/Uniq.total.emycut.mpileup.vcf.snp.fil $outdir/$case/$type/ASM/Uniq.total.sort.bam >$outdir/$case/$type/ASM/ASM.out\n";
			print ASM_Fil "$rscript --no-save --silent --slave <$Bin/ASM_snpfreq.R  --args $outdir/$case/$type/ASM/ASM.out $outdir/$case/$type/ASM/ASM.fil.out \n";

			###########################################################################################
			# rrbs_prepare.pl for each sample
			###########################################################################################
			# {{
			print RPRE "#!/bin/sh\n";
			print RPRE "perl $Bin/rrbs_prepare_script.pl $outdir/Partshell/$case/$type 30 $sgesel $bptsel $cpunum $Project $queue";
			close RPRE;
			# }}
		
			###########################################################################################
			# rrbs_asm.pl for each sample
			###########################################################################################
			# {{
			print RASM "#!/bin/sh\n";
			print RASM "perl $Bin/rrbs_asm_script.pl $outdir/Partshell/$case/$type 30 $sgesel $bptsel $cpunum $Project $queue";
			close RASM;
			# }}
		
			###########################################################################################
			# rrbs_dmr.pl for each sample
			###########################################################################################
			# {{
			print RDMR_1 "#!/bin/sh\n";
			print RDMR_1 "perl $Bin/rrbs_dmr_script_1.pl $outdir/Partshell/$case/ $type $normal 30 $sgesel $bptsel $cpunum $Project $queue";
			close RDMR_1;
			print RDMR_2 "#!/bin/sh\n";
			print RDMR_2 "perl $Bin/rrbs_dmr_script_2.pl $outdir/Partshell/$case/ $type $normal 30 $sgesel $bptsel $cpunum $Project $queue";
			close RDMR_2;
			# }}
		
			close ASM_MGT;
			close ASM_FILM;
			close ASM_SORT;
			close ASM_MERGE;
			close ASM_MPIP;
			close ASM_BCF;
			close ASM_SNP;
			close ASM_FILSNP;
			close ASM_ASM;
			close ASM_Fil;

			close PRE_AD;
			close PRE_DFQ;
			close PRE_Align;
			close PRE_RTN;
			close PRE_MG;
			close PRE_UNQ;
			close PRE_AFDV;

			close DMR_Uniq;
			close DMR_sort;
			close DMR_Uniqtotal;
			close DMR_Stat;
			close DMR_Pos;
			close DMR_Cat;
			close DMR_sortpos;
			close DMR_Qmap;
			close DMR_Cout;
			#close DMR_Errorate;
			close DMR_Getdmr;
			close DMR_Draw;
			close DMR_ANNO;

		}
		elsif($alignsoft=~/bsmap/){##for bsmap
			open PRE_AD,">$outdir/Partshell/$case/$type/RRBS-1-dealadaptor.sh" or die "can't create $outdir/Partshell/RRBS-1-dealadapter.sh\n";
			open PRE_BOWT, ">$outdir/Partshell/$case/$type/RRBS-2-bsmap.sh" or die "can't create $outdir/Partshell/RRBS-2-bsmap.sh\n";
            open PRE_Addread, ">$outdir/Partshell/$case/$type/RRBS-3-addread.sh" or die $!;
			open DMR_Dis, ">$outdir/Partshell/$case/$type/DMR-4-dis_pese.sh" or die $!;
			open DMR_Getmethy, ">$outdir/Partshell/$case/$type/DMR-5-Get_methy.sh" or die $!;
			open DMR_Sortpe, ">$outdir/Partshell/$case/$type/DMR-6-sortpealign.sh" or die $!;
			open DMR_Uniq, ">$outdir/Partshell/$case/$type/DMR-7-uniID.sh" or die $!;
			open DMR_Catu, ">$outdir/Partshell/$case/$type/DMR-8-CatuniqID.sh" or die $!;
			open DMR_Uniqtotal, ">$outdir/Partshell/$case/$type/DMR-9-UniLabel_combine.sh" or die $!;
			open DMR_Stat, ">$outdir/Partshell/$case/$type/DMR-10-Statistics.sh" or die $!;
			open DMR_Pos, ">$outdir/Partshell/$case/$type/DMR-11-PosOrd_combine.sh" or die $!;
			open DMR_MGT, ">$outdir/Partshell/$case/$type/DMR-12-cat.sh" or die $!;
			open DMR_sortpos, ">$outdir/Partshell/$case/$type/DMR-13-sort.sh" or die $!;
			open DMR_Qmap, ">$outdir/Partshell/$case/$type/DMR-14-qmap.sh" or die $!;
			open DMR_Cout, ">$outdir/Partshell/$case/$type/DMR-15-cout.sh" or die $!;
			open DMR_Enzyme, ">$outdir/Partshell/$case/$type/DMR-16-cout_40-220.sh" or die $!;
			open DMR_Getdmr, ">$outdir/Partshell/$case/$type/DMR-17-getdmr.sh" or die $!;
			open DMR_Draw, ">$outdir/Partshell/$case/$type/DMR-18-draw.sh" or die $!;
			open DMR_ANNO,">$outdir/Partshell/$case/$type/DMR-19-anno.sh" or die $!;
			
			open ASM_MGT,">$outdir/Partshell/$case/$type/ASM-4-mergetotal.sh" or die $!;
			open ASM_SNP, ">$outdir/Partshell/$case/$type/ASM-5-get_snp.sh" or die $!;
			open ASM_FILSNP , ">$outdir/Partshell/$case/$type/ASM-6-Filtersnp.sh" or die $!;
			open ASM_FILHET, ">$outdir/Partshell/$case/$type/ASM-7-Filterhet.sh" or die $!;
			open ASM_ASM, ">$outdir/Partshell/$case/$type/ASM-8-ASM.sh" or die $!;
			open ASM_Fil,">$outdir/Partshell/$case/$type/ASM-9-ASM_fil.sh" or die $!;
			
			my @samuniq;
			foreach my $lib(keys %{$sample{$case}{$type}}){
				foreach my $flow(keys %{$sample{$case}{$type}{$lib}}){
					my @fq=split /\t/,$sample{$case}{$type}{$lib}{$flow};
					system ("mkdir -p $outdir/$case/$type/$lib/$flow");
					
					if($fq[0]=~/PE/){
						my $adp=$adaptor;
						print PRE_AD "perl $Bin/Dealadaptor.pl $adp $fq[1] $fq[2] $outdir/$case/$type/$lib/$flow/t_clean1.fq.gz $outdir/$case/$type/$lib/$flow/a_clean2.fq.gz  >$outdir/$case/$type/$lib/$flow/adaptor.info \n";
						print PRE_BOWT "$alignsoft $alignval -d $rawref  -m $fq[-2] -x $fq[-1]  -a $outdir/$case/$type/$lib/$flow/t_clean1.fq.gz  -b $outdir/$case/$type/$lib/$flow/a_clean2.fq.gz -o $outdir/$case/$type/$lib/$flow/outfile.sam\n";
						print PRE_Addread "$picardaddread I=$outdir/$case/$type/$lib/$flow/outfile.sam O=$outdir/$case/$type/$lib/$flow/outfile.bam  LB=$lib PU=$flow SM=$case\-$type PL=$fq[-3]\n";
						
						print DMR_Dis "perl $Bin/bsmap_dist_pese.pl $outdir/$case/$type/$lib/$flow/outfile.bam $outdir/$case/$type/$lib/$flow/outfile.pe.sam $outdir/$case/$type/$lib/$flow/outfile.se.sam\n";
						print DMR_Getmethy "perl $Bin/Get_methy_bsmap.pl --mode pe $rawref $outdir/$case/$type/$lib/$flow/outfile.pe.sam $outdir/$case/$type/$lib/$flow/pe.align\nperl $Bin/Get_methy_bsmap.pl --mode pe $rawref $outdir/$case/$type/$lib/$flow/outfile.se.sam $outdir/$case/$type/$lib/$flow/se.align\n";
						print DMR_Sortpe "sort -s -k 1,1 -T $outdir/$case/$type/$lib/$flow $outdir/$case/$type/$lib/$flow/pe.align -o $outdir/$case/$type/$lib/$flow/pe.align\n";
						print DMR_Uniq "perl $Bin/Bsmap_UniLabel_combine.pl --mode pe $outdir/$case/$type/$lib/$flow/pe.align $outdir/$case/$type/$lib/$flow/pe.uniID 2>$outdir/$case/$type/$lib/$flow/pe.uniID.error\nperl $Bin/Bsmap_UniLabel_combine.pl --mode se $outdir/$case/$type/$lib/$flow/se.align $outdir/$case/$type/$lib/$flow/se.uniID 2>$outdir/$case/$type/$lib/$flow/se.uniID.error\n";
						print DMR_Catu "cat $outdir/$case/$type/$lib/$flow/pe.uniID $outdir/$case/$type/$lib/$flow/se.uniID >$outdir/$case/$type/$lib/$flow/all.uniID\n";
						print DMR_Uniqtotal "perl $Bin/Bsmap_UniLabel_combine.pl  --mode total $outdir/$case/$type/$lib/$flow/all.uniID  $outdir/$case/$type/$lib/$flow/withDupli.labelOrder.gz\n";
						print DMR_Stat "perl $Bin/Statistics.pl $emycut $outdir/$case/$type/$lib/$flow/withDupli.labelOrder.gz  $outdir/$case/$type/$lib/$flow/withDupli.labelOrder.enzyme.gz --out2 $outdir/$case/$type/$lib/$flow/withDupli.labelOrder.notinenzyme.gz --insert $insert1\:$insert2\n";
						
						print DMR_Pos "perl $Bin/PosOrd_combine.pl $outdir/$case/$type/$lib/$flow/withDupli.labelOrder.enzyme.gz  $outdir/$case/$type/$lib/$flow/posOrder.enzyme\n";
						push @samuniq,"I=$outdir/$case/$type/$lib/$flow/outfile.bam";
					}
					
					if($fq[0]=~/SE/){
						my $adp=$adaptor;
						print PRE_AD "perl $Bin/Dealadaptor.pl $adp $fq[1]  $outdir/$case/$type/$lib/$flow/t_clean1.fq.gz >$outdir/$case/$type/$lib/$flow/adaptor.info -mode se \n";
						print PRE_BOWT "$alignsoft $alignval -d $rawref  -a $outdir/$case/$type/$lib/$flow/t_clean1.fq.gz  -o $outdir/$case/$type/$lib/$flow/outfile.sam\n";
						print PRE_Addread "$picardaddread I=$outdir/$case/$type/$lib/$flow/outfile.sam O=$outdir/$case/$type/$lib/$flow/outfile.bam  LB=$lib PU=$flow SM=$case\-$type PL=$fq[-3]\n";
						print DMR_Getmethy "perl $Bin/Get_methy_bsmap.pl --mode se $rawref $outdir/$case/$type/$lib/$flow/outfile.bam $outdir/$case/$type/$lib/$flow/se.align.gz\n";
						
						print DMR_Uniq "perl $Bin/UniLabel_combine.pl --mode se $outdir/$case/$type/$lib/$flow/se.align.gz $outdir/$case/$type/$lib/$flow/se.uniID 2>$outdir/$case/$type/$lib/$flow/se.uniID.error\n";
						print DMR_Uniqtotal "perl $Bin/UniLabel_combine.pl  --mode total $outdir/$case/$type/$lib/$flow/se.uniID  $outdir/$case/$type/$lib/$flow/withDupli.labelOrder.gz\n";
						print DMR_Stat "perl $Bin/Statistics.pl $emycut $outdir/$case/$type/$lib/$flow/withDupli.labelOrder.gz  $outdir/$case/$type/$lib/$flow/withDupli.labelOrder.enzyme.gz --out2 $outdir/$case/$type/$lib/$flow/withDupli.labelOrder.notinenzyme.gz --insert $insert1\:$insert2\n";
						print DMR_Pos "perl $Bin/PosOrd_combine.pl $outdir/$case/$type/$lib/$flow/withDupli.labelOrder.enzyme.gz  $outdir/$case/$type/$lib/$flow/posOrder.enzyme\n";
						
						push @samuniq,"I=$outdir/$case/$type/$lib/$flow/outfile.bam";
					   
					}
				}
			}
			
			print  DMR_MGT  "cat $outdir/$case/$type/*/*/posOrder.enzyme >$outdir/$case/$type/qmap/posOrder.enzyme\n";
			print  DMR_sortpos "sort -s -k 1,1 -k 2,2n -k 3,3n -T $outdir/$case/$type/qmap/ $outdir/$case/$type/qmap/posOrder.enzyme -o $outdir/$case/$type/qmap/posOrder.enzyme \&\& perl $Bin/filter_enzyme.pl $blank $emycut $outdir/$case/$type/qmap/posOrder.enzyme --insert $insert1\:$insert2\n";
            foreach my $blkchr(@chrname){
                print DMR_Qmap "perl $Bin/Maplist_combine.pl $outdir/$case/$type/qmap/$blkchr\.gz $blank\n";
                
            }
            #print DMR_Qmap  "$Bin/Maplist_combine $outdir/$case/$type/qmap/posOrder.enzyme.filter $outdir/$case/$type/qmap/ $blank\n";
			
			my $cout_cov="$coutregion/cout_cov";
			foreach my $name(@chr){
				if(exists($hashchr{$name})){
					print DMR_Cout "perl $Bin/cout_9/Met_combine.pl $outdir/$case/$type/qmap/statistic\.$name\.gz $rawref 0.006 $outdir/$case/$type/cout/$name\.cout\.gz\n";
					print DMR_Enzyme "perl $Bin/cout_9/cut_cout_cov.pl $emycut $outdir/$case/$type/cout/$name\.cout\.gz $coutregion/$name\.cout.cov 4 $insert1 $insert2 $coutregion/$name\.cout.gz $cout_cov/$name\.cout.gz\n";
					
					if($type=~/Normal/){
						print DMR_Getdmr "echo \"I am normal.\"\n";
						print DMR_ANNO "echo \"I am normal.\"\n";
					} else{
						print DMR_Getdmr "perl $Bin/dmr_core.pl -i1 $outdir/$case/$normal/cout\_$insert1\-$insert2/$name\.cout.gz -i2 $outdir/$case/$type/cout\_$insert1\-$insert2/$name.cout.gz -o $outdir/$case/DMR/$type\_normal_$name\.dmr -model core\n";
						
						foreach my $elementfile(@element){
							print DMR_ANNO "perl $Bin/anno_dmr.pl $elementfile  $outdir/$case/DMR/$type\_normal_$name\.dmr \n";
						}
					}
				}
			}
			
			foreach my $elementone(@element){
				print DMR_Draw "perl $Bin/intersection_v2.pl $outdir/$case/$type/cout\_$insert1\-$insert2/ $elementone $emycut $rscript $insert1\,$insert2\n";
			}
			
			print ASM_MGT "$picardmergesam ".join ("\t",@samuniq)." O=$outdir/$case/$type/ASM/BSMAP.sort.bam  SO=coordinate TMP_DIR=$outdir/$case/$type/ASM \&\& $samtools index $outdir/$case/$type/ASM/BSMAP.sort.bam\n";
			print ASM_SNP "$bissnpsoft -T BisulfiteGenotyper -R $rawref -I $outdir/$case/$type/ASM/BSMAP.sort.bam -L $blank $dbsnp -vfn1 $outdir/$case/$type/ASM/BSMAP.CpGoutRG.vcf -vfn2 $outdir/$case/$type/ASM/BSMAP.SNPoutRG.vcf\n";
			print ASM_FILSNP "$bissnpsoft -T VCFpostprocess -R $rawref -oldVcf $outdir/$case/$type/ASM/BSMAP.SNPoutRG.vcf -newVcf $outdir/$case/$type/ASM/BSMAP.SNPoutRG.vcf.filter -snpVcf $outdir/$case/$type/ASM/BSMAP.SNPoutRG.vcf -o  $outdir/$case/$type/ASM/BSMAP.CpGoutRG.vcf.filter.summary.txt\n";
			print ASM_FILHET "perl $Bin/ASM_BSMAP_hetsnp.pl $outdir/$case/$type/ASM/BSMAP.SNPoutRG.vcf.filter >$outdir/$case/$type/ASM/BSMAP.SNPoutRG.vcf.filter.het\n";
			print ASM_ASM "perl  $Bin/ASM_get_CpG.pl $rawref $outdir/$case/$type/ASM/BSMAP.SNPoutRG.vcf.filter.het $outdir/$case/$type/ASM/BSMAP.sort.bam >$outdir/$case/$type/ASM/ASM.out\n";
			print ASM_Fil "$rscript --no-save --silent --slave <$Bin/ASM_snpfreq.R  --args $outdir/$case/$type/ASM/ASM.out $outdir/$case/$type/ASM/ASM.fil.out \n";
			
			
			
			
			###########################################################################################
			# rrbs_prepare.pl for each sample
			###########################################################################################
			# {{
			print RPRE "#!/bin/sh\n";
			print RPRE "perl $Bin/rrbs_prepare_script_BSMAP.pl $outdir/Partshell/$case/$type 30 $sgesel $bptsel $cpunum $Project $queue";
			close RPRE;
			# }}
		
			###########################################################################################
			# rrbs_asm.pl for each sample
			###########################################################################################
			# {{
			print RASM "#!/bin/sh\n";
			print RASM "perl $Bin/rrbs_asm_script_BSMAP.pl $outdir/Partshell/$case/$type 30 $sgesel $bptsel $cpunum $Project $queue";
			close RASM;
			# }}
		
			###########################################################################################
			# rrbs_dmr.pl for each sample
			###########################################################################################
			# {{
			print RDMR_1 "#!/bin/sh\n";
			print RDMR_1 "perl $Bin/rrbs_dmr_script_BSMAP_1.pl $outdir/Partshell/$case/ $type $normal 30 $sgesel $bptsel $cpunum $Project $queue";
			close RDMR_1;
			print RDMR_2 "#!/bin/sh\n";
			print RDMR_2 "perl $Bin/rrbs_dmr_script_BSMAP_2.pl $outdir/Partshell/$case/ $type $normal 30 $sgesel $bptsel $cpunum $Project $queue";
			close RDMR_2;
			# }}
			
			close ASM_MGT;
			close ASM_SNP;
			close ASM_FILSNP;
			close ASM_ASM;
			close ASM_FILHET;
			close ASM_Fil;
			
			close PRE_AD;
			close PRE_BOWT;
			close PRE_Addread;
			
			close DMR_Dis;
			close DMR_Getmethy;
			close DMR_Sortpe;
			close DMR_Uniq;
			close DMR_Uniqtotal;
			close DMR_Cat;
            close DMR_Catu;
			close DMR_Stat;
			close DMR_Pos;
			close DMR_Qmap;
			close DMR_MGT;
			close DMR_sortpos;
			close DMR_Cout;
			close DMR_Enzyme;
			close DMR_Getdmr;
			close DMR_Draw;
			close DMR_ANNO;    
		} # End of BSMAP
		elsif($alignsoft=~/bismark/){ ###for bismark
			open PRE_AD,">$outdir/Partshell/$case/$type/RRBS-1-dealadaptor.sh" or die "can't create $outdir/Partshell/RRBS-1-dealadapter.sh\n";
			open PRE_BOWT, ">$outdir/Partshell/$case/$type/RRBS-2-bsmap.sh" or die "can't create $outdir/Partshell/RRBS-2-bsmap.sh\n";
            open PRE_Addread, ">$outdir/Partshell/$case/$type/RRBS-3-addread.sh" or die $!;
			open DMR_Dis, ">$outdir/Partshell/$case/$type/DMR-4-dis_pese.sh" or die $!;
			open DMR_Getmethy, ">$outdir/Partshell/$case/$type/DMR-5-Get_methy.sh" or die $!;
			open DMR_Sortpe, ">$outdir/Partshell/$case/$type/DMR-6-sortpealign.sh" or die $!;
			open DMR_Uniq, ">$outdir/Partshell/$case/$type/DMR-7-uniID.sh" or die $!;
			open DMR_Catu, ">$outdir/Partshell/$case/$type/DMR-8-CatuniqID.sh" or die $!;
			open DMR_Uniqtotal, ">$outdir/Partshell/$case/$type/DMR-9-UniLabel_combine.sh" or die $!;
			open DMR_Stat, ">$outdir/Partshell/$case/$type/DMR-10-Statistics.sh" or die $!;
			open DMR_Pos, ">$outdir/Partshell/$case/$type/DMR-11-PosOrd_combine.sh" or die $!;
			open DMR_MGT, ">$outdir/Partshell/$case/$type/DMR-12-cat.sh" or die $!;
			open DMR_sortpos, ">$outdir/Partshell/$case/$type/DMR-13-sort.sh" or die $!;
			open DMR_Qmap, ">$outdir/Partshell/$case/$type/DMR-14-qmap.sh" or die $!;
			open DMR_Cout, ">$outdir/Partshell/$case/$type/DMR-15-cout.sh" or die $!;
			open DMR_Enzyme, ">$outdir/Partshell/$case/$type/DMR-16-cout_40-220.sh" or die $!;
			open DMR_Getdmr, ">$outdir/Partshell/$case/$type/DMR-17-getdmr.sh" or die $!;
			open DMR_Draw, ">$outdir/Partshell/$case/$type/DMR-18-draw.sh" or die $!;
			open DMR_ANNO,">$outdir/Partshell/$case/$type/DMR-19-anno.sh" or die $!;
			
			open ASM_MGT,">$outdir/Partshell/$case/$type/ASM-4-mergetotal.sh" or die $!;
			open ASM_SNP, ">$outdir/Partshell/$case/$type/ASM-5-get_snp.sh" or die $!;
			open ASM_FILSNP , ">$outdir/Partshell/$case/$type/ASM-6-Filtersnp.sh" or die $!;
			open ASM_FILHET, ">$outdir/Partshell/$case/$type/ASM-7-Filterhet.sh" or die $!;
			open ASM_ASM, ">$outdir/Partshell/$case/$type/ASM-8-ASM.sh" or die $!;
			open ASM_Fil,">$outdir/Partshell/$case/$type/ASM-9-ASM_fil.sh" or die $!;
			
			my @samuniq;
			my $refdir=dirname $rawref;
			foreach my $lib(keys %{$sample{$case}{$type}}){
				foreach my $flow(keys %{$sample{$case}{$type}{$lib}}){
					my @fq=split /\t/,$sample{$case}{$type}{$lib}{$flow};
					system ("mkdir -p $outdir/$case/$type/$lib/$flow");
					
					if($fq[0]=~/PE/){
						my $adp=$adaptor;
						print PRE_AD "perl $Bin/Dealadaptor.pl $adp $fq[1] $fq[2] $outdir/$case/$type/$lib/$flow/t_clean1.fq.gz $outdir/$case/$type/$lib/$flow/a_clean2.fq.gz  >$outdir/$case/$type/$lib/$flow/adaptor.info \n";
						print PRE_BOWT "$alignsoft $alignval -I $fq[-2] -X $fq[-1] --path_to_bowtie $Bin $refdir --temp_dir $outdir/$case/$type/$lib/$flow/  -o $outdir/$case/$type/$lib/$flow/ -1 $outdir/$case/$type/$lib/$flow/t_clean1.fq.gz  -2 $outdir/$case/$type/$lib/$flow/a_clean2.fq.gz \n";
						print PRE_Addread "$picardaddread I=$outdir/$case/$type/$lib/$flow/t_clean1.fq.gz_bismark_bt2_pe.sam O=$outdir/$case/$type/$lib/$flow/outfile.bam  LB=$lib PU=$flow SM=$case\-$type PL=$fq[-3]\n";
						
						print DMR_Dis "perl $Bin/bsmap_dist_pese.pl $outdir/$case/$type/$lib/$flow/outfile.bam $outdir/$case/$type/$lib/$flow/outfile.pe.sam $outdir/$case/$type/$lib/$flow/outfile.se.sam\n";
						print DMR_Getmethy "perl $Bin/Get_methy_bsmap.pl --mode pe $rawref $outdir/$case/$type/$lib/$flow/outfile.pe.sam $outdir/$case/$type/$lib/$flow/pe.align\nperl $Bin/Get_methy_bsmap.pl --mode pe $rawref $outdir/$case/$type/$lib/$flow/outfile.se.sam $outdir/$case/$type/$lib/$flow/se.align\n";
						print DMR_Sortpe "sort -s -k 1,1 -T $outdir/$case/$type/$lib/$flow $outdir/$case/$type/$lib/$flow/pe.align -o $outdir/$case/$type/$lib/$flow/pe.align\n";
						print DMR_Uniq "perl $Bin/Bsmap_UniLabel_combine.pl --mode pe $outdir/$case/$type/$lib/$flow/pe.align $outdir/$case/$type/$lib/$flow/pe.uniID 2>$outdir/$case/$type/$lib/$flow/pe.uniID.error\nperl $Bin/Bsmap_UniLabel_combine.pl --mode se $outdir/$case/$type/$lib/$flow/se.align $outdir/$case/$type/$lib/$flow/se.uniID 2>$outdir/$case/$type/$lib/$flow/se.uniID.error\n";
						print DMR_Catu "cat $outdir/$case/$type/$lib/$flow/pe.uniID $outdir/$case/$type/$lib/$flow/se.uniID >$outdir/$case/$type/$lib/$flow/all.uniID\n";
						print DMR_Uniqtotal "perl $Bin/Bsmap_UniLabel_combine.pl  --mode total $outdir/$case/$type/$lib/$flow/all.uniID  $outdir/$case/$type/$lib/$flow/withDupli.labelOrder.gz\n";
						print DMR_Stat "perl $Bin/Statistics.pl $emycut $outdir/$case/$type/$lib/$flow/withDupli.labelOrder.gz  $outdir/$case/$type/$lib/$flow/withDupli.labelOrder.enzyme.gz --out2 $outdir/$case/$type/$lib/$flow/withDupli.labelOrder.notinenzyme.gz --insert $insert1\:$insert2\n";
						
						print DMR_Pos "perl $Bin/PosOrd_combine.pl $outdir/$case/$type/$lib/$flow/withDupli.labelOrder.enzyme.gz  $outdir/$case/$type/$lib/$flow/posOrder.enzyme\n";
						push @samuniq,"I=$outdir/$case/$type/$lib/$flow/outfile.bam";
					}
					
					if($fq[0]=~/SE/){
						my $adp=$adaptor;
						print PRE_AD "perl $Bin/Dealadaptor.pl $adp $fq[1]  $outdir/$case/$type/$lib/$flow/t_clean1.fq.gz >$outdir/$case/$type/$lib/$flow/adaptor.info -mode se \n";
						print PRE_BOWT "$alignsoft $alignval --path_to_bowtie $Bin $refdir --temp_dir outdir/$case/$type/$lib/$flow/   -o $outdir/$case/$type/$lib/$flow/ $outdir/$case/$type/$lib/$flow/t_clean1.fq.gz \n";
						print PRE_Addread "$picardaddread I=$outdir/$case/$type/$lib/$flow/t_clean1.fq.gz_bismark_bt2.sam  O=$outdir/$case/$type/$lib/$flow/outfile.bam  LB=$lib PU=$flow SM=$case\-$type PL=$fq[-3]\n";
						print DMR_Getmethy "perl $Bin/Get_methy_bsmap.pl --mode se $rawref $outdir/$case/$type/$lib/$flow/outfile.bam $outdir/$case/$type/$lib/$flow/se.align.gz\n";
						
						print DMR_Uniq "perl $Bin/UniLabel_combine.pl --mode se $outdir/$case/$type/$lib/$flow/se.align.gz $outdir/$case/$type/$lib/$flow/se.uniID 2>$outdir/$case/$type/$lib/$flow/se.uniID.error\n";
						print DMR_Uniqtotal "perl $Bin/UniLabel_combine.pl  --mode total $outdir/$case/$type/$lib/$flow/se.uniID  $outdir/$case/$type/$lib/$flow/withDupli.labelOrder.gz\n";
						print DMR_Stat "perl $Bin/Statistics.pl $emycut $outdir/$case/$type/$lib/$flow/withDupli.labelOrder.gz  $outdir/$case/$type/$lib/$flow/withDupli.labelOrder.enzyme.gz --out2 $outdir/$case/$type/$lib/$flow/withDupli.labelOrder.notinenzyme.gz --insert $insert1\:$insert2\n";
						print DMR_Pos "perl $Bin/PosOrd_combine.pl $outdir/$case/$type/$lib/$flow/withDupli.labelOrder.enzyme.gz  $outdir/$case/$type/$lib/$flow/posOrder.enzyme\n";
						
						push @samuniq,"I=$outdir/$case/$type/$lib/$flow/outfile.bam";
					   
					}
				}
			}
			
			print  DMR_MGT  "cat $outdir/$case/$type/*/*/posOrder.enzyme >$outdir/$case/$type/qmap/posOrder.enzyme\n";
			print  DMR_sortpos "sort -s -k 1,1 -k 2,2n -k 3,3n -T $outdir/$case/$type/qmap/ $outdir/$case/$type/qmap/posOrder.enzyme -o $outdir/$case/$type/qmap/posOrder.enzyme \&\& perl $Bin/filter_enzyme.pl $blank $emycut $outdir/$case/$type/qmap/posOrder.enzyme --insert $insert1\:$insert2\n";
            foreach my $blkchr(@chrname){
                print DMR_Qmap "perl $Bin/Maplist_combine.pl $outdir/$case/$type/qmap/$blkchr\.gz $blank\n";
                
            }
            #print DMR_Qmap  "$Bin/Maplist_combine $outdir/$case/$type/qmap/posOrder.enzyme.filter $outdir/$case/$type/qmap/ $blank\n";
			
			my $cout_cov="$coutregion/cout_cov";
			foreach my $name(@chr){
				if(exists($hashchr{$name})){
					print DMR_Cout "perl $Bin/cout_9/Met_combine.pl $outdir/$case/$type/qmap/statistic\.$name\.gz $rawref 0.006 $outdir/$case/$type/cout/$name\.cout\.gz\n";
					print DMR_Enzyme "perl $Bin/cout_9/cut_cout_cov.pl $emycut $outdir/$case/$type/cout/$name\.cout\.gz $coutregion/$name\.cout.cov 4 $insert1 $insert2 $coutregion/$name\.cout.gz $cout_cov/$name\.cout.gz\n";
					
					if($type=~/Normal/){
						print DMR_Getdmr "echo \"I am normal.\"\n";
						print DMR_ANNO "echo \"I am normal.\"\n";
					} else{
						print DMR_Getdmr "perl $Bin/dmr_core.pl -i1 $outdir/$case/$normal/cout\_$insert1\-$insert2/$name\.cout.gz -i2 $outdir/$case/$type/cout\_$insert1\-$insert2/$name.cout.gz -o $outdir/$case/DMR/$type\_normal_$name\.dmr -model core\n";
						
						foreach my $elementfile(@element){
							print DMR_ANNO "perl $Bin/anno_dmr.pl $elementfile  $outdir/$case/DMR/$type\_normal_$name\.dmr \n";
						}
					}
				}
			}
			
			foreach my $elementone(@element){
				print DMR_Draw "perl $Bin/intersection_v2.pl $outdir/$case/$type/cout\_$insert1\-$insert2/ $elementone $emycut $rscript $insert1\,$insert2\n";
			}
			
			print ASM_MGT "$picardmergesam ".join ("\t",@samuniq)." O=$outdir/$case/$type/ASM/BSMAP.sort.bam  SO=coordinate TMP_DIR=$outdir/$case/$type/ASM \&\& $samtools index $outdir/$case/$type/ASM/BSMAP.sort.bam\n";
			print ASM_SNP "$bissnpsoft -T BisulfiteGenotyper -R $rawref -I $outdir/$case/$type/ASM/BSMAP.sort.bam -L $blank $dbsnp -vfn1 $outdir/$case/$type/ASM/BSMAP.CpGoutRG.vcf -vfn2 $outdir/$case/$type/ASM/BSMAP.SNPoutRG.vcf\n";
			print ASM_FILSNP "$bissnpsoft -T VCFpostprocess -R $rawref -oldVcf $outdir/$case/$type/ASM/BSMAP.SNPoutRG.vcf -newVcf $outdir/$case/$type/ASM/BSMAP.SNPoutRG.vcf.filter -snpVcf $outdir/$case/$type/ASM/BSMAP.SNPoutRG.vcf -o  $outdir/$case/$type/ASM/BSMAP.CpGoutRG.vcf.filter.summary.txt\n";
			print ASM_FILHET "perl $Bin/ASM_BSMAP_hetsnp.pl $outdir/$case/$type/ASM/BSMAP.SNPoutRG.vcf.filter >$outdir/$case/$type/ASM/BSMAP.SNPoutRG.vcf.filter.het\n";
			print ASM_ASM "perl  $Bin/ASM_get_CpG.pl $rawref $outdir/$case/$type/ASM/BSMAP.SNPoutRG.vcf.filter.het $outdir/$case/$type/ASM/BSMAP.sort.bam >$outdir/$case/$type/ASM/ASM.out\n";
			print ASM_Fil "$rscript --no-save --silent --slave <$Bin/ASM_snpfreq.R  --args $outdir/$case/$type/ASM/ASM.out $outdir/$case/$type/ASM/ASM.fil.out \n";
			
			
			
			
			###########################################################################################
			# rrbs_prepare.pl for each sample
			###########################################################################################
			# {{
			print RPRE "#!/bin/sh\n";
			print RPRE "perl $Bin/rrbs_prepare_script_BSMAP.pl $outdir/Partshell/$case/$type 30 $sgesel $bptsel $cpunum $Project $queue";
			close RPRE;
			# }}
		
			###########################################################################################
			# rrbs_asm.pl for each sample
			###########################################################################################
			# {{
			print RASM "#!/bin/sh\n";
			print RASM "perl $Bin/rrbs_asm_script_BSMAP.pl $outdir/Partshell/$case/$type 30 $sgesel $bptsel $cpunum $Project $queue";
			close RASM;
			# }}
		
			###########################################################################################
			# rrbs_dmr.pl for each sample
			###########################################################################################
			# {{
			print RDMR_1 "#!/bin/sh\n";
			print RDMR_1 "perl $Bin/rrbs_dmr_script_BSMAP_1.pl $outdir/Partshell/$case/ $type $normal 30 $sgesel $bptsel $cpunum $Project $queue";
			close RDMR_1;
			print RDMR_2 "#!/bin/sh\n";
			print RDMR_2 "perl $Bin/rrbs_dmr_script_BSMAP_2.pl $outdir/Partshell/$case/ $type $normal 30 $sgesel $bptsel $cpunum $Project $queue";
			close RDMR_2;
			# }}
			
			close ASM_MGT;
			close ASM_SNP;
			close ASM_FILSNP;
			close ASM_ASM;
			close ASM_FILHET;
			close ASM_Fil;
			
			close PRE_AD;
			close PRE_BOWT;
			close PRE_Addread;
			
			close DMR_Dis;
			close DMR_Getmethy;
			close DMR_Sortpe;
			close DMR_Uniq;
			close DMR_Uniqtotal;
			close DMR_Cat;
            close DMR_Catu;
			close DMR_Stat;
			close DMR_Pos;
			close DMR_Qmap;
			close DMR_MGT;
			close DMR_sortpos;
			close DMR_Cout;
			close DMR_Enzyme;
			close DMR_Getdmr;
			close DMR_Draw;
			close DMR_ANNO;    
	
		}

		############################################
		system("chmod +x $outdir/Partshell/$case/$type/*.sh");
    }
}

# Job schedule
if($sgesel == 1) {
	print PRE "foreach \(\@cmd){\n\$pm->start and next\;\nsystem(\$_)\;\$pm->finish\;\n}\n\$pm->wait_all_children\;\n";
	print DMR "foreach \(\@cmd_1){\n\$pm->start and next\;\nsystem(\$_)\;\$pm->finish\;\n}\n\$pm->wait_all_children\;\n";
	print DMR "foreach \(\@cmd_2){\n\$pm->start and next\;\nsystem(\$_)\;\$pm->finish\;\n}\n\$pm->wait_all_children\;\n";
	print ASM "foreach \(\@cmd){\n\$pm->start and next\;\nsystem(\$_)\;\$pm->finish\;\n}\n\$pm->wait_all_children\;\n";
}
else {
	print PRE "foreach \(\@cmd){\n\tsystem(\$_);\n}\n";
	print DMR "foreach \(\@cmd_1){\n\tsystem(\$_);\n}\n";
	print DMR "foreach \(\@cmd_2){\n\tsystem(\$_);\n}\n";
	print ASM "foreach \(\@cmd){\n\tsystem(\$_);\n}\n";
}

close PRE;
close DMR;
close ASM;

################################### 
# Final
################################### 
# {{
system("cp -f $Bin/RRBS_Kit_Run.sh $outdir/SMAP_Run.sh");
system("chmod +x $outdir/SMAP_Run.sh");
system("chmod +x $outdir/RRBS_report.sh");

print("#" x 80 . "\n");
print "Configuration process successfully finished. Please follow the steps below:\n";
print "Enter the output directory:\n";
print "\t\[$outdir\]\n";
print "Execute all the scripts in one step:\n";
print "\t   SMAP_Run.sh\n";
print "Execute the following scripts in order:\n";
print "\t1. RRBS_prepare.pl\n";
print "\t2. RRBS_dmr.pl\n";
print "\t3. RRBS_asm.pl\n";
print "\t4. RRBS_report.sh\n";
print("#" x 80 . "\n");
# }}


