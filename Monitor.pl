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
" . "#" x 80 . "\n" if($HELP || !($CONFIG && &abs_path($CONFIG)) || !($OUT_DIR));

if($OUT_DIR){
	`mkdir -p $OUT_DIR` unless(-d $OUT_DIR);
}

my $lib=abs_path($CONFIG);
my $outdir=abs_path($OUT_DIR);

my ($Bin,$rawref,$refc2t,$refg2a,$adaptor,$rscript,$emycut,$blank,$mpileupregion,$dbsnp,$region);
my ($alignsoft,$javasoft,$samtools,$bcftools,$bissnpsoft,$vcfutils,$picardmergesam,$picardsort,$picardaddread,$insert1,$insert2);
my ($Mode,$Bpt,$SgeProj,$SgeQueue,$SlurmPart);
my @element;
my %sample;
my $Target;
my $refdir;
my %samplese;

unless(-d "$outdir"){
	`mkdir $outdir`;
}
	
open LIB,$lib or die "Unable to open configure file!\n";
while(<LIB>){
	chomp;
	next if(/\#/);
	
######################### Variables
	#Working mode : multicore, sge, slurm
	if(/^Mode/){
		my @tmp=split /\=/;
        $tmp[1]=~s/^\s+//;
        $Mode=$tmp[1];
    }
	
	#Breakpoint switch : on, off
	if(/^Bpt/){
		my @tmp=split /\=/;
        $tmp[1]=~s/^\s+//;
        $Bpt=$tmp[1];
    }
	
	#Sge project (qsub command arguments "-p")
	if(/^SgeProj/){
		my @tmp=split /\=/;
		$tmp[1]=~s/^\s+//g;
		if($tmp[1] eq "") {
			$SgeProj="";
		}
		else {
			$SgeProj=$tmp[1];
		}
	}

	#Sge queue (qsub command arguments "-q")
	if(/^SgeQueue/){
		my @tmp=split /\=/;
		$tmp[1]=~s/^\s+//g;
		if($tmp[1] eq "") {
			$SgeQueue="";
		} 
		else {
			$SgeQueue=$tmp[1];	
		}
	}
	
	#Slurm partition (sbatch command arguments "-p")
	if(/^SlurmPart/){
		my @tmp=split /\=/;
		$tmp[1]=~s/^\s+//g;
		if($tmp[1] eq "") {
			$SlurmPart="";
		} 
		else {
			$SlurmPart=$tmp[1];	
		}
	}
	
######################### Softwares
	if(/^Alignsoft/){
		my @tmp=split /\=/;
		$alignsoft=$tmp[1];
		$alignsoft=~s/\s+$//;
		$alignsoft=~s/^\s+//;
	}

	if(/^Javasoft/){
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
    if(/^Region/){
        my @region=split /\=/;
        $region[1]=~s/\s+//g;
        $region=$region[1];
		$region=abs_path($region);
		die "Invalid region file!\n" unless($region && -e $region);
    }
	
	if(/^Adapter/){
		my @tmp=split /\=/;
		$adaptor=$tmp[1];
		$adaptor=~s/\s+//g;
		if($adaptor ne "NULL"){
			$adaptor=abs_path($adaptor);
			die "Error! Invalid adaptor file!\n" unless($adaptor && -e $adaptor);
		}
	}

	if(/^Reference/){
		my @reference=split /\=/;
		$reference[1]=~s/\s+//;
		$rawref=$reference[1];	
		$refdir=dirname $rawref;	
	}
	
    if(/^Target/){
		my @targetnumber=split /\=/;
		$targetnumber[1]=~s/\s+//g;
		$Target=$targetnumber[1];
		$Target=abs_path($Target);
		die "Invalid target file 1!\n" unless($Target && -e $Target);
    }
	
    if(/^Annodir/){
        my @direlement=split /\=/;
       	$direlement[1]=~s/\s+//; 
        my $direlement=abs_path($direlement[1]);
        @element=glob("$direlement/*");
    }
	
	if(/^Sample/){
		my @sp=split /\=/;
		$sp[1]=~s/^\s+//;
		my @saminfo=split /\s+/,$sp[1];
		$saminfo[5] = abs_path($saminfo[5]);
		$saminfo[6] = abs_path($saminfo[6]);
		if($saminfo[2] eq "PE"){
            die "Error! Invalid sample path of $saminfo[5]!\n" unless(-e $saminfo[5]);
            die "Error! Invalid sample path of $saminfo[6]!\n" unless(-e $saminfo[6]);
			$sample{$saminfo[0]}{$saminfo[1]}{$saminfo[3]}{$saminfo[4]}="$saminfo[2]\t$saminfo[5]\t$saminfo[6]\t$saminfo[7]";
		}
		elsif($saminfo[2] eq "SE"){
            die "Error! Invalid sample path of $saminfo[5]!\n" unless(-e $saminfo[5]);
			$sample{$saminfo[0]}{$saminfo[1]}{$saminfo[3]}{$saminfo[4]}="$saminfo[2]\t$saminfo[5]";
		}	
	}
}
close LIB;

$dbsnp ||=" ";
my $nthread;
if($Mode eq 'multicore') {
	$nthread = 2;
}
else {
	$nthread = 4;
}
print "Working at $Mode mode.\n";
my $shPath = "#!/bin/sh\n";

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
if($alignsoft=~/bsmap/) {
	$alignsoft="$Bin/bsmap";
}
elsif($alignsoft=~/bismark/i){
	$alignsoft="$Bin/bismark/bismark";
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
# others
$samtools=$Bin."/samtools";
$bcftools=$Bin."/bcftools";
$vcfutils=$Bin."/vcfutils.pl";
###############################
# Common Data
###############################
system ("mkdir -p $outdir/Partshell/Parcore") unless (-d "$outdir/Partshell/Parcore");
my $parDir = "$outdir/Partshell/Parcore";

open BS_1, ">$parDir/BS_1.sh" or die $!;
open BS_2, ">$parDir/BS_2.sh" or die $!;
open BS_3, ">$parDir/BS_3.sh" or die $!;
open BS_4, ">$parDir/BS_4.sh" or die $!;
open BS_5, ">$parDir/BS_5.sh" or die $!;
open BS_6, ">$parDir/BS_6.sh" or die $!;
open BS_7, ">$parDir/BS_7.sh" or die $!;
open BS_8, ">$parDir/BS_8.sh" or die $!;          
open BS_9, ">$parDir/BS_9.sh" or die $!;

# Top shell
open PAR_SH, ">$outdir/RRBS_Run.sh";		
print PAR_SH "perl $Bin/parExe.pl -outdir $outdir -chkint 30 -mode $Mode -bpt $Bpt -nthread $nthread -sgeproj $SgeProj -sgequeue $SgeQueue -slurmpart $SlurmPart";
	
foreach my $case(keys %sample){
	system ("mkdir -p $outdir/Partshell/$case") unless (-d "$outdir/Partshell/$case");
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
		system ("mkdir -p $outdir/$case/$type/ASM/") unless (-d "$outdir/$case/$type/ASM/");
		system ("mkdir -p $outdir/Partshell/$case/$type") unless (-d "$outdir/Partshell/$case/$type");
        system ("mkdir -p $outdir/$case/$type/DMR") unless (-d "$outdir/$case/$type/DMR");
        
        open PRE_AD, ">$outdir/Partshell/$case/$type/BS-1-dealadaptor.sh" or die $!;
        open PRE_BOWT, ">$outdir/Partshell/$case/$type/BS-2-bsmap.sh" or die $!;
        open PRE_Addread, ">$outdir/Partshell/$case/$type/BS-3-addread.sh" or die $!;
		print BS_1 "$outdir/Partshell/$case/$type/BS-1-dealadaptor.sh\n";
		print BS_2 "$outdir/Partshell/$case/$type/BS-2-bsmap.sh\n";
		print BS_3 "$outdir/Partshell/$case/$type/BS-3-addread.sh\n";
            
        # Merge
        open BS_MGT, ">$outdir/Partshell/$case/$type/BS-4-mergetotal.sh" or die $!;
        open FASTQ, ">$outdir/Partshell/$case/$type/BS-5-fastqc.sh" or die $!;
        open BSSNP, ">$outdir/Partshell/$case/$type/BS-6-bssnper.sh" or die $!;
		print BS_4 "$outdir/Partshell/$case/$type/BS-4-mergetotal.sh\n";
        print BS_5 "$outdir/Partshell/$case/$type/BS-5-fastqc.sh\n";
        print BS_6 "$outdir/Partshell/$case/$type/BS-6-bssnper.sh\n";
        
		# ASM
        open BS_ASM, ">$outdir/Partshell/$case/$type/ASM-7-callasm.sh" or die $!;
        open ASM_FIL, ">$outdir/Partshell/$case/$type/ASM-8-asmfil.sh" or die $!;
		print BS_7 "$outdir/Partshell/$case/$type/ASM-7-callasm.sh\n";
        print BS_8 "$outdir/Partshell/$case/$type/ASM-8-asmfil.sh\n";
            
        # DMR
        open DMR_call, ">$outdir/Partshell/$case/$type/DMR-7-calldmr.sh" or die $!;
        open DMR_ANNO, ">$outdir/Partshell/$case/$type/DMR-8-anno.sh" or die $!;
        open DMR_summary, ">$outdir/Partshell/$case/$type/DMR-9-summary.sh" or die $!;
		print BS_7 "$outdir/Partshell/$case/$type/DMR-7-calldmr.sh\n";
        print BS_8 "$outdir/Partshell/$case/$type/DMR-8-anno.sh\n";
        print BS_9 "$outdir/Partshell/$case/$type/DMR-9-summary.sh\n";
        
        
        
	
        my @samuniq;
        foreach my $lib(keys %{$sample{$case}{$type}}){
            foreach my $flow(keys %{$sample{$case}{$type}{$lib}}){
                my @fq=split /\t/,$sample{$case}{$type}{$lib}{$flow};
                system ("mkdir -p $outdir/$case/$type/$lib/$flow");
					
                if($fq[0]=~/PE/){
					my $adp=$adaptor;
					# hellbelly
					if($adaptor eq "NULL"){
						print PRE_AD "ln -s $fq[1] $outdir/$case/$type/$lib/$flow/t_clean1.fq\n";
						print PRE_AD "ln -s $fq[2] $outdir/$case/$type/$lib/$flow/a_clean2.fq\n";
					}
					else{
						print PRE_AD "perl $Bin/Dealadaptor.pl $adp $fq[1] $fq[2] $outdir/$case/$type/$lib/$flow/t_clean1.fq $outdir/$case/$type/$lib/$flow/a_clean2.fq >$outdir/$case/$type/$lib/$flow/adaptor.info \n";
					}
					if($alignsoft =~ /bsmap/i){
						print PRE_BOWT "$alignsoft $alignval -d $rawref  -a $outdir/$case/$type/$lib/$flow/t_clean1.fq  -b $outdir/$case/$type/$lib/$flow/a_clean2.fq |$Bin/samtools view -bS - > $outdir/$case/$type/$lib/$flow/outfile.bam\n";
                        print PRE_Addread "$picardaddread I=$outdir/$case/$type/$lib/$flow/outfile.bam O=$outdir/$case/$type/$lib/$flow/outfileadd.bam  LB=$lib PU=$flow SM=$case\-$type PL=$fq[-1]\n";
                        
                    }
					elsif($alignsoft =~ /bismark/i){
						print PRE_BOWT "$alignsoft $alignval --path_to_bowtie $Bin $refdir --temp_dir $outdir/$case/$type/$lib/$flow/  -o $outdir/$case/$type/$lib/$flow/ -1 $outdir/$case/$type/$lib/$flow/t_clean1.fq  -2 $outdir/$case/$type/$lib/$flow/a_clean2.fq \n";
                        print PRE_Addread "$picardaddread I=$outdir/$case/$type/$lib/$flow/t_clean1.fq_bismark_bt2_pe.sam  O=$outdir/$case/$type/$lib/$flow/outfileadd.bam  LB=$lib PU=$flow SM=$case\-$type PL=$fq[-1]\n";
                    }
					
                    push @samuniq,"I=$outdir/$case/$type/$lib/$flow/outfileadd.bam";
				}
					
				if($fq[0]=~/SE/){
					my $adp=$adaptor;
					if($adaptor eq "NULL"){
						print PRE_AD "ln -s $fq[1] $outdir/$case/$type/$lib/$flow/t_clean1.fq\n";
					}
					else{
						print PRE_AD "perl $Bin/Dealadaptor.pl $adp $fq[1]  $outdir/$case/$type/$lib/$flow/t_clean1.fq >$outdir/$case/$type/$lib/$flow/adaptor.info -mode se \n";
                    }    
                    if($alignsoft =~ /bsmap/i){
					print PRE_BOWT "$alignsoft $alignval -d $rawref  -a $outdir/$case/$type/$lib/$flow/t_clean1.fq |$Bin/samtools view -bS - > $outdir/$case/$type/$lib/$flow/outfile.bam\n";
                        print PRE_Addread "$picardaddread I=$outdir/$case/$type/$lib/$flow/outfile.bam O=$outdir/$case/$type/$lib/$flow/outfileadd.bam  LB=$lib PU=$flow SM=$case\-$type PL=$fq[-1]\n";
                    }
					elsif($alignsoft =~ /bismark/i){
						print PRE_BOWT "$alignsoft $alignval  --path_to_bowtie $Bin $refdir --temp_dir $outdir/$case/$type/$lib/$flow/  -o $outdir/$case/$type/$lib/$flow/  $outdir/$case/$type/$lib/$flow/t_clean1.fq\n";
                        print PRE_Addread "$picardaddread I=$outdir/$case/$type/$lib/$flow/t_clean1.fq_bismark_bt2.sam O=$outdir/$case/$type/$lib/$flow/outfileadd.bam  LB=$lib PU=$flow SM=$case\-$type PL=$fq[-1]\n";
                        
                    }
					
					push @samuniq,"I=$outdir/$case/$type/$lib/$flow/outfileadd.bam";
				}
			}
		}
		
		##first 3 steps over
        #Merge
		print BS_MGT "$picardmergesam ".join ("\t",@samuniq)." O=$outdir/$case/$type/ASM/BSMAP.sort.bam  SO=coordinate TMP_DIR=$outdir/$case/$type/ASM \&\& $samtools index $outdir/$case/$type/ASM/BSMAP.sort.bam\n";
		print FASTQ "$Bin/FastQC/fastqc   $outdir/$case/$type/ASM/BSMAP.sort.bam\n";
		print BSSNP "perl $Bin/BS-Snper.pl --fa $rawref --input $outdir/$case/$type/ASM/BSMAP.sort.bam --output $outdir/$case/$type/ASM/snp.candidate --methcg $outdir/$case/$type/ASM/meth.cg --methchg $outdir/$case/$type/ASM/meth.chg --methchh $outdir/$case/$type/ASM/meth.chh --minhetfreq 0.1 --minhomfreq 0.85 --minquali 15 --mincover 10 --maxcover 1000 --minread2 2 --errorate 0.02 --mapvalue 20 >$outdir/$case/$type/ASM/SNP.out\n";
            
		###ASM
		print BS_ASM "perl $Bin/ASM_get_CpG.pl $rawref $outdir/$case/$type/ASM/SNP.out $outdir/$case/$type/ASM/BSMAP.sort.bam >$outdir/$case/$type/ASM/ASM.out\n";
		print ASM_FIL "$rscript --no-save --silent --slave < $Bin/ASM_snpfreq.R  --args $outdir/$case/$type/ASM/ASM.out $outdir/$case/$type/ASM/ASM.fil.out\n";
        
		##DMR
		system("mkdir $outdir/$case/$type/DMR") unless(-e "$outdir/$case/$type/DMR");
        if($type ne $normal){
            print DMR_call "$Bin/dmr $region $outdir/$case/$normal/ASM/meth.cg  $outdir/$case/$type/ASM/meth.cg $outdir/$case/$type/DMR/dmr.out\n";
            foreach my $elementone(@element){
                print DMR_ANNO "perl $Bin/anno_dmr.pl $elementone  $outdir/$case/$type/DMR/dmr.out\n";
            }
        }
		print DMR_summary "perl $Bin/CpG_QC.pl $Target $outdir/$case/$type/ASM/meth.cg $rscript\n";
		system("chmod +x $outdir/Partshell/$case/$type/*.sh");
	}
}
system("chmod -R 755 $Bin");
close PAR_SH;
close BS_1;
close BS_2;
close BS_3;
close BS_4;
close BS_5;
close BS_6;
close BS_7;
close BS_8;
close BS_9;


=head
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
=cut


