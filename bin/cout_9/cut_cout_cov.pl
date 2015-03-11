#  input :< region sorted by chr and startpose >< snp sorted by chr and pose><coverage and depth.out>
#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use PerlIO::gzip;
# perl /ifshk7/BC_TUMOR/PROJECT/qbe130308_gaosj_Methylation/software/RRBS_Kit/bin//cout_9/cut_cout_cov.pl /ifshk7/BC_TUMOR/PROJECT/qbe130308_gaosj_Methylation/software/RRBS_Kit/example/hg19.ccgg.fragment.bed.frag_40-220.example /ifshk7/BC_TUMOR/PROJECT/qbe130308_gaosj_Methylation/software/RRBS_Kit/Testout/K1/MB/cout/chr17.cout.gz /ifshk7/BC_TUMOR/PROJECT/qbe130308_gaosj_Methylation/software/RRBS_Kit/Testout/K1/MB/cout_40-220/chr17.cout.cov 4 40 220 /ifshk7/BC_TUMOR/PROJECT/qbe130308_gaosj_Methylation/software/RRBS_Kit/Testout/K1/MB/cout_40-220/chr17.cout.gz /ifshk7/BC_TUMOR/PROJECT/qbe130308_gaosj_Methylation/software/RRBS_Kit/Testout/K1/MB/cout_40-220/cout_cov/chr17.cout.gz

die "perl $0 <region file> <cout> <out file> <dep> <start> <end> < cut.cout><cov.cout>" if(@ARGV<8);
my $region = shift @ARGV;
my $hapmap = shift @ARGV;
my $out = shift @ARGV;
my $depfilter = shift @ARGV;
my $start  = shift @ARGV;
my $end = shift @ARGV;
my $cout_cut = shift @ARGV;
my $cov_cout=shift @ARGV;
my %regionData;
my %chrIndex;

my @tmp;
my $regStart;
my $regEnd;
my $name = basename $hapmap;
my $chrName;
#chrX.cout
my $chr;
if($name=~/(chr\S+)\.\d+\.cout/){
	 $chr=$1;
}elsif($name=~/(chr\S+)\.cout/){
	$chr=$1;
}

if($region=~/gz$/){
	open (IN,"<:gzip","$region") or die "Cannot open $region\n";
}else{

	open (IN,"$region") or die "Cannot open $region\n";
}
while (<IN>){
	chomp;
	@tmp = split;
	$chrName = $tmp[0];
	last unless($chr);
	next unless ($chrName eq $chr);
	#$tmp[1] =~ m/(\d+)\-(\d+)/;
	#$regStart = $1;
	#$regEnd = $2;
	$regStart = $tmp[1]+2;
	$regEnd = $tmp[2];
	next unless ($tmp[3] >= $start && $tmp[3] <= $end);
	push (@{$regionData{$chrName}},[$regStart,$regEnd]);
	$chrIndex{$chrName} = 0;
}
close (IN);


if($hapmap=~/gz$/){
	open (FILE,"<:gzip","$hapmap") or die "Cannot open $hapmap\n";
}else{	
	open (FILE,"$hapmap") or die "Cannot open $hapmap\n";
}	

#chr8    2       +       CHH     CAA     0       0       0       0       1

my %hash;
$hash{"C"}{"total"}=0;
$hash{"CG"}{"total"}=0;
$hash{"CHG"}{"total"}=0;
$hash{"CHH"}{"total"}=0;
$hash{"C"}{"cov"}=0;
$hash{"CG"}{"cov"}=0;
$hash{"CHG"}{"cov"}=0;
$hash{"CHH"}{"cov"}=0;

my %hash2;
my %hash22;
$hash2{"C"}{"cov"}=0;
$hash2{"CG"}{"cov"}=0;
$hash2{"CHG"}{"cov"}=0;
$hash2{"CHH"}{"cov"}=0;

$hash2{"C"}{"total"}=0;
$hash2{"CG"}{"total"}=0;
$hash2{"CHG"}{"total"}=0;
$hash2{"CHH"}{"total"}=0;

$hash22{"C"}{"cov"}=0;
$hash22{"CG"}{"cov"}=0;
$hash22{"CHG"}{"cov"}=0;
$hash22{"CHH"}{"cov"}=0;

my $site=0;
open Cout , ">:gzip","$cout_cut" or die "$!";
open Cov_cout,">:gzip","$cov_cout"or die $!;


my %emycover;
while (<FILE>){#FILE cout
	chomp;
	@tmp = split;
	$site = $tmp[1];
	my $dep;
	$hash{$tmp[3]}{"total"}++;
	$dep = $tmp[6]+$tmp[7];
	my $level;
	if($dep>=$depfilter){
		$hash{$tmp[3]}{"cov"}++;
	}
	next if (!exists $regionData{$chr});
	for (my $i = $chrIndex{$chr}; $i < @{$regionData{$chr}}; $i++){
		if ($site < $regionData{$chr}->[$i][0]){ $chrIndex{$chr} = $i;last;} #next point
		next if ($site > $regionData{$chr}->[$i][1]);#next region
		my $region="$chr\t$regionData{$chr}->[$i][0]\t$regionData{$chr}->[$i][1]";
		$hash2{$tmp[3]}{"total"}++;
		if($dep>0){
			$emycover{"$region"}=1;		
			$hash2{$tmp[3]}{"cov"}++;
		}
		if($dep>=$depfilter){
			$hash22{$tmp[3]}{"cov"}++;
		}
		print Cout "$_\n";
		$level=$tmp[6]/$dep*100 if $dep>0;
		$level="-" if $dep==0;
 	   	print Cov_cout "$_\t$level\n";	
		$chrIndex{$chr} = $i;
		last;
	}
}
close (FILE);
close(Cout);
close (Cov_cout);

my $cover=0;
foreach my $rname(keys %emycover){
	$cover+=$emycover{$rname};
}
my $allregion; my $emycoverate;

if($chr){
 $allregion=@{$regionData{$chr}};
 $emycoverate="$cover\/$allregion";
}
open OUT , ">$out" or die "$!";
$hash{"C"}{"total"}=$hash{"CG"}{"total"}+$hash{"CHG"}{"total"}+$hash{"CHH"}{"total"};
$hash{"C"}{"cov"}=$hash{"CG"}{"cov"}+$hash{"CHG"}{"cov"}+$hash{"CHH"}{"cov"};
$hash2{"C"}{"total"}=$hash2{"CG"}{"total"}+$hash2{"CHG"}{"total"}+$hash2{"CHH"}{"total"};
$hash2{"C"}{"cov"}=$hash2{"CG"}{"cov"}+$hash2{"CHG"}{"cov"}+$hash2{"CHH"}{"cov"};
$hash22{"C"}{"cov"}=$hash22{"CG"}{"cov"}+$hash22{"CHG"}{"cov"}+$hash22{"CHH"}{"cov"};
#print OUT "chr\tC\tCG\tCHG\tCHG\tCcov\tCGcov\tCHGcov\tCHHcov\tCregion\tCGregion\tCHGregion\tCHHregion\tCcov\tCGcov\tCHGcov\tCHHcov\tCcov\tCGcov(depfilter)\tCHGcov(depfilter)\tCHHcov(depfilter)\n";
print OUT "Chr\tEmycoverate\tTotal_C\tTotal_CG\tTotal_CHG\tTotal_CHH\tCov_C\tCov_CG\tCov_CHG\tCov_CHH\tTotal2_C\tTotal2_CG\tTotal2_CHG\tTotal2_CHH\tCov2_C\tCov2_CG\tCov2_CHG\tCov2_CHH\tCovf_C\tCovf_CG\tCovf_CHG\tCovf_CHH\n";
if($chr){
print OUT "$chr\t$emycoverate\t".$hash{"C"}{"total"}."\t".$hash{"CG"}{"total"}."\t".$hash{"CHG"}{"total"}."\t".$hash{"CHH"}{"total"}."\t".$hash{"C"}{"cov"}."\t".$hash{"CG"}{"cov"}."\t".$hash{"CHG"}{"cov"}."\t".$hash{"CHH"}{"cov"}."\t".$hash2{"C"}{"total"}."\t".$hash2{"CG"}{"total"}."\t".$hash2{"CHG"}{"total"}."\t".$hash2{"CHH"}{"total"}."\t".$hash2{"C"}{"cov"}."\t".$hash2{"CG"}{"cov"}."\t".$hash2{"CHG"}{"cov"}."\t".$hash2{"CHH"}{"cov"}."\t".$hash22{"C"}{"cov"}."\t".$hash22{"CG"}{"cov"}."\t".$hash22{"CHG"}{"cov"}."\t".$hash22{"CHH"}{"cov"}."\n";
}
close OUT;

