use strict;
use File::Basename qw(basename dirname);
die "perl $0 lamda errorout\n" unless(@ARGV==2);
my $lamda=shift;
my $error=shift;
if($lamda=~/\.gz$/){
open IN,"gunzip -c $lamda|" or die "no $lamda\n";
}else{
open IN,$lamda or die "no $lamda\n";
}
open OUT,">$error" or die "cant create $error\n";
my $coutdir=dirname($lamda);
my @chr=glob("$coutdir/chr2*.cout*");
print @chr;
my $methy;my $total;
while(<IN>){
	chomp;
	my @tmp=split;
	my $tt=$tmp[6]+$tmp[7];
	$methy+=$tmp[6];
	$total+=$tt;
	
}
close IN;

my $total2;
if($methy>0){
	my $erate=$methy/$total;
	$erate=1-$erate;
	print OUT "Error_Cnot2T_rate\t$methy\t$total\t$erate\n";
}


foreach my $chr(@chr){
	if($chr=~/\.gz$/){
		open IN,"gunzip -c $chr|" or die "no $lamda\n";
	}else{
		open IN,$chr or die "no $lamda\n";
	}
	while(<IN>){
		chomp;
		my @tmp=split;
		my $tt=$tmp[6]+$tmp[7];
		unless($tmp[3] eq 'CG' ){
			$methy+=$tmp[6];
			$total2+=$tt;
		}
	}
	close IN;
}
if($total2==0){
	print STDERR "This-work-is-completed\n";
	#die;
}else{
	my $erate=$methy/$total2;
	print OUT "CHH_methyrate\t$methy\t$total2\t$erate\n";
}							


