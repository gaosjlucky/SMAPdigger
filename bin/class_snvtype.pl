use strict;
use File::Basename qw/basename dirname/;
use Cwd 'abs_path';
use FindBin qw($Bin);
use lib "$Bin/lib";
#to class somatic germline and unkown type snp
die "perl $0 cancersnv norm_c2t_mpileup norm_g2a_mpileup out_germline out_somatic out_unknown >threenumber.out\n" unless(@ARGV==6);
my ($cancersnv,$c2t_mp,$g2a_mp,$out_germ,$out_somatic,$out_unknow) = @ARGV;

open IN,$cancersnv or die "no ASM.out\n";
open CMP,"$Bin/bcftools view $c2t_mp|" or die $!;
open GMP,"$Bin/bcftools view $g2a_mp|" or die $!;
open O_GER,">$out_germ" or die "can't create outgermsnv\n";
open O_SOM,">$out_somatic" or die "can't create out_somatic\n";
open O_UN,">$out_unknow" or die "can't create out_unknow";


my %hash;my %gcount;my %scount;
#chr10_100017811 100017749       A       G       A       19      9       6 
while(<IN>){
	chop;
	my @a=split;
	my @b=split /\_/,$a[0];
	$hash{$b[0]}{$b[1]}=$_;
	$gcount{$b[0]}{$b[1]}=0;
	$scount{$b[0]}{$b[1]}=0;
}
close IN;
#chr10	101287837	.	G	X	0	.	DP=1;I16=0,1,0,0,37,1369,0,0,42,1764,0,0,3,9,0,0;QS=1.000000,0.000000,0.000000,0.000000	PL:DP:S0,3,37:1:0
while(<CMP>){
	chop;
	my @a=split;
	if(exists($hash{$a[0]}{$a[1]})){
		if($a[7]=~/DP\=(\d+)\;/){
                        next if($1<9);
                }
		if($a[4] =~/[ACGT]/){
			$gcount{$a[0]}{$a[1]}++;	
		}else{
			$scount{$a[0]}{$a[1]}++;
		}	
	}
	
}
close CMP;

while(<GMP>){
        chop;
        my @a=split;
        if(exists($hash{$a[0]}{$a[1]})){
		if($a[7]=~/DP\=(\d+)\;/){
			next if($1<9);	
		}
		if($a[4] =~/[ACGT]/){
                        $gcount{$a[0]}{$a[1]}++;
                }else{
                        $scount{$a[0]}{$a[1]}++;
                }
        }
        
}
close GMP;
my $germline=0;
my $somatic=0;
my $unknow=0;

foreach my $chr(sort keys %hash){
	foreach my $pos(sort {$a<=>$b} keys %{$hash{$chr}}){
		my $coverage=$gcount{$chr}{$pos}+$scount{$chr}{$pos};
		if($gcount{$chr}{$pos}>0){
			$germline++;
			print O_GER $hash{$chr}{$pos}."\n";
		}elsif($scount{$chr}{$pos}>1){
			print O_SOM $hash{$chr}{$pos}."\n";
			$somatic++;
		}else{
			$unknow++;
			print O_UN $hash{$chr}{$pos}."\n";
		}
	}
}

print "germline\tsomatic\tunknow\n$germline\t$somatic\t$unknow\n";

