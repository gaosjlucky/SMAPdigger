use strict;
use PerlIO::gzip;

die "perl $0 infile outfile\n" unless(@ARGV==2);
my $i=shift;
my $out =shift;

if($i=~/\.gz$/){
	open IN,"gunzip -c $i|" or die "no infile\n";
}else{
	open IN,$i or die "no infile\n";
}
if($out=~/\.gz/){
	open OUT,">:gzip",$out  or die "can't creat outfile\n";
}else{
	open OUT,">$out" or die "can't creat outfile\n";
}
while(<IN>){
	chomp;
	my @a=split /\t/;
	$a[0]=~s/W$//;
	$a[0]=~s/C$//;
	if(@a>12){
		next unless( /\:\+\:len/ && /\:\-\:len/ );
		
		$a[3]="p";
		print OUT join("\t",@a)."\n";
	} elsif(@a<12){
		$a[3]="s";

		if($a[4] =~ /\+\:len(\d+)/){
			$a[2]=$a[1]+$1-1;
		}elsif($a[4]=~/\-\:len(\d+)/){
			$a[1]=$a[2]-$1+1;
		}
		print OUT join("\t",@a)."\n";
	}
}
close IN;		
#chr1W   177210298       177210404       105     FCC01P5ACXX:8:1101:10000:102007#ACTTGAAT/1:rep1:+:len90:W
#chr17H  034095448       034095537       s       >FCC01P5ACXX:8:1101:10004:131444#CGATGTAT/2:rep1:+:len90:Cs     chr17   34095448        CAAAATAACTCTTACTAAAAT
#chr6H   131117709       131117814       p       >FCC01P5ACXX:8:1101:10004:136670#CGATGTAT/1:rep1:+:len90:Wp     chr6    131117709       CGGAAGGTTTTTTAGT
	
