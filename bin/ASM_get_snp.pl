use strict;
my $i=shift;
open IN,$i or die "no vcffile\n";
#chr1    763080  .       G       A,T,X   128     .       DP=128;VDB=0.0000;AF1=0.5;CI95=0.5,0.5;DP4=68,0,60,0;MQ=42;FQ=131;PV4=1,1,1,1
while(<IN>){
	chomp;
	my @a=split;
	my $alle=$a[4]=~tr/,/,/;
	#next if($alle>1);
	#next if($a[4]=~/\,/);
	next if($a[5]<10);
#	next if($a[3] eq 'C' && $a[4]=~/T/);
#	next if($a[3] eq 'G' && $a[4]=~/A/);
	
	if(/DP4\=(\S+)\;/){
		my @num=split /\,/,$1;
		my $ref1num=$num[0];
		my $ref2num=$num[1];
		my $var1num=$num[2];
		my $var2num=$num[3];
		#next if($num[2]==0);
		#next if($num[3]==0);
		my $refnum=$ref1num+$ref2num;
		my $varnum=$var1num+$var2num;
		my $total=$refnum+$varnum;
		next if($total<10);
		my $varfreq=$varnum/$total;
		next if($varfreq<0.2);
		next if($varfreq>0.7);
		#next if($ref1num == 0 && $var1num==0);
		#next if($ref2num==0 && $var2num==0);	
		print $_."\n";
	}	


}
close IN;

