use strict;
my $snp=shift;
#chr1    762827  .       T       A       133     .       DP=12;VDB=0.0000;AF1=0.5;CI95=0.5,0.5;DP4=3,0,0,9;MQ=42;FQ=18.1;PV4=0.0045,1,1
open SNP,$snp or die "no snp file\n";
my @tmp;
 $tmp[0]=<SNP>;
chomp $tmp[0];
my $cout=0;
while(<SNP>){
	chomp;
	$tmp[1]=$_;
	my @line1=split /\t/,$tmp[0];
	my @line2=split /\t/,$tmp[1];
	if($line2[1]-$line1[1]<=5){
		$tmp[0]=$tmp[1];
		$cout++;
		next;
	}elsif($line2[1]-$line1[1]<=50){
		if($line1[3] eq "T" && $line1[4] eq "C" ||($line1[3]="A" && $line1[4] eq "G")){
			$tmp[0]=$tmp[1];
			$cout++;
			next;
		}else{
			if($cout==0){
				print $tmp[0]."\n";
				$tmp[0]=$tmp[1];
			}else{
				$cout=0;
				$tmp[0]=$tmp[1];
			}
			
		}
	}else{
		if($cout==0){
			print $tmp[0]."\n";
			$tmp[0]=$tmp[1];
		}
		if($cout>=1){
			#print $tmp[1]."\n";
			#$line1=$line2;
			$cout=0;
			$tmp[0]=$tmp[1];
		}		
	}
		

}
my @tmp1=split /\t/,$tmp[0];
my @tmp2=split /\t/,$tmp[2];
if($tmp2[1]-$tmp1[1]>10){
	print "$tmp[1]\n";
}

close SNP;
