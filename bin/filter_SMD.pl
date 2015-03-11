use strict;
my $i=shift;
open IN, $i or die $!;
while(<IN>){
	chomp;
	my @a=split;
	my $totalnormal=$a[6]+$a[7];
	my $totalcancer=$a[9]+$a[10];
	next if($totalnormal==0);
	next if($totalcancer==0);
    my $ratenormal=$a[6]/$totalnormal;
    my $ratecancer=$a[9]/$totalcancer;
	my $diff=$ratecancer-$ratenormal;
	if($diff >0.15 && $a[12]<0.000001){
		print "Up\t$_\n";	
	}
	if($diff<-0.15 && $a[12] <0.000001){
		print "Down\t$_\n";
	}
}
close IN;
