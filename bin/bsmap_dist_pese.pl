use strict;
use FindBin '$Bin';
die "perl $0 insam pe se \n" unless(@ARGV==3);
my $i=shift;
my $pe=shift;
my $se=shift;

if($i=~/.bam$/){
	open IN,"$Bin/samtools view -h $i|" or die $!;
}else{

	open IN,$i or die $!;
}
open PE,">$pe" or die "can't create filepe\n";
open SE,">$se" or die "can't create filese\n";

=head
 p=0x1 (paired), P=0x2 (properly paired), u=0x4 (unmapped),
     U=0x8 (mate unmapped), r=0x10 (reverse), R=0x20 (mate reverse)
     1=0x40 (first), 2=0x80 (second), s=0x100 (not primary), 
     f=0x200 (failure) and d=0x400 (duplicate)

=cut

while(<IN>){
	chomp;
	my @a=split;
	if(/^@/){
		print SE $_."\n";
		print PE $_."\n";
		next;
	}	
	if($a[1]&0x2){
		print PE $_."\n";
	}else{
		print SE $_."\n";
	}

}
close IN;
