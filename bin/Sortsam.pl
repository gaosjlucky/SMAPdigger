use strict;
my $indir=shift;

open C2T, "$indir/C2T.sam" or die $!;
open G2A, "$indir/G2A.sam" or die $!;
open H1, ">$indir/c2t.head" or die $!;
open H2, ">$indir/g2a.head" or die $!;
open C2TS, ">$indir/c2t.sort.sam" or die $!;
open G2AS, ">$indir/g2a.sort.sam" or die $!;

while(<C2T>){
	chomp;
	if(/^@/){
		print H1 $_."\n";
	}else{
		print C2TS $_."\n";
	}
}

while(<G2A>){
	chomp;
	if(/^@/){
		print H2 $_."\n";
	}else{
		print G2AS $_."\n";
	}
}

close C2T;
close C2TS;
close G2A;
close G2AS;
close H1;
close H2;

system ("sort -s -k 1,1 -T $indir $indir/c2t.sort.sam -o $indir/c2t.sort.sam");
system ("sort -s -k 1,1 -T $indir $indir/g2a.sort.sam -o $indir/g2a.sort.sam");

system("cat $indir/c2t.head $indir/c2t.sort.sam >$indir/C2T.sam");
system("cat $indir/g2a.head $indir/g2a.sort.sam >$indir/G2A.sam");
system("rm -f $indir/c2t.head $indir/g2a.head $indir/c2t.sort.sam $indir/g2a.sort.sam" );
