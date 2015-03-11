#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory
use File::Basename qw(basename dirname);
use PerlIO::gzip;
use FindBin '$Bin';
use lib "$Bin/lib/compress";
use Compress::Raw::Zlib;
use Compress::Zlib;


my ($outdir,$mode , $Help);
GetOptions(
        "o:s"=>\$outdir,
	"mode:s"=>\$mode,
        "help"=>\$Help
);
#print  "$outdir\n";die;
$mode ||="pe";

print $mode."\n";
if($mode=~/pe/){
die "perl $0 fq1 fq2 -o outdir  \n" if(@ARGV<2  || !$outdir || $Help);
my $read1=shift;
my $read2=shift;
if($read1=~/\.gz$/){
        open RD1, "gunzip -c $read1|" or die "no read1";
}else{
        open (RD1,$read1)||die "no read1\n";
}
if($read2=~/\.gz$/){
        open RD2, "gunzip -c $read2|" or die "no read1";
}else{
        open (RD2,$read2)||die "no read2\n";
}

`mkdir -p $outdir` unless(-d "$outdir");
#open $outchr,">:gzip","$outdir/Methy.qmap.$blank[0].gz" or die
open CTOUT, ">:gzip","$outdir/ori.t.gz" or die $!;
open GAOUT , ">:gzip","$outdir/ori.a.gz" or die $!;
open RAW, ">$outdir/v1.fq" or die $!;

while(<RD1>){
	chomp;
	my $label1=$_;
	my $base1=<RD1>;
	my $symbol1=<RD1>;
	my $quality1=<RD1>;
	chomp ($base1,$symbol1,$quality1);
	
	my $label2=<RD2>;
	my $base2=<RD2>;
        my $symbol2=<RD2>;
        my $quality2=<RD2>;
        chomp ($label2,$base2,$symbol2,$quality2);
	#$rawfq1_2gz->gzwrite($read) ||die $!;
	$label1=~s/^@//;
	$label2=~s/^@//;
	my $read1="$label1\t$base1\t$quality1";
	my $read2="$label2\t$base2\t$quality2";

	print RAW "$read1\t$read2\n";
	$base1=~tr/Cc/Tt/;
	$base2=~tr/Gg/Aa/;
	print CTOUT "\@$label1\n$base1\n$symbol1\n$quality1\n";
	print GAOUT "\@$label2\n$base2\n$symbol2\n$quality2\n";
}
close RD1;
close RD2;
}

if($mode=~/se/){
 die "perl $0 fq1 -o outdir  \n" if(@ARGV<1  || !$outdir || $Help);
my $read1=shift;
if($read1=~/\.gz$/){
        open RD1, "gunzip -c $read1|" or die "no read1";
}else{
        open (RD1,$read1)||die "no read1\n";
}

`mkdir -p $outdir` unless(-d "$outdir");
#open $outchr,">:gzip","$outdir/Methy.qmap.$blank[0].gz" or die
open CTOUT, ">$outdir/ori.t" or die $!;
open RAW, ">$outdir/v1.fq" or die $!;
while(<RD1>){
        chomp;
        my $label1=$_;
        my $base1=<RD1>;
        my $symbol1=<RD1>;
        my $quality1=<RD1>;
        chomp ($base1,$symbol1,$quality1);
        $label1=~s/^@//;
        my $read1="$label1\t$base1\t$quality1";
        print RAW "$read1\n";
        $base1=~tr/Cc/Tt/;
        print CTOUT "\@$label1\n$base1\n$symbol1\n$quality1\n";		
}
close RD1;
close CTOUT;
close RAW;
}

