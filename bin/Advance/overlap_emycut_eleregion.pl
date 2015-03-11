#!/usr/bin/perl
use strict;
use File::Basename;
use warnings;
use PerlIO::gzip;
my $usage =<<"USAGE";
Command: perl $0 <sorted region file> <gene region file> <outdir>
USAGE
my ($pos_file, $region_file,$outdir) = @ARGV;
die $usage unless $pos_file && $region_file && $outdir;
my %region;
my %index;
my $basepos=basename($pos_file);
my $baseele=basename($region_file);
my @chrname=split /\./,$basepos;
my @elename=split /\./,$baseele;
my $overlap="$chrname[0]\.$elename[0]\.overlap";
open OUT,">$outdir/$overlap" or die $!;

read_region($region_file, \%region, \%index);
parse_file($pos_file, \%region, \%index);

sub read_region{
	my ($file, $ref, $index_ref) = @_;
	my $fh;
	if($file=~/\.gz$/){
		open  $fh,"<:gzip","$file" or die "no $file\n";	
	}else{
		open  $fh,"<$file" or die $!;
	}
	#open  $fh,"<:gzip","$file" if $file=~/\.gz$/ ;
	
	while (<$fh>){
		chomp;
		my @tmp = split;
		#push @{$ref->{$tmp[0]}}, [$tmp[1], $tmp[2],$tmp[3],$tmp[4],$tmp[5]];
		push @{$ref->{$tmp[0]}}, [$tmp[1], $tmp[2],$tmp[3]];
	}
	close $fh;

	for my $chr (keys %$ref){
		@{$ref->{$chr}} = sort {$a->[0] <=> $b->[0]} @{$ref->{$chr}};
		$index_ref->{$chr} = 0;
	}
}
sub parse_file{
    my ($file, $ref, $index_ref) = @_;
    my $fh;
    open  $fh, $file unless $file=~/\.gz$/;
    open  $fh,"<:gzip","$file" if $file=~/\.gz$/ ;
	while (<$fh>){ #chr10   50006   +       CHH     CCT     0       0       0       1
		chomp;
		my @tmp = split;
		next unless defined $index_ref->{$tmp[0]} ;	 ## out of region			
		my $flag = 1;
		my $isOverlap = 0;
		my $i;
		for( $i = $index_ref->{$tmp[0]}; $i < @{$ref->{$tmp[0]}}; $i++){
			if ($tmp[1] < $ref->{$tmp[0]}[$i][0]){$index_ref->{$tmp[0]} = $i;last;}
			next if $tmp[1] > $ref->{$tmp[0]}[$i][1];
			if ($flag){
				$index_ref->{$tmp[0]} = $i;
				$flag = 0;
			}
			$isOverlap = 1;
			last;
		}

		if ($isOverlap){
			print OUT "$_\t$ref->{$tmp[0]}[$i][0]\t$ref->{$tmp[0]}[$i][1]\t$ref->{$tmp[0]}[$i][2]\n";

		}
	}
}

