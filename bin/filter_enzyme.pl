#!/usr/bin/perl -w
use strict;
use PerlIO::gzip;
use File::Basename;
use Getopt::Long;

print join("\n",@ARGV)."\n";
die "perl $0 blank hs.cut.emy  posOrder  --insert [40:220]" if (@ARGV<3);
my $blank=shift;
my $emycut=shift;
my $infile=shift;
my $outdir=dirname $infile;
my ($fqlen,$insert,$Help);

GetOptions(
        "fqlen:i"=>\$fqlen,
        "insert:s"=>\$insert,
        "help"=>\$Help,
);
print join("\t",@ARGV)."\n";
die "perl $0  hs.cut.emy  posOrder  --insert [40:220]" if ($Help);
die "perl $0  hs.cut.emy  posOrder  --insert [40:220]" unless ($insert);
#perl /nas/RD_12A/gaoshengjie/software/Methylation/RRBS_Kit/rrbs_kit_pure/bin/Statistics.pl /nas/RD_12A/gaoshengjie/software/Methylation/RRBS_Kit/rrbs_kit_pure/data/hg19.fragment.bed.frag.40-220 /nas/RD_12A/gaoshengjie/software/Methylation/RRBS_Kit/rrbs_kit_pure/qsub/K1/MB/HUMkfoHADDDAAPEI-2/FCC01P5ACXX_L8/withDupli.labelOrder.gz  --insert 40:220

my $length= $fqlen || -1;
#$insert || ="40:220";
my ($short,$long)=split (/\:/,$insert);

if($emycut=~/gz$/){
	open EMC,"<:gzip",$emycut  or die "$!" ;
}else{
	open EMC,$emycut  or die "$!";
}

my %hash;
 my %region;
my %chr;
my %blank;
my $silo=0;
open BK,$blank or die $!;
while(<BK>){
	chomp;
	my @a=split;
	$blank{$a[0]}=1;
}

while(<EMC>){
	chomp;	
	my @a=split;
	next unless(exists($blank{$a[0]}));
	$hash{$a[0]}{$a[1]}="$a[3]\tstart";
	$hash{$a[0]}{$a[2]}="$a[3]\tend";
	open $chr{$a[0]},">:gzip", "$outdir/$a[0]\.gz" unless(exists($chr{$a[0]}));
	if($a[3]>=$short && $a[3]<=$long){
		$region{$a[0]}{$a[1]}="$a[3]\tstart";
		$region{$a[0]}{$a[2]}="$a[3]\tend";
		$silo+=2;
	}
}
close EMC;

if($infile=~/gz$/){
	open IN,"<:gzip",$infile  || die "$!" ;
}else{
	open IN,$infile || die "$!";
}

my $pairtwo=0;
my $pairone=0;
my $pairzero=0;
my $singleone=0;
my $singlezero=0;
open NUM,">$outdir/mapnum.out";
my $headchr;
while(<IN>){
	next unless(/\w/);
	chomp;
	my @a=split /\t/;
	next unless(exists($chr{$a[0]}));
	next if($_=~/chr(\w+)_random/ || $_=~/chr(\w+)_hap/);
	my $lenth1=length($a[7]);
	my $lenth2;
	if($a[4]=~/len(\d+)/){
        	next unless($lenth1 == $1);
	}
	if($a[3] eq "p"){
        	next if($a[5] ne $a[12]);
        	$lenth2=length($a[14]);
        	if($a[11]=~/len(\d+)/){
                	next unless($lenth2 == $1);
       		}

        	my $dist=abs($a[2]-$a[1]);
        	next if($dist>1000);
		&pair($_,$length,\*OUT);
		$headchr=$chr{$a[0]};
		print $headchr $_."\n";
	}
	elsif($a[3] eq "s"){
		&single($_,$length,\*OUT);
		$headchr=$chr{$a[0]};
		print $headchr $_."\n";
	}
}
close IN;

my %emy;
my %emypos;
sub pair{
	my ($line,$len,$o1)=@_;
	
	my $score=0;
	my @a=split /\t/,$line;
	if($len==-1){
		if($a[4]=~/\:len(\d+)\:/) {$len=$1;}
	}

#chr1    7217260 7217331 p       FCC01P5ACXX:8:1207:17713:200479#CGATGTAT/1:rep1:-:len71:C:score12       chr1    7217260 TCAAACGTAATAACAAACACCTATAATCCCAACGACTAAAAAAACTAAAACAACAAAATAACGTAAACCCG  |###|m|##|##||##|||||#||#|||||#|m|||####|#|||#|##||#||#|||##|m|#|||||m ccgggcgtggtggcaggcacctgtagtcccagcgactggggagactgaggcagcagaatggcgtgaacccg C@EEED>FFFFHHFHEHEECJIJGIGDGEHGDHFGJJIJIIJHHIIJIHFHHGIJJJJFHFHHFFFFFC@@ FCC01P5ACXX:8:1207:17713:200479#CGATGTAT/2:rep1:+:len71:C:score30       chr1    7217261                                                                       A                                                                       # cgggcgtggtggcaggcacctgtagtcccagcgactggggagactgaggcagcagaatggcgtgaacccgg BC@FFADFGHHHHJGIJJJJJIJJJIIGIGJJJIJGIJIIJJIGIJGGIIHHIHHHECEFFF>AECCDDDB


	if($a[4]=~/\:\+/){
		my $axis1=$a[6];
		my $end=$axis1+5;
		foreach my $pos($axis1..$end){
	 		if (exists $hash{$a[5]}{$pos}){
				$emy{$a[5]}{$pos}=$hash{$a[5]}{$pos};
				$score++;
				last;
			}
		}
		
		my $end1=$axis1+$len-2;
		my $end2=$axis1+$len-3;	
########################readone	
		if(exists($hash{$a[5]}{$end1})){#last one
			my @reg=split /\t/,$hash{$a[5]}{$end1};
			if($len>=$reg[0] && $reg[1] eq 'end'){
				my @methy=split //,$a[8];
				$methy[-1]=' ';
				$a[8]=join("",@methy);
			}
		}
		if(exists($hash{$a[5]}{$end2})){#last two
			my @reg=split /\t/,$hash{$a[5]}{$end2};
                        if($len>=$reg[0] && $reg[1] eq 'end'){
				my @methy=split //,$a[8];
                                $methy[-1]=" ";
				$methy[-2]=" ";
                                $a[8]=join("",@methy);
                        }
		}						
		my $axis2=$a[13]+$len-1;
		my $start=$axis2-5;
		my $emyend=$axis2-2;
#########################readtwo
		if(exists($hash{$a[5]}{$emyend})){
			my @regionlenth=split /\t/,$hash{$a[5]}{$emyend};
                        my $type=$regionlenth[1];
                        if($type eq 'end'){#last two
				my @methy=split //,$a[15];
                                $methy[-1]=" ";
                                $methy[-2]=" ";
                                $a[15]=join("",@methy);
                        }
		}
		my $emyend2=$axis2-1;
		if(exists($hash{$a[5]}{$emyend2})){ #last one
                        my @regionlenth=split /\t/,$hash{$a[5]}{$emyend2};
                        my $type=$regionlenth[1];
                        if($type eq 'end'){
                                my @methy=split //,$a[15];
                                $methy[-1]=" ";
                                $a[15]=join("",@methy);
                        }
                }

		foreach my $pos($start..$axis2){
			if (exists $hash{$a[5]}{$pos}){
				$emy{$a[5]}{$pos}=$hash{$a[5]}{$pos};
				$score++;
				last;
			}
		}
	}
	elsif($a[4]=~/\:\-\:/){
		my $axis1=$a[6]+$len-1;
		my $start=$axis1-5;
		foreach my $pos($start..$axis1){
                        if( exists $hash{$a[5]}{$pos}){
				$emy{$a[5]}{$pos}=$hash{$a[5]}{$pos};
                                $score++;
                                last;
                        }
                }
########################################	read1
		if(exists($hash{$a[5]}{$a[6]})){
			my @reg=split /\t/,$hash{$a[5]}{$a[6]};
			if($len>=$reg[0] && $reg[1] eq 'start'){
                                my @methy=split //,$a[8];
                                $methy[0]=" ";
                                $methy[1]=" ";
                                $a[8]=join("",@methy);
                        }
		
		}
		my $onebase=$a[6]-1;
		if(exists($hash{$a[5]}{$onebase})){
			my @reg=split /\t/,$hash{$a[5]}{$onebase};
			if($len>=$reg[0] && $reg[1] eq 'start'){
				my @methy=split //,$a[8];
                                $methy[0]=" ";
				$a[8]=join("",@methy);
			}			
		}

				
		
######################################  read2 		
		my $axis2=$a[13];
		my $end=$axis2+5;
		my $startone=$a[13]-1;

		if(exists($hash{$a[5]}{$a[13]})){
			my @reg2=split /\t/,$hash{$a[5]}{$a[13]};
                        my $type=$reg2[1];
                        if($type eq 'start'){
                                my @methy=split //,$a[15];
                                $methy[0]=" ";
				$methy[1]=" ";
                                $a[15]=join("",@methy);
                        }
		}
		if(exists($hash{$a[5]}{$startone})){
                        my @reg2=split /\t/,$hash{$a[5]}{$startone};
                        my $type=$reg2[1];
                        if($type eq 'start'){
                                my @methy=split //,$a[15];
                                $methy[0]=" ";
                                $a[15]=join("",@methy);
                        }
                }

		
		foreach my $pos($axis2..$end){
			if (exists $hash{$a[5]}{$pos}){
				$emy{$a[5]}{$pos}=$hash{$a[5]}{$pos};
				$score++;
                                last;
			}
		}
	}
	

	$line=join("\t",@a);	
	if($score==2){
		$pairtwo++;
		#print $$o1 "$line\n"; 
	}
	elsif($score==1){
		$pairone++;
		#print $$o1 "$line\n";
	}
	elsif($score==0){
		$pairzero++;
		#print $$o1 "$line\n";
	}
}

sub single{
	my ($line,$len,$o1)=@_;
	my $score=0;
	my @a=split /\t/,$line;
	if($len==-1){
                $a[4]=~/\:len(\d+)\:/; $len=$1;
        }

	if($a[4]=~/\+/){ 
		my $axis1=$a[6] ;
		my $end=$axis1+5;
		foreach my $pos($axis1..$end){
			if (exists $hash{$a[5]}{$pos}){
				$emy{$a[5]}{$pos}=$hash{$a[5]}{$pos};
				$score++;
				last;
			}
		}
	}	
	elsif($a[4]=~/\-/){
		my $axis2=$a[6]+$len-1;
		my $start=$axis2-5;
		foreach my $pos($start..$axis2){
                        if (exists $hash{$a[5]}{$pos}){
				$emy{$a[5]}{$pos}=$hash{$a[5]}{$pos};
                                $score++;
                                last;
                        }
                }

	}
	if($score==1){
		$singleone++;
		#print $$o1 "$line\n"
	}
	elsif($score==0){
		$singlezero++;
		#print $$o1 "$line\n"

	}
}
my $frag;
foreach my $chr(keys %emy){
	foreach my $pos(keys %{$emy{$chr}}){
		$frag++;
	}
}

my $rate_emy="$frag\/$silo";


print NUM "total\tenzyme\tpairtwo(readsnum)\tpairone(PEnum)\tpairzero(PEnum)\tsingleone\tsinglezero\tration_enzyme\n";
my $total= $pairtwo*2+$pairone*2+$pairzero*2+$singleone+$singlezero;
my $enzyme = $pairtwo*2+$pairone+$singleone;
my $pair = $pairtwo*2;
print NUM "$total\t$enzyme\t$pair\t$pairone\t$pairzero\t$singleone\t$singlezero\t$rate_emy\n";

