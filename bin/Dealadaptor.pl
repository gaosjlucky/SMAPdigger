#!/usr/bin/perl


=head1 Name

adaptor.pl

=head1 Description
clean adaptor

=head1 Version

  Author: Gaosj, gaoshengjie@genomics.org.cn
  Version: 1.4,  Date: 2010-11-7 Modify: 2011-9-15
		

=head1 Comman-line Option
 
 --n <num>	the max number of N in read
 --seed <num>	the min legth matched in adaptor
 --lowqrate <num>	the max rate of low q base
=head1 Usage

 perl $0  adaptor read1 read2 out1 out2 >adaptor.info option

=head1 Exmple

perl $0 adaptor.fa XXX.fq1.gz XXX.fq2.gz cleanfq1.gz cleanfq2.gz -n 10 -seed 10 -lowqrate 0.5 -mode pe

=cut



use strict;
use Data::Dumper;
use File::Basename qw(basename dirname);
use File::Path;
use FindBin '$Bin';
#use lib "$Bin/lib/compress";
#use Compress::Raw::Zlib;
#use Compress::Zlib;
use Getopt::Long;

my %hash;
my ($nnum,$seed,$lowqrate,$Help,$mode,$quli_limit);
my ($totalread,$cleanread,$totalbase,$cleanbase);


GetOptions(
        "n:i"=>\$nnum,
	"seed:i"=>\$seed,
	"mode:s"=>\$mode,
	"lowqrate:f"=>\$lowqrate,
	"qlmit:i"=>\$quli_limit,
	"h!"=>\$Help,
	
);

$nnum ||=10;
$seed ||=10;
$lowqrate ||=0.5;
$mode ||="pe";
$quli_limit ||=70;
#print "$nnum\t$seed\t$lowqrate\n";
#die "perl $0 adaptor read1 read2 out1 out2 option\n -n \[10\] -seed \[10\] \n"  unless($adaptor && $read1 && $read2 && $clean1 && $clean2 );
if($mode=~/pe/){
my $adaptor=shift;
my $read1=shift;
my $read2=shift;
my $clean1=shift;
my $clean2=shift;
die `pod2text $0` unless($adaptor && $read1 && $read2 && $clean1 && $clean2 );
die `pod2text $0` if ($Help);

open (AD,$adaptor)||die "no adptor\n";
while(<AD>){
	chomp;
	next if(/\>/ );
	next unless(/\w/);
	my $seq=$_;
	my $i;
	my $len=length($seq);
	my $end=$len-$seed;
	my $reverse=reverse($seq);
	$reverse=~tr/ACGT/TGCA/;
	for($i=0;$i<$end;$i++){
		my $seedseq=substr($seq,$i,$seed);
		$hash{$seedseq}=1;
		my $reverseq=substr($reverse,$i,$seed);
		$hash{$reverseq}=1;
	}
}
close AD;
print  STDERR "adaptor read done\n";
if($read1=~/\.gz$/){
	die "Can't read read1" if(!-e $read1);
	open RD1, "gunzip -c $read1|" or die "no read1";
}else{
	open (RD1,$read1)||die "no read1\n";
}
if($read2=~/\.gz$/){
	die "Can't read read2" if(!-e $read2);
	open RD2, "gunzip -c $read2|" or die "no read1";
}else{
	open (RD2,$read2)||die "no read2\n";
}
open (CLN1,">$clean1") || die "can not creat clean1\n";
open (CLN2,">$clean2") || die "cant creat clean2\n";
#$fq1gz = gzopen("$pathname\_1.fq.gz","rb") or die "ERROR: cann't ungzip fq: $!";
#$fq2gz = gzopen("$pathname\_2.fq.gz","rb") or die "ERROR: cann't ungzip fq: $!";

#my $clean1gz = gzopen("$clean1","wb") or die $!;
#my $clean2gz = gzopen("$clean2","wb") or die $!;



while(<RD1>){
	chomp;
	$totalread+=2;
	my $read1name=$_;
	my $seq1=<RD1>;
	chomp $seq1;
	my $sym1=<RD1>;
	chomp $sym1;
	my $qual1=<RD1>;
	chomp $qual1;

	my $read2name=<RD2>;
	chomp $read2name;
	my $seq2=<RD2>;
        chomp $seq2;
	my $sym2=<RD2>;
        chomp $sym2;
        my $qual2=<RD2>;
        chomp $qual2;
	my $location1=0;
	my $seed1=0;
	my $location2=0;
	my $seed2=0;
	my $first1=0;
	my $first2=0;	

	my $n1=$seq1=~s/N/N/g;
	my $n2=$seq2=~s/N/N/g;
	#print $n1;die;
	if($n1>10 || $n2>10){
		#print "$read1name\_$read2name\tNtomany\n";
		$cleanread++;
		next;
	}

	my @qu1=split //,$qual1;
	my @qu2=split //,$qual2;
	my $lowbase1;my $lowbase2;
	foreach my $baseq(@qu1){
		if(ord($baseq)<$quli_limit){
			$lowbase1++;
		}
	}
	
	foreach my $baseq(@qu2){
                if(ord($baseq)<$quli_limit){
                        $lowbase2++;
                }
        }	
	my $lengthread=length($seq1);
	$totalbase+=$lengthread*2;
	my $lowrate1=$lowbase1/$lengthread; 
	my $lowrate2=$lowbase2/$lengthread;
	if($lowrate1>$lowqrate || $lowrate2>$lowqrate){
		#print "$read1name\_$read2name\tLowqtoomany\n";
		$cleanread++;
		next;		
	}
		
	my $end=length($seq1)-$seed;
	my @pos1;my @pos2;
	for(my $i=0;$i<$end;$i++){
		my $seedseq1=substr($seq1,$i,$seed);
		my $seedseq2=substr($seq2,$i,$seed);
		
		
		if(exists($hash{$seedseq1})){
			
			$location1=$i;
			$seed1++;
			push @pos1,$i;
			
		}
		if(exists($hash{$seedseq2})){
                        $location2=$i;
                        $seed2++;
			push @pos2,$i;
                }
	}
	my $adplen1=$pos1[-1]-$pos1[0]+1+$seed;
	my $adplen2=$pos2[-1]-$pos2[0]+1+$seed;

	my $startlen1=$pos1[0]+1;
	my $startlen2=$pos2[0]+1;
	
	my $endlen1=$lengthread-$pos1[-1]-$seed;
	my $endlen2=$lengthread-$pos2[-1]-$seed;
	
	my $endpos1=$pos1[-1]+$seed;
	my $endpos2=$pos2[-1]+$seed;

	my $startbase1=substr($seq1,0,$startlen1);
	my $startbase2=substr($seq2,0,$startlen2);

	my $endbase1=substr($seq1,$endpos1,$endlen1);
	my $endbase2=substr($seq2,$endpos2,$endlen2);

	
	my $startqual1=substr($qual1,0,$startlen1);
	my $startqual2=substr($qual2,0,$startlen2);
	
        my $endqual1=substr($qual1,$endpos1,$endlen1);
        my $endqual2=substr($qual2,$endpos2,$endlen2);	
 	
	if($location1!=0 && $pos1[0]<31 && $endlen1<31){
		#print "$read1name\t$pos1[0]\t$adplen1\ttwoadaptor\n";
		$cleanread++;
		next;
	}
	if($location2!=0 && $pos2[0]<31 && $endlen2<31){
		#print "$read2name\t$pos2[0]\t$adplen2\ttwoadaptor\n";
		$cleanread++;
		next;
	}

	if($location1==0){
		#print CLN1 "$read1name\n$seq1\n$sym1\n$qual1\n";
		my $read1="$read1name\n$seq1\n$sym1\n$qual1\n";
		#$clean1gz->gzwrite($read1) || die $!;
		print CLN1 $read1 || die $!;
		$cleanbase+=length($seq1);
	}elsif($startlen1>=31){
		my $read1="$read1name\n$startbase1\n$sym1\n$startqual1\n";
		#print CLN1 "$read1name\n$startbase1\n$sym1\n$startqual1\n";
		#$clean1gz->gzwrite($read1) || die $!;
		print CLN1 $read1 || die $!;
		$cleanbase+=length($startbase1);
		#print "$read1name\t$pos1[0]\t$adplen1\n";#start from 0
	}elsif($startlen1<31){
		#print CLN1 "$read1name\n$endbase1\n$sym1\n$endqual1\n";
		my $read1="$read1name\n$endbase1\n$sym1\n$endqual1\n";
		#$clean1gz->gzwrite($read1) || die $!;
		print CLN1 $read1 || die $!;
		$cleanbase+=length($endbase1);
		#print "$read1name\t$pos1[0]\t$adplen1\n";#start fron 0
	}

	if($location2==0){
		#print CLN2 "$read2name\n$seq2\n$sym2\n$qual2\n";
		my $read2="$read2name\n$seq2\n$sym2\n$qual2\n";
		#$clean2gz->gzwrite($read2) || die $!;
		print CLN2 $read2 || die $!;
		$cleanbase+=length($seq2);
	}elsif($startlen2>=31){
                #print CLN2 "$read2name\n$startbase2\n$sym1\n$startqual2\n";
		my $read2="$read2name\n$startbase2\n$sym1\n$startqual2\n";
		#$clean2gz->gzwrite($read2) || die $!;		
		print CLN2 $read2 || die $!;
		$cleanbase+=length($startbase2);
		#print "$read2name\t$pos2[0]\t$adplen2\n";
        }elsif($startlen2<31){
                #print CLN2 "$read2name\n$endbase2\n$sym2\n$endqual2\n";
		my $read2="$read2name\n$endbase2\n$sym2\n$endqual2\n";
		#$clean2gz->gzwrite($read2) || die $!;
		print CLN2 $read2 || die $!;
		$cleanbase+=length($endbase2);
		#print "$read2name\t$pos2[0]\t$adplen2\n";

        }
}

close CLN1;
close CLN2;
close RD1;
close RD2;
$cleanread=$totalread-$cleanread;
print "Totalread\tTotalbase\tCleanread\tCleanbase\n$totalread\t$totalbase\t$cleanread\t$cleanbase\n";
}
##################################################################

if($mode=~/se/){
my $adaptor=shift;
my $read1=shift;
my $clean1=shift;
die `pod2text $0` unless($adaptor && $read1  && $clean1 );
die `pod2text $0` if ($Help);
open (AD,$adaptor)||die "no adptor\n";
while(<AD>){
        chomp;
        next if(/\>/ );
        next unless(/\w/);
        my $seq=$_;
        my $i;
        my $len=length($seq);
        my $end=$len-$seed;
        my $reverse=reverse($seq);
        $reverse=~tr/ACGT/TGCA/;
        for($i=0;$i<$end;$i++){
                my $seedseq=substr($seq,$i,$seed);
                $hash{$seedseq}=1;
                my $reverseq=substr($reverse,$i,$seed);
                $hash{$reverseq}=1;
        }
}
close AD;
print  STDERR "adaptor read done\n";

if($read1=~/\.gz$/){
        open RD1, "gunzip -c $read1|" or die "no read1";
}else{
        open (RD1,$read1)||die "no read1\n";
}
#my $clean1gz = gzopen("$clean1","wb") or die $!;
while(<RD1>){
        chomp;
	$totalread++;
        my $read1name=$_;
        my $seq1=<RD1>;
        chomp $seq1;
        my $sym1=<RD1>;
        chomp $sym1;
        my $qual1=<RD1>;
        chomp $qual1;
        my $location1=0;
        my $seed1=0;
        my $first1=0;
	my $n1=$seq1=~s/N/N/g;
	if($n1>10){
                #print "$read1name\_$read2name\tNtomany\n";
                $cleanread++;
		next;
        }

        my @qu1=split //,$qual1;
        my $lowbase1;
        foreach my $baseq(@qu1){
                if(ord($baseq)<$quli_limit){
                        $lowbase1++;
                }
        }
        my $lengthread=length($seq1);
	$totalbase+=$lengthread;	
        my $lowrate1=$lowbase1/$lengthread;
	if($lowrate1>$lowqrate){
                #print "$read1name\_$read2name\tLowqtoomany\n";
                $cleanread++;
		next;
        }
	my $end=length($seq1)-$seed;
	my @pos1;
	for(my $i=0;$i<$end;$i++){
                my $seedseq1=substr($seq1,$i,$seed);
                if(exists($hash{$seedseq1})){
                        $location1=$i;
                        $seed1++;
                        push @pos1,$i;
                }
        }
        my $adplen1=$pos1[-1]-$pos1[0]+1+$seed;
	my $startlen1=$pos1[0]+1;
        my $endlen1=$lengthread-$pos1[-1]-$seed;
        my $endpos1=$pos1[-1]+$seed;
        my $startbase1=substr($seq1,0,$startlen1);
        my $endbase1=substr($seq1,$endpos1,$endlen1);
	
	my $startqual1=substr($qual1,0,$startlen1);
        my $endqual1=substr($qual1,$endpos1,$endlen1);

        if($location1!=0 && $pos1[0]<31){
		$cleanread++;
                next;
        }
	
	if($location1==0){
                my $read1="$read1name\n$seq1\n$sym1\n$qual1\n";
                #$clean1gz->gzwrite($read1) || die $!;
		print CLN1 $read1 || die $!;
		$cleanbase+=length($seq1);
        }elsif($startlen1>=31){
                my $read1="$read1name\n$startbase1\n$sym1\n$startqual1\n";
                #$clean1gz->gzwrite($read1) || die $!;
		print CLN1 $read1 || die $!;
		$cleanbase+=length($startbase1);
        }elsif($startlen1<31){
                my $read1="$read1name\n$endbase1\n$sym1\n$endqual1\n";
                #$clean1gz->gzwrite($read1) || die $!;
		print CLN1 $read1 || die $!;
		$cleanbase+=length($startbase1);
        }
}
close CLN1;
close RD1;
$cleanread=$totalread-$cleanread;
print "Totalread\tTotalbase\tCleanread\tCleanbase\n$totalread\t$totalbase\t$cleanread\t$cleanbase\n";
}






