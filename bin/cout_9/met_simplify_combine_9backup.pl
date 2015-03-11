use strict;
use FindBin '$Bin';
use File::Path;  ## function " mkpath" and "rmtree" deal with directory
use File::Basename qw(basename dirname);
use lib qw($Bin/lib);
use PerlIO::gzip;
use Math::CDF qw/pbinom/;
die "perl $0 qmap errorate count\n" unless(@ARGV==3);
my $qmap=shift;
my $rate=shift;
my $cout=shift;
#my $rate=shift;
if($qmap=~/\.gz$/){
	open IN,"gunzip -c $qmap|" or die $!;
}else{
	open IN,$qmap or die "cant open $qmap\n";
}
if($cout=~/\.gz$/){
	open OUT,">:gzip","$cout" or die $!; 
}else{
	open OUT,">$cout" or die $!;
}


#A       0       0       6       0       0       0       0       0       0       0       eeeb_b
#chrM    35      -       CHG     CCG     1       1       200     8
#chrM    31      +       CHH     CAC     1       1       94      6
#chrM    33      +       CG      CGG     1       0       97      98
#chrM    34      -       CG      CGT     1       2       197     8

my $chr=<IN>;
chomp $chr;
$chr=~s/^\>//;
my $i=0;
my  ($start,$middle,$end);
$start=<IN>;
chomp $start;
$middle=<IN>;
chomp $middle;

while(<IN>){
	$i++;
	my @start=split /\t/,$start;
	my @middle=split /\t/,$middle;
	chomp $_;
	$end=$_;
	my @end=split /\t/,$_;

	if($start[0]=~/[c|C]/){
		my $copyall=$start[1]+$start[2]+$start[3]+$start[4]+$start[9];
		my $all=$start[1]+$start[2]+$start[3]+$start[4]+$start[5]+$start[7]+$start[6]+$start[8];
		my $copynum;
		my $n=$start[1]+$start[2]+$start[3];
		my $unmethy1=$start[2]+$start[3];
		my $thre;
		if($start[2]==0){
			$thre=0.01*$start[1];
		}else{		
			$thre=0.01*$start[1]/$start[2];
		}
		my $bin=0;
		
		if($all==0){
			$copynum=0;
			$bin=0;
		}else{
			$copynum=$copyall/$all;
			my $smallest=1-pbinom($n-1,$n,$rate);
			if($smallest>=$thre){
				$bin=$n+1;
			}else{
			for my $k(1..$n){
				my $tmp=pbinom($k,$n,$rate) - pbinom($k-1,$n,$rate);
				if($tmp<$thre){
					$bin=$k;
					last;
				}		
			}
			}	
			
		}	
		if($middle[0] =~/[G|g]/){
			print OUT "$chr\t$i\t\+\tCG\t$start[0]$middle[0]$end[0]\t$copynum\t$start[1]\t$unmethy1\t$bin\n";
		}elsif($end[0]=~/[G|g]/){
			print OUT"$chr\t$i\t\+\tCHG\t$start[0]$middle[0]$end[0]\t$copynum\t$start[1]\t$unmethy1\t$bin\n";
		}else{
			print OUT "$chr\t$i\t\+\tCHH\t$start[0]$middle[0]$end[0]\t$copynum\t$start[1]\t$unmethy1\t$bin\n";		
		}
		$start=$middle;
		$middle=$end;
		
	}

	if($end[0]=~/[g|G]/){
		my $j=$i+2;
		my $copyall=$end[1]+$end[2]+$end[3]+$end[4]+$end[9];
                my $all=$end[1]+$end[2]+$end[3]+$end[4]+$end[5]+$end[7]+$end[6]+$end[8];
                my $copynum;
		my $n=$end[1]+$end[2]+$end[3];
		my $unmethy2=$end[2]+$end[3];
		my $thre;
		if($end[2]==0){
			$thre=0.01*$end[1];
		}else{
			$thre=0.01*$end[1]/$end[2];
		}
		my $bin=0;

                if($all==0){
                        $copynum=0;
                }else{
                        $copynum=$copyall/$all;
			my $smallest=1-pbinom($n-1,$n,$rate);
                        if($smallest>=$thre){
                                $bin=$n+1;
                        }else{
                        for my $k(1..$n){
                                my $tmp=pbinom($k,$n,$rate) - pbinom($k-1,$n,$rate);
                                if($tmp<$thre){
                                        $bin=$k;
                                        last;
                                }
                        }
                        }			
                }

		
		if($middle[0] =~/[c|C]/){
			print OUT "$chr\t$j\t\-\tCG\t$end[0]$middle[0]$start[0]\t$copynum\t$end[1]\t$unmethy2\t$bin\n";
		}elsif($start[0]=~/[c|C]/){
			print OUT "$chr\t$j\t\-\tCHG\t$end[0]$middle[0]$start[0]\t$copynum\t$end[1]\t$unmethy2\t$bin\n";
		}else{
			print OUT "$chr\t$j\t\-\tCHH\t$end[0]$middle[0]$start[0]\t$copynum\t$end[1]\t$unmethy2\t$bin\n";
		}	
		$start=$middle;
		$middle=$start;		
	}
	#else{
	#	$start=$middle;$middle=$end;
	#}

	unless($start[0]=~/[c|C]/ || $end[0]=~/[[G|g]]/){
		$start=$middle;
		$middle=$end;	
	}
}
close IN;
close OUT;
	
