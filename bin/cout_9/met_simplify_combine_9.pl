use strict;
use FindBin '$Bin';
use File::Path;  ## function " mkpath" and "rmtree" deal with directory
use File::Basename qw(basename dirname);
use lib qw($Bin/lib);
use PerlIO::gzip;
use Math::CDF qw/pbinom/;
die "perl $0 qmap ref errorate count\n" unless(@ARGV==4);
my $qmap=shift;
my $ref=shift;
my $rate=shift;
my $cout=shift;
#my $rate=shift;
if($qmap=~/\.gz$/){
	open IN,"gunzip -c $qmap|" or die $!;
}else{
	open IN,$qmap or die "cant open $qmap\n";
}
if($ref =~/\.gz/){
        open REF,"gunzip -c $ref|" or die $!;
}else{
        open REF,$ref or die $!;
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
my %hash;
my $seq;
while(<REF>){
        chomp;
        if(/\>(chr\S+)/){
                if($chr eq $1){
                        $/=">";
                         $seq=<REF>;
                        chomp $seq;
                        #print $seq."\n";
                        $seq=~s/\s//g;
                        $/="\n";
                        print STDERR "$chr read done!\n";
                        last;

                }
        }
}
close REF;


my $i=0;
my  ($start,$middle,$end,$startbase,$middlebase,$endbase);
$start=<IN>;
chomp $start;
$middle=<IN>;
chomp $middle;
$startbase=substr($seq,$i,1);
$middlebase=substr($seq,$i+1,1);

while(<IN>){
	$i++;
	my @start=split /\t/,$start;
	my @middle=split /\t/,$middle;
	chomp $_;
	$end=$_;
	my @end=split /\t/,$_;
	$endbase=substr($seq,$i+1,1);

	if($startbase=~/[c|C]/){
		my $copyall=$start[1]+$start[2]+$start[3]+$start[4]+$start[9];
		my $copyall2=$middle[1]+$middle[2]+$middle[3]+$middle[4]+$middle[9];
		my $all=$start[1]+$start[2]+$start[3]+$start[4]+$start[5]+$start[7]+$start[6]+$start[8];
		my $all2=$middle[1]+$middle[2]+$middle[3]+$middle[4]+$middle[5]+$middle[7]+$middle[6]+$middle[8];
		my $copynum;my $copynum2;
		my $n=$start[1]+$start[2]+$start[3];
		my $n2=$middle[1]+$middle[2]+$middle[3];
		my $unmethy1=$start[2]+$start[3];
		my $unmethy2=$middle[2]+$middle[3];
		my $thre;my $thre2;
		if($start[2]==0){
			$thre=0.01*$start[1];
		}else{		
			$thre=0.01*$start[1]/$start[2];
		}
		if($middle[2]==0){
                        $thre2=0.01*$middle[1];
                }else{
                        $thre2=0.01*$middle[1]/$middle[2];
                }
	
		my $bin=0;my $bin2=0;
		
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
		if($all2==0){
                        $copynum2=0;
                        $bin2=0;
                }else{
                        $copynum2=$copyall2/$all2;
                        my $smallest=1-pbinom($n2-1,$n2,$rate);
                        if($smallest>=$thre2){
                                $bin2=$n+1;
                        }else{
                        for my $k(1..$n){
                                my $tmp=pbinom($k,$n,$rate) - pbinom($k-1,$n,$rate);
                                if($tmp<$thre2){
                                        $bin2=$k;
                                        last;
                                }
                        }
                        }

                }
		
		my $mpos=$i+1;
		if($middlebase =~/[G|g]/){
			print OUT "$chr\t$i\t\+\tCG\t$startbase$middlebase\t$copynum\t$start[1]\t$unmethy1\t$bin\n";
			print OUT "$chr\t$mpos\t\-\tCG\t$middlebase$startbase\t$copynum2\t$middle[1]\t$unmethy2\t$bin2\n";
		}elsif($endbase=~/[G|g]/){
			print OUT"$chr\t$i\t\+\tCHG\t$startbase$middlebase$endbase\t$copynum\t$start[1]\t$unmethy1\t$bin\n";
		}else{
			print OUT "$chr\t$i\t\+\tCHH\t$startbase$middlebase$endbase\t$copynum\t$start[1]\t$unmethy1\t$bin\n";		
		}
		
	}

	if($endbase=~/[g|G]/){
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

		
		if($middlebase =~/[c|C]/){			
			#print OUT "$chr\t$j\t\-\tCG\t$endbase$middlebase$startbase\t$copynum\t$end[1]\t$unmethy2\t$bin\n";
		}elsif($startbase=~/[c|C]/){
			print OUT "$chr\t$j\t\-\tCHG\t$endbase$middlebase$startbase\t$copynum\t$end[1]\t$unmethy2\t$bin\n";
		}else{
			print OUT "$chr\t$j\t\-\tCHH\t$endbase$middlebase$startbase\t$copynum\t$end[1]\t$unmethy2\t$bin\n";
		}	
	}

	$start=$middle;
	$middle=$end;	
	$startbase=$middlebase;
	$middlebase=$endbase;
}
close IN;
close OUT;
	
