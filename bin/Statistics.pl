#!/usr/bin/perl -w
use strict;
use PerlIO::gzip;
use File::Basename;
use Getopt::Long;
=head1 Name
	 statistic.pl 
=head1 Usage
	perl statistic.pl.pl  hs.cut.emy  posOrder  enzyme_output --out2 notinenzyme_output --fqlen  fq_readlength[default=posOrder:len(*)] --insert [38:230]
=cut

die `pod2text $0` if (@ARGV<3);
my $emycut=shift;
my $infile=shift;
my $outfile=shift;

my ($fqlen,$insert,$out2,$Help);

GetOptions(
        "fqlen:i"=>\$fqlen,
        "insert:s"=>\$insert,
	#"enzyme"=>\$enzyme,
	"out2:s"=>\$out2,
        "help"=>\$Help
);

die `perl statistic.pl.pl  hs.cut.emy  posOrder  enzyme_output  --out2 notinenzyme --fqlen  fq_readlength[default=posOrder:len(*)] --insert [38:230]` if ($Help);

my $length= $fqlen || -1;
$insert= $insert || "40:300";
my ($short,$long)=split (/\:/,$insert);

if($emycut=~/gz$/){
	open EMC,"<:gzip",$emycut  or die "$!" ;
}else{
	open EMC,$emycut  or die "$!";
}


my %hash; my %region;
my $silo=0;
while(<EMC>){
	chomp;
	my @a=split;
	$hash{$a[0]}{$a[1]}="$a[3]\tstart";
	$hash{$a[0]}{$a[2]}="$a[3]\tend";
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

if($outfile=~/gz$/){
	open OUT,">:gzip",$outfile  || die "$!" ;
}else{
	open OUT,">$outfile" || die "$!";
}

if($out2){
if($out2=~/gz$/){
        open OUT2,">:gzip",$out2  || die "$!" ;
}else{
        open OUT2,">$out2" || die "$!";
}
}

my $outdir = dirname $outfile;
my $out6="$outdir/enzymeNUMs.out";
my $pairtwo=0;
my $pairone=0;
my $pairzero=0;
my $singleone=0;
my $singlezero=0;
#open PT,">$out1";open PO,">$out2";open PZ,">$out3";open SO,">$out4";open SZ,">$out5";open NUM,">$out6";
open NUM,">$out6";

while(<IN>){
	chomp;
	my @a=split /\t/;
	if(@a>11){
#		&pair($_,$length,\*PT,\*PO,\*PZ);
		&pair($_,$length,\*OUT);
	}
	elsif(@a<12){
#		&single($_,$length,\*SO,\*SZ);
		&single($_,$length,\*OUT);
	}
}
close IN;

my %emy;
my %emypos;
sub pair{
#	my ($line,$len,$o1,$o2,$o3)=@_;
	my ($line,$len,$o1)=@_;
	my $score=0;
	my @a=split /\t/,$line;

	if($len==-1){
		if($a[4]=~/\:len(\d+)\:/) {$len=$1;}
	}
#chr5W   1513080 1513260 180     FCC01P5ACXX:8:1101:10000:6391#TTAGGCAT/1:rep1:+:len90:W:score42 chr5    1513080 CGGTATTTTGGTGTGGGTTTGCGGGAAAAGGGCGTGGGGTGGAGGTGTAAGTTTAGTTGAGAAATTAGGAGTAGAGAGTGTAGAGTTTTA      m||*||**|||||||||**||m||||||||||m||||||||||||||*|||*|*||*||||||||*||||||||||||||*||||****|      CGGCATCCTGGTGTGGGCCTGCGGGAAAAGGGCGTGGGGTGGAGGTGCAAGCTCAGCTGAGAAATCAGGAGTAGAGAGTGCAGAGCCCCA      CCCDDFFFHHHDFEGIJHIJJJJJJHGJJJJJJJFHIJH5@CABD5=>CDDDEEEDDDDDDDCCDDDDDDD?CDDCDDACACDEDCEDDD      FCC01P5ACXX:8:1101:10000:6391#TTAGGCAT/2:rep1:-:len90:W:score42 chr5    1513171 GGAGGGAGGGTTTAGGTAGGGTGGTTATGTGAGTATTATAGTCGATGTTTTAGTTTTTGTTTTAAGATGATCGATGGTAGGTAGGTGTTG      ||||||||||***|||*|||||||*||||||||*||*|*||*m||||****||||**||*||*|||||||*m|||||*|||*|||||**|      GGAGGGAGGGCCCAGGCAGGGTGGCTATGTGAGCATCACAGCCGATGCCCCAGTTCCTGCTTCAAGATGACCGATGGCAGGCAGGTGCCG      DDDDDDDDCDDEEEEEFFFFHHHHJJJIJJJIJJJJJJHIJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJHHHHHFFFFFCCC

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
#chr16C  1811744 1811894 90      FCC01P5ACXX:8:1101:10000:100563#TTAGGCAT/1:rep1:-:len90:C:score42       chr16   1811805 CGCTAAACGACCACCATTCACCATTCACCATTCACACACAATATATCCTACTATTCACACTCGAAACATAAAAATACGAAATACCTCCCG      |m||###|m#||||||||||||||||||||||||||||||#|||||||||||||||||||||m###|||####||#|m###|#||||||m      CGCTGGGCGGCCACCATTCACCATTCACCATTCACACACAGTATATCCTACTATTCACACTCGGGGCATGGGGATGCGGGGTGCCTCCCG      DDDDDDDFFEHHHHEJIIIHJJJIIHCJJJIHHIIGJIJJIGIJIEIIF?IIJJIIHHGJIJJJGHJJIJJJJIHFJHHHHFAFFFFCBB      FCC01P5ACXX:8:1101:10000:100563#TTAGGCAT/2:rep1:+:len90:C:score42       chr16   1811744 CAAATAACCACTATTCATACGCATTCACACGTAAAACAACCATCATTCCCACACATTCACA                                   |###|#||||||#|||||||m|||||||||m|####||#||||||||||||||||||||||                                   CGGGTGACCACTGTTCATACGCATTCACACGTGGGGCAGCCATCATTCCCACACATTCACACGCTGGGCGGCCACCATTCACCATTCACC      CCCFFFFFHHHHHJJJJJJJJJJJJJJJJJJHIJJJJJJJJJJJJJJJJJHHIJJJJJJJJJIJHGHFFFFDDDDDDDDEEDDDDDFEDD
	elsif($a[4]=~/\-/){
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
		select($o1);
		print "$line\n";
	}
	elsif($score==1){
		$pairone++;
		select($o1);
		print "$line\n";
	#	select($o2);
	}
	elsif($score==0){
		$pairzero++;
		if($out2){
			print OUT2 "$line\n";
		}
	#	select($o3);
		
	}
}

sub single{
#	my ($line,$len,$o1,$o2)=@_;
	my ($line,$len,$o1)=@_;
	my $score=0;
	my @a=split /\t/,$line;
	if($len==-1){
                $a[4]=~/\:len(\d+)\:/; $len=$1;
        }

#	$a[4]=~/\:len(\w+)\:/;my $len=$1;
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
		select($o1);
		print "$line\n";
	}
	elsif($score==0){
		$singlezero++;
		if($out2){
			print OUT2 "$line\n";
		}

	}
#	print "$line\n";
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

