#perl /nas/RD_12A/gaoshengjie/Database/Human/RRBS/Pipeline/bin/bowtie/ASM_cat.pl /nas/RD_12A/gaoshengjie/Database/Human/RRBS/Pipeline/bin/newasm_stimu//8x/Normal/ASM/C2T.emycut.mpileup.vcf /nas/RD_12A/gaoshengjie/Database/Human/RRBS/Pipeline/bin/newasm_stimu//8x/Normal/ASM/G2A.emycut.mpileup.vcf Uniq.total.emycut.mpileup.vcf
use strict;
die "perl $0 ref c2t.vcf g2a.vcf cat.vcf" unless(@ARGV==4);
my $reference=shift;
my $c2t=shift;
my $g2a=shift;
my $out=shift;
my %Ref;
open CT, $c2t or die "no c2t.vcf\n";
open GA, $g2a or die "no g2a.vcf\n";
open OUT,">$out" or die "can't create outfile\n";

&Getref($reference);
my %hash;my %hash2;
while(<CT>){
	chomp;
	my @a=split;
	$a[3]=substr($Ref{$a[0]},$a[1]-1,1);
	if(/\;MQ\=(\d+)\;/){
		my $mq=$1;
		if($a[3] eq 'G' && $a[4] =~/A/){
		 if($mq>=37 && $a[5]>100 ){
			#print OUT $_."\n";
			$hash{$a[0]}{$a[1]}=join("\t",@a);
		  }
		}
		elsif($1>30 && $a[5]>20){
			$hash{$a[0]}{$a[1]}=join("\t",@a);
		}
	}
}
close CT;
#chr1    778855  .       T       A,X     77      .       DP=8;VDB=0.0356;AF1=0.5;CI95=0.5,0.5;DP4=2,2,2,2

while(<GA>){
        chomp;
        my @a=split;
	$a[3]=substr($Ref{$a[0]},$a[1]-1,1);
	my @add=split /\;/,$a[7];
        if(/\;MQ\=(\d+)\;/){
                my $mq=$1;
                if($a[3] eq 'C' && $a[4] =~/T/){
                 if($mq>=37 && $a[5]>100 ){
			if(exists($hash{$a[0]}{$a[1]})){
				delete $hash{$a[0]}{$a[1]};
			}else{
				$hash{$a[0]}{$a[1]}=join("\t",@a);
			}
		 }
                }
                elsif($mq>30 && $a[5]>20){
		  if(exists($hash{$a[0]}{$a[1]})){
			#delete $hash{$a[0]}{$a[1]};
                        my @b=split /\t/,$hash{$a[0]}{$a[1]};
			my @change=split /\;/,$b[7];
			my ($ctdp,$gadp,$ctdp1,$ctdp2,$ctdp3,$ctdp4,$gadp1,$gadp2,$gadp3,$gadp4);			
			if($a[7]=~/DP\=(\d+)\;/){
				$gadp=$1;	
			}
			if($a[7]=~/DP4\=(\d+)\,(\d+)\,(\d+)\,(\d+)/){
				$gadp1=$1;$gadp2=$2;$gadp3=$3;$gadp4=$4;
			}
			if($b[7]=~/DP\=(\d+)\;/){
				$ctdp=$1;
                        }
                        if($b[7]=~/DP4\=(\d+)\,(\d+)\,(\d+)\,(\d+)/){
				$ctdp1=$1;$ctdp2=$2;$ctdp3=$3;$ctdp4=$4;
                        }
	
			my $dp=$ctdp + $gadp;
			my $dp1=$ctdp1 + $gadp1;
			my $dp2=$ctdp2 + $gadp2;
			my $dp3=$ctdp3 + $gadp3;
			my $dp4=$ctdp4 + $gadp4;
			$change[0]="DP\=$dp";
			$change[4]="DP4\=$dp1\,$dp2\,$dp3\,$dp4";		
			$b[7]=join("\;",@change);
			$a[7]=join("\;",@change);
            		my $depthab=$a[5]+$b[5];
			$a[5]=$depthab;
			$b[5]=$depthab;
			if($a[4] eq $b[4]){
			  #$hash{$a[0]}{$a[1]}=join("\t",@b);
			  if($a[3] eq 'A'){
				$hash{$a[0]}{$a[1]}=join("\t",@b);
			  }elsif($b[3] eq 'T'){
				$hash{$a[0]}{$a[1]}=join("\t",@a);	
			  }else{
				$hash{$a[0]}{$a[1]}=join("\t",@b);	
			  }
			}else{
				delete ($hash{$a[0]}{$a[1]});
				#print STDERR "warning tribase\t$a[3]\t$b[3]\t$a[4]\t$b[4]\n$hash{$a[0]}{$a[1]}\n";
			}
#=cut		
		  }else{
			$hash{$a[0]}{$a[1]}=$_;
		  }
                }
        }
}
close GA;

foreach my $chr(sort keys %hash){
	foreach my $pos(sort{$a<=>$b} keys %{$hash{$chr}}){
		print OUT $hash{$chr}{$pos}."\n";
=head
		my @info=split /\t/,$hash{$chr}{$pos};
		if($info[3] eq 'C' && $info[4] =~/T/){
			$hash2{$chr}{$pos}=$hash{$chr}{$pos};	
		}elsif($info[3] eq 'G' && $info[4] =~/A/){
			$hash2{$chr}{$pos}=$hash{$chr}{$pos};
		}
=cut
	}

}

sub Getref
{
        my $ref=shift;
        open RF,$ref or die "need reference\n";
        my $chr;
        while(<RF>){
                chop;
                if(/^>(chr\w+)/){
                        $chr=$1;
                        print STDERR "$chr\n";
                }else{
                        $_=~s/\s+//g;
                        $Ref{$chr}.=$_;
                }
        }
        close RF;
}

