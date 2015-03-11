use strict;
use FindBin '$Bin';
use lib "$Bin/../lib";
use Math::CDF qw/pchisq/;
use Text::NSP::Measures::2D::Fisher::twotailed qw/calculateStatistic/;
#use lib '/ifs5/PC_HUMAN_US/USER/zhouquan/work/dmr/Statistics-TTest-1.1.0/lib/perl5/site_perl/5.8.8';
use Statistics::PointEstimation;
use Statistics::TTest;
my $debug = 0;
my $nodebug = 0;
die "Usage: $0 <normal.cout> <cancer.cout> <targetregion>  <cover_threshold> <output>\n" if @ARGV != 5;
open NC, ( $ARGV[0] =~ /.gz$/ ? "gzip -cd $ARGV[0] |" : $ARGV[0] ) or die $!;
open TC, ( $ARGV[1] =~ /.gz$/ ? "gzip -cd $ARGV[1] |" : $ARGV[1] ) or die $!;
open TR, ( $ARGV[2] =~ /.gz$/ ? "gzip -cd $ARGV[1] |" : $ARGV[2] ) or die $!;
open OP,    ">$ARGV[4]"         or die $!;
my $thre_cover=$ARGV[3];
my %hash;
#chr21   9437444 +       CG      CG      1       28      57      3
my $chr;
while(<NC>){
	chomp;
	my  @norm=split;
	my $tc=<TC>;
	chomp $tc;
	my @case=split /\t/,$tc;
	die unless($case[1] eq $norm[1]);
	next unless($norm[5]==1);
	if($norm[3] eq 'CG'){
		$chr=$norm[0];
		my $total_norm=$norm[6]+$norm[7];
		my $total_case=$case[6]+$case[7];
		my ($rate_norm,$rate_case);
		if($total_norm==0){
			$rate_norm='#';
		}else{
			$rate_norm=$norm[6]/$total_norm;
		}
		if($total_case==0){
			$rate_case='#';	
		}else{
			$rate_case=$case[6]/$total_case;
		}
		$hash{$norm[0]}{$norm[1]}="$total_norm\t$norm[6]\t$total_case\t$case[6]";
	}
	#$chr=$a[0];		
}

close NC;
close TC;

while(<TR>){
        chomp;
        my @a=split;
	next if($a[0] ne $chr);
	my $total_CpG;
	my $norm_cover;
	my $case_cover;
	my $overlap_cover;
	my $cpvalue;
	my $tpvalue;
	my $methy_normal;
	my $methy_case;
	my $total_normal=0;
	my $total_case=0;
	my @arrnorm;
	my @arrtumor;
	foreach my $pos($a[1]..$a[2]){
		if(exists($hash{$chr}{$pos})){
			$total_CpG++;
			my @info=split /\t/,$hash{$chr}{$pos};
			if($info[0]>=$thre_cover){
				$norm_cover++;		
			}
			if($info[2]>=$thre_cover){
                                $case_cover++;
                        }
			if($info[0]>=$thre_cover && $info[2]>=$thre_cover){
                                $overlap_cover++;
				$total_normal+=$info[0];
				$total_case+=$info[2];
				$methy_normal+=$info[1];
				$methy_case+=$info[3];
				push @arrnorm, ($info[1]/$info[0]);	
				push @arrtumor,($info[3]/$info[2]);	
                        }		
		}	

	}
	next if($total_normal==0);
	next if($total_case==0);
	my $normrate= $methy_normal/$total_normal;
	my $caserate = $methy_case/$total_case;
	my $cpvalue = sprintf( "%e", &chisqTest( $methy_case, $total_case-$methy_case, $methy_normal, $total_normal-$methy_normal ) );	
        #my $region="$a[1]\-$a[2]";
        #$hash{$a[0]}{$a[1]}=$a[3];
	my $tTestPValue;
	if ( @arrtumor >= 2 ){ 
            {   
                my $tTest = new Statistics::TTest;
                $tTest->set_significance(99);
                $tTest->load_data( \@arrtumor, \@arrnorm );
                $tpvalue = sprintf( "%e", $tTest-> {t_prob} );
            }   
	}
	print OP join("\t",@a[0..3])."\t"."$overlap_cover\/$norm_cover\/$case_cover\/$total_CpG\t$caserate\t$normrate\t$cpvalue\t$tpvalue\n";
	
}
close TR;

close OP;

sub chisqTest()
{
    my ( $a, $b, $c, $d ) = @_;
    if ( $a == 0 && $b == 0 && $c == 0 && $d == 0 )
    {
        return 1;
    }
    my $tt  = $a + $b + $c + $d;
    my $np1 = $a + $c;
    my $np2 = $b + $d;
    my $n1p = $a + $b;
    my $n2p = $c + $d;
    if ( $np1 < 5 || $np2 < 5 || $n1p < 5 || $n2p < 5 || $tt < 40 )
    {
        return &Text::NSP::Measures::2D::Fisher::twotailed::calculateStatistic(
                   n11 => $a,
                   n1p => $n1p,
                   np1 => $np1,
                   npp => $tt
               );
    }
    my $squareValue = ( abs( $a * $d - $b * $c ) - $tt / 2 )**2 * $tt / ( $np1 * $np2 * $n1p * $n2p );
    return 1 - &Math::CDF::pchisq( $squareValue, 1 );
}
