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
die "Usage: $0 <normal.cout> <cancer.cout> <output>\n" if @ARGV != 3;
open NC, ( $ARGV[0] =~ /.gz$/ ? "gzip -cd $ARGV[0] |" : $ARGV[0] ) or die $!;
open TC, ( $ARGV[1] =~ /.gz$/ ? "gzip -cd $ARGV[1] |" : $ARGV[1] ) or die $!;
open OP,    ">$ARGV[2]"         or die $!;
#open THROW, ">./throw_pos.list" or die $!;

# each nucleosome as a methylation unit, size = 250bp #
my $nucleosomeSize   = 250;
# significance of difference threshold #
my $pThreshold       = 0.05;
my $countCpgNum      = 0;

my ( @nCpgInfo, @tCpgInfo );
my ( $blockStart, $blockEnd ) = ( 0, 0 );
my @arrRes;
my $outFlag  = 1;
my $forsaken = "";

while ( my $line = <TC> )
{
    if ( ( $outFlag % 100000 ) == 0 )
    {
        print STDERR "-- Processed $outFlag lines --\n";
    }
    ++$outFlag;
    my @t = split /\s+/, $line;
    $line = <NC>;
    my @n = split /\s+/, $line;

    ($t[3] ne "CG") and next;

    push @tCpgInfo, "$t[1],$t[6],$t[7]";
    push @nCpgInfo, "$n[1],$n[6],$n[7]";

    $line = <TC>;
    my @t2 = split /\s+/, $line;
    $line = <NC>;
    my @n2 = split /\s+/, $line;

    push @tCpgInfo, "$t2[1],$t2[6],$t2[7]";
    push @nCpgInfo, "$n2[1],$n2[6],$n2[7]";

    ++$countCpgNum;

    if ( $countCpgNum >= 5 )
    {
        my $maxPosNum = 2 * $countCpgNum;
        my ( $tmn, $tunmn, $nmn, $nunmn ) = ( 0, 0, 0, 0 );
        for ( my $index = 2 * $countCpgNum - 10 ; $index != 2 * $countCpgNum ; ++$index )
        {
            my ( $pos, $mn, $unmn ) = split /,/, $tCpgInfo[$index];
            if ( $index >= 2 )
            {
# calculate the distance between current pos and previous pos #
                my ($post) = ( split /,/, $tCpgInfo[ $index - 2 ] )[0];
                if ( $pos - $post > $nucleosomeSize )
                {
                    my $previousIdx = $index - 1;
                    if ( $previousIdx >= 10 )
                    {
                        $maxPosNum = $previousIdx + 1;
                        last;
                    }
                    else
                    {
                        for (my $index = 0; $index <= $previousIdx; ++$index)
                        {
                            shift @tCpgInfo;
                            shift @nCpgInfo;
                        }
                        $countCpgNum = @tCpgInfo / 2;
                        if ( ( @tCpgInfo % 2 ) != 0 )
                        {
                            print "--error, check. Odd number postions number.-- @ line " . __LINE__ and die;
                        }
                        goto LABLE_END;
                    }
                }
            }
            $tmn   += $mn;
            $tunmn += $unmn;
        }

        for ( my $index = 2 * $countCpgNum - 10 ; $index != $maxPosNum ; ++$index )
        {
            my ( $pos, $mn, $unmn ) = split /,/, $nCpgInfo[$index];
            $nmn   += $mn;
            $nunmn += $unmn;
        }

		$blockStart = ( split /,/, $tCpgInfo[0] )[0];
        $blockEnd   = ( split /,/, $tCpgInfo[$maxPosNum - 1] )[0];
        my $pValue = sprintf( "%e", &chisqTest( $tmn, $tunmn, $nmn, $nunmn ) );
		
		my ( $tumorMethyRatio, $normalMethyRatio ) = ( 0, 0 );
        $tumorMethyRatio  = sprintf( "%.8f", $tmn / ( $tmn + $tunmn ) ) if ( ( $tmn + $tunmn ) != 0 );
        $normalMethyRatio = sprintf( "%.8f", $nmn / ( $nmn + $nunmn ) ) if ( ( $nmn + $nunmn ) != 0 );
        $forsaken = "$t[0]\t$blockStart\t$blockEnd\t---/---/---/---\t";
        $forsaken .= "$tmn/$tunmn\t$tumorMethyRatio\t$nmn/$nunmn\t$normalMethyRatio\t$pValue\t----\n";

		
		print @tCpgInfo . "-- @ line ". __LINE__ ."\n" if ($debug);
        print STDERR __LINE__ . "----\n" if ($nodebug);
        if ( $pValue > $pThreshold)
        {
            my %hashCountTmp;
            my ( @arrTumor, @arrNormal );
            my ( $chiTm, $chiTum, $chiNm, $chiNum ) = ( 0, 0, 0, 0 );
            my ( $methyInTumor, $methyInNormal, $shareNum ) = ( 0, 0, 0 );
			print @tCpgInfo . "-- @ line ". __LINE__ ."\n" if ($debug);
            for ( my $index = 0 ; $index <= $maxPosNum - 2 ; $index += 2 )
            {
                my ( $pos,  $mn,  $unmn )  = split /,/, $tCpgInfo[$index];
                my ( $pos2, $mn2, $unmn2 ) = split /,/, $tCpgInfo[ $index + 1 ];

                if ( $index >= 2 )
                {
# calculate the distance between current pos and previous pos #
                    my ($post) = ( split /,/, $tCpgInfo[ $index - 2 ] )[0];
                    if ( $pos - $post > $nucleosomeSize )
                    {
                        my $previousIdx = $index - 2;
                        if ( $previousIdx >= 10 )
                        {
                            $maxPosNum = $previousIdx + 2;
                            last;
                        }
                        else
                        {							
                            for (my $index = 0; $index <= $previousIdx + 1; ++$index)
                            {
                                shift @tCpgInfo;
                                shift @nCpgInfo;
                            }							
                            $countCpgNum = @tCpgInfo / 2;
                            if ( ( @tCpgInfo % 2 ) != 0 )
                            {
                                print "--error, check. Odd number postions number.-- @ line " . __LINE__ and die;
                            }
                            goto LABLE_END;
                        }
                    }
                }

                $hashCountTmp {$pos} = 1;
                $chiTm  += ( $mn + $mn2 );
                $chiTum += ( $unmn + $unmn2 );
                ++$methyInTumor if ( ( $mn + $unmn + $mn2 + $unmn2 ) != 0 );
                push @arrTumor, ($mn + $mn2)/( $mn + $unmn + $mn2 + $unmn2 );
            }
			print @tCpgInfo . "-- @ line ". __LINE__ ."\n" if ($debug);
            for ( my $index = 0 ; $index <= $maxPosNum - 2 ; $index += 2 )
            {
                my ( $pos,  $mn,  $unmn )  = split /,/, $nCpgInfo[$index];
                my ( $pos2, $mn2, $unmn2 ) = split /,/, $nCpgInfo[ $index + 1 ];
                $chiNm  += ( $mn + $mn2 );
                $chiNum += ( $unmn + $unmn2 );
                ++$shareNum if ( exists( $hashCountTmp {$pos} ) );
                ++$methyInNormal if ( ( $mn + $unmn + $mn2 + $unmn2 ) != 0 );
                push @arrNormal, ($mn + $mn2)/( $mn + $unmn + $mn2 + $unmn2 );
            }

            my $pValue2 = sprintf( "%e", &chisqTest( $chiTm, $chiTum, $chiNm, $chiNum ) );
            my $tTestPValue = "----";

            if ( @arrTumor >= 2 )
            {
                my $tTest = new Statistics::TTest;
                $tTest->set_significance(99);
                $tTest->load_data( \@arrTumor, \@arrNormal );
                $tTestPValue = sprintf( "%e", $tTest-> {t_prob} );
            }
			my $cutNum = $maxPosNum;
            %hashCountTmp = ();
            $maxPosNum /= 2;
			
			$blockStart = ( split /,/, $tCpgInfo[0] )[0];
            $blockEnd   = ( split /,/, $tCpgInfo[$cutNum - 1] )[0];
			
            my ( $tumorMethyRatio, $normalMethyRatio ) = ( 0, 0 );
            $tumorMethyRatio  = sprintf( "%.8f", $chiTm / ( $chiTm + $chiTum ) ) if ( ( $chiTm + $chiTum ) != 0 );
            $normalMethyRatio = sprintf( "%.8f", $chiNm / ( $chiNm + $chiNum ) ) if ( ( $chiNm + $chiNum ) != 0 );
            $forsaken = "$t[0]\t$blockStart\t$blockEnd\t$shareNum/$methyInNormal/$methyInTumor/$maxPosNum\t";
            $forsaken .= "$chiTm/$chiTum\t$tumorMethyRatio\t$chiNm/$chiNum\t$normalMethyRatio\t$pValue2\t$tTestPValue\n";

            if ( $countCpgNum == 5 )
            {
                shift @tCpgInfo;
                shift @tCpgInfo;
                shift @nCpgInfo;
                shift @nCpgInfo;
                $countCpgNum = 4;
                goto LABLE_END;
            }

            if ($pValue2 <= $pThreshold && ($tumorMethyRatio >= 2 * $normalMethyRatio || 2 * $tumorMethyRatio <= $normalMethyRatio))
            #if ($pValue2 <= $pThreshold)
			#if (($tumorMethyRatio >= 2 * $normalMethyRatio || 2 * $tumorMethyRatio <= $normalMethyRatio))
			{
                push @arrRes, $forsaken;
            }
			print @tCpgInfo . "-- @ line ". __LINE__ ."\n" if ($debug);
            for (my $index = 0; $index != $cutNum; ++$index)
            {
                shift @tCpgInfo;
                shift @nCpgInfo;
            }	
			print @tCpgInfo . "-- @ line ". __LINE__ ."\n" if ($debug);
        }
LABLE_END:
        $countCpgNum = @tCpgInfo / 2;
		print @tCpgInfo . "-- @ line ". __LINE__ ."\n" if ($debug);
    }
}

if ( ( @arrRes != 0 and $forsaken ne $arrRes[$#arrRes] ) ) {
    push @arrRes, $forsaken;
}

print OP "#chr\tpos_start\tpos_end\tshared_num/methy_in_normal/methy_in_tumor/pos_num\ttumor_methy_num/tumor_NONE_methy_num\ttumor_methy_ratio\tnormal_methy_num/normal_NONE_methy_num\tnormal_methy_ratio\tchisq.test_p-value\tt.test_p-value\tannoation\n";

foreach (@arrRes)
{
    print OP $_;
}

close NC;
close TC;
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
