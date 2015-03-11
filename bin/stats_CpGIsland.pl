use strict;
use File::Basename qw/basename dirname/;
use Cwd 'abs_path';
use FindBin qw($Bin);
use lib "$Bin/lib";
use Math::CDF qw/pchisq/;
use Text::NSP::Measures::2D::Fisher::twotailed qw/calculateStatistic/;

die "Usage: $0 <normal_cout-dir> <tumor_cout-dir> <CGI.bed>\n" if @ARGV < 1;

my ( $d_ncout, $d_tcout, $f_cgi ) = @ARGV;

# default PE90 #
$f_cgi ||= "/nas/RD_12A/gaoshengjie/Database/Human/RRBS/bin/prepare/hg19/element/CpGIsland.bed.sort";

#check parameters#
my $abs_d_ncout = abs_path($d_ncout);
my @all_ncout   = glob "$abs_d_ncout/chr*.cout.gz";
my $abs_d_tcout = abs_path($d_tcout);
my @all_tcout   = glob "$abs_d_tcout/chr*.cout.gz";

my @all_chr =
  qw/chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM/;

my $flag_count = 1;

# title #
my @title =
  qw/chr start-pos end-pos total-CpG normal-CpG N-methy-CpG N-unmet-CpG N-methy-ratio tumor-CpG T-methy-CpG T-unmet-CpG T-methy-ratio shared-CpG chisq_test-p_value/;
print "#\t", join( "\t", @title ), "\n";

foreach my $cur_chr (@all_chr)
{

    # process CGI file #
    my %hash_CGI;
    my %hash_CGI_map;
    my %hash_CGI_pos;
    open I_CGI, $f_cgi or die $!;
    while (<I_CGI>)
    {
        chomp;
        my ( $chr, $sp, $ep, $seq ) = split;
        next if ( $cur_chr ne $chr );

        $hash_CGI{$flag_count}{"CpG_total"} = [ 0, 0 ];    # total CpG number, [normal, tumor]

        $hash_CGI{$flag_count}{"CpG_rtotal"} = [ 0, 0 ];   # total real CpG number, [normal, tumor]

        $hash_CGI{$flag_count}{"CpG_nmethy"} = [ 0, 0 ];   # [methy CpG number, reads number]
        $hash_CGI{$flag_count}{"CpG_tmethy"} = [ 0, 0 ];   # [methy CpG number, reads number]

        $hash_CGI{$flag_count}{"CpG_nunmethy"} = [ 0, 0 ]; # [un-methy CpG number, reads number]
        $hash_CGI{$flag_count}{"CpG_tunmethy"} = [ 0, 0 ]; # [un-methy CpG number, reads number]
        $hash_CGI{$flag_count}{"CpG_share"}    = [ 0, 0 ]; # non zero and shared CpG number, [theoretical, real]

        #		$hash_CGI{$flag_count}{"CpG_1x"} = 0;		# ...
        #		$hash_CGI{$flag_count}{"CpG_4x"} = 0;
        #		$hash_CGI{$flag_count}{"CpG_10x"} = 0;
        #		$hash_CGI{$flag_count}{"CpG_30x"} = 0;
        #   	$hash_CGI{$flag_count}{"CpG_50x"} = 0;

        $hash_CGI_map{$flag_count} = "$sp\t$ep";

        for ( $sp .. $ep )
        {
            $hash_CGI_pos{$_} = $flag_count;
        }

        ++$flag_count;
    }
    close I_CGI;

    my $f_ncout = "$abs_d_ncout/$cur_chr.cout.gz";
    my $f_tcout = "$abs_d_tcout/$cur_chr.cout.gz";

    die "-- $f_ncout does not exists. --\n" if not -e $f_ncout;
    die "-- $f_tcout does not exists. --\n" if not -e $f_tcout;

    open NCOUT, "<:gzip(autopop)", $f_ncout or die $!;
    my %hash_share;
    while (<NCOUT>)
    {
        chomp;
        my ( $pos, $m, $unm ) = (split)[ 1, 6, 7 ];
        next if ( !exists( $hash_CGI_pos{$pos} ) );

        my $count = $hash_CGI_pos{$pos};
        $hash_share{$pos} = 1 if $m != 0 or $unm != 0;

        #print STDERR "$pos, $m, $unm\t--$count--\n";

        ++$hash_CGI{$count}{"CpG_total"}->[0];
        ++$hash_CGI{$count}{"CpG_rtotal"}->[0] if $m != 0 or $unm != 0;

        ++$hash_CGI{$count}{"CpG_nmethy"}->[0] if $m != 0;
        $hash_CGI{$count}{"CpG_nmethy"}->[1] += $m;

        ++$hash_CGI{$count}{"CpG_nunmethy"}->[0] if $unm != 0;
        $hash_CGI{$count}{"CpG_nunmethy"}->[1] += $unm;
    }
    close NCOUT;

    open TCOUT, "<:gzip(autopop)", $f_tcout or die $!;
    while (<TCOUT>)
    {
        chomp;
        my ( $pos, $m, $unm ) = (split)[ 1, 6, 7 ];
        next if ( !exists( $hash_CGI_pos{$pos} ) );
        my $count = $hash_CGI_pos{$pos};

        ++$hash_CGI{$count}{"CpG_total"}->[1];
        ++$hash_CGI{$count}{"CpG_rtotal"}->[1] if $m != 0 or $unm != 0;

        ++$hash_CGI{$count}{"CpG_tmethy"}->[0] if $m != 0;
        $hash_CGI{$count}{"CpG_tmethy"}->[1] += $m;

        ++$hash_CGI{$count}{"CpG_tunmethy"}->[0] if $unm != 0;
        $hash_CGI{$count}{"CpG_tunmethy"}->[1] += $unm;

        if ( exists( $hash_share{$pos} ) )
        {
            ++$hash_CGI{$count}{"CpG_share"}->[0];
            if ( $m != 0 or $unm != 0 )
            {
                ++$hash_CGI{$count}{"CpG_share"}->[1];
            }
        }
    }
    close TCOUT;

    foreach my $count ( sort { $a <=> $b } keys %hash_CGI )
    {
        my $pair = $hash_CGI_map{$count};

        print "$cur_chr\t$pair\t";
        print $hash_CGI{$count}{"CpG_total"}->[0], "\t";

        print $hash_CGI{$count}{"CpG_rtotal"}->[0],   "\t";
        print $hash_CGI{$count}{"CpG_nmethy"}->[1],   "\t";
        print $hash_CGI{$count}{"CpG_nunmethy"}->[1], "\t";
		if ($hash_CGI{$count}{"CpG_nunmethy"}->[1] + $hash_CGI{$count}{"CpG_nmethy"}->[1] != 0){
			print sprintf("%.3f\t", $hash_CGI{$count}{"CpG_nmethy"}->[1] / ($hash_CGI{$count}{"CpG_nunmethy"}->[1] + $hash_CGI{$count}{"CpG_nmethy"}->[1]));
		}
		else{
			print "0.000\t";
		}

        print $hash_CGI{$count}{"CpG_rtotal"}->[1],   "\t";
        print $hash_CGI{$count}{"CpG_tmethy"}->[1],   "\t";
        print $hash_CGI{$count}{"CpG_tunmethy"}->[1], "\t";
		if ($hash_CGI{$count}{"CpG_tmethy"}->[1] + $hash_CGI{$count}{"CpG_tunmethy"}->[1] != 0){
			print sprintf("%.3f\t", $hash_CGI{$count}{"CpG_tmethy"}->[1] / ($hash_CGI{$count}{"CpG_tmethy"}->[1] + $hash_CGI{$count}{"CpG_tunmethy"}->[1]));
		}
		else{
			print "0.000\t";
		}

        print $hash_CGI{$count}{"CpG_share"}->[1], "\t";

        print sprintf(
            "%e\n",
            &chisqTest(
                $hash_CGI{$count}{"CpG_nmethy"}->[1], $hash_CGI{$count}{"CpG_nunmethy"}->[1],
                $hash_CGI{$count}{"CpG_tmethy"}->[1], $hash_CGI{$count}{"CpG_tunmethy"}->[1]
            )
        );
    }
}

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

