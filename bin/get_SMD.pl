use strict;
use File::Basename qw/basename dirname/;
use Cwd 'abs_path';
use FindBin qw($Bin);
use lib "$Bin/lib";
use Math::CDF qw/pchisq/;
use Text::NSP::Measures::2D::Fisher::twotailed qw/calculateStatistic/;
die "perl $0 normal_cout cancer_cout >SMD.out" unless(@ARGV==2);
print "Chr\tPos\tStrand\tType\tRef\tCopy_number\t#normal_reads_of_methy\t#normal_read_of_unmethy\tnormal_threshhold\t#cancer_reads_of_methy\t#cancer_read_of_unmethy\tcancer_threshhold\tChisq_pvalue\n";
my ($normalfile,$cancerfile)=@ARGV;
if($normalfile=~/\.gz$/){
	open NM,"gunzip -c $normalfile|" or die $!;
}else{
	open NM,$normalfile or die $!;
}
if($cancerfile=~/\.gz$/){
	open CE,"gunzip -c $cancerfile|" or die $!;
}else{
	open CE,$cancerfile or die $!;
}
while(<NM>){
	chomp;
	my $cout_normal=$_;
	my $cout_cancer=<CE>;
	next unless($cout_normal=~/\w/);
	chomp $cout_cancer;
	my @normal=split /\s+/,$cout_normal;
	my @cancer=split /\s+/,$cout_cancer;
	die unless($normal[1] eq $cancer[1]);
	my $pvalue=&chisqTest(
                $normal[6], $normal[7],
                $cancer[6], $cancer[7]
            );
	#print $pvalue."\#";
	#next if($pvalue>0.1);

	print "$cout_normal\t".join("\t",@cancer[6..8])."\t";
	 print sprintf(
            "%e\n",
            &chisqTest(
                $normal[6], $normal[7],
                $cancer[6], $cancer[7]
            )
	)
	
}
close NM;
close CE;

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
 
