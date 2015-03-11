use strict;
use IO::Handle;
use Cwd;
use Getopt::Long;
use File::Path;  
use File::Basename qw(basename dirname);

use threads;
use threads::shared;
use Parallel::ForkManager;

my $usage =<<"USAGE";
 name:    $0
          this perl pro is written for recognizing the Restrict Enzyme Digestion site on the reference sequence,
          and it will print out the recognized fragments, coordinate, and cutting coordinate
 usage:   perl $0
	       -B    path/Bowtie2
	       -R    Insert region(40-220)
               -I    xxx.fa; the sequences you want to find Digestion Site within.
               -P    Restrict Enzyme Digestion site; Enzyme Digestion site here use the literal IUPAC
                     consensus sequence, The definition of degenerate characters are as follows:
                        M = A or C
                        R = A or G
                        W = A or T
                        S = C or G
                        Y = C or T
                        K = G or T
                        V = not T
                        H = not G
                        D = not C
                        B = not A
                        N = any base
                      Example: the digestion site of SmlI is that: CTYRAG
                               the digestion site of Sau96I is that: GGNCC
                      -----> if you use one more Restrict Enzymes,you must join Restrict Enzyme Digestion fragnments with "-",
                             for example:
                                 CTYRAG-GGNCC
                                 AAAAAA-GGGGG-TTTTTT-CCCCC
               -S     the cutting site related the first base in digestion framnment, default is 0,
                      for example:  1,  |CCCCC  is 0
                                    2,  C|CCCC  is 1
                                    3,  CC|CCC  is 2
                                    ...
                      -----> if you use one more Restrict Enzymes,you must join cutting site with "-",
                             for example:
                                 1-3
                                 4-2-0
               -O     output files;

 example:     perl $0 -B path/Bowtie2 -R 40-220 -I hg18.fa -P YNCGNR -S 3 -O cut.site
              perl $0 -B path/Bowtie2 -R 40-220 -I XXX.fa -P YNCGNR-GCWGC -S 3-2 -O results.out
 author:      2826
 date:        2011.07.06
 modify author: 548
 date:	      2014.04.24
USAGE

my ($help,$in,$bowtie2,$site,$left,$right,$out,$region);
GetOptions(
    "help"=>\$help,
    "B=s"=>\$bowtie2,
    "I=s" =>\$in,
    "P=s" =>\$site,
    "S=s" =>\$left,
    "O:s" =>\$out,
    "R:s" =>\$region
);
die $usage unless ($in && $site && $left && $bowtie2 && $out);
die $usage if ($help);
$out ||= "./cut.site";
$region ||="40-220";
my @reg_cut=split /\-/,$region;
my @Ps=split("-",$site);
my @Ss=split("-",$left);
my $base_fa = basename($in);
my $p_len;
my $time=`date`;
open OUT,">$out" or die "Cannot open : $out\n";
my $dir=dirname $in;
mkpath "$dir/C2T";
mkpath "$dir/G2A";
my $build_bowtie="$bowtie2/bowtie2-build";
if($#Ps != $#Ss){
    die "reset your -P or -S parameters!\n"
}
else{
    print STDERR "Starting cut at -> $time";
    for (my $j=0;$j<=$#Ps;$j++){
        $site = uc($Ps[$j]);
        my $pattern = convert($site);
        $p_len = length $site;
        print STDERR "Starting cut $base_fa with $site  -------------------->\n";
        open IN,"<$in" or die "Cannot open : $in\n";
        my $total_site;
        $/ = ">"; <IN>; $/ = "\n";
        while (<IN>){
            chomp;
            my $chr = $_;
            $chr =~ s/^(\w+).*$/$1/;
            $/ = ">";
            my $seq = <IN>;
            $/ = "\n";
            chomp($seq);
            $seq =~ s/\s+//g;
            $seq =~ s/>$//;
            $seq = uc($seq);
			my $total = cut($chr,\$seq,$Ss[$j],$pattern);
            print STDERR "$chr\tcontains $site:\t $total\n";
            $total_site += $total;
        }
        print STDERR "Total\tcontains $site:\t $total_site\n";
        close IN;
    }
}
close OUT;
###################
`sort -k1,1 -k5,5g -u -T ./ $out -o $out`;
###################
$time=`date`;
print STDERR "--------------------> Finish Cutting $in with $site\n";
print STDERR "--------------------> Finish Cutting at  -> $time";


open REF, $in or die $!;
my ($chrID, %seq);
open C2T,">$dir/C2T/c2t.trans.fa" or die $!;
open G2A,">$dir/G2A/g2a.trans.fa" or die $!;
while(<REF>){
	chomp;
	if(/>(\S+)/){
		print C2T $_."\n";
		print G2A $_."\n";
		$chrID = $1;	
	}else{	
		my $line1=$_;
		my $line2=$_;
		$seq{$chrID} .= $_;
		$line1=~tr/Cc/Tt/ ;
		$line2=~tr/Gg/Aa/ ;
		print C2T $line1."\n";
		print G2A $line2."\n";
	}
	
}
close REF;
close C2T;
close G2A;
open CUT,"<$out" or die "can not open file: $out\n";
open OUTA,">$out.frag" or die "can not open file: $out.frag\n";
my ($sta,$n,$len,$gene,$ends);
while (<CUT>){
    chomp;
    my @line = (split(/\t/,$_))[0,4];
    if($.==1 || $line[0] ne $n){
        $len = $line[1];
        $gene = substr($seq{$line[0]}, 0, $len);
        print OUTA join "\t",$line[0],1,$line[1],$len,$gene,"\n";
    }
    elsif($line[0] eq $n){
        $len = $line[1] - $sta;
        $ends = $line[1];
        $gene = substr($seq{$line[0]}, $sta, $len);
        print OUTA join "\t",$line[0],$sta + 1,$ends,$len,$gene,"\n";
    }
    $sta = $line[1];
    $n = $line[0];
}
close CUT;
close OUTA;
####################
open FG,"<$out.frag" or die $!;
open OUTM,">$out.mpileup" or die $!;
open OUTR,">$out.frag.$region" or die $!;
my $sinsert=$reg_cut[0];
my $binsert=$reg_cut[1];
while(<FG>){
	chomp;
	my @a=split;
	if($a[3]>=$sinsert && $a[3]<=$binsert){
		print OUTR $_."\n";
		foreach my $pose($a[1]..$a[2]){
			print OUTM "$a[0]\t$pose\n";
		}
	}	
	
}
close FG;
close OUTM;
###################
sub cut{
    my ($chr,$seq,$ss,$pat) = @_;
    my $num;
    my $regex = eval {qr/$pat/};
    print STDERR "check your pattern! $@" if $@;
    while ($$seq=~m/$regex/ig){
        $num++;
        my $seqs = $1;
        my $end = pos ($$seq);
        my $star = $end - $p_len + 1;
        my $st = $ss + $star - 1;
        # $ed = $end + $right;
        print OUT join "\t",$chr,$seqs,$star,$end,$st,"\n";
    }
    return $num;
}
###################
sub convert{
    my $s = shift;
    my %Iupac = (
        'A' => 'A',     'C' => 'C',     'G' => 'G',     'T' => 'T',
        'M' => '[AC]',  'R' => '[AG]',  'W' => '[AT]',  'S' => '[CG]',
        'Y' => '[CT]',  'K' => '[GT]',  'V' => '[ACG]', 'H' => '[ACT]',
        'D' => '[AGT]', 'B' => '[CGT]', 'N' => '[ACGTN]',
    );
    #print "$s\n";
    my @a = split(//,$s);
    my $par;
    foreach my $b (@a){                
        $par .= $Iupac{$b} ? "$Iupac{$b}" : '[ACGTN]';
    }
    $par = '('.$par.')';
    return $par;
}
###################
$time=`date`;
print STDERR "--------------------> Finish CUT OUTput at $time\n";
print STDERR "$build_bowtie $dir/C2T/c2t.trans.fa $dir/C2T/c2t.trans.fa\n ";
###################
my @cmd;
push @cmd, "$build_bowtie $dir/C2T/c2t.trans.fa $dir/C2T/c2t.trans.fa";
push @cmd, "$build_bowtie $dir/G2A/g2a.trans.fa $dir/G2A/g2a.trans.fa";
# Dispatch
my $pm = new Parallel::ForkManager(2);
for(my $ii = 0; $ii < 2; $ii++) {
	print "$cmd[$ii]\n";
	$pm->start and next;
	system("$cmd[$ii]");
	$pm->finish;
}
#system ("$build_bowtie",  "$dir/C2T/c2t.trans.fa","$dir/C2T/c2t.trans.fa");
#system ("$build_bowtie",  "$dir/G2A/g2a.trans.fa","$dir/G2A/g2a.trans.fa");
$pm->wait_all_children;
$time=`date`;
###################
print STDERR "--------------------> Finish build_bowtie2 at  -> $time \n";

