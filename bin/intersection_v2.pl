use strict;
use File::Basename qw/basename dirname/;
use Cwd 'abs_path';
use FindBin qw($Bin);

die "Usage: $0 <cout-dir> <target.bed> <CCGG.bed> <R_exec_path> <fragment-region>\n" if @ARGV != 5;

my ( $d_cout, $f_target, $f_ccgg, $R_exec, $region_size ) = @ARGV;

# default frag size #
# $region_size = "40,220";
my ($r_sp, $r_ep) = (split /,/, $region_size)[0,1];

#"/ifs5/PC_HUMAN_US/USER/zhouquan/work/methy_statas_table/CpGIsland.bed.sort.seq"
#"/ifs5/PC_HUMAN_US/USER/zhouquan/work/methy_statas_table/Upstream2k.bed.sort.seq"
#"/nas/RD_12A/gaoshengjie/Database/Human/RRBS/bin/prepare/hg19/Emy_cut/hg19.ccgg.fragment.bed.frag"
#"/opt/blc/genome/biosoft/R/bin/R"

#check parameters#
my $abs_d_cout = abs_path($d_cout);
my @all_cout   = glob "$abs_d_cout/chr*.cout.gz";
my @all_chr;
foreach my $chr(@all_cout){
    my $name=basename($chr);
    $name=~s/\.cout.gz//g;
    push @all_chr,$name;
    
}

#my @all_chr =
# qw/chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM/;

my ( $target_inter_pos_num, $ccgg_target_share_num) = ( 0, 0 );

my %hash_stats;
my %hash_share_pos;

my $target_base = basename $f_target;
$target_base =~ s/\..*//;

foreach my $cur_chr (@all_chr) {
    my %hash_target;
    open I_target, $f_target or die $!;
    while (<I_target>) {
        chomp;
        my ( $chr, $sp, $ep, $anno, $seq ) = split;
        next if ( $cur_chr ne $chr );
        $seq =~ tr/atgcn/ATGCN/;
        my $pos = &seq_split( $seq, $sp );

        for (@$pos) {
            #$hash_target{$_} = "$sp,$ep";
            $hash_target{$_} = '';
        }
    }
    close I_target;

    my %hash_ccgg;
    open I_CCGG, $f_ccgg or die $!;
    while (<I_CCGG>) {
        chomp;
        my ( $chr, $sp, $ep, $size, $seq ) = split;
        next if ( $cur_chr ne $chr );
        next if ( $size > $r_ep or $size < $r_sp );
        $seq =~ tr/atgcn/ATGCN/;
        my $pos = &seq_split( $seq, $sp );
        for (@$pos) {
            #$hash_ccgg{$_} = "$sp,$ep";
            $hash_ccgg{$_} = '';
        }
    }
    close I_CCGG;

    foreach my $key_pos ( sort keys %hash_ccgg ) {
        if ( exists( $hash_target{$key_pos} ) ) {
            ++$target_inter_pos_num;
            $hash_share_pos{"$cur_chr,$key_pos"} = '';
        }
    }

    $ccgg_target_share_num += keys %hash_target;
    
    %hash_target      = ();
    %hash_ccgg     = ();

    #print STDERR "---start buffer---\n";
    #sleep(180);

    print STDERR "-- $cur_chr process done. --\n";
}

foreach my $file (@all_cout) {
    open FILE, "gzip -cd $file | " or die $!;
    while (<FILE>) {
        #chomp;
        my ( $chr, $pos, $num ) = (split)[ 0, 1, 7 ];
        if ( exists( $hash_share_pos{"$chr,$pos"} ) ) {
            ++$hash_stats{"a1x"}  if ( $num > 0 );
            ++$hash_stats{"b4x"}  if ( $num > 3 );
            ++$hash_stats{"c10x"} if ( $num > 9 );
            ++$hash_stats{"d20x"} if ( $num > 19 );

            #            ++$hash_stats{"e30x"} if ( $num > 29 );
            #            ++$hash_stats{"f50x"} if ( $num > 49 );
        }       
    }
    close FILE;
}

print STDERR "----------------------------------------------------------------\n";

open STAT_OUT, ">$abs_d_cout/$target_base.overlap" or die $!;
print STAT_OUT "$target_base";
my $denominator = $ccgg_target_share_num;
print STAT_OUT "\ttheoretically\t";
print STAT_OUT $target_inter_pos_num / $denominator;

foreach my $type ( sort keys %hash_stats ) {
    print STAT_OUT "\t", $hash_stats{$type} / $denominator;      
}
print STAT_OUT "\n";

open STAT_R, ">$abs_d_cout/$target_base.R" or die $!;
my $rscript = " 
# Rscript #
pdf(\"$abs_d_cout/$target_base.pdf\", height=8, width=12)
par(font.lab=1,font.axis=1,cex.lab=1.9,cex.main=3.5,cex.axis=1.5,mar=c(6,6,4,0.5),mgp=c(3,1,0))
data<-read.table(\"$abs_d_cout/$target_base.overlap\", head=F)
b<-t(data[3:6])
h<-as.matrix(b)
barplot(h*100,yaxt=\"n\",ylim=c(0,70),ylab=\"Coverage of Overlaped-CpG(%)\",legend.text = c(\"Theoretical\",\">=1X\",\">=4X\",\">=10X\",\">=20X\"),args.legend = list(cex=1.5),beside = TRUE,col=c(rainbow(7)),axes=T)
mtext(\"$target_base\",line=2,side=1,adj=0.5)
axis(side=2,labels=T,at=c(seq(0,70,10)),lwd=0.8)
dev.off()
";

print STAT_R $rscript;

system(
#"$R_exec CMD BATCH $abs_d_cout/$target_base.R && rm -rf $abs_d_cout/$target_base.overlap $abs_d_cout/$target_base.R*"
"$R_exec CMD BATCH $abs_d_cout/$target_base.R"
);

# input (seq, sp), return (positions).
sub seq_split {
    my @bases = split //, $_[0];
    my @pos;
    my $real_sp = $_[1];
    my $len     = length $_[0];
    for ( my $idx = 0 ; $idx < $len - 1 ; ++$idx ) {
        if ( $bases[$idx] eq 'C' and $bases[ $idx + 1 ] eq 'G' ) {
            push @pos, $idx + $real_sp;
        }
    }
    return ( \@pos );
}
