use strict;
use File::Basename qw/basename dirname/;
use Cwd 'abs_path';
use FindBin qw($Bin);


die "perl $0 target CG(CHG,CHH) /nas/RD_12A/gaoshengjie/software/R-3.0.0/bin/R >table.out \n" if(@ARGV<2);
my $target=shift;
my $file=shift;
my $R_exec=shift;

my $abs_d_cout = abs_path($file);
$abs_d_cout = dirname $abs_d_cout;
my $target_base = basename $file;
$target_base =~ s/\..*//;
#open TABLE,">$abs_d_cout/coverage.table" or die "can't create the table\n";
#print TABLE "baba";


my %hash;my %thero;
open TG,$target or die $!;
while(<TG>){
	chomp;
	next if(/^#/);
	my @a=split;
	#push @{$hash{$a[0]}},[$a[1],$a[2]];	
	$hash{$a[0]}=$a[2];
	$thero{$a[0]}=$a[5];
}
close TG;

#chr1    1465980 CG      12      14      66      .       .       .
if($file=~/\.gz$/){
	open FILE, "gzip -cd $file | " or die $!;
}else{
	open FILE,$file or die $!;
}
my %hash_stats;
while (<FILE>) {
	chomp;
       # my ( $chr, $pos, $num ) = (split)[ 0, 1, 7 ];
        #if ( exists( $hash_share_pos{"$chr,$pos"} ) ) {
	my @a=split;
	my $num=$a[4]+$a[7];
	#print $num."\n";
        ++$hash_stats{"a1x"}  if ( $num > 0 );
        ++$hash_stats{"b4x"}  if ( $num > 3 );
        ++$hash_stats{"c10x"} if ( $num > 9 );
        ++$hash_stats{"d20x"} if ( $num > 19 );
            #++$hash_stats{"e30x"} if ( $num > 29 );
            #++$hash_stats{"f50x"} if ( $num > 49 );
        #}
}
close FILE;

#print "$file\.overlap\n";die;
open STAT_OUT, ">$file\.overlap" or die $!; 
print STAT_OUT "$file";
my $denominator = $hash{"CG"};
print STAT_OUT "\ttheoretically\t";
print STAT_OUT $thero{"CG"}/$hash{"CG"};

foreach my $type ( sort keys %hash_stats ) { 
    print STAT_OUT "\t", $hash_stats{$type} / $denominator;    
}
print STAT_OUT "\n";


open STAT_R, ">$file.R" or die $!;
my $rscript = " 
# Rscript #
pdf(\"$file.pdf\", height=8, width=12)
par(font.lab=1,font.axis=1,cex.lab=1.9,cex.main=3.5,cex.axis=1.5,mar=c(6,6,4,0.5),mgp=c(3,1,0))
data<-read.table(\"$file.overlap\", head=F)
b<-t(data[3:7])
h<-as.matrix(b)
barplot(h*100,yaxt=\"n\",ylim=c(0,60),ylab=\"Coverage of Overlaped-CpG(%)\",legend.text = c(\"Theoretical\",\">=1X\",\">=4X\",\">=10X\",\">=20X\"),args.legend = list(cex=1.2),beside = TRUE,col=c(rainbow(7)),axes=T)
mtext(\"$target_base\",line=2,side=1,adj=0.5)
axis(side=2,labels=T,at=c(seq(0,70,10)),lwd=0.8)
dev.off()
";


print STAT_R $rscript;
system(
"$R_exec CMD BATCH $file.R"
# && rm -rf $file.overlap $file.R*"
#"$R_exec CMD BATCH $abs_d_cout/$target_base.R"
);




