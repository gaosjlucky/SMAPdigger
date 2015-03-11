use strict;
use Cwd 'abs_path';
use File::Basename;
use FindBin qw($Bin);

#my $abs_path = abs_path($file);
die "Usage: $0 <asm.input> <output> <Rpath>\n" if @ARGV != 3;
my ($input, $output,$rpath) = @ARGV;
system("rm $output") if -e $output;
my $outdir=abs_path($output);
$outdir=dirname($outdir);
system("$rpath --no-save --silent --slave < $Bin/chisq.stat.R --args $input $output");
print "Chisq.test Done.\n";

open FILE, $output or die $!;
my %hash_stats;
my %pos;
while(<FILE>){
        my ($chrpos,$ref, $alt) = (split /\s+/)[0,2,3];
	if(exists($pos{$chrpos})){
        	next;
	}else{
		 ++ $hash_stats{"$ref>$alt"};
		$pos{$chrpos}=1;
	}
}
close FILE;
#my $outdir=abs_path($output);
#$outdir=dirname($outdir);
#print $outdir ;die;

open OUT, ">$outdir/R_draw_temp_file_used_once_and_never_seen.data" or die $!;
my $title = "";
my $items = "";
foreach my $key_type ("A>C","A>G","A>T","C>A","C>G","C>T","G>A","G>C","G>T","T>A","T>C","T>G"){
        #next if $key_type eq "\"C->T\"" or $key_type eq "G\"->A\"";
        #$title .= "$key_type,";
        $items .= "$hash_stats{$key_type}\t";
}
chop $title;
print OUT $items."\n";
close OUT;

my $rscript =<<RS;

pdf(\"$outdir/stats_results.pdf\", width=8, height=6)
data = read.table(\"$outdir/R_draw_temp_file_used_once_and_never_seen.data\", header=F, as.is=F)
max_num=max(data)
RS
#$rscript .= "colnames(data) <- c($title)\n";

$rscript .=<<RAS;
max_num=ceiling(max_num / 1000) * 1000
barplot(as.matrix(data), xlab=\"Type of Mutations\", ylab=\"No. of Mutations\", col=\"lightblue\", ylim=c(0, max_num),name=c(\"A>C\",\"A>G\",\"A>T\",\"C>A\",\"C>G\",\"C>T\",\"G>A\",\"G>C\",\"G>T\",\"T>A\",\"T>\C",\"T>G\"))
dev.off()
RAS

open ORS,">$outdir/frequency.R";
print ORS $rscript;
close ORS;
print "-"x75,"\n",$rpath,"\n","-"x75,"\n";
system("$rpath --no-save --silent --slave < $outdir/frequency.R");
system("rm $outdir/R_draw_temp_file_used_once_and_never_seen.data");

