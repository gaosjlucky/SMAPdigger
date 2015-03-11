#!/usr/bin/perl -w
use strict;
use File::Basename;

die"perl $0 <in_dir> <out_dir>
example sampledir:<in_dir>/sample_name/sample_type/cout_40-220/chr1.cout.cov
"if(@ARGV != 1);

my $dir=shift;
my $out_dir;
$out_dir = shift;
$out_dir ||= ".";

my @in=<$dir/*/*/*/chr1.cout.cov>;

open OUT_cov,">$out_dir/Total-cov.xls";
foreach my $in(@in){
		my $sample_type = basename(dirname(dirname $in));
		my $sample = basename(dirname(dirname(dirname $in)));
		my $dir_name =dirname $in;
		my %hash;my @head;my %hash_t;my %hash_t0;my %hash_t1;my $max=0;my $total_stat=0;my $total_line;
		open OUT_qc,">$out_dir/$sample.$sample_type.stat.xls";
		foreach my $i(1 .. 22,"X","Y","M"){
				my $chr = "chr$i";
				open IN,"$dir_name/$chr.cout.cov";
				my @b = split/\s+/,<IN>;	
				while(<IN>){
						chomp;
						my @c=split/\t/,$_;
						my $n = @c - 1;
						my @temp;
						foreach my $t(1 .. $n){
							if(exists $hash{$b[$t]}){	
								$hash{$b[$t]} = $hash{$b[$t]} + $c[$t];
							}else{
								$hash{$b[$t]} = $c[$t];		
							}
							push @temp,$b[$t];
						}
						@head=@temp;
				}
				close IN;
				open INGZ,"gzip -cd $dir_name/$chr.cout.gz|";
				while(<INGZ>){
						chomp;
						my @c=split;
						$total_line++ if(/^chr/);
						my $n=$c[5]+$c[6];
						$hash_t0{$n}++;
						$max=$n if($n > $max);
						$total_stat=$total_stat + $n;
				}
				close INGZ;
		}
		my $mean_depth=$total_stat/$total_line;
		print OUT_qc "Total Digest Region: $total_line\n";
		print OUT_qc "Mean Depth: $mean_depth\n";
		foreach my $i1(0 .. $max){
				print OUT_qc "$i1 X: $hash_t0{$i1}\n"if(exists $hash_t0{$i1});
				foreach my $kk($i1 .. $max){
						$hash_t1{$i1}=$hash_t1{$i1}+$hash_t0{$kk} if(exists $hash_t0{$kk});		
				}
		}
		my $temp_i=0;
		my $temp_total = $hash_t1{$temp_i};
		open OUTTemp,">$out_dir/$sample-$sample_type.cov";
		open OUTTempR,">$out_dir/$sample-$sample_type.cov.R";
		print OUTTemp"Depth(MeanDepth:$mean_depth.X)\tDepth_rate\n";
		foreach my $k0(0 .. $max){
				my $ratio = $hash_t1{$k0}/$temp_total;
				print OUTTemp "$k0\t$ratio\n";
		}
		close OUTTemp;
		print OUTTempR"pdf(\"$out_dir/$sample-$sample_type.cov.cumulate.depth.pdf\")
n =read.table(\"$sample-$sample_type.cov\",head=T)	
xb = c(\"Depth\",paste(\"MeanDepth:$mean_depth.X\"))
plot(n[,1],n[,2],type = \"n\", xlab = xb, ylab = \"rate\",ylim=c(0,1),pch=\".\",col=\"red\")
lines(n[,1],n[,2],col=c(1,1,0,alpha=0.5))
points(n[,1],n[,2],ylim=c(0,1),type=\"p\",pch=\".\",col=\"red\")
";

		close OUTTempR;
		`Rscript $out_dir/$sample-$sample_type.cov.R`;

		print OUT_cov join("\t",@head)."\n";
		print OUT_cov "$sample\t$sample_type\t";
		foreach my $k3(@head){
				print OUT_cov "$hash{$k3}\t";
		}
		print OUT_cov "\n";

}
