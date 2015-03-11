use strict;
use warnings;
use FindBin '$Bin';
use lib "$Bin/lib";
use Math::CDF qw/pchisq/;
use Text::NSP::Measures::2D::Fisher::twotailed qw/calculateStatistic/;
use Statistics::PointEstimation;
use Statistics::TTest;
use Getopt::Long;

my ($cout1, $cout2, $out, $model, $interval, $cpg_num, $help);

GetOptions(
	"i1:s"=>\$cout1,
	"i2:s"=>\$cout2,
	"o:s"=>\$out,
	"model:s"=>\$model,
	"it:i"=>\$interval,
	"cg:i"=>\$cpg_num,
	"h!"=>\$help,
);

$help =<<qq;
	
	****************************************************
	Author: meijunpu\@genomics.org.cn; 2014-05-14
	****************************************************

	perl $0 -i1 tumor.cout -i2 normal.cout -o chr*.CpG -model [core|tend] [-it (280) -cg (10) -h]

	-i1:	tumor.cout[.gz]
	-i2:	normal.cout[.gz]
	-o:	output
	-model:	tend or core
	-it:	interval of adjacent cpg [280bp]
	-cg:	number of cp (2 * number of cpg)[10]
	-h:	help
qq

$interval ||= 280;
$cpg_num ||= 10;

die "$help\n" if(!$cout1 or !$cout2 or !$out or !$model or ($model ne "core" && $model ne "tend"));

open CU1,($cout1 =~ /\.gz/)?"gunzip -c $cout1|":$cout1 or die $!;
open CU2,($cout2 =~ /\.gz/)?"gunzip -c $cout2|":$cout2 or die $!;
open OU,">$out" or die $!;

print OU "#chr\tpos_start\tpos_end\tshared_num/methy_in_normal/methy_in_tumor/pos_num\ttumor_methy_num/tumor_NONE_methy_num\ttumor_methy_ratio\tnormal_methy_num/normal_NONE_methy_num\tnormal_methy_ratio\tchisq.test_p-value\tt.test_p-value\n";
my $ncpg = 0;
my (@ratio_tum, @ratio_nor);
my ($t_m_n, $t_un_m_n, $n_m_n, $n_un_m_n) = (0,0,0,0);
my ($share, $in_normal, $in_tumor, $pos_num) = (0, 0, 0, 0);
my ($chiq_Pvalue, $t_Pvalue);
my ($start, $end);
my @m;
my %first;

###  tumor  normal  ###
while(my $line1 = <CU1>){
	my $line2 = <CU2>;
	chomp ($line1, $line2);
	@m = split /\t/,(join "\t",$line1, $line2);
	next if(($m[3] !~ /CG/ && $m[3] !~ /GC/) or ($m[12] !~ /CG/ && $m[12] !~ /GC/) );

	if(!$ncpg){	### 起始，$ncpg 为0，直接赋值
		############    recover all info	###########
		$ncpg = 1;
		@ratio_tum = ();
		@ratio_nor = ();
		%first = ();
		for my $i(0..$#m){  $first{$ncpg}->[$i] = $m[$i];}
		($start, $end) = ($m[1], $m[1]);
		############    recover all info    ###########
	}elsif($ncpg < $cpg_num){	### $ncpg < $cpg_num，只看距离是否符合要求来决定是否延伸cpg区域
		if($ncpg == 1){	### $ncpg 为1 时，看第一个cpg的甲基化是否为0，如果为0则初始化
			if($m[1] - $end != 1){	### $ncpg 为奇数时，看2个cg是否距离相差1；不为1 则重置hash及数组
				$ncpg = 1;
				@ratio_tum = ();
				@ratio_nor = ();
				%first = ();
				for my $i(0..$#m){  $first{$ncpg}->[$i] = $m[$i];}
				($start, $end) = ($m[1], $m[1]);
			}else{	### 2 个 cg距离相差1，看tumor及normal中是否含有甲基化reads；如果没有甲基化，则ncpg归0；有甲基化则延伸，获得第一个cpg
				my $z_sum = $m[6] + $m[7] + $m[15] + $m[16] + $first{1}->[6] + $first{1}->[7] + $first{1}->[15] + $first{1}->[16];
#				print STDERR "$.\t$ncpg\t$m[1]\t$z_sum\n";
				if (!$z_sum){	### 第一个cpg无甲基化reads，直接剔除
					$ncpg = 0;
				}else{
					$ncpg ++;	### ncpg为偶数，更新数据
					for my $i(0..$#m){  $first{$ncpg}->[$i] = $m[$i];}
					my ($ratio_tum, $ratio_nor) = &m_ratio($ncpg, \%first);
					push @ratio_tum, $ratio_tum;
					push @ratio_nor, $ratio_nor;
					$end = $m[1];
#					print STDERR "$.\t$ncpg\t$m[1]\t$end\t$z_sum\n";
				}
			}
		}else{	### ncpg > 1；根据ncpg的奇偶性来判断是否延伸
#			print STDERR "$.\t$ncpg\t$m[1]\t$end\n";
			if($ncpg % 2 == 0){	### ncpg为偶数，则看相邻的2个cpg的距离是否小于interval
				if($m[1] - $end > $interval){	### ncpg为偶数，相邻2个cpg的距离大于interval，ncpg重置为1，起始位点重置
					############    recover all info    ###########
					$ncpg = 1;
					@ratio_tum = ();
					@ratio_nor = ();
					%first = ();
					for my $i(0..$#m){	$first{$ncpg}->[$i] = $m[$i];}
					($start, $end) = ($m[1], $m[1]);
					############    recover all info    ###########
				}else{	### ncpg为偶数，相邻2个cpg距离小于interval，ncpg++
					$ncpg ++;	### ncpg为奇数，不更新数组
					for my $i(0..$#m){  $first{$ncpg}->[$i] = $m[$i];}
					$end = $m[1];
				}
			}else{	### ncpg为奇数，即cpg不完整，看相邻2个cg距离是否为1
				if($m[1] - $end != 1){	### ncpg为奇数，cpg不完整。相邻2个cg距离不为1，ncpg重置为1，起始位点重置
					$ncpg = 1;
					@ratio_tum = ();
					@ratio_nor = ();
					%first = ();
					for my $i(0..$#m){	$first{$ncpg}->[$i] = $m[$i];}
					($start, $end) = ($m[1], $m[1]);
				}else{	### ncpg为奇数，cpg不完整，相邻2个cg距离为1，ncpg ++，cpg ++
					$ncpg ++;	### ncpg为偶数，更新数组
					for my $i(0..$#m){  $first{$ncpg}->[$i] = $m[$i];}
					my ($ratio_tum, $ratio_nor) = &m_ratio($ncpg, \%first);
					push @ratio_tum, $ratio_tum;
					push @ratio_nor, $ratio_nor;
					$end = $m[1];
				}
			}
		}
	}elsif($ncpg == $cpg_num){	### ncpg == cpg_num，根据model的情况来判断是否进行额外检验，以决定是否延伸。
		if($model eq "tend"){	### model为tend，即不需进行额外检验，只需判断前后2个cpg的距离是否小于interval来决定是否延伸cpg区域。
			if($m[1] - $end > $interval){	### model 为tend，ncpg == cpg_num，接下来的一个cpg距离过远，直接输出dmr区，重置ncpg，重置起始位点
				my $tTest = new Statistics::TTest;
				$tTest->set_significance(99);
				$tTest->load_data( \@ratio_tum, \@ratio_nor);
				$t_Pvalue = sprintf("%e", $tTest->{t_prob});
				############    T test    ############
				($share, $in_normal, $in_tumor, $pos_num, $t_m_n, $t_un_m_n, $n_m_n, $n_un_m_n) = &cal($ncpg, \%first);
				############	chiq test      ############
				$chiq_Pvalue = sprintf("%e", &chisqTest( $t_m_n, $t_un_m_n, $n_m_n, $n_un_m_n ));
				############	T test    ############
				my $info = join "\t",$m[0],$start, $end, (join "/",$share, $in_normal, $in_tumor, $pos_num), (join "/", $t_m_n, $t_un_m_n), ($t_m_n + $t_un_m_n)?(sprintf("%.8f", $t_m_n / ($t_m_n + $t_un_m_n))):0, (join "/", $n_m_n, $n_un_m_n), ($n_m_n + $n_un_m_n)?(sprintf("%.8f", $n_m_n / ($n_m_n + $n_un_m_n))):0, $chiq_Pvalue, $t_Pvalue;
				print OU $info,"\n";
				for my $j(1..$ncpg){my $mt = $first{$j}->[0];for my $i(1..$#m)
					{$mt .= "\t".$first{$j}->[$i];}#print STDERR "$mt\n";
					}# print STDERR "129\n";
				#######  重置ncpg、起始位点、哈希、数组
				$ncpg = 1;
				@ratio_tum = ();
				@ratio_nor = ();
				%first = ();
				for my $i(0..$#m)
				{	$first{$ncpg}->[$i] = $m[$i];}
				($start, $end) = ($m[1], $m[1]);
			}else{	### model 为tend，ncpg == cpg_num，接下来的cpg在interval范围内，延伸cpg区域。
				$ncpg ++;	### ncpg为奇数，不更新数组
				for my $i(0..$#m)
					{  $first{$ncpg}->[$i] = $m[$i];}
				$end = $m[1];
			}
		}else{	### model 为 core，需检测 1. T检验及卡方检验p值 2. tumor及normal甲基化率的差异
			my $tTest = new Statistics::TTest;
			$tTest->set_significance(99);
			$tTest->load_data( \@ratio_tum, \@ratio_nor);
			$t_Pvalue = sprintf("%e", $tTest->{t_prob});
			############    T test    ############
			($share, $in_normal, $in_tumor, $pos_num, $t_m_n, $t_un_m_n, $n_m_n, $n_un_m_n) = &cal($ncpg, \%first);
			############	chiq test      ############
			$chiq_Pvalue = sprintf("%e", &chisqTest( $t_m_n, $t_un_m_n, $n_m_n, $n_un_m_n ));
			############    methylation difference test    #########
			my ($T_rat, $N_rat) = (&rat($t_m_n, $t_un_m_n,), &rat($n_m_n, $n_un_m_n));
			my $T_N = 0;
			if($T_rat - $N_rat > 0.1 or $T_rat - $N_rat < -0.1){
				$T_N = 1;
			}
			############    methylation difference test    ############
			if($t_Pvalue > 0.05 or $chiq_Pvalue > 0.05 or !$T_N){	#### 未通过检验，剔除第一个cpg
				############     remove the first cpg    ###########
				shift @ratio_tum;
				shift @ratio_nor;
				$start = $first{3}->[1];
				for my $j(3..$ncpg){
					for my $i(0..$#m)
					{$first{$j-2}->[$i] = $first{$j}->[$i];	}
				}
				$ncpg -= 2;
				############     remove the first cpg    ###########
				if($m[1] - $end <= $interval){	########  第一个cpg被剔除，只剩4个；判断下一个cg的距离
					$ncpg ++;	### ncpg为奇数，不更新数组
					$end = $m[1];
					for my $i(0..$#m)
					{  $first{$ncpg}->[$i] = $m[$i];}
				}else{	################  do not extend, recover all info； 第4个cpg与下一个cpg距离大于interval，重置ncpg、起始位点、哈希、数组
					$ncpg = 1;
					@ratio_tum = ();
					@ratio_nor = ();
					%first = ();
					for my $i(0..$#m) 
					{	$first{$ncpg}->[$i] = $m[$i];}
					($start, $end) = ($m[1], $m[1]);
				}
			}else{	#### 通过检验，已有5个完整cpg，进行检验，判断是否延伸cpg区域
				if($m[1] - $end > $interval){	#########	do not extend; print out if the last element > 0；判断第5个和第6个cpg之间距离大于interval，不延伸，直接输出，重置ncpg、起始位点、哈希、数组
					my $info = join "\t",$m[0],$start, $end, (join "/",$share, $in_normal, $in_tumor, $pos_num), (join "/", $t_m_n, $t_un_m_n), ($t_m_n + $t_un_m_n)?(sprintf("%.8f", $t_m_n / ($t_m_n + $t_un_m_n))):0, (join "/", $n_m_n, $n_un_m_n), ($n_m_n + $n_un_m_n)?(sprintf("%.8f", $n_m_n / ($n_m_n + $n_un_m_n))):0, $chiq_Pvalue, $t_Pvalue;
					print OU $info,"\n";
					for my $j(1..$ncpg){	my $mt = $first{$j}->[0];for my $i(1..$#m)
						{$mt .= "\t".$first{$j}->[$i];}#print STDERR "$mt\n";
						}# print STDERR "193\n";
					#######  重置ncpg、起始位点、哈希、数组
					$ncpg = 1;
					@ratio_tum = ();
					@ratio_nor = ();
					%first = ();
					for my $i(0..$#m)
					{	$first{$ncpg}->[$i] = $m[$i];}
					($start, $end) = ($m[1], $m[1]);
				}else{	### 第5个和第6个cpg之间距离小于interval，延伸，获得第6个cpg的上半部分。
				###############   tmp extend    ##########
					$ncpg ++;	### ncpg为奇数，不更新数组
					for my $i(0..$#m)
					{  $first{$ncpg}->[$i] = $m[$i];}
					$end = $m[1];
				}
			}
		}
	}else{	### 已经得到至少5个cpg;如ncpg为奇数，根据model来判断是否延伸；如ncpg为偶数，根据interval来判断是否延伸
		if($ncpg % 2 == 0){	### 得到至少6个完整cpg，根据interval来判断是否延伸
			if($m[1] - $end > $interval){	### 相邻cpg距离大于interval，不延伸，直接输出，重置ncpg、起始位点、哈希、数组
				my $tTest = new Statistics::TTest;
				$tTest->set_significance(99);
				$tTest->load_data( \@ratio_tum, \@ratio_nor);
				$t_Pvalue = sprintf("%e", $tTest->{t_prob});
				($share, $in_normal, $in_tumor, $pos_num, $t_m_n, $t_un_m_n, $n_m_n, $n_un_m_n) = &cal($ncpg, \%first);
				$chiq_Pvalue = sprintf("%e", &chisqTest( $t_m_n, $t_un_m_n, $n_m_n, $n_un_m_n ));
				if($model eq "tend"){
					my $info = join "\t",$m[0],$start, $end, (join "/",$share, $in_normal, $in_tumor, $pos_num), (join "/", $t_m_n, $t_un_m_n), ($t_m_n + $t_un_m_n)?(sprintf("%.8f", $t_m_n / ($t_m_n + $t_un_m_n))):0, (join "/", $n_m_n, $n_un_m_n), ($n_m_n + $n_un_m_n)?(sprintf("%.8f", $n_m_n / ($n_m_n + $n_un_m_n))):0, $chiq_Pvalue, $t_Pvalue;
					print OU $info,"\n";
					for my $j(1..$ncpg){	my $mt = $first{$j}->[0];for my $i(1..$#m)
						{$mt .= "\t".$first{$j}->[$i];}#print STDERR "$mt\n";
						}# print STDERR "193\n";
				}else{
					my ($T_rat, $N_rat) = (&rat($t_m_n, $t_un_m_n,), &rat($n_m_n, $n_un_m_n));
					my $T_N = 0;
					if($T_rat - $N_rat > 0.1 or $T_rat - $N_rat < -0.1){
						$T_N = 1;
					}
					if($t_Pvalue > 0.05 or $chiq_Pvalue > 0.05 or !$T_N){
						######  do not print
					}else{
						my $info = join "\t",$m[0],$start, $end, (join "/",$share, $in_normal, $in_tumor, $pos_num), (join "/", $t_m_n, $t_un_m_n), ($t_m_n + $t_un_m_n)?(sprintf("%.8f", $t_m_n / ($t_m_n + $t_un_m_n))):0, (join "/", $n_m_n, $n_un_m_n), ($n_m_n + $n_un_m_n)?(sprintf("%.8f", $n_m_n / ($n_m_n + $n_un_m_n))):0, $chiq_Pvalue, $t_Pvalue;
						print OU $info,"\n";
						for my $j(1..$ncpg){    my $mt = $first{$j}->[0];for my $i(1..$#m)
							{$mt .= "\t".$first{$j}->[$i];}#print STDERR "$mt\n";
						}# print STDERR "252\n";
					}
				}
				#######  重置ncpg、起始位点、哈希、数组
				$ncpg = 1;
				@ratio_tum = ();
				@ratio_nor = ();
				%first = ();
				for my $i(0..$#m)
				{	$first{$ncpg}->[$i] = $m[$i];}
				($start, $end) = ($m[1], $m[1]);
			}else{	### 相邻cpg距离小于interval，直接延伸
				$ncpg ++;	### ncpg为奇数，不更新数组
				for my $i(0..$#m)
				{  $first{$ncpg}->[$i] = $m[$i];}
				$end = $m[1];
			}
		}else{	### 得到至少5个cpg，一个cpg的起点，1. 判断与下一个cg距离是否为1 2. 根据model来判断是否延伸； ncpg为奇数
			if($m[1] - $end > 1){	### cpg不完整，剔除最末端的一个cg，保留前面完整的cpg，输出;重置ncpg，重置起始位点
				$ncpg -= 1;	### ncpg减去1变为偶数，重新计算P值
				my $tTest = new Statistics::TTest;
				$tTest->set_significance(99);
				$tTest->load_data( \@ratio_tum, \@ratio_nor);
				$t_Pvalue = sprintf("%e", $tTest->{t_prob});
				############    T test    ############
				($share, $in_normal, $in_tumor, $pos_num, $t_m_n, $t_un_m_n, $n_m_n, $n_un_m_n) = &cal($ncpg, \%first);
				############	chiq test      ############
				$chiq_Pvalue = sprintf("%e", &chisqTest( $t_m_n, $t_un_m_n, $n_m_n, $n_un_m_n ));
				############	T test    ############
				my $info = join "\t",$m[0],$start, $end, (join "/",$share, $in_normal, $in_tumor, $pos_num), (join "/", $t_m_n, $t_un_m_n), ($t_m_n + $t_un_m_n)?(sprintf("%.8f", $t_m_n / ($t_m_n + $t_un_m_n))):0, (join "/", $n_m_n, $n_un_m_n), ($n_m_n + $n_un_m_n)?(sprintf("%.8f", $n_m_n / ($n_m_n + $n_un_m_n))):0, $chiq_Pvalue, $t_Pvalue;
				print OU $info,"\n";
				for my $j(1..$ncpg){my $mt = $first{$j}->[0];for my $i(1..$#m)
					{$mt .= "\t".$first{$j}->[$i];}#print STDERR "$mt\n";
				}# print STDERR "129\n";
				#######  重置ncpg、起始位点、哈希、数组
				$ncpg = 1;
				@ratio_tum = ();
				@ratio_nor = ();
				%first = ();
				for my $i(0..$#m)
				{	$first{$ncpg}->[$i] = $m[$i];}
				($start, $end) = ($m[1], $m[1]);
			}else{	### cpg 完整，根据model来判断是否延伸
				if($model eq "tend"){	### model为tend，不需进行检验，直接延伸
					$ncpg ++;	### ncpg为偶数，更新数组
					for my $i(0..$#m){  $first{$ncpg}->[$i] = $m[$i];}
					my ($ratio_tum, $ratio_nor) = &m_ratio($ncpg, \%first);
					push @ratio_tum, $ratio_tum;
					push @ratio_nor, $ratio_nor;
					$end = $m[1];
				}else{	### model 为 core，进行检验以判断是否延伸cpg区域 1. 与前面的core进行比较，是否趋势相同  2. p值、甲基化差异是否更大？
					### 此时，ncpg为奇数； 先暂时延伸dmr区域，检验之后，来判断是否进行延伸。
					######## 得到前面的dmr区域的p值
					my $tTest = new Statistics::TTest;
					$tTest->set_significance(99);
					$tTest->load_data( \@ratio_tum, \@ratio_nor);
					my $t_Pvalue = sprintf("%e", $tTest->{t_prob});
					($share, $in_normal, $in_tumor, $pos_num, $t_m_n, $t_un_m_n, $n_m_n, $n_un_m_n) = &cal($ncpg-1, \%first);
					my $chiq_Pvalue = sprintf("%e", &chisqTest( $t_m_n, $t_un_m_n, $n_m_n, $n_un_m_n ));
					######## 得到前面的dmr区域的p值
					######## 先暂时延伸dmr区域，得到新的P值、甲基化差异情况，判断是否达标以决定是否进行延伸
					$ncpg ++;	### ncpg 为偶数，更新数组
					for my $i(0..$#m)
					{  $first{$ncpg}->[$i] = $m[$i];}
					$end = $m[1];
					my ($ratio_tum, $ratio_nor) = &m_ratio($ncpg, \%first);
					push @ratio_tum, $ratio_tum;
					push @ratio_nor, $ratio_nor;

					my $o_tTest = new Statistics::TTest;
					$o_tTest->set_significance(99);
					$o_tTest->load_data( \@ratio_tum, \@ratio_nor);
					my $o_t_Pvalue = sprintf("%e", $tTest->{t_prob});
					my ($o_share, $o_in_normal, $o_in_tumor, $o_pos_num, $o_t_m_n, $o_t_un_m_n, $o_n_m_n, $o_n_un_m_n) = &cal($ncpg, \%first);
					my $o_chiq_Pvalue = sprintf("%e", &chisqTest( $o_t_m_n, $o_t_un_m_n, $o_n_m_n, $o_n_un_m_n ));
					my ($T_rat, $N_rat) = (&rat($t_m_n, $t_un_m_n,), &rat($n_m_n, $n_un_m_n));
					my $T_N = 0;
					if($T_rat - $N_rat > 0.1 or $T_rat - $N_rat < -0.1){
						$T_N = 1;
					}

					if($o_chiq_Pvalue > $chiq_Pvalue or $o_t_Pvalue > $t_Pvalue or !$T_N){	#########   do not extend; print out if the last element > 0 && $ncpg >= 5    ############
					##### 没通过检验，不延伸dmr区域，只输出前面的原始dmr区域；保留临时dmr区域里面的最后一个cpg
					#########  recover, pop the last element
						$end = $first{$ncpg-2}->[1];
						my $info = join "\t",$m[0],$start, $end, (join "/",$share, $in_normal, $in_tumor, $pos_num), (join "/", $t_m_n, $t_un_m_n), ($t_m_n + $t_un_m_n)?(sprintf("%.8f", $t_m_n / ($t_m_n + $t_un_m_n))):0, (join "/", $n_m_n, $n_un_m_n), ($n_m_n + $n_un_m_n)?(sprintf("%.8f", $n_m_n / ($n_m_n + $n_un_m_n))):0, $chiq_Pvalue, $t_Pvalue;
						print OU $info,"\n";
						for my $j(1..$ncpg-2){my $mt = $first{$j}->[0];for my $i(1..$#m){$mt .= "\t".$first{$j}->[$i];}#print STDERR "$mt\n";
						}# print STDERR "297\n";
						#########   do not extend; print out if the last element > 0 && $ncpg >= 5    ############
						#########   保留临时dmr区域最后一个cpg，更新ncpg、起始位点、数组、哈希
						for my $j($ncpg-1,$ncpg) {
							for my $i(0..$#m) {
								$first{$j-$ncpg+2}->[$i] = $first{$j}->[$i];
							}
						}
						my ($ratio_tum_1, $ratio_nor_1) = ($ratio_tum[-1],$ratio_nor[-1]);
						@ratio_tum = ();
						@ratio_nor = ();
						push @ratio_tum, $ratio_tum_1;
						push @ratio_nor, $ratio_nor_1;
						($start, $end) = ($first{1}->[1], $first{2}->[1]);
						$ncpg = 2;
					}else{
						########  extend dmr successfully  #########
					}
				}
			}
		}
	}
	my $zt = join "\t",@m;
	#print "$ncpg\t$zt\n";
}
close CU1;
close CU2;

if($ncpg >= $cpg_num){
	($share, $in_normal, $in_tumor, $pos_num, $t_m_n, $t_un_m_n, $n_m_n, $n_un_m_n) = &cal($ncpg, \%first);
	my $tTest = new Statistics::TTest;
	$tTest->set_significance(99);
	$tTest->load_data( \@ratio_tum, \@ratio_nor);
	$t_Pvalue = sprintf("%e", $tTest->{t_prob});
	$chiq_Pvalue = sprintf("%e", &chisqTest( $t_m_n, $t_un_m_n, $n_m_n, $n_un_m_n ));
	my ($T_rat, $N_rat) = (&rat($t_m_n, $t_un_m_n,), &rat($n_m_n, $n_un_m_n));
	my $T_N = 0;
	if($T_rat - $N_rat > 0.1 or $T_rat - $N_rat < -0.1){
		$T_N = 1;
	}
	next if($t_Pvalue > 0.05 or $chiq_Pvalue > 0.05 or !$T_N);
	my $info = join "\t",$m[0],$start, $end, (join "/",$share, $in_normal, $in_tumor, $pos_num), (join "/", $t_m_n, $t_un_m_n), ($t_m_n + $t_un_m_n)?(sprintf("%.8f", $t_m_n / ($t_m_n + $t_un_m_n))):0, (join "/", $n_m_n, $n_un_m_n), ($n_m_n + $n_un_m_n)?(sprintf("%.8f", $n_m_n / ($n_m_n + $n_un_m_n))):0, $chiq_Pvalue, $t_Pvalue;
	print OU $info,"\n";
	for my $j(1..$ncpg){my $mt = $first{$j}->[0];for my $i(1..$#m){$mt .= "\t".$first{$j}->[$i];}print STDERR "$mt\n";} print STDERR "346\n";
}
close OU;

sub cal(){
	my ($n,$first) = @_;
	if($n % 2 == 1){
		print STDERR "error\t$.\n";
	}
	my ($a,$b,$c,$d,$e,$f,$g,$h) = (0,0,0,0,0,0,0,0);
	my $m = int ($n / 2);
	my %hash;
	for my $i(1..$m){
		$hash{$i}->[6] = ${$first}{2*$i-1}->[6] + ${$first}{2*$i}->[6];
		$hash{$i}->[7] = ${$first}{2*$i-1}->[7] + ${$first}{2*$i}->[7];
		$hash{$i}->[15] = ${$first}{2*$i-1}->[15] + ${$first}{2*$i}->[15];
		$hash{$i}->[16] = ${$first}{2*$i-1}->[16] + ${$first}{2*$i}->[16];
	}
	for my $j(1..$m){
		$d ++;
		$c ++ if($hash{$j}->[6]);
		$b ++ if($hash{$j}->[15]);
		$a ++ if($hash{$j}->[6] && $hash{$j}->[15]);

		$e += $hash{$j}->[6];
		$f += $hash{$j}->[7];
		$g += $hash{$j}->[15];
		$h += $hash{$j}->[16];
	}
	return ($a, $b, $c, $d, $e, $f, $g, $h);
}

sub m_num_cal(){
	my ($a, $b, $c, $d, $e, $f, $g, $h, @m) = @_;
	$d ++;
	$c ++ if($m[6]);
	$b ++ if($m[15]);
	$a ++ if($m[6] && $m[15]);

	$e += $m[6];
	$f += $m[7];
	$g += $m[15];
	$h += $m[16];

	my $p_z = join "\t", $a, $b, $c, $d, $e, $f, $g, $h;
#	print STDERR $p_z,"\t",$m[1],"\n";
	return ($a, $b, $c, $d, $e, $f, $g, $h);
}

sub m_num_move(){
	my ($a, $b, $c, $d, $e, $f, $g, $h, @m) = @_;
	$d --;
	$c -- if($m[6]);
	$b -- if($m[15]);
	$a -- if($m[6] && $m[15]);

	$e -= $m[6];
	$f -= $m[7];
	$g -= $m[15];
	$h -= $m[16];

#	my $im = join "\t",@m;	print STDERR $im,"\n";
	return ($a, $b, $c, $d, $e, $f, $g, $h);
}

sub m_ratio(){
	my ($n, $hash) = @_;
	my $T_m = ${$hash}{$n-1}->[6] + ${$hash}{$n}->[6];
	my $T_um = ${$hash}{$n-1}->[7] + ${$hash}{$n}->[7];
	my $N_m = ${$hash}{$n-1}->[15] + ${$hash}{$n}->[15];
	my $N_um = ${$hash}{$n-1}->[16] + ${$hash}{$n}->[16];
	my ($rt, $rn) = (&rat($T_m, $T_um), &rat($N_m, $N_um));
	return ($rt, $rn);
}

sub rat(){
	my ($m, $um) = @_;
	my $rt;
	if(!$m && !$um){
		$rt = -1;
	}else{
		$rt = $m / ($m + $um);
	}
	return $rt;
}

sub chisqTest(){
	my ( $a, $b, $c, $d ) = @_;
	if ( $a == 0 && $b == 0 && $c == 0 && $d == 0 ){
		return 1;
	}
	my $tt  = $a + $b + $c + $d;
	my $np1 = $a + $c;
	my $np2 = $b + $d;
	my $n1p = $a + $b;
	my $n2p = $c + $d;
	if ( $np1 < 5 || $np2 < 5 || $n1p < 5 || $n2p < 5 || $tt < 40 ){
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
