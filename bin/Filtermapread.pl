use strict;
use FindBin '$Bin';
die "Usage: $0 <*.sam>\n" if @ARGV != 1;
my ($file_sam) = @ARGV;

my @same_name_reads;
my @buffer;
#FCC076MACXX:1:1101:10000:118419#CGATGTAT/1      97      chr22   48885963        42      70M     =       48885962
if($file_sam =~/\.bam$/){
        open SAM, "$Bin/samtools view -h  $file_sam|" or die $!;
}else{
        open SAM, $file_sam or die $!;
}
my $first=<SAM>;
chomp $first;
my @end;
while(<SAM>){
	chomp;
	my $second=$_;
	my @first=split /\t/,$first;
	my @second=split /\t/,$second;
	if($first=~/\@/){
		print $first."\n";
		$first=$second;
		next;	
	}else{
	   if($first[0] eq $second[0]){
		@end=();	
		if(!($first[1] & 0x8) && !($first[1] & 0x4)){
			print $first."\n";
			$first=$second;
		}elsif(!($second[1] & 0x8) && !($second[1] & 0x4)){
			print $second."\n";
			$first=$second;
		}else{
			die "some error in $first\t$second\n";
		}
	   }else{
		$end[0]=$first;
		$end[1]=$second;
		print $first."\n";
		$first=$second;
	   }
	}		
}
close SAM;

if(@end >0){
	print $end[1]."\n";
}
