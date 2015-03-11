#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin '$Bin';
use File::Path;  ## function " mkpath" and "rmtree" deal with directory
use File::Basename qw(basename dirname);
use lib "$Bin/lib/compress";
use Compress::Raw::Zlib;
use Compress::Zlib;
my ($mode,$Help);
GetOptions(
        "mode:s"=>\$mode,
        "help"=>\$Help
);

&usage() if( !$mode || @ARGV!=3 || $Help);
$mode ||="pe";

sub usage{
        die "perl $0  --mode <se|pe> Ref sam out.align\n" ;
}


my %Ref;

my $reference=shift;
my $soap=shift;
my $out=shift;
#my $newsam=shift;
my $pout;
if($out=~/\.gz/){
	$pout = gzopen("$out","wb") or die $!;	
}else{
	open $pout,">$out" or die $!;
}

#open SM,">$newsam" or die $!;

if($soap=~/\.gz$/){
	open SOP,"gunzip -c $soap|" ||die "no soap\n";
}elsif($soap=~/\.bam$/){
	open SOP,"$Bin/samtools view $soap|" or die "can't open bam\n";
}else{
	open SOP,$soap ||die "no soap \n";
}

&Getref($reference);

while(<SOP>){	
	chomp;
	if(/^@/){
		#print SM "$_\n";
		next;
	}
	my @line=split;
	my $ab=$line[1];
	#my $strand=$line[6];
	my $word;
	my $fqbase;
	if($mode=~/pe/){	
	if(($ab & 0x40  && !($ab & 0x10)) || ($ab & 0x80 &&  $ab & 0x10)){
		$word='W';
	}
	if(($ab & 0x40 && $ab & 0x10) || ($ab &0x80  && !($ab & 0x10)) ){
		$word='C';
	}
	}elsif($mode=~/se/){
		if(!($ab&0x10) ){
			$word='W';
		}
		if($ab&0x10){
			$word='C';	
		}
	}
	$fqbase=$line[9];
	my $lenth;     
	my $cigarforlenth=$line[5];
        if($line[5]=~/(\d+)M/){
                $lenth=$1;
        }else{
           while ($cigarforlenth =~ /^(\d+)(\D)/){
                my ($len,$flag) = ($1,$2);
                if($flag eq 'M'){
                        $lenth += $len;
                }
                elsif($flag eq 'D'){
                	$lenth += $len;
		}
                elsif($flag eq 'I'){
                }
                elsif($flag eq 'S'){
                }
                $cigarforlenth =~ s/^\d+\D//;
          }
	}	
	my $readlenth=length($line[9]);
	my @methylation;
	my $refbase;
	if($lenth ==  $readlenth){
		$refbase=&Substr($line[2],$line[3],$lenth);	
		@methylation=&Methy($fqbase,$refbase,$word);
	}else{
		$refbase=&Substrmd($line[2],$line[3],$line[5]);
		my $fqbasemodify=&Modify_Rseq_Based_on_Cigar($line[5],$line[9]);
		@methylation=&Methy($fqbasemodify,$refbase,$word);
	}
	my @head;
	my @score;
	$head[0]=0;
	$score[0]=1;	
	if($methylation[0] eq ' '){ $score[0]=-3;}
	for (my $i=1;$i<@methylation;$i++){
		if($methylation[$i] eq ' '){
			$score[$i] = $score[$i-1]-3;
			$head[$i] = $head[$i-1];
		}else{
			if($score[$i-1]>=0){
				$score[$i] = $score[$i-1]+1;
				$head[$i] = $head[$i-1];
			}else{
				$score[$i] = 1;
				$head[$i] = $i;#start pos;
			}
		}
	}
	my $k=0; #last postition for best map;del end;
	for(my $i=1,my $j=$score[0];$i<@methylation;$i++){
		if($score[$i]>$j){
			$j=$score[$i];
			$k=$i;
		}
	}

	if($head[$k]>0){
		for(my $start=0;$start<$head[$k];$start++){
			$methylation[$start]=' ';
		}
	}
	if($k<@methylation){
		for(my $end=$k+1;$end<@methylation;$end++){
			$methylation[$end]=' ';
		}			
	}
	
	my $methylation=(join "", @methylation);
	my $space=$methylation=~tr/ / /;
#make a newsam for calling snp
	
	#if($space<=2){
	#print SM $_."\n";
	#}
	my $repeat; my $strand;
	if($ab & 0x10){
		$strand = '-';	
	}else{
		$strand = '+';
	}
	if($line[1] & 0x100){
		$repeat=2;
	}else{
		$repeat=1;
	}
	
	my $align="$line[0]\:rep$repeat\:$strand\:len$lenth\:$word\:score$line[4]\t$line[2]\t$line[3]\t$line[9]\t$methylation\t$refbase\t$line[10]\n";
	if($out=~/\.gz$/){
		$pout->gzwrite($align) ||die $!;	
	}else{
		print $pout "$line[0]\:rep$repeat\:$strand\:len$lenth\:$word\:score$line[4]\t$line[2]\t$line[3]\t$line[9]\t$methylation\t$refbase\t$line[10]\n";
	}
	
}
close SOP;

if($out=~/\.gz/){
	$pout->gzclose();
	
}else{
	close  $pout;
}

sub Methy
{
	my $fqbase=shift;
	my $refbase=shift;
	my $word=shift;
	my @Query=split //,$fqbase;	
	my @Ref=split //, $refbase;
	my @align;

	if($word eq 'W'){
		for (my $i=0;$i<@Query;$i++){
			$align[$i]=' ';
			$Ref[$i]=uc($Ref[$i]);
			$Query[$i]=uc($Query[$i]);
#base eq ref
			if($Ref[$i] eq $Query[$i]){
				if($Query[$i] eq 'C'){
					$align[$i]='m';
				}else{	
					$align[$i]='|';
				}
			}elsif($Query[$i] eq 'T' && $Ref[$i] eq 'C'){
				$align[$i]='*';
			}
#base ne ref
			if($Query[$i] eq 'G' && ($Ref[$i] eq 'R' || $Ref[$i] eq 'K' || $Ref[$i] eq 'S' || $Ref[$i] eq 'B' || $Ref[$i] eq 'D' || $Ref[$i] eq 'V')){
				$align[$i]='|';
			}
			if($Query[$i] eq 'A' && ($Ref[$i] eq 'R' || $Ref[$i] eq 'M' || $Ref[$i] eq 'W' || $Ref[$i] eq 'H' || $Ref[$i] eq 'D' || $Ref[$i] eq 'V')){
                                $align[$i]='|';
                        }
			if($Query[$i] eq 'T' && ($Ref[$i] eq 'Y' || $Ref[$i] eq 'K' || $Ref[$i] eq 'W' || $Ref[$i] eq 'H' || $Ref[$i] eq 'D' || $Ref[$i] eq 'B')){
                                $align[$i]='|';
                        }
			if($Query[$i] eq 'T' && ($Ref[$i] eq 'M' || $Ref[$i] eq 'S' || $Ref[$i] eq 'V' )){
                                $align[$i]='*';
                        }
			if($Query[$i] eq 'C' && ($Ref[$i] eq 'Y' || $Ref[$i] eq 'M' || $Ref[$i] eq 'S' || $Ref[$i] eq 'B' || $Ref[$i] eq 'H' || $Ref[$i] eq 'V')){
                                $align[$i]='m';
                        }			
				
		}
		
	}
	else{
		for (my $i=0;$i<@Query;$i++){
			$align[$i]=' ';
                    $Ref[$i]=uc($Ref[$i]);
                    $Query[$i]=uc($Query[$i]);
#G pos
                    if($Ref[$i] eq $Query[$i]){
                        if($Query[$i] eq 'G'){
                                        $align[$i]='m';
                        }else{
                                        $align[$i]='|';
                        }
                     }elsif($Query[$i] eq 'A' && $Ref[$i] eq 'G'){
                                $align[$i]='#';
                     }
	
#nonG pos
                        if($Query[$i] eq 'T' && ($Ref[$i] eq 'Y' || $Ref[$i] eq 'K' || $Ref[$i] eq 'W' || $Ref[$i] eq 'B' || $Ref[$i] eq 'D' || $Ref[$i]  eq 'H')){
                                $align[$i]='|';
                        }
                        if($Query[$i] eq 'C' && ($Ref[$i] eq 'Y' || $Ref[$i] eq 'M' || $Ref[$i] eq 'S' || $Ref[$i] eq 'B' || $Ref[$i] eq 'H' || $Ref[$i] eq 'V')){
                                $align[$i]='|';
                        }
                        if($Query[$i] eq 'A' && ($Ref[$i] eq 'R' || $Ref[$i] eq 'M' || $Ref[$i] eq 'W' || $Ref[$i] eq 'D' || $Ref[$i] eq 'H' || $Ref[$i] eq 'V')){
                                $align[$i]='|';
                        }
                        if($Query[$i] eq 'A' && ($Ref[$i] eq 'K' || $Ref[$i] eq 'S' || $Ref[$i] eq 'B' )){
                                $align[$i]='*';
                        }
                        if($Query[$i] eq 'G' && ($Ref[$i] eq 'R' || $Ref[$i] eq 'K' || $Ref[$i] eq 'S' || $Ref[$i] eq 'B' || $Ref[$i] eq 'D' || $Ref[$i] eq 'V')){
                                $align[$i]='m';
                        }

                }
	}

	#my $meth =(join "", @align);
	return @align;
}

sub Getref
{
	my $ref=shift;
	open RF,$ref or die "need reference\n";
	my $chr;
	while(<RF>){
		chop;
		if(/^>(chr\w+)/){
			$chr=$1;
			print STDERR "$chr\n";
		}else{
			$_=~s/\s+//g;
			$Ref{$chr}.=$_;
		}
	}
	close RF;	
}

sub Substr
{
	my $chr=shift;
	my $start=shift;
	my $length=shift;
	$start--;
	my $refbase=substr($Ref{$chr},$start,$length);
	return $refbase;
}

sub Substrmd
{
	my $chr=shift;
	my $start=shift;
	my $cigar=shift;
	$start--;
	my $anchor_pos = $start;
	my $mod_rseq = '';
	
	while ($cigar =~ /^(\d+)(\D)/){
		my ($len,$flag) = ($1,$2);
		if($flag eq 'M'){
			$mod_rseq .= uc(substr($Ref{$chr},$anchor_pos,$len)); # upper case letter
			$anchor_pos += $len;
		}
		elsif($flag eq 'D'){
			$mod_rseq .= uc(substr($Ref{$chr},$anchor_pos,$len)); # upper case letter
			$anchor_pos += $len;
			
		}
		elsif($flag eq 'I'){
			#$mod_rseq .= 'I' x $len; # upper case letter
			
		}
		elsif($flag eq 'S'){
			$anchor_pos += $len;
		}
		$cigar =~ s/^\d+\D//;
	}
	return $mod_rseq;
	
}

#-------- modify rseq based on cigar string
sub Modify_Rseq_Based_on_Cigar{
	my ($cigar,$rseq) = @_;
	my $anchor_pos = 0;
	my $mod_rseq = '';
	while ($cigar =~ /^(\d+)(\D)/){
		my ($len,$flag) = ($1,$2);
		if($flag eq 'M'){
			$mod_rseq .= uc(substr($rseq,$anchor_pos,$len)); # upper case letter
			$anchor_pos += $len;
		}
		elsif($flag eq 'D'){
			$mod_rseq .= 'N' x $len; # upper case letter
		}
		elsif($flag eq 'I'){
			#$mod_rseq .= lc(substr($rseq,$anchor_pos,$len)); # lower case letter
			$anchor_pos += $len;
		}
		elsif($flag eq 'S'){
			$anchor_pos += $len;
		}
		$cigar =~ s/^\d+\D//;
	}
	return $mod_rseq;
}



