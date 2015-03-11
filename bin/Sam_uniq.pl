#!/bin/env perl
use strict;
#use warnings;
use Data::Dumper;
use File::Basename;
use Getopt::Long;
use PerlIO::gzip;
use File::Basename qw(basename dirname);
use FindBin '$Bin';
my ($mode , $Help);
GetOptions(
        "mode:s"=>\$mode,
        "help"=>\$Help
);

&usage() if(@ARGV<2 || $Help || ! $mode);

my $in=shift;
my $out=shift;

if($mode=~/se/){
	&seID(\$in,\$out);
}
if($mode=~/pe/){
	&peID(\$in,\$out);
}


if($mode=~/to/){
	&totalID(\$in,\$out);
}


sub seID

{
	my $input=shift;
	my $output=shift;
	if($$input=~/.bam$/){
		open IN,"$Bin/samtools view  -h $$input|" ||die "no inputse\n";
	}else{	
		open IN,$$input ||die "no se_input"		
	}
	open OUT,">$$output" || die "can't creat output\n";

	my $lastline=<IN>;
	chomp $lastline;
	my @array;
	push (@array,$lastline);
#FCC076MACXX:1:1101:10000:127569#TGACCAAT/1      83      chr19   11516938        42      90M     =       11516872        -156    TGGAGAGAGATAGGGGGGAGTTTAGGTGAGGGAGGGAATGGGGTGGAAATTTGATATTTAGGGTGAGATTTGGGTTTAGGGAGGAGGGTT      [Qaa^Vaababbb^[^a___b`_`a_Xabbdddcb`bbbbb`_Wa[Uaabba`Uaaa_aaabb]bbbbb^[RUHhhhgggggccccc___      AS:i:-11        XN:i:0  XM:i:2  XO:i:0  XG:i:0  NM:i:2  MD:Z:40C32A16   YS:i:0  YT:Z:CP
	while(<IN>) {
		chomp;
		if(/^\@/){
			print OUT $_."\n";
			next;
		}
		my $line=$_;
		my @b=split /\t/,$line;
		my @a=split /\t/,$lastline;
		my $readname;
		#my ($readname1,$readname2,$rep1,$rep2,$strand1,$strand2,$len1,$len2,$word1,$word2);
		if($a[0] eq $b[0]){	
			push @array,$line;
		}else{
			#print OUT $array[0]."\n";
			if(@array==1 && !($a[1]&0x100)){
				
				print OUT "$lastline\n";
			}else{
				&getseline(\@array,\*OUT);
			}
			@array=();
			push(@array, $line);
		}
			
		$lastline=$line;	
	}
	close IN;

	if( @array>0 ){
		my @a=split /\t/, $array[0];
		unless($a[1]&0x100){
		print OUT $array[0]."\n";
		}
	}
	close OUT;
}


sub peID
{
	my $input=shift;
	my $output=shift;
	if($$input=~/.bam$/){
                open IN,"$Bin/samtools view -h  $$input|" ||die "no inputse\n";
        }else{
                open IN,$$input ||die "no se_input"
        }
	open OUT,">$$output" || die "can't create output\n";
	while(<IN>){
		chomp;
		if(/^\@/){
			print OUT $_."\n";
			next;
		}
		
		my $line1=$_;
		my $line2=<IN>;
		my $line3;
		my $line4;
		chomp $line2;
		
		
		my @line1=split /\t/,$line1;
		my @line2=split /\t/,$line2;	
		my @line3;
		my @line4;
		
		my ($name1,$name2,$name3,$name4);
		if($line1[0] ne $line2[0]){
			if($line1[0]=~/(\S+)(\d)$/){
				$name1=$1;
			}
			if($line1[0]=~/(\S+)(\d)$/){
				$name2=$1;
			}
			#################################
			if($name1 eq $name2){
				&getpeline(\@line1,\@line2,$name1,$name2,\*OUT);
			}else{
				print STDERR "$name1 and $name2 are not pe\n";
			}
		}
=head		
		else{
			$line3=<IN>;
			$line4=<IN>;
			chomp($line3,$line4);
			@line3=split (/\t/,$line3);
			@line4=split (/\t/,$line4);
				
			if($line3[0]=~/^(\S+\/)(\d)$/){
             	        	$name3=$1;;
               		}
                	if($line4[0]=~/^(\S+\/)(\d)$/){
                 		$name4=$1;;
                	}
			######################################################################			
			if(($line1[4]+$line3[4])==($line2[4]+$line4[4])){
				print STDERR "$line1[0]  $line3[0]  this pair cant match correctly\n";
			}elsif(($line1[4]==255 || $line3[4]==255) && ($line2[4]!=255 && $line4[4]!=255)){
				&getpeline(\@line2,\@line4,$name2,$name4,\*OUT);
			}elsif(($line2[4]==255 || $line4[4]==255) && ($line1[4]!=255 && $line3[4]!=255)){
				&getpeline(\@line1,\@line3,$name1,$name3,\*OUT);
			}elsif($line1[4]==255 || $line2[4]==255 || $line3[4]==255 || $line4[4]==255){
				print STDERR "$name1 this pair cant match correctly\n";
			}elsif(($line1[4]+$line3[4])>($line2[4]+$line4[4])){
				&getpeline(\@line1,\@line3,$name1,$name3,\*OUT);
			}elsif(($line1[4]+$line3[4])<($line2[4]+$line4[4])){
				&getpeline(\@line2,\@line4,$name2,$name4,\*OUT);
			}
				
		}
=cut
	}
	close IN;
	close OUT;	
}


#0x4 fragment unmapped
#0x8 next fragment in the template unmapped

sub totalID
{
	my $input=shift;
        my $output=shift;
	if($$input=~/\.bam$/){
		open IN,"$Bin/samtools view -h $$input|" ||die "no input\n";
	}else{
        	open IN,$$input ||die "no input se\n";
	}
	if($$output=~/\.gz$/){
		open OUT,">:gzip","$$output"  || die "$!" ;	
	}else{
        	open OUT,">$$output" || die "can't creat output\n";
	}
	my $lastline=<IN>;
        chomp $lastline;
        my @array;
        push (@array,$lastline);
	while(<IN>){
		chomp;
		if(/^\@/){
			print OUT $_."\n";
			next;
		}
		my $line=$_;
                my @a=split /\t/,$line;
                my @b=split /\t/,$lastline;
                if($a[0] eq $b[0]){
                        push @array,$line;
                }else{
                        if(@array==1){
				my @tmp==split /\t/,$array[0];
				unless($tmp[5]=~/I/ or $tmp[5]=~/D/){		
                                	print OUT $array[0]."\n";
				}
                        }else{
                                &gettotalline(\@array,\*OUT);
                        }
                        @array=();
                        push(@array, $line);
                }

                $lastline=$line;

        }
        close IN;
        if(@array==1){
		my @tmp=split /\t/,$array[0];
		unless($tmp[1] &0x100){
		unless($tmp[5]=~/I/ or $tmp[5]=~/D/){
                	print OUT $array[0]."\n";
		}
		}
        }
        close OUT;
}


#########################sub_get
sub getpeline
{

	my ($lineI,$lineII,$readname,$readname2,$out)=@_;
	
	unless( $$lineI[1]&0x100 || $$lineII[1]&0x100){
		print $out join ("\t",@$lineI)."\n".join("\t",@$lineII)."\n";
	}

=head	
	#select($out);
	if($wordI eq $wordII){
		
	}else{
		print STDERR "$wordI\t $wordII\t $readname \t $readname2\tthis pair has not same word\n";
	}
=cut
			
		
}

sub getseline
{
	my ($setwo,$out)=@_;
	my @a=split /\t/,$$setwo[0];
	my @b=split /\t/,$$setwo[1];
	if(@{$setwo} != 2){
		print STDERR "something is wrong in ".join("\t",@{$setwo})."\n";
	}elsif(!($a[1]&0x100) && !($b[1]&0x100)){
################################################################ deal same se	
		if($a[4]==$b[4]){  
			print STDERR "$a[0] $b[0]  this read cant match correctly\n";
		}elsif($a[4]==255 ){
			print $out "$$setwo[1]\n";
		}elsif($b[4]==255){
			print $out "$$setwo[0]\n";
		}elsif($a[4] > $b[4]){
			print $out "$$setwo[0]\n";
		}elsif($a[4] < $b[4]){
			print $out "$$setwo[1]\n";
		}				
	}
}

sub gettotalline
{
	my ($setwo,$out)=@_;
	my @a=split /\t/,$$setwo[0];
	my @b=split /\t/,$$setwo[1];
	if(@{$setwo}==2){
		if( !($a[1] & 0x8) && ($b[1] & 0x8)){
			my @tmp=split /\t/, $$setwo[0];
			unless($tmp[5]=~/I|D/){
				print $out $$setwo[0]."\n";
			}	
                }elsif(!($b[1] & 0x8) && ($a[1] & 0x8)){
			my @tmp=split /\t/,$$setwo[1];
			unless($tmp[5]=~/I|D/){
				print $out $$setwo[1]."\n";	
			}
                }else{	
			die "some error in read $a[0] and $b[0] \n";
                }
        }else{
		die "some error in read $a[0] \n";
	}
	
}

sub usage
{

die "perl $0 --mode <pe/se/total> sam outsam\n";

}


