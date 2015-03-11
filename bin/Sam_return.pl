#!/usr/bin/perl;
use strict;

use Getopt::Long;
use PerlIO::gzip;
use File::Basename qw(basename dirname);
use FindBin '$Bin';
my $mode;
my $Help;
my $type;
GetOptions(
        "mode:s"=>\$mode,
	"type:s"=>\$type,
        "help"=>\$Help
);

&usage() if( !$mode || @ARGV!=4 || $Help);
$type ||="pe";

sub usage{
	die "perl $0  --mode <W|C> sam fq soappe soapse\n" ;
}
my $sam=shift;
my $fq=shift;
my $pe=shift;
my $se=shift;

open SAM,$sam or die "no samfile\n";
if($fq=~/\.gz$/){  
#FCB01Y3ACXX:8:1101:2369:2194#ACAGTGAT/1 TGGGTTGGTTTTTGTTGCGTGAGGAGTGTTGTTTTTGTTGTTTTTTTTGTTGA   ___ecceegggggiiihiighfhidgafghhigiiihiiiiiiiiige[ad`]   FCB01Y3ACXX:8:1101:2369:2194#ACAGTGAT/2 CAACAAAAAAAACAACAAAAACAACACTCCTCACGCAACAAAAACCAACCCAA   bbbeeeeeggggghihiiiii_cgfgfggiihihghfgfhhiigfgdddcbbc
	open FQ, "gunzip -c $fq|" or die "no fq\n";
}else{
	open FQ,  $fq or die "no fq\n";
}
open PE,">$pe" or die "cant create pe\n";
open SE,">$se" or die "cant create se\n";


if($type=~/pe/){
while(<SAM>){
	chomp;
	if(/^\@/){
		print PE "$_\n";
		print SE "$_\n";
		next;
	}
#FCC076MACXX:1:1101:1247:2319#TGACCAAT   83      chr1    201898265       2       60M     =       201897879       -446    TTATTACCCAAACTAAAATACAATAACATAATCTCAACTCACTACAACCTCTACCTCCCA    hgebeaihiggehhhhfaaeeg_fhhgghhffgagebffggfdfhgfee^e`Yccaeb__    AS:i:-24        XS:i:0  XN:i:0  XM:i:4  XO:i:0  XG:i:0  NM:i:4  MD:Z:0C9T23T23T1        YS:i:-6 YT:Z:CP
	my $read1=$_;
	my $read2=<SAM>;
	chomp $read2;
	my @read1=split /\t/,$read1;
        my @read2=split /\t/,$read2;
	
	my $fqpair=<FQ>;
	next if(($read1[4]<30 || $read1[4]==255) && ($read2[4]<30 || $read2[4]==255));
	#next if(!($read1[1] & 0x2) && !($read1[1] & 0x4) && !($read1[1] & 0x8));
	chomp $fqpair;
	my @fqp=split /\t/,$fqpair;
	
    if($mode eq 'W'){	
	if($fqp[0] eq "$read1[0]\/1"){		
		if($read1[1] & 0x40){
			$read1[9]=$fqp[1];
			$read2[9]=$fqp[4];
			$read1[0]=$fqp[0];	
			$read2[0]=$fqp[3];
			$read2[9]=reverse ($read2[9]);
			$read2[9]=~tr/ACGTacgt/TGCAtgca/;
			
			if(!(($read1[1] & 0x4) || ($read1[1] & 0x8)) ){
	                        if(!($read1[1] & 0x10) && ($read2[1] & 0x10)){
				  if(!($read1[1] & 0x100) && !($read2[1] & 0x100)){	
                               		print PE join ("\t",@read1)."\n".join ("\t",@read2)."\n";
				  }elsif(!($read1[1] & 0x100) && ($read2[1] & 0x100)){
				  	print SE join ("\t",@read1)."\n";
				  }elsif(($read1[1] & 0x100) && !($read2[1] & 0x100)){
					print SE join ("\t",@read2)."\n";	
				  }
        	                }
               		}

               		if(($read1[1] & 0x8) && !($read1[1] & 0x4) && !($read1[1] & 0x10) && !($read1[1] & 0x100)){
                        	print SE join ("\t",@read1)."\n";
                	}if(($read2[1] & 0x8) && !($read2[1] & 0x4) && ($read2[1] & 0x10) && !($read2[1] & 0x100)){
                        	print SE join ("\t",@read2)."\n";
                	}
			
						
		}else{
			$read1[9]=$fqp[4];
                        $read2[9]=$fqp[1];
                        $read1[0]=$fqp[3];
                        $read2[0]=$fqp[0];
			$read1[9]=reverse ($read1[9]);
                        $read1[9]=~tr/ACGTacgt/TGCAtgca/;
			if(!(($read1[1] & 0x4) || ($read1[1] & 0x8)) ){
                                if(!($read2[1] & 0x10) && ($read1[1] & 0x10)){
				  if(!($read1[1] & 0x100) && !($read2[1] & 0x100)){
                                        print PE join ("\t",@read1)."\n".join ("\t",@read2)."\n";
                                  }elsif(!($read1[1] & 0x100) && ($read2[1] & 0x100)){
                                        print SE join ("\t",@read1)."\n";
                                  }elsif(($read1[1] & 0x100) && !($read2[1] & 0x100)){
                                        print SE join ("\t",@read2)."\n";
                                  }
                                }
                        }

                        if(($read1[1] & 0x8) && !($read1[1] & 0x4) && ($read1[1] & 0x10) && !($read1[1] & 0x100)){
                                print SE join ("\t",@read1)."\n";
                        }
			if(($read2[1] & 0x8) && !($read2[1] & 0x4) && !($read2[1] & 0x10) && !($read2[1] & 0x100)){
                                print SE join ("\t",@read2)."\n";
                        }						
		}

	}else{
		die "some erroin $fqpair\n$read1[0]\t$read2[0]\n";
	
	}
    }
    if($mode eq 'C'){
	if($fqp[0] eq "$read1[0]\/1"){

                if($read1[1] & 0x40){
                        $read1[9]=$fqp[1];
                        $read2[9]=$fqp[4];
                        $read1[0]=$fqp[0];
                        $read2[0]=$fqp[3];
			$read1[9]=reverse ($read1[9]);
                        $read1[9]=~tr/ACGTacgt/TGCAtgca/;
			
			
			if(!(($read1[1] & 0x4) || ($read1[1] & 0x8))){
	                        if(($read1[1] & 0x10) && !($read2[1] & 0x10)){
				  if(!($read1[1] & 0x100) && !($read2[1] & 0x100)){
                                        print PE join ("\t",@read1)."\n".join ("\t",@read2)."\n";
                                  }elsif(!($read1[1] & 0x100) && ($read2[1] & 0x100)){
                                        print SE join ("\t",@read1)."\n";
                                  }elsif(($read1[1] & 0x100) && !($read2[1] & 0x100)){
                                        print SE join ("\t",@read2)."\n";
                                  }
				}	
                	}

               		if(($read1[1] & 0x8) && !($read1[1] & 0x4) && ($read1[1] & 0x10) && !($read1[1] & 0x100)){
                        	print SE join ("\t",@read1)."\n";
                	}if(($read1[1] & 0x4) && !($read1[1] & 0x8) && !($read2[1] & 0x10) && !($read2[1] & 0x100)){
                        	print SE join ("\t",@read2)."\n";
                	}

                }else{
                        $read1[9]=$fqp[4];
                        $read2[9]=$fqp[1];
                        $read1[0]=$fqp[3];
                        $read2[0]=$fqp[0];
			$read2[9]=reverse ($read2[9]);
                        $read2[9]=~tr/ACGTacgt/TGCAtgca/;

			if(!(($read1[1] & 0x4) || ($read1[1] & 0x8))){
                                if(($read2[1] & 0x10) && !($read1[1] & 0x10)){
				  if(!($read1[1] & 0x100) && !($read2[1] & 0x100)){
                                        print PE join ("\t",@read1)."\n".join ("\t",@read2)."\n";
                                  }elsif(!($read1[1] & 0x100) && ($read2[1] & 0x100)){
                                        print SE join ("\t",@read1)."\n";
                                  }elsif(($read1[1] & 0x100) && !($read2[1] & 0x100)){
                                        print SE join ("\t",@read2)."\n";
                                  }
                                }
                        }

                        if(($read1[1] & 0x8) && !($read1[1] & 0x4) && !($read1[1] & 0x10) && !($read1[1] & 0x100)){
                                print SE join ("\t",@read1)."\n";
                        }if(($read1[1] & 0x4) && !($read1[1] & 0x8) && ($read2[1] & 0x10) && !($read2[1] & 0x100)){
                                print SE join ("\t",@read2)."\n";
                        }

		
                }
        }else{
                die "some erroin $fqpair\n";

        }


    }
	
}
close SAM;

}


if($type=~/se/){
while(<SAM>){
        chomp;
        if(/^\@/){
                print PE "$_\n";
                print SE "$_\n";
                next;
        }
#FCC076MACXX:1:1101:1247:2319#TGACCAAT   83      chr1    201898265       2       60M     =       201897879       -446    TTATTACCCAAACTAAAATACAATAACATAATCTCAACTCACTACAACCTCTACCTCCCA    hgebeaihiggehhhhfaaeeg_fhhgghhffgagebffggfdfhgfee^e`Yccaeb__    AS:i:-24        XS:i:0  XN:i:0  XM:i:4  XO:i:0  XG:i:0  NM:i:4  MD:Z:0C9T23T23T1        YS:i:-6 YT:Z:CP
        my $read1=$_;
        my @read1=split /\t/,$read1;

        my $fqpair=<FQ>;
        next if($read1[4]<30 || $read1[4]==255);
        #next if(!($read1[1] & 0x2) && !($read1[1] & 0x4) && !($read1[1] & 0x8));
        chomp $fqpair;
        my @fqp=split /\t/,$fqpair;

	if($mode eq 'W'){
        if($fqp[0] eq "$read1[0]"){
                if(!($read1[1] & 0x10) && !($read1[1] & 0x4) && !($read1[1] & 0x100)){
                        $read1[9]=$fqp[1];
                        $read1[0]=$fqp[0];
			print SE join ("\t",@read1)."\n";
##########################################################################3
			
                }
        }else{
                die "some erroin $fqpair\n$read1[0]\n";

        }
    }
#####################################
    if($mode eq 'C'){
        if($fqp[0] eq "$read1[0]"){
                if($read1[1] & 0x10 && !($read1[1] & 0x4) && !($read1[1] & 0x100) ){
                        $read1[9]=$fqp[1];
                        $read1[0]=$fqp[0];
                        $read1[9]=reverse ($read1[9]);
                        $read1[9]=~tr/ACGTacgt/TGCAtgca/;
			print SE join ("\t",@read1)."\n";
		}
        }else{
                die "some erroin $fqpair\n";

        }
    }

}
close SAM;
}
