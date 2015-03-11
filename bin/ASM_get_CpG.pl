#author:gaoshengjie
#modify in 2013/1/30
#modify2 in 2013/12/05
#modify3 in 2015/01/14
#email: gaoshengjie@genomics.cn

use strict;
use Cwd qw/abs_path/;
use File::Basename qw/basename dirname/;
use FindBin '$Bin';
die "perl $0 ref snp bam >out 2>snp.info\n" unless(@ARGV==3);
my $ref=shift;
my $snp=shift;
my $bam=shift;
my $outdir=dirname($bam);
my %Ref;
&Getref($ref);
open SNP,$snp or die "no snp\n";
if($bam=~/\.bam$/){
	open BAM,"$Bin/samtools view $bam|" or die "no bam\n";
}else{
	open BAM,$bam or die "on bam\n";	
}
my %snp_cpg;
my %snppos;
my @tmp;
#chr1    846338  .       A       G,X     140     .       DP=65

my $rep;my $threebase;
while(<SNP>){
	chomp;
	my @a=split;
	$a[4]=~s/\,X//g;
	my @basetype=split /\,/,$a[4];
	if(@basetype>2){
		$threebase++;
		next;
	}
	my $Ref=&Substr($a[0],$a[1],1);
	if($Ref=~/[acgt]/){
		$rep++;
		$Ref=uc $Ref;
		#print "$ref\t$a[3]\t$a[4]\n" if($ref ne $a[3]);
		#next;
	}
	if(@basetype==1){		
		$snppos{$a[0]}{$a[1]}="$Ref\t$a[4]\t$a[3]";
	}else{
		if($Ref eq $a[3]){
			$threebase++;
			next;
		}
		elsif($Ref eq $basetype[0]){
			$snppos{$a[0]}{$a[1]}="$Ref\t$basetype[1]\t$a[3]";
		}
		elsif($Ref eq $basetype[1]){
			$snppos{$a[0]}{$a[1]}="$Ref\t$basetype[0]\t$a[3]";
		}else{
			$threebase++;
		}
	}
}
close SNP;

print STDERR "Snp read done\n";
my (%cpg_ref, %cpg_snp);	
my %hash;

my %SNP_TA;
open NUMRD,">$outdir/mapreadnum.out" or die $!;
my $mapnumread=0;
while(<BAM>){
	chomp;
	my @bamline=split;
	next if($bamline[4]<30);
    	$mapnumread++;
	my $lenth=length($bamline[9]);
	my $end=$bamline[3]+$lenth-1;
	foreach my $pos($bamline[3]..$end){
		if(exists($snppos{$bamline[2]}{$pos})){
			my @base=split /\s+/,$snppos{$bamline[2]}{$pos};
                        my $snp_basepos=$pos-$bamline[3];
			my $snpref=$base[0];
                        my $snpbase=substr($bamline[9],$snp_basepos,1);

			if($snpref eq "A"){ #ref is A
				if($base[1] eq "T"){
					if($snpbase eq "T"){
						&GetASM($_,$pos);
					}
					if($snpbase eq "A"){
						&RefASM($_,$pos);
					}	
				}
				if($base[1] eq "C"){
					if($snpbase eq "A"){
						&RefASM($_,$pos);
					}
					if($snpbase eq "C"){
						&GetASM($_,$pos);
					}	
					if($snpbase eq "T" && ($bamline[1] & 0x40) && !($bamline[1] & 0x10) ){##note Watson srand 
						&GetASM($_,$pos);
					}
					if($snpbase eq "T" && ($bamline[1] & 0x80) && ($bamline[1] & 0x10) ){##note Watson srand 
                                                &GetASM($_,$pos);
                                        }  	
				}
				if($base[1] eq "G"){#watson strand
					if($snpbase eq "G" && ($bamline[1] & 0x40) && !($bamline[1] & 0x10) ){
						&GetASM($_,$pos);	
					}
					if($snpbase eq "G" && ($bamline[1] & 0x80) && ($bamline[1] & 0x10) ){
						&GetASM($_,$pos);
					}
					if($snpbase eq "A" && ($bamline[1] & 0x40) && !($bamline[1] & 0x10) ){
                                                &RefASM($_,$pos);    
                                        }   
                                        if($snpbase eq "A" && ($bamline[1] & 0x80) && ($bamline[1] & 0x10) ){
                                                &RefASM($_,$pos);
                                        }	

				}
			}
			
			if($snpref eq "T"){ #ref is T
                                if($base[1] eq "A"){
                                        if($snpbase eq "T"){
                                                &RefASM($_,$pos);
                                        }
                                        if($snpbase eq "A"){
                                                &GetASM($_,$pos);
                                        }
                                }
                                if($base[1] eq "C"){
                                        if($snpbase eq "T" && ($bamline[1] & 0x40) && ($bamline[1] & 0x10) ){
                                                &RefASM($_,$pos);
                                        }
                                        if($snpbase eq "C"){
                                                &GetASM($_,$pos);
                                        }
                                        if($snpbase eq "T" && ($bamline[1] & 0x80) && !($bamline[1] & 0x10)){##note Cricksrand 
                                                &RefASM($_,$pos);
                                        }
                                }
                                if($base[1] eq "G"){#watson strand
                                        if($snpbase eq "A" && ($bamline[1] & 0x40) && ($bamline[1] & 0x10) ){
                                                &GetASM($_,$pos);
                                        }
                                        if($snpbase eq "A" && ($bamline[1] & 0x80) && !($bamline[1] & 0x10)){
                                                &GetASM($_,$pos);
                                        }
                                        if($snpbase eq "T" && ($bamline[1] & 0x40) && !($bamline[1] & 0x10) ){
                                                &RefASM($_,$pos);
                                        }
                                        if($snpbase eq "G"){
                                                &GetASM($_,$pos);
                                        }

                                }
                        }
			if($snpref eq "C"){ #ref is C
                                if($base[1] eq "A"){
                                        if($snpbase eq "T" && ($bamline[1] & 0x40) && !($bamline[1] & 0x10)){
                                                &RefASM($_,$pos);
                                        }
					if($snpbase eq "T" && ($bamline[1] & 0x80) && ($bamline[1] & 0x10) ){
                                                &RefASM($_,$pos);
                                        }
                                        if($snpbase eq "A"){
                                                &GetASM($_,$pos);
                                        }
                                }
                                if($base[1] eq "T"){
                                        if($snpbase eq "T" && ($bamline[1] & 0x40) && ($bamline[1] & 0x10)){
                                                &GetASM($_,$pos);
                                        }
                                        if($snpbase eq "C"){
                                                &RefASM($_,$pos);
                                        }
                                        if($snpbase eq "T" && ($bamline[1] & 0x80) && !($bamline[1] & 0x10)){##note Cricksrand 
                                                &GetASM($_,$pos);
                                        }
                                }
                                if($base[1] eq "G"){#watson strand
                                        if($snpbase eq "A" && ($bamline[1] & 0x40) && ($bamline[1] & 0x10) ){
                                                &GetASM($_,$pos);
                                        }
                                        if($snpbase eq "A" && ($bamline[1] & 0x80) && !($bamline[1] & 0x10)){
                                                &GetASM($_,$pos);
                                        }
                                        if($snpbase eq "C"  ){
                                                &RefASM($_,$pos);
                                        }
                                        if($snpbase eq "G"){
                                                &GetASM($_,$pos);
                                        }
					if($snpbase eq "T" && ($bamline[1] & 0x40) && !($bamline[1] & 0x10)){
						&RefASM($_,$pos);
					}
					if($snpbase eq "T" && ($bamline[1] & 0x80) && ($bamline[1] & 0x10) ){
                                                &RefASM($_,$pos);
                                        }

                                }
                        }
			if($snpref eq "G"){ #ref is G
                                if($base[1] eq "A"){
                                        if($snpbase eq "G"){
                                                &RefASM($_,$pos);
                                        }
                                        if($snpbase eq "A" && ($bamline[1] & 0x40) && !($bamline[1] & 0x10)){
                                                &GetASM($_,$pos);
                                        }
					if($snpbase eq "A" && ($bamline[1] & 0x80) && ($bamline[1] & 0x10) ){
                                                &GetASM($_,$pos);
                                        }
                                }
                                if($base[1] eq "C"){#G>C
                                        if($snpbase eq "T" && ($bamline[1] & 0x40) && !($bamline[1] & 0x10)){
                                                &GetASM($_,$pos);
                                        }
					if($snpbase eq "T" && ($bamline[1] & 0x80) && ($bamline[1] & 0x10) ){
                                                &GetASM($_,$pos);
                                        }
                                        if($snpbase eq "C"){
                                                &GetASM($_,$pos);
                                        }
                                        if($snpbase eq "A" && ($bamline[1] & 0x80) && !($bamline[1] & 0x10)){##note Cricksrand 
                                                &RefASM($_,$pos);
                                        }
					if($snpbase eq "A" && ($bamline[1] & 0x40) && ($bamline[1] & 0x10)){##note Cricksrand 
                                                &RefASM($_,$pos);
                                        }
					if($snpbase eq "G"){
						&RefASM($_,$pos);
					}
                                }
                                if($base[1] eq "T"){# G>T
                                        if($snpbase eq "A" && ($bamline[1] & 0x40) && ($bamline[1] & 0x10) ){
                                                &RefASM($_,$pos);
                                        }
                                        if($snpbase eq "A" && ($bamline[1] & 0x80) && !($bamline[1] & 0x10)){
                                                &RefASM($_,$pos);
                                        }
                                        if($snpbase eq "T" ){
                                                &GetASM($_,$pos);
                                        }
                                        if($snpbase eq "G"){
                                                &RefASM($_,$pos);
                                        }

                                }
                        }

		}		
	}			
}
close BAM;
print NUMRD "$mapnumread\n";
print "SNPpos\tCpGpos\tPosRef\tPosvar\tTransRef\tRefmeth\tRefunmethy\tSNPmeth\tSNPunmethy\n";

foreach my $posinfo(sort keys %hash){
	my @info=split /\t/,$posinfo;
	my @snpinfo=split /\_/,$info[0];
	unless(exists ($cpg_ref{$posinfo}{'Methy'})){
		$cpg_ref{$posinfo}{'Methy'}=0;
	}	
	unless(exists ($cpg_ref{$posinfo}{'Unmethy'})){
                $cpg_ref{$posinfo}{'Unmethy'}=0;
        }
	unless(exists ($cpg_snp{$posinfo}{'Methy'})){
                $cpg_snp{$posinfo}{'Methy'}=0;
        }
	unless(exists ($cpg_snp{$posinfo}{'Unmethy'})){
                $cpg_snp{$posinfo}{'Unmethy'}=0;
        }
	next if($cpg_snp{$posinfo}{'Unmethy'}==0 && $cpg_snp{$posinfo}{'Methy'}==0);	
	my $snpnumber=$cpg_snp{$posinfo}{'Methy'}+$cpg_snp{$posinfo}{'Unmethy'};
	my $refnumber=$cpg_ref{$posinfo}{'Methy'}+$cpg_ref{$posinfo}{'Unmethy'};
	my $snprate=$snpnumber/($refnumber+$snpnumber);
	next if($snpnumber<=3 || $snprate<=0.1);	
	print "$posinfo\t$snppos{$snpinfo[0]}{$snpinfo[1]}\t".$cpg_ref{$posinfo}{'Methy'}."\t$cpg_ref{$posinfo}{'Unmethy'}\t$cpg_snp{$posinfo}{'Methy'}\t$cpg_snp{$posinfo}{'Unmethy'}\n";
}


#############
sub RefASM
{
 my ($line, $pos)=@_;
 my @bamline=split /\s+/,$line;
 my $lenth=length $bamline[9];
 my $refaim=&Substr($bamline[2],$bamline[3],$lenth);
 $refaim = uc($refaim);
  for (my $i=0;$i<$lenth-2;$i++){
  my $refcpg=substr($refaim,$i,2);
  if($refcpg eq "CG"){
   my $cpgpos=$bamline[3]+$i;
   my $base1cpg=substr($bamline[9],$i,1);
   my $base2cpg=substr($bamline[9],$i+1,1);
   my $cpginfo="$bamline[2]\_$pos\t$cpgpos";
   $hash{$cpginfo}++;
   if($base1cpg eq "C" && ($bamline[1] & 0x40 && !($bamline[1] & 0x10))){
    $cpg_ref{$cpginfo}{'Methy'}++;
   }

   if($base1cpg eq "C" && ($bamline[1] & 0x80 && ($bamline[1] & 0x10))) {
    $cpg_ref{$cpginfo}{'Methy'}++;
   }   

   if($base2cpg eq "G" && ($bamline[1] & 0x40 && ($bamline[1] & 0x10)) ){
    $cpg_ref{$cpginfo}{'Methy'}++;
   }
   if($base2cpg eq "G" && ($bamline[1] & 0x80 && !($bamline[1] & 0x10))){
    $cpg_ref{$cpginfo}{'Methy'}++;
   }

   if($base1cpg eq "T" && ($bamline[1] & 0x40 && !($bamline[1] & 0x10))){
    $cpg_ref{$cpginfo}{'Unmethy'}++;
   }
   if($base1cpg eq "T" && ($bamline[1] & 0x80 && ($bamline[1] & 0x10))){
    $cpg_ref{$cpginfo}{'Unmethy'}++
   }

   if($base2cpg eq "A" && ($bamline[1] & 0x40 && ($bamline[1] & 0x10))){
     $cpg_ref{$cpginfo}{'Unmethy'}++;
   }
   
   if($base2cpg eq "A" && ($bamline[1] & 0x80 && !($bamline[1] & 0x10))){
     $cpg_ref{$cpginfo}{'Unmethy'}++;
    }
  }
 }
}

sub GetASM
{
 my ($line,$pos)=@_;
 my @bamline=split /\s+/,$line;
 my $lenth=length $bamline[9];
 my $refaim=&Substr($bamline[2],$bamline[3],$lenth);
 $refaim = uc($refaim);
 for (my $i=0;$i<$lenth-2;$i++){
  my $refcpg=substr($refaim,$i,2);
  if($refcpg eq "CG"){
   my $cpgpos=$bamline[3]+$i;
   my $base1cpg=substr($bamline[9],$i,1);
   my $base2cpg=substr($bamline[9],$i+1,1);
   my $cpginfo="$bamline[2]\_$pos\t$cpgpos";
   $hash{$cpginfo}++;
   if($base1cpg eq "C" && ($bamline[1] & 0x40 && !($bamline[1] & 0x10))){
    $cpg_snp{$cpginfo}{'Methy'}++;
   }
   if($base1cpg eq "C" && ($bamline[1] & 0x80 && ($bamline[1] & 0x10))) {
    $cpg_snp{$cpginfo}{'Methy'}++;
   }
   
   if($base2cpg eq "G" && ($bamline[1] & 0x40 && ($bamline[1] & 0x10))){
    $cpg_snp{$cpginfo}{'Methy'}++;
   }
   if($base2cpg eq "G" && ($bamline[1] & 0x80 && !($bamline[1] & 0x10))){
    $cpg_snp{$cpginfo}{'Methy'}++;
   }  
   if($base1cpg eq "T" && ($bamline[1] & 0x40 && !($bamline[1] & 0x10))){
    $cpg_snp{$cpginfo}{'Unmethy'}++;
   }
   if($base1cpg eq "T" && ($bamline[1] & 0x80 && ($bamline[1] & 0x10))){
    $cpg_snp{$cpginfo}{'Unmethy'}++;
   }  

   if($base2cpg eq "A" && ($bamline[1] & 0x40 && ($bamline[1] & 0x10))){
     $cpg_snp{$cpginfo}{'Unmethy'}++;
   }
   if($base2cpg eq "A" && ($bamline[1] & 0x80 && !($bamline[1] & 0x10))){
     $cpg_snp{$cpginfo}{'Unmethy'}++;
   } 
  }
 }
}

sub Getref
{
        my $refe=shift;
	if($refe=~/\.gz$/){
		open REF,"gunzip -c $refe|" or die "no ref\n";
	}else{
       		open REF,$refe or die "no ref\n";
	}
        my $chr;
        while(<REF>){
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
			$mod_rseq .= lc(substr($rseq,$anchor_pos,$len)); # lower case letter
			$anchor_pos += $len;
		}
		elsif($flag eq 'S'){
			$anchor_pos += $len;
		}
		$cigar =~ s/^\d+\D//;
	}
	return $mod_rseq;
}

