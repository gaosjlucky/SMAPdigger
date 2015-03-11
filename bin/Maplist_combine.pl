use strict;
use warnings;
use PerlIO::gzip;
use File::Basename qw(basename dirname);
die "perl $0 methymap blank\n" if(@ARGV!=2);
my $methy=shift;
my $blank=shift;
my $outdir=dirname($methy);
#my $blankdir=dirname($blank);
my $namebase=basename $methy;
if($methy=~/\.gz$/){
	open MT,"gunzip -c $methy|" or die "no methy_filter file\n";
}else{
	open MT,$methy or die "no methy_filter file\n";
}

open OUT, ">:gzip", "$outdir/statistic.$namebase" or die $!;

open BK,$blank or die "no blank file\n";
my %hash;my %head;
my $flag;

my @tmp2;my @tmp;
my $hd;
while(<BK>){
	chomp;
	next unless(/\w/);
	my @blank=split;
	$hash{$blank[0]}="$blank[1]\t$blank[2]";
}
close BK;
my $chr=$namebase;
$chr=~s/\.gz//;
my @chrlen=split /\s+/,$hash{$chr};
print OUT "\>$chr\n";

for(my $j=0;$j<=$chrlen[1];$j+=1000000){
    my @region;
    my $jend;
    if($chrlen[1]>$j+999999){
        $jend=$j+999999;
    }else{
        $jend=$chrlen[1];
    }
    foreach my $k($j..$jend){
        my $line="N\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t";
        push @region,$line;
    }
    if(@tmp2>0){#@tmp2 is former end map;
        for (my $i=0;$i<@tmp2;$i++){
            if($tmp2[$i][1]>=$j && $tmp2[$i][1]<=$jend){
                &count(\@{$tmp2[$i]},\@region,$j,$jend,\@tmp);
                shift @tmp2;
            }            
        }
        @tmp=();
    }
    while(<MT>){
        chomp;
        next unless(/\w+/);
        my @map=split /\t/;
        $map[0]=~s/W$//; $map[0]=~s/C$//;$map[0]=~s/H$//;
        
        if(exists($hash{$map[0]})){
            if($map[1]>$jend){
                push @tmp2,\@map;
                last;
            }elsif( $j<=$map[1] && $map[1]<= $jend){
                &count(\@map,\@region,$j,$jend,\@tmp2);
            }
        }else{
            die "no this type chr in blank file\n";	
        }
    }
    print OUT join ("\n",@region)."\n";
}
close MT;

sub pecount
{
	my $mp=shift;
	my $region=shift;
	my $start=shift;
	my $end=shift;
	my $tmp=shift;
	my @map=@{$mp};
	if(@map>0){
        	unless($map[1]=~/\d+/){
            	#print @map."\n";
           	 print join ("\t",@map)."\n";
            	die;
        	}
     	 }
    
    if($map[1] >= $start && $map[2]<=$end ){
	for my $i($map[1]..$map[2]){
	    my $pos=$i-$start;
	    my @matrix=split /\t/,$$region[$pos-1];
            if($map[3] eq 'p'){
                my ($maptype1,$maptype2,$ref1,$ref2,$quli1,$quli2,$lenth1,$lenth2);
                if($map[4]=~/\:len(\d+)\:/){
                    $lenth1=$1;
                }
                if($map[11]=~/\:len(\d+)\:/){
                    $lenth2=$1;
                }
                next if(($i>$map[6]+$lenth1-1 && $i< $map[13]) || ($i>$map[13]+$lenth2-1 && $i < $map[6]) );
                
                if($i>=$map[6] && $i<=$map[6]+$lenth1-1 && $i>=$map[13] &&  $i<=$map[13] +$lenth2-1){
                    my $pos1=$i-$map[6]+1;
                    my $pos2=$i-$map[13]+1;
                    $maptype1=substr($map[8],$pos1-1,1);
                    $ref1=substr($map[9],$pos1-1,1);
                    $quli1=substr($map[10],$pos1-1,1);
                    $maptype2=substr($map[15],$pos2-1,1);
                    $ref2=substr($map[16],$pos2-1,1);
                    $quli2=substr($map[17],$pos2-1,1);
                    #rep1:+:len90:W:score
		
                    if($map[4]=~/\:rep1\:[+|-]\:len(\d+)\:(\w)\:score/){
                        my $type=$2;
                        if($maptype1 eq 'm' ){
                            $matrix[1]++;
                        }
                        if($maptype1 eq '#' && $type eq 'C') {
                            $matrix[2]++;
                        }
                        if($maptype1 eq '*' && $type eq 'W'){
                            $matrix[2]++;
                        }
                        if($maptype1 eq '|' ){
                            $matrix[3]++;
                        }
                        if($maptype1 eq ' ' && $maptype2 eq ' '  ){
                            $matrix[4]++;
                        }
                    }elsif($map[4]=~/\:rep(\d+)\:[+|-]\:len(\d+)\:(\w)\:score/){
                        $matrix[9]+=$1;
                        my $type=$3;
                        if($maptype1 eq 'm' ){
                            $matrix[5]++;
                        }
                        if($maptype1 eq '#' && $type eq 'C'){                            $matrix[6]++;
                        }
                        if($maptype1 eq '*' && $type eq 'W'){                                $matrix[6]++;
                        }
                        if($maptype1 eq '|'  ){
                            $matrix[7]++;
                        }
                        if($maptype1 eq " " && $maptype2 eq " "){
                            $matrix[8]++;
                        }
                    }
                    #################################
                    if($map[11]=~/\:rep1\:[+|-]\:len(\d+)\:(\w)\:score/){
                        my $type=$2;
                        if($maptype2 eq 'm' ){
                            $matrix[1]++;
                        }
                        if($maptype2 eq '#' && $type eq 'C'){
                            $matrix[2]++;
                        }
                        if($maptype2 eq '*' && $type eq 'W'){
                            $matrix[2]++;
                        }
                        if($maptype2 eq '|' ){
                            $matrix[3]++;
                        }
                        if($maptype1 eq ' ' && $maptype2 eq ' '  ){
                            
                        }
                    }elsif($map[11]=~/\:rep(\d+)\:[+|-]\:len(\d+)\:(\w)\:score/){
                        $matrix[9]+=$1;
                        my $type=$3;
                        if($maptype2 eq 'm' ){
                            $matrix[5]++;
                        }
                        if($maptype2 eq '#' &&  $type eq 'C' ){
                            $matrix[6]++;
                        }
                        if($maptype2 eq '*' &&  $type eq 'W' ){
                            $matrix[6]++;
                        }
                        if($maptype2 eq '|'  ){
                            $matrix[7]++;
                        }
                        if($maptype1 eq ' ' && $maptype2 eq ' '){
                        }
                    }
                    ###################################
                }else{
                    if($i>=$map[6] && $i<=$map[6]+$lenth1-1){
                        my $pos1=$i-$map[6]+1;
                        $maptype1=substr($map[8],$pos1-1,1);
                        $ref1=substr($map[9],$pos1-1,1);
                        $quli1=substr($map[10],$pos1-1,1);
                        if($map[4]=~/\:rep1\:[+|-]\:len(\d+)\:(\w)\:score/){
                            my $type=$2;
                            if($maptype1 eq 'm' ){
                                $matrix[1]++;
                            }
                            if($maptype1 eq '#' && $type eq 'C'){
                                $matrix[2]++;
                            }
                            if($maptype1 eq '*' && $type eq 'W'){
                                $matrix[2]++;
                            }
                            if($maptype1 eq '|' ){
                                $matrix[3]++;
                            }
                            if($maptype1 eq ' '  ){
                                $matrix[4]++;
                            }
                        }elsif($map[4]=~/\:rep(\d+)\:[+|-]\:len(\d+)\:(\w)\:score/){
                            $matrix[9]+=$1;
                            my $type=$3;
                            if($maptype1 eq 'm' ){
                                $matrix[5]++;
                            }
                            if($maptype1 eq '#' && $type eq 'C'){
                                $matrix[6]++;
                            }
                            if($maptype1 eq '*' && $type eq 'W'){
                                $matrix[6]++;
                            }
                            if($maptype1 eq '|'  ){
                                $matrix[7]++;
                            }
                            if($maptype1 eq " " ){
                                $matrix[8]++;
                            }
                        }
                        
                    }elsif($i>=$map[13] && $i<=$map[13]+$lenth2-1){
                        my $pos2=$i-$map[13]+1;
                        $maptype2=substr($map[15],$pos2-1,1);
                        $ref2=substr($map[16],$pos2-1,1);
                        $quli2=substr($map[17],$pos2-1,1);
                        if($map[11]=~/\:rep1\:[+|-]\:len(\d+)\:(\w)\:score/){
                            my $type=$2;
                            if($maptype2 eq 'm' ){
                                $matrix[1]++;
                            }
                            if($maptype2 eq '#' && $type eq 'C') {
                                $matrix[2]++;
                            }
                            if($maptype2 eq '*' && $type eq 'W') {
                                $matrix[2]++;
                            }
                            if($maptype2 eq '|' ){
                                $matrix[3]++;
                            }
                            if($maptype2 eq ' '  ){
                                $matrix[4]++;
                            }
                        }elsif($map[11]=~/\:rep(\d+)\:[+|-]\:len(\d+)\:(\w)\:score/){
                            $matrix[9]+=$1;
                            my $type=$3;
                            if($maptype2 eq 'm' ){
                                $matrix[5]++;
                            }
                            if($maptype2 eq '#' && $type eq 'C' ){
                                $matrix[6]++;
                            }
                            if($maptype2 eq '*' && $type eq 'W' ){
                                $matrix[6]++;
                            }
                            if($maptype2 eq '|'  ){
                                $matrix[7]++;
                            }
                            if($maptype2 eq ' '){
                                $matrix[8]++;
                            }
                        }
                        
                    }else{
                        print STDERR "$i in $map[0] have some error\n";
						
                    }
                    
                }
            }else{
		die $!;	
	    }
            $$region[$pos-1]=join("\t",@matrix);
        }
    }else{
        die "no in between $start and $end \n";
    }
  
}
sub secount
{
    my $mp=shift;
	my $region=shift;
	my $start=shift;
	my $end=shift;
	my $tmp=shift;
	my @map=@{$mp};
	if(@map>0){
        unless($map[1]=~/\d+/){
            print join ("\t",@map)."\n";
            die;
        }
    }
    
    for my $i($map[1]..$map[2]){
        my $pos=$i-$start;
        my @matrix=split /\t/,$$region[$pos-1];
        my $spos=$i-$map[1];
        my $maptype=substr($map[8],$spos,1);
        my $ref=substr($map[9],$spos,1);
        my $quli=substr($map[10],$spos,1);
        $matrix[0]=$ref;
        if($map[4]=~/\:rep1\:[+|-]\:len(\d+)\:(\w)\:score/){
            my $type=$2;
            if($maptype eq 'm' ){
                $matrix[1]++;
            }
            if($maptype eq '#'  && $type eq 'C'){
                $matrix[2]++
            }
            if($maptype eq '*'  && $type eq 'W'){
                $matrix[2]++
            }
            if($maptype eq '|' ){
                $matrix[3]++
            }
            if($maptype eq ' '  ){
                $matrix[4]++
            }
        }elsif($map[4]=~/\:rep(\d+)\:[+|-]\:len(\d+)\:(\w)\:score/){
            $matrix[9]+=$1;
            my $type=$3;
            if($maptype eq 'm' ){
                $matrix[5]++;
            }
            if($maptype eq '#' && $type eq 'C' ){
                $matrix[6]++;
            }
            if($maptype eq '*' && $type eq 'W' ){
                $matrix[6]++;
            }
            if($maptype eq '|'  ){
                $matrix[7]++;
            }
            if($maptype eq ' '){
                $matrix[8]++;
            }
        }
        $$region[$pos-1]=join("\t",@matrix);
    }
}

sub count
{
	my $mp=shift;
	my $region=shift;
	my $start=shift;
	my $end=shift;
	my $tmp=shift;
	my @map=@{$mp};
	if(@map>0){
	 unless($map[1]=~/\d+/){
		#print @map."\n";
		print join ("\t",@map)."\n";
		die;
	 }
	}
		
   if($map[1] >= $start && $map[2]<=$end ){
        if($map[3] eq 'p'){
            & pecount(\@{$mp},\@{$region},$start,$end,$tmp);
        }
        if($map[3] eq 's'){
            & secount(\@{$mp},\@{$region},$start,$end,$tmp);
        }
        
    }elsif($map[1]<=$end && $map[2]>$end){#transport the region end
	if($map[3] eq 's'){ #se read
		my @inmap=@map;
		my @outmap=@map;
		my $inlen=$end-$map[1]+1;
		$inmap[8]=substr($map[8],0,$inlen);
		my $outlen=$map[2]-$end;
		if($map[4]=~/^(\S+\:len)(\d+)(\:\w)/){
			$inmap[4]="$1$inlen$3";	
			$outmap[4]="$1$outlen$3";		
		}	
			
		$inmap[2]=$end;
		$inmap[7]=substr($map[7],0,$inlen);
		$inmap[8]=substr($map[8],0,$inlen);
		$inmap[9]=substr($map[9],0,$inlen);
		$inmap[10]=substr($map[10],0,$inlen);
		$outmap[1]=$end+1;		
		$outmap[6]=$end+1;
		$outmap[7]=substr($map[7],$inlen,$outlen);
		$outmap[8]=substr($map[8],$inlen,$outlen);
		$outmap[9]=substr($map[9],$inlen,$outlen);
		$outmap[10]=substr($map[10],$inlen,$outlen);
		&secount(\@inmap,\@{$region},$start,$end,\@{$tmp});	#repair need supplementary
		push @{$tmp},\@outmap;
	}

	
	if($map[3] eq 'p'){ #pe reads; turn to this step change pe to 2 se
		my $len1;my $len2;
		if($map[4]=~/^(\S+\:len)(\d+)(\:\w)/){
		  		$len1=$2;
           	}
		if($map[11]=~/^(\S+\:len)(\d+)(\:\w)/){
                	$len2=$2;
            	}

	    	my $seend1=$map[6]+$len1-1;
		my $seend2=$map[13]+$len2-1;
            	my @pe2se1=split /\t/,"$map[0]\t$map[6]\t$seend1\ts\t$map[4]\t$map[5]\t$map[6]\t$map[7]\t$map[8]\t$map[9]\t$map[10]";
            	my @pe2se2=split /\t/,"$map[0]\t$map[13]\t$seend2\ts\t$map[11]\t$map[12]\t$map[13]\t$map[14]\t$map[15]\t$map[16]\t$map[17]";
            
            
            
	    	if($map[6]>$seend2 || $map[13] >$seend1){##ok no overlap pe change to 2se
                	my @pe2se1=split /\t/,"$map[0]\t$map[6]\t$seend1\ts\t$map[4]\t$map[5]\t$map[6]\t$map[7]\t$map[8]\t$map[9]\t$map[10]";
                	my @pe2se2=split /\t/,"$map[0]\t$map[13]\t$seend2\ts\t$map[11]\t$map[12]\t$map[13]\t$map[14]\t$map[15]\t$map[16]\t$map[17]";
                	&count(\@pe2se1,\@{$region},$start,$end,\@{$tmp});
                	&count(\@pe2se2,\@{$region},$start,$end,\@{$tmp});
					
		 }elsif($map[6]>=$map[13] && $seend1<=$seend2){#ok #read2 covers read1
                	my @pe2se2=split /\t/,"$map[0]\t$map[13]\t$seend2\ts\t$map[11]\t$map[12]\t$map[13]\t$map[14]\t$map[15]\t$map[16]\t$map[17]";
                &secount(\@pe2se2,\@{$region},$start,$end,\@{$tmp});
           	 }elsif($map[6]<=$map[13] && $seend1>=$seend2){#ok #read1 covers read2
                	my @pe2se1=split /\t/,"$map[0]\t$map[6]\t$seend1\ts\t$map[4]\t$map[5]\t$map[6]\t$map[7]\t$map[8]\t$map[9]\t$map[10]";
                	&secount(\@pe2se1,\@{$region},$start,$end,\@{$tmp});
            	}else{#overlap
		   if($map[6]<=$end && $map[13]<=$end ){#one
                    if($seend1>=$end && $seend2>=$end){ #two reads through endpoint
			my $inlen1=$end-$map[6]+1;
			my $inmap1=substr($map[8],0,$inlen1);
			my $inquli1=substr($map[10],0,$inlen1);
			my $inref1=substr($map[9],0,$inlen1);
			my $inlen2=$end-$map[13]+1;					
                   	my $inmap2=substr($map[15],0,$inlen2);
                    	my $inquli2=substr($map[17],0,$inlen2);
                    	my $inref2=substr($map[16],0,$inlen2);
				
			my $outlen1=$seend1-$end;
			my $outmap1=substr($map[8],$inlen1,$outlen1);
			my $outquli1=substr($map[10],$inlen1,$outlen1);
			my $outref1=substr($map[9],$inlen1,$outlen1);					
                    	my $outlen2=$seend2-$end;
                    	my $outmap2=substr($map[15],$inlen2,$outlen2);
                   	my $outquli2=substr($map[17],$inlen2,$outlen2);
                    	my $outref2=substr($map[16],$inlen2,$outlen2);
                    	my $outstart=$end+1;
					
			my $inname1;my $outname1;my $inname2;my $outname2;
					
			if($map[4]=~/^(\S+\:len)(\d+)(\:\w)/){
                        	$inname1="$1$inlen1$3";
                        	$outname1="$1$outlen1$3";
                   	 }
			if($map[11]=~/^(\S+\:len)(\d+)(\:\w)/){
                        	$inname2="$1$inlen2$3";
                        	$outname2="$1$outlen2$3";
                   	 }
			my @inmapt=split /\t/,"$map[0]\t$map[1]\t$end\tp\t$inname1\t$map[5]\t$map[6]\t$map[7]\t$inmap1\t$inref1\t$inquli1\t$inname2\t$map[12]\t$map[13]\t$map[14]\t$inmap2\t$inref2\t$inquli2";
			my @outmapt=split /\t/,"$map[0]\t$outstart\t$map[2]\tp\t$outname1\t$map[5]\t$outstart\t$map[7]\t$outmap1\t$outref1\t$outquli1\t$outname2\t$map[12]\t$outstart\t$map[14]\t$outmap2\t$outref2\t$outquli2";
			&pecount(\@inmapt,\@{$region},$start,$end,\@{$tmp});
			push @{$tmp},\@outmapt;	

                    }
                    if($seend1 <=$end && $seend2>$end ){#2read2 through endpoint
			my $inlen2=$end-$map[13]+1;
                 	my $inmap2=substr($map[15],0,$inlen2);
                    	my $inquli2=substr($map[17],0,$inlen2);
                    	my $inref2=substr($map[16],0,$inlen2);		
			my $outlen2=$seend2-$end;
                    	my $outmap2=substr($map[15],$inlen2,$outlen2);
                    	my $outquli2=substr($map[17],$inlen2,$outlen2);
                    	my $outref2=substr($map[16],$inlen2,$outlen2);
                    	my $outstart=$end+1;
                    	my $inname2;my $outname2;
####################################################################################################################################
                    if($map[11]=~/^(\S+\:len)(\d+)(\:\w)/){
                          $inname2="$1$inlen2$3";
                           $outname2="$1$outlen2$3";
                     }
		    my @inmapt=split /\t/,"$map[0]\t$map[1]\t$end\tp\t$map[4]\t$map[5]\t$map[6]\t$map[7]\t$map[8]\t$map[9]\t$map[10]\t$inname2\t$map[12]\t$map[13]\t$map[14]\t$inmap2\t$inref2\t$inquli2";
		    my @outmapt=split /\t/,"$map[0]\t$outstart\t$map[2]\ts\t$outname2\t$map[5]\t$outstart\t$map[7]\t$outmap2\t$outref2\t$outquli2";
		    &pecount(\@inmapt,\@{$region},$start,$end,\@{$tmp});
                    push @{$tmp},\@outmapt;
                    }
                    if($seend2<=$end && $seend1>$end){#3 read1 through endpoint
			my $inlen1=$end-$map[6]+1;
                    	my $inmap1=substr($map[8],0,$inlen1);
                    	my $inquli1=substr($map[10],0,$inlen1);
                    	my $inref1=substr($map[9],0,$inlen1);
			my $outlen1=$seend1-$end;
                    	my $outmap1=substr($map[8],$inlen1,$outlen1);
                    	my $outquli1=substr($map[10],$inlen1,$outlen1);
                    	my $outref1=substr($map[9],$inlen1,$outlen1);
			my $outstart=$end+1;
			my $inname1;my $outname1;
			if($map[4]=~/^(\S+\:len)(\d+)(\:\w)/){
                            $inname1="$1$inlen1$3";
                            $outname1="$1$outlen1$3";
                    	}
					
			my @inmapt=split /\t/,"$map[0]\t$map[1]\t$end\tp\t$inname1\t$map[5]\t$map[6]\t$map[7]\t$inmap1\t$inref1\t$inquli1\t$map[11]\t$map[12]\t$map[13]\t$map[14]\t$map[15]\t$map[16]\t$map[17]";
                    my @outmapt=split /\t/,"$map[0]\t$outstart\t$map[2]\ts\t$outname1\t$map[12]\t$outstart\t$map[13]\t$outmap1\t$outref1\t$outquli1";
			&pecount(\@inmapt,\@{$region},$start,$end,\@{$tmp});
                	push @{$tmp},\@outmapt;
                    }
                }
#############################################################################

                if($map[6]<=$end && $map[13]>$end){#two ,read1 through endpoint
			my $inlen1=$end-$map[6]+1;
                	my $inmap1=substr($map[8],0,$inlen1);
                	my $inquli1=substr($map[10],0,$inlen1);
                	my $inref1=substr($map[9],0,$inlen1);
			my $outlen1=$seend1-$end;
                	my $outmap1=substr($map[8],$inlen1,$outlen1);
                	my $outquli1=substr($map[10],$inlen1,$outlen1);
                	my $outref1=substr($map[9],$inlen1,$outlen1);
			my $outstart=$end+1;
			my ($inname1,$outname1);
			if($map[4]=~/^(\S+\:len)(\d+)(\:\w)/){
                        	$inname1="$1$inlen1$3";
                        	$outname1="$1$outlen1$3";
               		}
			my @inmapt=split /\t/,"$map[0]\t$map[1]\t$end\ts\t$inname1\t$map[5]\t$map[6]\t$map[7]\t$inmap1\t$inref1\t$inquli1";
			my @outmap=split /\t/,"$map[0]\t$outstart\t$map[2]\tp\t$outname1\t$map[5]\t$outstart\t$map[7]\t$outmap1\t$outref1\t$outquli1\t$map[11]\t$map[12]\t$map[13]\t$map[14]\t$map[15]\t$map[16]\t$map[17]";
			&secount(\@inmapt,\@{$region},$start,$end,\@{$tmp});
                	push @{$tmp},\@outmap;	
                }
                if($map[6]>$end && $map[13]<=$end){#three, read2 through endpoint
                	my $inlen2=$end-$map[13]+1;
                	my $inmap2=substr($map[15],0,$inlen2);
                	my $inquli2=substr($map[17],0,$inlen2);
                	my $inref2=substr($map[16],0,$inlen2);
                	my $outlen2=$seend2-$end;
                	my $outmap2=substr($map[15],$inlen2,$outlen2);
                	my $outquli2=substr($map[17],$inlen2,$outlen2);
                	my $outref2=substr($map[16],$inlen2,$outlen2);
                	my $outstart=$end+1;
                	my $inname2;my $outname2;
                	if($map[11]=~/^(\S+\:len)(\d+)(\:\w)/){
                    		$inname2="$1$inlen2$3";
                   		$outname2="$1$outlen2$3";
               		}
			my @inmapt=split /\t/,"$map[0]\t$map[1]\t$end\ts\t$inname2\t$map[12]\t$map[13]\t$map[14]\t$inmap2\t$inref2\t$inquli2";
			my @outmapt=split /\t/,"$map[0]\t$outstart\t$map[2]\tp\t$map[4]\t$map[5]\t$map[6]\t$map[7]\t$map[8]\t$map[9]\t$map[10]\t$outname2\t$map[12]\t$outstart\t$map[14]\t$outmap2\t$outref2\t$outquli2";	
			&secount(\@inmapt,\@{$region},$start,$end,\@{$tmp});
                	push @{$tmp},\@outmapt;	
		}
         }
	}	
     }elsif($map[1]>$end){
		push @tmp2,\@map;
		#print STDERR   "some erro in ".join("\t",@map)."\t$map[0]\t$map[1]\t$map[2]\t$end \n";
     }
	
}


