use strict;
die "perl $0 C2T.sam G2A.sam >countread.out\n" unless(@ARGV==2);
my $c2t=shift;
my $g2a=shift;
open CT,$c2t or die "No c2t sam\n";
open GA,$g2a or die "No g2a sam\n";
my $totalreads=0;
my $mapreads=0;
my $uniqmapreads=0;
while(<CT>){
	my $ctline1=$_;
	my $ctline2=<CT>;
	my $galine1=<GA>;
	my $galine2=<GA>;
	chomp ($ctline1,$ctline2,$galine2,$galine1);
	next if(/^@/);
	my @ct1=split /\t/,$ctline1;
	my @ct2=split /\t/,$ctline2;
	my @ga1=split /\t/,$galine1;
	my @ga2=split /\t/,$galine2;
	die "Line's error\n" if($ct1[0] ne $ga1[0]);
	$totalreads+=2;
	next if($ct1[2] eq '*' && $ct2[2] eq '*' && $ga1[2] eq '*' && $ga2[2] eq '*');	
	my $ctmapqulity=$ct1[4]+$ct2[4];
	my $gamapqulity=$ga1[4]+$ga2[4];
	next if($ctmapqulity == $gamapqulity);		
	if($ctmapqulity > $gamapqulity){
		if($ct1[4]>0){
			$mapreads++;
		}
		if($ct2[4]>0){
			$mapreads++;
		}
		if($ct1[4]>30){
                        $uniqmapreads++;
                }
                if($ct2[4]>30){
                        $uniqmapreads++;
                }
	}else{
		if($ga1[4]>0){
                        $mapreads++;
                }
                if($ga2[4]>0){
                        $mapreads++;
                }
                if($ga1[4]>30){
                        $uniqmapreads++;
                }
                if($ga2[4]>30){
                        $uniqmapreads++;
                }

	}	

}
close CT;
close GA;
print "Totalreads\tMapreads\tUniqMapreads\n$totalreads\t$mapreads\t$uniqmapreads\n";
