#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use File::Basename;
my %hash;
my %sample;
die "perl $0 patient.smd .." unless(@ARGV>1);
foreach my $file(@ARGV){
	my @pat=split /\//,$file;
	open $hash{$pat[-3]},$file or die $!;
}

my @file=sort keys %hash;
#Chr     Pos     Strand  Type    Ref     Copy_number     #normal_reads_of_methy  #normal_read_of_unmethy normal_threshhold       #cancer_reads_of_methy  #cance

print "Chr\tPos\tStrand\tType\tRef\tCopy_number";
foreach my $sample(@file){
	print "\t$sample";
	my $pp=$hash{$sample};
	<$pp>;
		
}
print "\n";
my $pt=$hash{$file[0]};

while(<$pt>){
	my %hash2;
	my @line=split /\t/,$_;
	#push @{$hash2{$patient}},$_;
	print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\,$line[7]\,$line[9]\,$line[10]";

	for(my $i=1;$i<@file;$i++){
		my $pt2=$hash{$file[$i]};
		my $line2=<$pt2>;
		chomp $line2;
		my @aa=split /\t/,$line2;
		die "$aa[1] ne $line[1]" if($aa[1] ne $line[1]);
		print "\t$aa[6],$aa[7],$aa[9],$aa[10]";
		#push @{$hash2{$patient2}},$line2;
	}
	print "\n";				
		
}	

for(my $k=0;$k<@file;$k++){
	close $hash{$file[$k]};

}			


	
