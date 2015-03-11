use strict;
my $i=shift;
open IN,$i or die $!;
while(<IN>){
next if($_=~/chr(\w+)_random/ || $_=~/chr(\w+)_hap/);
print $_;
}
close IN;
