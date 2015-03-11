#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw/$Bin/;
use Cwd qw/abs_path/;
use threads;
use threads::shared;
use Parallel::ForkManager;

use lib "$FindBin::Bin/lib";
use Parcore::parManager qw/parmanager/;

my ($shell,$zpath,$zt,$zw,$zb,$zc,$zp,$zq,@zfile);
my ($format_time,$sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst);
$zpath = shift @ARGV;
$zt = shift @ARGV;	# Check interval
$zw = shift @ARGV;	# With/Without SGE
$zb = shift @ARGV;	# With/Without Breakpoint
$zc = shift @ARGV;	# Core number
$zp = shift @ARGV;	# SGE project
$zq = shift @ARGV;	# SGE queue
chomp($zpath);

# -s: script file
# -o: output directory
# -m: memory required
# -d: disk required
# -t: check interval
# -w: with/without sge
# -c: cpu core number
# -q: sge queue
# -p: sge project

# Step 8 #
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 8 start at $format_time.\n";
$shell = "$zpath/ASM-8-mergetotal.sh";
die "Step 8 error!" unless !parmanager("-s:$shell","-o:$zpath","-m:3G","-d:50G","-t:$zt","-w:$zw","-b:$zb","-c:$zc","-q:$zq","-p:$zp");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 8 finish at $format_time.\n";
print "#" x 45 . "\n";
# Step 9 #
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 9 start at $format_time.\n";
$shell = "$zpath/ASM-9-filtermapread.sh";
die "Step 9 error!" unless !parmanager("-s:$shell","-o:$zpath","-m:200M","-d:50G","-t:$zt","-w:$zw","-b:$zb","-c:$zc","-q:$zq","-p:$zp");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 9 finish at $format_time.\n";
print "#" x 45 . "\n";
# Step 10 #
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 10 start at $format_time.\n";
$shell = "$zpath/ASM-10-sort.sh";
die "Step 10 error!" unless !parmanager("-s:$shell","-o:$zpath","-m:2G","-d:50G","-t:$zt","-w:$zw","-b:$zb","-c:$zc","-q:$zq","-p:$zp");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 10 finish at $format_time.\n";
print "#" x 45 . "\n";
# Step 11 #
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 11 start at $format_time.\n";
$shell = "$zpath/ASM-11-mergesam.sh";
die "Step 11 error!" unless !parmanager("-s:$shell","-o:$zpath","-m:3G","-d:50G","-t:$zt","-w:$zw","-b:$zb","-c:$zc","-q:$zq","-p:$zp");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 11 finish at $format_time.\n";
print "#" x 45 . "\n";
# Step 12 #
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 12 start at $format_time.\n";
$shell = "$zpath/ASM-12-mpileup.sh";
die "Step 12 error!" unless !parmanager("-s:$shell","-o:$zpath","-m:200M","-d:50G","-t:$zt","-w:$zw","-b:$zb","-c:$zc","-q:$zq","-p:$zp");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 12 finish at $format_time.\n";
print "#" x 45 . "\n";
# Step 13 #
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 13 start at $format_time.\n";
$shell = "$zpath/ASM-13-bcf.sh";
die "Step 13 error!" unless !parmanager("-s:$shell","-o:$zpath","-m:200M","-d:50G","-t:$zt","-w:$zw","-b:$zb","-c:$zc","-q:$zq","-p:$zp");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 13 finish at $format_time.\n";
print "#" x 45 . "\n";
# Step 14 #
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 14 start at $format_time.\n";
$shell = "$zpath/ASM-14-cat.sh";
die "Step 14 error!" unless !parmanager("-s:$shell","-o:$zpath","-m:200M","-d:50G","-t:$zt","-w:$zw","-b:$zb","-c:$zc","-q:$zq","-p:$zp");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 14 finish at $format_time.\n";
print "#" x 45 . "\n";
# Step 15 #
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 15 start at $format_time.\n";
$shell = "$zpath/ASM-15-get_snp.sh";
die "Step 15 error!" unless !parmanager("-s:$shell","-o:$zpath","-m:200M","-d:50G","-t:$zt","-w:$zw","-b:$zb","-c:$zc","-q:$zq","-p:$zp");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 15 finish at $format_time.\n";
print "#" x 45 . "\n";
# Step 16 #
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 16 start at $format_time.\n";
$shell = "$zpath/ASM-16-Filtersnp.sh";
die "Step 16 error!" unless !parmanager("-s:$shell","-o:$zpath","-m:200M","-d:50G","-t:$zt","-w:$zw","-b:$zb","-c:$zc","-q:$zq","-p:$zp");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 16 finish at $format_time.\n";
print "#" x 45 . "\n";
# Step 17 #
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 17 start at $format_time.\n";
$shell = "$zpath/ASM-17-ASM.sh";
die "Step 17 error!" unless !parmanager("-s:$shell","-o:$zpath","-m:3G","-d:50G","-t:$zt","-w:$zw","-b:$zb","-c:$zc","-q:$zq","-p:$zp");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 17 finish at $format_time.\n";
print "#" x 45 . "\n";
# Step 18 #
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 18 start at $format_time.\n";
$shell = "$zpath/ASM-18-ASM_fil.sh";
die "Step 18 error!" unless !parmanager("-s:$shell","-o:$zpath","-m:200M","-d:50G","-t:$zt","-w:$zw","-b:$zb","-c:$zc","-q:$zq","-p:$zp");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 18 finish at $format_time.\n";
print "#" x 45 . "\n";