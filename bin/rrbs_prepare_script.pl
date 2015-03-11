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

my ($shell,$zpath,$zt,$zw,$zb,$zc,$zp,$zq,$zi,$zf,@zfile);
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

# Step0 #
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 0 start at $format_time.\n";
$shell = "$zpath/RRBS-0-dealadptor.sh";
die "Step 0 error!" unless !parmanager("-s:$shell","-o:$zpath","-m:200M","-d:50G","-t:$zt","-w:$zw","-b:$zb","-c:$zc","-q:$zq","-p:$zp");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 0 finish at $format_time.\n";
print "#" x 45 . "\n";
# Step1 #
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 1 start at $format_time.\n";
$shell = "$zpath/RRBS-1-dealfq.sh";
die "Step 1 error!" unless !parmanager("-s:$shell","-o:$zpath","-m:200M","-d:50G","-t:$zt","-w:$zw","-b:$zb","-c:$zc","-q:$zq","-p:$zp");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 1 finish at $format_time.\n";
print "#" x 45 . "\n";
# Step2 #
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 2 start at $format_time.\n";
$shell = "$zpath/RRBS-2-align.sh";
die "Step 2 error!" unless !parmanager("-s:$shell","-o:$zpath","-m:5G","-d:50G","-t:$zt","-w:$zw","-b:$zb","-c:$zc","-q:$zq","-p:$zp");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 2 finish at $format_time.\n";
print "#" x 45 . "\n";
# Step3 #
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 3 start at $format_time.\n";
$shell = "$zpath/RRBS-3-sort.sh";
die "Step 3 error!" unless !parmanager("-s:$shell","-o:$zpath","-m:2G","-d:50G","-t:$zt","-w:$zw","-b:$zb","-c:$zc","-q:$zq","-p:$zp");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 3 finish at $format_time.\n";
print "#" x 45 . "\n";
# Step4 #
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 4 start at $format_time.\n";
$shell = "$zpath/RRBS-4-sam_return.sh";
die "Step 4 error!" unless !parmanager("-s:$shell","-o:$zpath","-m:200M","-d:50G","-t:$zt","-w:$zw","-b:$zb","-c:$zc","-q:$zq","-p:$zp");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 4 finish at $format_time.\n";
print "#" x 45 . "\n";
# Step5 #
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 5 start at $format_time.\n";
$shell = "$zpath/RRBS-5-mergerepective.sh";
die "Step 5 error!" unless !parmanager("-s:$shell","-o:$zpath","-m:2G","-d:50G","-t:$zt","-w:$zw","-b:$zb","-c:$zc","-q:$zq","-p:$zp");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 5 finish at $format_time.\n";
print "#" x 45 . "\n";
# Step6 #
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 6 start at $format_time.\n";
$shell = "$zpath/RRBS-6-sam_uniq.sh";
die "Step 6 error!" unless !parmanager("-s:$shell","-o:$zpath","-m:200M","-d:50G","-t:$zt","-w:$zw","-b:$zb","-c:$zc","-q:$zq","-p:$zp");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 6 finish at $format_time.\n";
print "#" x 45 . "\n";
# Step7 #
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 7 start at $format_time.\n";
$shell = "$zpath/RRBS-7-sam2align.sh";
die "Step 7 error!" unless !parmanager("-s:$shell","-o:$zpath","-m:3G","-d:50G","-t:$zt","-w:$zw","-b:$zb","-c:$zc","-q:$zq","-p:$zp");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 7 finish at $format_time.\n";
print "#" x 45 . "\n";
