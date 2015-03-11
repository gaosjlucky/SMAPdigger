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

# Step1 #
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 1 start at $format_time.\n";
$shell = "$zpath/RRBS-1-dealadaptor.sh";
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
$shell = "$zpath/RRBS-2-bsmap.sh";
die "Step 2 error!" unless !parmanager("-s:$shell","-o:$zpath","-m:10G","-d:50G","-t:$zt","-w:$zw","-b:$zb","-c:$zc","-q:$zq","-p:$zp");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 2 finish at $format_time.\n";
print "#" x 45 . "\n";
# Step3 #
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 3 start at $format_time.\n";
$shell = "$zpath/RRBS-3-addread.sh";
die "Step 3 error!" unless !parmanager("-s:$shell","-o:$zpath","-m:3G","-d:50G","-t:$zt","-w:$zw","-b:$zb","-c:$zc","-q:$zq","-p:$zp");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 3 finish at $format_time.\n";
print "#" x 45 . "\n";
