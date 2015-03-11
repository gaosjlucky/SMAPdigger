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

my ($shell,$ipath,$tpath,$npath,$type,$normal,$zt,$zw,$zb,$zc,$zp,$zq,@zfile);
my ($format_time,$sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst);
$ipath = shift @ARGV;
$type = shift @ARGV;
$normal = shift @ARGV;
$zt = shift @ARGV;	# Check interval
$zw = shift @ARGV;	# With/Without SGE
$zb = shift @ARGV;	# With/Without Breakpoint
$zc = shift @ARGV;	# Core number
$zp = shift @ARGV;	# SGE project
$zq = shift @ARGV;	# SGE queue
chomp($ipath);
$tpath = $ipath . $type;
$npath = $ipath . $normal;
# -s: script file
# -o: output directory
# -m: memory required
# -d: disk required
# -t: check interval
# -w: with/without sge
# -c: cpu core number
# -q: sge queue
# -p: sge project

# Step 17 #
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 17 start at $format_time.\n";
$shell = "$tpath/DMR-17-getdmr.sh";
die "Step 17 error!" unless !parmanager("-s:$shell","-o:$tpath","-m:200M","-d:50G","-t:$zt","-w:$zw","-b:$zb","-c:$zc","-q:$zq","-p:$zp");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 17 finish at $format_time.\n";
print "#" x 45 . "\n";
# Step 18 #
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 18 start at $format_time.\n";
$shell = "$tpath/DMR-18-draw.sh";
die "Step 18 error!" unless !parmanager("-s:$shell","-o:$tpath","-m:200M","-d:50G","-t:$zt","-w:$zw","-b:$zb","-c:$zc","-q:$zq","-p:$zp");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 18 finish at $format_time.\n";
print "#" x 45 . "\n";
# Step 19 #
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 19 start at $format_time.\n";
$shell = "$tpath/DMR-19-anno.sh";
die "Step 19 error!" unless !parmanager("-s:$shell","-o:$tpath","-m:200M","-d:50G","-t:$zt","-w:$zw","-b:$zb","-c:$zc","-q:$zq","-p:$zp");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 19 finish at $format_time.\n";
print "#" x 45 . "\n";