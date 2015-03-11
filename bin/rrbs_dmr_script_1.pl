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

# Step 8 #
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 8 start at $format_time.\n";
$shell = "$tpath/DMR-8-Uniq.sh";
die "Step 8 error!" unless !parmanager("-s:$shell","-o:$tpath","-m:200M","-d:50G","-t:$zt","-w:$zw","-b:$zb","-c:$zc","-q:$zq","-p:$zp");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 8 finish at $format_time.\n";
print "#" x 45 . "\n";
# Step 9 #
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 9 start at $format_time.\n";
$shell = "$tpath/DMR-9-sort.sh";
die "Step 9 error!" unless !parmanager("-s:$shell","-o:$tpath","-m:2G","-d:50G","-t:$zt","-w:$zw","-b:$zb","-c:$zc","-q:$zq","-p:$zp");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 9 finish at $format_time.\n";
print "#" x 45 . "\n";
# Step 10 #
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 10 start at $format_time.\n";
$shell = "$tpath/DMR-10-Uniqtotal.sh";
die "Step 10 error!" unless !parmanager("-s:$shell","-o:$tpath","-m:200M","-d:50G","-t:$zt","-w:$zw","-b:$zb","-c:$zc","-q:$zq","-p:$zp");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 10 finish at $format_time.\n";
print "#" x 45 . "\n";
# Step 11 #
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 11 start at $format_time.\n";
$shell = "$tpath/DMR-11-Stat.sh";
die "Step 11 error!" unless !parmanager("-s:$shell","-o:$tpath","-m:3G","-d:50G","-t:$zt","-w:$zw","-b:$zb","-c:$zc","-q:$zq","-p:$zp");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 11 finish at $format_time.\n";
print "#" x 45 . "\n";
# Step 12 #
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 12 start at $format_time.\n";
$shell = "$tpath/DMR-12-posord.sh";
die "Step 12 error!" unless !parmanager("-s:$shell","-o:$tpath","-m:200M","-d:50G","-t:$zt","-w:$zw","-b:$zb","-c:$zc","-q:$zq","-p:$zp");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 12 finish at $format_time.\n";
print "#" x 45 . "\n";
# Step 13 #
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 13 start at $format_time.\n";
$shell = "$tpath/DMR-13-cat.sh";
die "Step 13 error!" unless !parmanager("-s:$shell","-o:$tpath","-m:200M","-d:50G","-t:$zt","-w:$zw","-b:$zb","-c:$zc","-q:$zq","-p:$zp");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 13 finish at $format_time.\n";
print "#" x 45 . "\n";
# Step 14 #
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 14 start at $format_time.\n";
$shell = "$tpath/DMR-14-sortpos.sh";
die "Step 14 error!" unless !parmanager("-s:$shell","-o:$tpath","-m:2G","-d:50G","-t:$zt","-w:$zw","-b:$zb","-c:$zc","-q:$zq","-p:$zp");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 14 finish at $format_time.\n";
print "#" x 45 . "\n";
# Step 15 #
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 15 start at $format_time.\n";
$shell = "$tpath/DMR-15-qmap.sh";
die "Step 15 error!" unless !parmanager("-s:$shell","-o:$tpath","-m:3G","-d:50G","-t:$zt","-w:$zw","-b:$zb","-c:$zc","-q:$zq","-p:$zp");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 15 finish at $format_time.\n";
print "#" x 45 . "\n";
# Step 16 #
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 16 start at $format_time.\n";
$shell = "$tpath/DMR-16-cout.sh";
die "Step 16 error!" unless !parmanager("-s:$shell","-o:$tpath","-m:200M","-d:50G","-t:$zt","-w:$zw","-b:$zb","-c:$zc","-q:$zq","-p:$zp");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step 16 finish at $format_time.\n";
print "#" x 45 . "\n";
