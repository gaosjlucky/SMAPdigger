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

my($format_time,$sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst);	
my($mode,$bpt,$shell,$outdir,$memsz,$disksz,$chkint,$nthread,$sgequeue,$sgeproj,$slurmpart);
GetOptions(
	"-mode:s"		=> \$mode,
	"-bpt:s"		=> \$bpt,
	"-outdir:s"		=> \$outdir,
	"-chkint:s"		=> \$chkint,
	"-nthread:s"	=> \$nthread,
	"-sgequeue:s"	=> \$sgequeue,
	"-sgeproj:s"	=> \$sgeproj,
	"-slurmpart:s"	=> \$slurmpart
);

my $parDir = "$outdir/Partshell/Parcore";
my $stepId = 0;
# Step 1 ############################################
$stepId = 1;
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step $stepId start at $format_time.\n";
$shell = "$parDir/BS_$stepId.sh";
$memsz = "200M";
$disksz = "10G";
die "Step $stepId error!" unless !parmanager("-mode:$mode","-bpt:$bpt","-shell:$shell","-outdir:$outdir","-memsz:$memsz","-disksz:$disksz","-chkint:$chkint","-nthread:$nthread","-sgequeue:$sgequeue","-sgeproj:$sgeproj","-slurmpart:$slurmpart");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step $stepId finish at $format_time.\n";
print "#" x 45 . "\n";
# Step 2 ############################################
$stepId = 2;
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step $stepId start at $format_time.\n";
$shell = "$parDir/BS_$stepId.sh";
$memsz = "10G";
$disksz = "10G";
die "Step $stepId error!" unless !parmanager("-mode:$mode","-bpt:$bpt","-shell:$shell","-outdir:$outdir","-memsz:$memsz","-disksz:$disksz","-chkint:$chkint","-nthread:$nthread","-sgequeue:$sgequeue","-sgeproj:$sgeproj","-slurmpart:$slurmpart");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step $stepId finish at $format_time.\n";
print "#" x 45 . "\n";
# Step 3 ############################################
$stepId = 3;
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step $stepId start at $format_time.\n";
$shell = "$parDir/BS_$stepId.sh";
$memsz = "1G";
$disksz = "10G";
die "Step $stepId error!" unless !parmanager("-mode:$mode","-bpt:$bpt","-shell:$shell","-outdir:$outdir","-memsz:$memsz","-disksz:$disksz","-chkint:$chkint","-nthread:$nthread","-sgequeue:$sgequeue","-sgeproj:$sgeproj","-slurmpart:$slurmpart");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step $stepId finish at $format_time.\n";
print "#" x 45 . "\n";
# Step 4 ############################################
$stepId = 4;
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step $stepId start at $format_time.\n";
$shell = "$parDir/BS_$stepId.sh";
$memsz = "2G";
$disksz = "10G";
die "Step $stepId error!" unless !parmanager("-mode:$mode","-bpt:$bpt","-shell:$shell","-outdir:$outdir","-memsz:$memsz","-disksz:$disksz","-chkint:$chkint","-nthread:$nthread","-sgequeue:$sgequeue","-sgeproj:$sgeproj","-slurmpart:$slurmpart");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step $stepId finish at $format_time.\n";
print "#" x 45 . "\n";
# Step 5 ############################################
$stepId = 5;
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step $stepId start at $format_time.\n";
$shell = "$parDir/BS_$stepId.sh";
$memsz = "1G";
$disksz = "10G";
die "Step $stepId error!" unless !parmanager("-mode:$mode","-bpt:$bpt","-shell:$shell","-outdir:$outdir","-memsz:$memsz","-disksz:$disksz","-chkint:$chkint","-nthread:$nthread","-sgequeue:$sgequeue","-sgeproj:$sgeproj","-slurmpart:$slurmpart");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step $stepId finish at $format_time.\n";
print "#" x 45 . "\n";
# Step 6 ############################################
$stepId = 6;
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step $stepId start at $format_time.\n";
$shell = "$parDir/BS_$stepId.sh";
$memsz = "20G";
$disksz = "10G";
die "Step $stepId error!" unless !parmanager("-mode:$mode","-bpt:$bpt","-shell:$shell","-outdir:$outdir","-memsz:$memsz","-disksz:$disksz","-chkint:$chkint","-nthread:$nthread","-sgequeue:$sgequeue","-sgeproj:$sgeproj","-slurmpart:$slurmpart");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step $stepId finish at $format_time.\n";
print "#" x 45 . "\n";
# Step 7 ############################################
$stepId = 7;
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step $stepId start at $format_time.\n";
$shell = "$parDir/BS_$stepId.sh";
$memsz = "4G";
$disksz = "10G";
die "Step $stepId error!" unless !parmanager("-mode:$mode","-bpt:$bpt","-shell:$shell","-outdir:$outdir","-memsz:$memsz","-disksz:$disksz","-chkint:$chkint","-nthread:$nthread","-sgequeue:$sgequeue","-sgeproj:$sgeproj","-slurmpart:$slurmpart");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step $stepId finish at $format_time.\n";
print "#" x 45 . "\n";
# Step 8 ############################################
$stepId = 8;
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step $stepId start at $format_time.\n";
$shell = "$parDir/BS_$stepId.sh";
$memsz = "1G";
$disksz = "10G";
die "Step $stepId error!" unless !parmanager("-mode:$mode","-bpt:$bpt","-shell:$shell","-outdir:$outdir","-memsz:$memsz","-disksz:$disksz","-chkint:$chkint","-nthread:$nthread","-sgequeue:$sgequeue","-sgeproj:$sgeproj","-slurmpart:$slurmpart");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step $stepId finish at $format_time.\n";
print "#" x 45 . "\n";
# Step 9 ############################################
$stepId = 9;
print "#" x 45 . "\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step $stepId start at $format_time.\n";
$shell = "$parDir/BS_$stepId.sh";
$memsz = "200M";
$disksz = "10G";
die "Step $stepId error!" unless !parmanager("-mode:$mode","-bpt:$bpt","-shell:$shell","-outdir:$outdir","-memsz:$memsz","-disksz:$disksz","-chkint:$chkint","-nthread:$nthread","-sgequeue:$sgequeue","-sgeproj:$sgeproj","-slurmpart:$slurmpart");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time());
$format_time = sprintf("%d-%d-%d %d:%d:%d", $year+1990, $mon+1, $mday, $hour, $min, $sec);
print "Step $stepId finish at $format_time.\n";
print "#" x 45 . "\n";



