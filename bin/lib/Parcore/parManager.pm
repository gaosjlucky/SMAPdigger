=head1 NAME
Parcore::parManager - A simple parallel processing fork manager for multicore CPU
=cut

package Parcore::parManager;
use strict;
use warnings;
use File::Basename qw/basename dirname/;
use FindBin qw/$Script $RealBin/;
use Cwd qw/abs_path/;
use threads;
use threads::shared;
use Parallel::ForkManager;
require Exporter;

our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw(parmanager);

sub parmanager {
	####################################################################
	# Argument parsing
	####################################################################
	# {{
	my %arg_hash;
	&parmanager_arg(\%arg_hash, @_);
	#################################
	# sbatch -p test -N 1 --output=/xxx/test.log -D /xxx/ /xxx/test.sh
	# -x: slurm partition (-p)
	# Switch
	my $mode 		= 	$arg_hash{'mode'};		# multicore, sge, slurm
	my $bpt 		= 	$arg_hash{'bpt'};		# on: without breakpoint, off: with breakpoint
	# Basic arg
	my $shell		= 	$arg_hash{'shell'};		# shell script
	my $outdir   	= 	$arg_hash{'outdir'};	# output directory
	my $memsz 		= 	$arg_hash{'memsz'};		# memory size requirement
	my $disksz 		= 	$arg_hash{'disksz'};	# disk size requirement
	my $chkint 		= 	$arg_hash{'chkint'};	# check interval (resource & state check)
	my $nthread 	= 	$arg_hash{'nthread'};	# thread number
	# For sge mode
	my $sgequeue 	= 	$arg_hash{'sgequeue'};	# sge queue
	my $sgeproj 	= 	$arg_hash{'sgeproj'};	# sge project
	# For slurm env
	my $slurmpart 	= 	$arg_hash{'slurmpart'};	# slurm parition
	#################################
	if(!(-e $outdir)) {
		system("mkdir -p $outdir");
	}
	# }}
	
	####################################################################
	# Shell command splitting
	####################################################################
	# {{
	my ($FIA,$FOA,$bufA,$errFile);
	my $shPath = dirname $shell;
	# Execute state shell & log 
	my $logDir = "$shell.tmp";
	# Clear history record
	if($bpt eq 'off') {
		system("rm -rf $logDir");
	}
	# Make new log directory
	if(!(-d $logDir)) {
		system("mkdir -p $logDir");
	}
	# Command split
	my @cmd;
	open (SHS,"<$shell") || die "fail $shell: $!\n";
	while(<SHS>){
		# Each lines in the shell script could be parallelly executed
		next if(/^\s+$/ || /^\#/);
		s/\s+$//; 
		chomp;
		push @cmd, $_;
	}
	close SHS;
	# Save to shell scripts
	my $shCmdCnt = @cmd;
	my @shCmdFin; 
	my $shCmdLeft = $shCmdCnt;
	my $shCmdOk = "This-work-is-completed";
	# Init
	for(my $i = 0; $i < $shCmdCnt; $i++) {
		$shCmdFin[$i] = 0;
	}
	# Re-check
	if($bpt eq 'on') {
		for(my $i = 0; $i < $shCmdCnt; $i++) {
			# Check execute log
			$errFile = "$logDir/step$i.sh.z.log";
			if(-e $errFile) {
				open $FIA,"<$errFile";
				$bufA = <$FIA>;
				chomp($bufA);
				close $FIA;
				if($bufA eq $shCmdOk) {
					$shCmdFin[$i] = 1;
					$shCmdLeft--;
				}
			}
		}
	}
	# Generate script
	for(my $i = 0; $i < $shCmdCnt; $i++) {	
		if($shCmdFin[$i] == 0) {
			$errFile = "$logDir/step$i.sh.z.log";
			open $FOA,">$logDir/step$i.sh";
			print $FOA "#!/bin/sh\n$cmd[$i] && perl -e 'print STDERR \"$shCmdOk\"' >&$errFile";
			close $FOA;
			system("chmod +x $logDir/step$i.sh");
		}
	}
	# }}
	
	####################################################################
	# Shell command execute
	####################################################################
	# {{
	my ($cmdLine,$diskAvail,$memAvail,$pm,$tag);
	if($mode eq 'multicore') {
		$pm = new Parallel::ForkManager($nthread);
	}
	for(my $i = 0; $i < $shCmdCnt; $i++) {
		if($shCmdFin[$i] == 1) {
			next;
		}
		###############################
		# Check disk resources
		###############################
		# {{
		for(my $j = 0; $j < 30; $j++) {
			$tag = 0;
			$diskAvail = &chk_disk($outdir,0);
			if(&compare_space($diskAvail, $disksz) eq 'NA'){ 
				print "Not enought disk space for $outdir. $disksz needed, $diskAvail available.\n";
				print "Retry in $chkint seconds.";
				sleep($chkint);
				next;
			}
			else {
				$tag = 1;
				last;
			}
		}
		
		if($tag == 0) {
			die "Not enought disk space for $outdir. $disksz needed, $diskAvail available. Exit.\n";
		}
		# }}
		
		###############################
		# Check local memory resources
		# For multicore mode only
		###############################
		# {{
		if($mode eq 'multicore') {
			for(my $j = 0; $j < 30; $j++) {
				$tag = 0;
				$memAvail = &chk_mem(0);
				$memAvail .= 'K';
				if(&compare_space($memAvail, $memsz) eq 'NA'){ 
					print "Not enought memory space. $memsz needed, $memAvail available.\n";
					print "Retry in $chkint seconds.";
					sleep($chkint);
					next;
				}
				else {
					$tag = 1;
					last;
				}
			}
			
			if($tag == 0) {
				die "Not enought memory space. $memsz needed, $memAvail available. Exit.\n";
			}
		}
		# }}
		
		###############################
		# Work dispatch
		###############################
		# {{
		if($mode eq 'multicore') {
			# Load commands from script file
			open $FIA,"<$logDir/step$i.sh";
			$cmdLine = <$FIA>;
			chomp($cmdLine);
			close $FIA;
			# Dispatch
			$pm->start and next;
			system("$cmdLine");
			$pm->finish;
		}
		elsif($mode eq 'sge') {
			# Submit task on SGE 
			my $SGELOG;
			open $SGELOG,">$logDir/step$i.log";
			print $SGELOG "qsub -o $logDir -e $logDir -q $sgequeue -P $sgeproj -l vf=$memsz $logDir/step$i.sh";
			close $SGELOG;
			system("qsub -o $logDir -e $logDir -q $sgequeue -P $sgeproj -l vf=$memsz $logDir/step$i.sh");
		}
		elsif($mode eq 'slurm') {
			# Convert mem size to MB
			my $slurmMem = &transform_space($memsz);
			$slurmMem = $slurmMem * 1024;	
			# Submit task on SLURM
			system("sbatch -p $slurmpart -N 1 --mem=$slurmMem -D $logDir $logDir/step$i.sh");
		}
		# }}
	}
	
	####################################################################
	# Work collection
	####################################################################
	# {{
	my $totalSec = 0;
	if($mode eq 'multicore') {
		# For multicore mode
		$pm->wait_all_children;
		sleep(2);
		for(my $i = 0; $i < $shCmdCnt; $i++) {
			if($shCmdFin[$i] == 1) {
				next;
			}
			else {
				# Check execute log
				$errFile = "$logDir/step$i.sh.z.log";
				if(-e $errFile) {
					open $FIA,"<$errFile";
					$bufA = <$FIA>;
					chomp($bufA);
					if($bufA eq $shCmdOk) {
						$shCmdFin[$i] = 1;
						$shCmdLeft--;
					}
					close $FIA;
				}
			}
		}
	}
	else {
		# For sge and slurm mode
		while($shCmdLeft > 0 && $totalSec <= 90000) {
			for(my $i = 0; $i < $shCmdCnt; $i++) {
				if($shCmdFin[$i] == 1) {
					next;
				}
				else {
					# Check execute log
					$errFile = "$logDir/step$i.sh.z.log";
					if(-e $errFile) {
						open $FIA,"<$errFile";
						$bufA = <$FIA>;
						chomp($bufA);
						if($bufA eq $shCmdOk) {
							$shCmdFin[$i] = 1;
							$shCmdLeft--;
						}
						close $FIA;
					}
				}
			}
			$totalSec = $totalSec + $chkint;
			sleep($chkint);
		}
	}
	# }}
	
	# Success
	if($shCmdLeft == 0) {
		return 0;
	}
	else {
		return -1;
	}
}

#------- deal the para of routine jobguard ---------
sub parmanager_arg {
	my $phash = shift;
	# -s: shell script file
	# -m: min free memory space 
	# -d: output directory 
	# -t: resource check interval 
	# -f: min free disk space
	# -p: max thread number
	foreach my $para (@_) {
		my ($sign,$value) = ($para =~ /^\s*\-+([^\s:]+)\s*[:=]\s*(\S+)\s*$/);
		$$phash{$sign} = $value; 
	}
}

#------ check the free disk space -------
sub chk_disk{
	my ($disk,$chkHeader) = @_[0,1];
	my ($i,$avail,$df_info,@df_info_array);
	my $time = 0;
	
	DF: 
	{
		$df_info = `df -h $disk`;
		$time++;
		
		if($chkHeader == 1) {
			# Needed ONLY when the system language is ENGLISH!
			if($df_info !~ /Filesystem.+Size.+Used.+Avail/i){
				die "Cannot get the disk $disk space info!\n" if($time == 3);
				redo DF;
			}
		}
		
		@df_info_array = split(/\n+/,$df_info);
		shift @df_info_array;	# Remove header
		$df_info = join('', @df_info_array);		
		$df_info =~ s/^\s+//g;
		$df_info =~ s/\s+$//g;
		$avail = (split /\s+/,$df_info)[3];
		if($avail !~ /\d/) {
			die "Command resolved. But cannot get the numerical space info of disk $disk!\n" if($time == 3);
			redo DF;
		}
	}
	
	return $avail;
}

#------ check the free disk space -------
sub chk_mem{
	my $chkHeader = $_[0];
	my $time = 0;
	my ($avail,$free_info,@free_info_array);
	
	FREE: 
	{
		$free_info = `free`;
		$time++;
		
		if($chkHeader == 1) {
			# Needed ONLY when the system language is ENGLISH!
			die "Cannot get the memory space info!\n" if($time == 3);
			redo FREE;
		}
		
		@free_info_array = split(/\n+/,$free_info);
		shift @free_info_array;	# Remove header
		$free_info = join('', @free_info_array);		
		$free_info =~ s/^\s+//g;
		$free_info =~ s/\s+$//g;
		$avail = (split /\s+/,$free_info)[3];

		if($avail !~ /\d/) {
			die "Cannot get the avail space info of memory!\n" if($time == 3);
			redo FREE;
		}
	}
	return $avail;
}

#------ compare the spaces -------
sub compare_space{
	my ($want_bigger,$want_smaller) = @_[0,1];
	$want_bigger = &transform_space($want_bigger);
	$want_smaller = &transform_space($want_smaller);
	if($want_bigger >= $want_smaller){
		return 'OK';
	}
	else{
		return 'NA';
	}
}

#------ transform space based on the unit --------
sub transform_space{
	my ($old) = $_[0];
	my ($num, $unit) = ($old =~ /^([\d\.]+)(.+)?$/);
	$unit ||= 'k'; ## default
	my (%multiple, @multiple);
	push @multiple,(1024**($_)) for (-2 .. 1);
	@multiple{qw/T t G g M m K k/} = @multiple[3,3,2,2,1,1,0,0];
	
	if(exists($multiple{$unit})){
		$num *= $multiple{$unit};
	}
	else{
		die "Cannot distinguish the unit $unit from $old!\n";
	}
	
	# Convert to size with GB unit
	return $num;
}

1;
