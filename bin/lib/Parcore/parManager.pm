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
	# -s: script file
	# -o: output directory
	# -m: memory required
	# -d: disk required
	# -t: check interval
	# -w: with/without sge
	# -b: breakpoint off/on
	# -c: cpu core number
	# -q: sge queue
	# -p: sge project
	# Common arg
	my $shScript = $arg_hash{'s'};		# shell script
	my $outputDir = $arg_hash{'o'};		# output directory
	my $memSpace = $arg_hash{'m'};		# memory requirement
	my $diskSpace = $arg_hash{'d'};		# disk requirement
	my $chkInterval = $arg_hash{'t'};	# check interval (resource & state check)
	# Switch
	my $withSge = $arg_hash{'w'};		# 0: without SGE, 1: with SGE
	my $withBpt = $arg_hash{'b'};		# 0: without breakpoint, 1: with breakpoint
	# For multi-core env	
	my $threadNum = $arg_hash{'c'};		# max process number
	# For sge env
	my $queueName = $arg_hash{'q'};		# sge queue name
	my $projectName = $arg_hash{'p'};	# sge project name
	####################
	if(!(-e $outputDir)) {
		system("mkdir -p $outputDir");
	}
	# }}
	
	####################################################################
	# Shell command splitting
	####################################################################
	# {{
	my ($FIA,$FOA,$bufA,$errFile);
	my $shPath = dirname $shScript;
	# Execute state shell & log 
	system("mkdir -p $shScript.tmp");
	# Command split
	my @cmd;
	open (SHS,"<$shScript") || die "fail $shScript: $!\n";
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
	if($withBpt == 1) {
		for(my $i = 0; $i < $shCmdCnt; $i++) {
			# Check execute log
			$errFile = "$shScript.tmp/step$i.sh.z.log";
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
			$errFile = "$shScript.tmp/step$i.sh.z.log";
			open $FOA,">$shScript.tmp/step$i.sh";
			print $FOA "$cmd[$i] && perl -e 'print STDERR \"$shCmdOk\"' >&$errFile";
			close $FOA;
			system("chmod +x $shScript.tmp/step$i.sh");
		}
	}
	# }}
	
	####################################################################
	# Shell command execute
	####################################################################
	# {{
	my ($cmdLine,$diskAvail,$memAvail,$pm,$tag);
	if($withSge == 0) {
		$pm = new Parallel::ForkManager($threadNum);
	}
	for(my $i = 0; $i < $shCmdCnt; $i++) {
		if($shCmdFin[$i] == 1) {
			next;
		}
		###############################
		# Resource check
		###############################
		# {{
		if($withSge == 0) {
			for(my $j = 0; $j < 20; $j++) {
				$tag = 1;
				$diskAvail = &chk_disk($outputDir,0);
				if(&compare_space($diskAvail,$diskSpace) eq 'NA'){ 
					print "Not enought disk space for $outputDir. $diskSpace needed, $diskAvail available.\n";
					print "Retry in $chkInterval seconds.";
					$tag = 2;
					sleep($chkInterval);
					next;
				}
				$memAvail = &chk_mem(0);
				$memAvail .= 'K';
				if(&compare_space($memAvail,$memSpace) eq 'NA'){ 
					print "Not enought memory space. $memSpace needed, $memAvail available.\n";
					print "Retry in $chkInterval seconds.";
					$tag = 3;
					sleep($chkInterval);
					next;
				}
			}
			
			if($tag == 2) {
				die "Not enought disk space for $outputDir. $diskSpace needed, $diskAvail available. Exit.\n";
			}	
			elsif($tag == 3) {
				die "Not enought memory space. $memSpace needed, $memAvail available. Exit.\n";
			}
		}
		# }}
		
		###############################
		# Work dispatch
		###############################
		# {{
		if($withSge == 0) {
			# Load from script file
			open $FIA,"<$shScript.tmp/step$i.sh";
			$cmdLine = <$FIA>;
			chomp($cmdLine);
			close $FIA;
			# Dispatch
			$pm->start and next;
			system("$cmdLine");
			$pm->finish;
		}
		else {
			system("qsub -o $shScript.tmp -e $shScript.tmp -q $queueName -P $projectName -l vf=$memSpace $shScript.tmp/step$i.sh");
		}
		# }}
	}
	
	####################################################################
	# Work collection
	####################################################################
	# {{
	if($withSge == 0) {
		$pm->wait_all_children;
		for(my $i = 0; $i < $shCmdCnt; $i++) {
			if($shCmdFin[$i] == 1) {
				next;
			}
			else {
				# Check execute log
				$errFile = "$shScript.tmp/step$i.sh.z.log";
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
		while($shCmdLeft > 0) {
			for(my $i = 0; $i < $shCmdCnt; $i++) {
				if($shCmdFin[$i] == 1) {
					next;
				}
				else {
					# Check execute log
					$errFile = "$shScript.tmp/step$i.sh.z.log";
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
			sleep($chkInterval);
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
	my ($num,$unit) = ($old =~ /^([\d\.]+)(.+)?$/);
	$unit ||= 'k'; ## default
	my (%multiple,@multiple);
	push @multiple,(1024**($_)) for (-2 .. 1);
	@multiple{qw/T t G g M m K k/} = @multiple[3,3,2,2,1,1,0,0];
	if(exists($multiple{$unit})){
		$num *= $multiple{$unit};
	}
	else{
		die "Cannot distinguish the unit $unit from $old!\n";
	}
	return $num;
}

1;
