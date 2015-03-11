#! /usr/bin/perl -w
use strict;
use Cwd qw/abs_path/;
use File::Basename qw/basename dirname/;
use FindBin qw/$RealScript $RealBin/;
use Getopt::Long;

die "Usage: perl $0 <binPath> <outputPath> <tableName>\n" if @ARGV != 3;
my($binPath, $outputPath, $tableName) = @ARGV;
my ($out_dir) = $outputPath;
my $abs_out_dir = abs_path($out_dir);

my (@fileList, %hash_dmr_pair);
my ($file, $base, $suffix, $f_f_dir, $f_f_dir_base, $father_dir, $dir, $f_dir, $f_dir_base);
my ($pair_key, $cout_dir, $key, $pair_files, $file_for_header, $base_key);

die "Invalid output path!\n" unless -d $abs_out_dir;

if (-d "$abs_out_dir/report") {
    print "Clear old report files...";
    system("rm -rf $abs_out_dir/report/*");
    print "Completed.\n";
}

system("mkdir -p $abs_out_dir/report/");

&traverse_dir($abs_out_dir, \@fileList);

foreach $file (@fileList) {
    $base       = basename $file;
    $base       =~ /\w+\.(.*)/;
    $suffix     = $1;
    $father_dir = dirname $file;
    $dir        = basename $father_dir;
    $f_dir      = dirname $father_dir;
    $f_dir_base = basename $f_dir;
	
    if($suffix eq "pdf" or $suffix eq "out.pdf") {
     	$f_f_dir      = dirname $f_dir;
     	$f_f_dir_base = basename $f_f_dir;
     	system("cp -fr $file $abs_out_dir/report/$f_f_dir_base\_$f_dir_base\_$base");
    }

    if ($dir eq "ASM" and ($suffix eq "fil.out")) {
        $f_f_dir      = dirname $f_dir;
        $f_f_dir_base = basename $f_f_dir;
		system("cp -fr $file $abs_out_dir/report/$f_f_dir_base\_$f_dir_base\_$base");
    }
    elsif ($dir eq "DMR" and $suffix eq "dmr") {
		next if $file !~ /_chr.*.dmr/;
		$pair_key = $file;
		$pair_key =~ s/_chr.*//;
        push @{$hash_dmr_pair{$pair_key}}, $file;
	}
	else {
	}
}

foreach $key (sort keys %hash_dmr_pair) {
    $father_dir = dirname $key;
    $dir        = basename $father_dir;
    $f_dir      = dirname $father_dir;
    $f_dir_base = basename $f_dir;

    $pair_files = join(" ", @ {$hash_dmr_pair{$key}});
    $file_for_header = ${$hash_dmr_pair{$key}}[0];
    $base_key = basename $key;

    system("head -1 $file_for_header > $abs_out_dir/report/$f_dir_base\_DMR\_.tmp.dmr.header && cat $pair_files | grep -v '#chr' > $abs_out_dir/report/$f_dir_base\_DMR\_$base_key.dmr.tmp && cat $abs_out_dir/report/$f_dir_base\_DMR\_.tmp.dmr.header $abs_out_dir/report/$f_dir_base\_DMR\_$base_key.dmr.tmp > $abs_out_dir/report/$f_dir_base\_DMR\_$base_key.dmr && rm -f $abs_out_dir/report/$f_dir_base\_DMR\_$base_key.dmr.tmp && rm -f $abs_out_dir/report/$f_dir_base\_DMR\_.tmp.dmr.header");
}


##################################################################################################
# Generate table
##################################################################################################

$tableName = "$abs_out_dir/report/$tableName";

open OUT, ">$tableName" or die "Unable to open table file!\n";

my @basic = <$outputPath/*/*/*/*/adaptor.info>;
my @file= <$outputPath/*/*/qmap/mapnum.out>;
my @cout_file = <$outputPath/*/*/cout_*-*/chr*.cout.cov>;

my %hash;
my ($l, $r);
foreach my $k (@cout_file) {
    ( $l, $r ) = ( $1, $2 ) if ( $k =~ /cout_(\d+)-(\d+)\// );
}
my ( $totalbase, $totalread, $cleanbase,$cleanread,$mapread, $rateofmap, $mapreadintarget,$rateintarget, $total_target_region,$covertargetnum,$rateofcover ) = ( "# of Totalbase", "# of Totalread", "# of Cleanbase", "# of Cleanread", "# of mapped reads", "Rate of mapping reads", "# of mapped reads in enzyme target region", "Rate of mapped reads in enzyme target region", "# of total enzyme $l-$r target regions", "# of covered enzyme $l-$r target regions","Rate of covered enzyme $l-$r target regions");

################## get adaptor.info  @a1   #####################
foreach my $fq (@basic) {
    my $name1 = basename( dirname( dirname( dirname( dirname $fq) ) ) );
    my $name2 = basename( dirname( dirname( dirname $fq) ) );
    my $name="$name1\-$name2";
    open IN, "$fq" || die "no in adaptor.info in $! \n";
    while (<IN>) {
        chomp;
	next unless(/\w/);
	my @head = split;
	my $line=<IN>;
	chomp $line;
        my @b = split /\s+/,$line;
        
        $hash{$name}{ $totalread } += $b[0];
	$hash{$name}{$totalbase}+=$b[1];
	$hash{$name}{$cleanread}+=$b[2];
	$hash{$name}{$cleanbase}+=$b[3];
	
        
    }
    close IN;
}




################## get all elements  #####################
foreach my $a (@file) {
    my $name1 = basename( dirname( dirname( dirname $a) ) );
    my $name2 = basename( dirname( dirname $a) );
    my $name="$name1\-$name2" ;
    open IN, $a;
    while (<IN>) {
        chomp;
	next unless(/\w/);
	my $head=$_;
	my $line=<IN>;
	chomp $line;
        my @b = split /\s+/,$line;
        $hash{$name}{$mapread}=$b[0];
	$hash{$name}{$mapreadintarget}=$b[1];
	die "some error in total" unless($hash{$name}{$totalread});
	$hash{$name}{$rateofmap}=$b[0]/$hash{$name}{$totalread};
	$hash{$name}{$rateintarget}=$b[1]/$b[0];
	
    }
    close IN;
}



################## get enzyme_region_coverage_rate  #####################

foreach my $cov (@cout_file) {
    my $name1 = basename( dirname( dirname( dirname $cov) ) );
    my $name2 = basename( dirname( dirname $cov) );
    my $name="$name1\-$name2";
    open IN, "$cov" || die "no cov\n";
    while (<IN>) {
        chomp;
        next if (/Emycoverate/);
        my @a = split;
        my @b = split /\//, $a[1];
	
        $hash{$name}{$total_target_region}     += $b[1];
        $hash{$name}{$covertargetnum} += $b[0];
    }
    close IN;
}
foreach my $k ( keys %hash ) {
        my $tmp = $hash{$k}{$covertargetnum} / $hash{$k}{$total_target_region};
        $hash{$k}{$rateofcover} = $tmp;
}

###################  out print #######################################


my @hh=( "# of Totalbase", "# of Totalread", "# of Cleanbase", "# of Cleanread", "# of mapped reads", "Rate of mapping reads", "# of mapped reads in enzyme target region", "Rate of mapped reads in enzyme target region", "# of total enzyme $l-$r target regions", "# of covered enzyme $l-$r target regions","Rate of covered enzyme $l-$r target regions");

my @sample=keys %hash;
my $head="\t".join("\t",@sample)."\n";
print OUT $head;
foreach my $type(@hh){
	print OUT "$type";
	foreach my $case(@sample){
		print OUT "\t$hash{$case}{$type}";
	
	}
	print OUT "\n";

}


close OUT;

##################################################################################################
# Generate frontpage
##################################################################################################
system("cp -f $binPath/reportLib/report.html $abs_out_dir/report.html");
system("mkdir -p $abs_out_dir/image/");
system("cp -rf $binPath/reportLib/image/* $abs_out_dir/image/");
# }}

##################################################################################################
# Descend only level directories deep. Dir must be end of a '/'. Call by traverse_dir("/dir/", \@arr) ##
##################################################################################################
sub traverse_dir() {
    # avoid stack overflow, set the maximum traverse number. 
    my ( $cur_dir, $arr ) = @_;

    if ( !opendir( CURDIR, $cur_dir ) ) {
        printf( "open dir %s failed.\n", $cur_dir );
    }
	
    my @all_cur_dir = readdir(CURDIR);
    closedir(CURDIR);

    # drop current dir separator
    $cur_dir =~ s/\/$//;
    foreach my $sub_file (@all_cur_dir) {
        next if ( $sub_file =~ /^\./ );

        if ( -d "$cur_dir/$sub_file" ) {
            &traverse_dir( "$cur_dir/$sub_file", \@$arr );
        }
        else {
            push @$arr, "$cur_dir/$sub_file";
        }
    }
}


