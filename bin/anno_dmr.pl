use strict;
use File::Basename qw/dirname basename/;

die "Usage: $0 <element> [ <dmr_file_list> ... ]\n" if @ARGV < 2;
my ( $db, @dmr_list ) = @ARGV;

my @db_all;

## databases dir ##
#$db ||= "/nas/RD_12A/gaoshengjie/Database/Human/RRBS/bin/prepare/hg19/element/";
push @db_all, $db;
## interval size ##
my $MAX_INTERVAL = 1e6;

my %hash_db;
my @db_used;
foreach my $db_type (@db_all) {
    open DB, $db_type;
    my $db_name = basename $db_type;
    $db_name =~ s/\..*//g;
#    print $db_name;die;
    push @db_used, $db_name;
    while (<DB>) {
        chomp;
        my ( $chr, $sp, $ep, $info ) = split /\s+/;

        my $int_trunc = int( $sp / $MAX_INTERVAL );
        my $lower     = $int_trunc * $MAX_INTERVAL + 1;
        ++$int_trunc;
        my $upper = $int_trunc * $MAX_INTERVAL;

        if ( $ep > $upper ) {

            # save [sp, upper)
            push @{ $hash_db{$db_name}{$chr}{"$lower-$upper"} }, "$sp,$ep,$info";

            my $int_trunc = int( $ep / $MAX_INTERVAL );
            my $lower     = $int_trunc * $MAX_INTERVAL + 1;
            ++$int_trunc;
            my $upper = $int_trunc * $MAX_INTERVAL;

            # save [upper, ep];
            push @{ $hash_db{$db_name}{$chr}{"$lower-$upper"} }, "$sp,$ep,$info";
        }
        else {
            push @{ $hash_db{$db_name}{$chr}{"$lower-$upper"} }, "$sp,$ep,$info";
        }
    }
    close DB;
}



foreach my $dmr_file (@dmr_list) {
    open FILE, $dmr_file or die $!;
    #print $dmr_file; die;
    my %hash_fd;
    foreach my $db_name (@db_used) {
        my $FD;
        open $FD, ">$dmr_file.$db_name" or die $!;
	#print "$dmr_file.$db_name\n";die;
        $hash_fd{$db_name} = $FD;
    }

    while ( my $line = <FILE> ) {
        chomp $line;
        my ( $chr, $sp, $ep ) = ( split /\s+/, $line )[ 0, 1, 2 ];

        my $int_trunc = int( $sp / $MAX_INTERVAL );
        my $lower     = $int_trunc * $MAX_INTERVAL + 1;
        ++$int_trunc;
        my $upper = $int_trunc * $MAX_INTERVAL;

        foreach my $db ( sort keys %hash_db ) {
            select $hash_fd{$db};
            $| = 1;
            print $line;
            goto LABLE_END if ( !exists( $hash_db{$db}{$chr} ) );
            goto LABLE_END if ( !exists( $hash_db{$db}{$chr}{"$lower-$upper"} ) );

            foreach my $key_pair ( @{ $hash_db{$db}{$chr}{"$lower-$upper"} } ) {
                my ( $key_sp, $key_ep, $info ) = split /,/, $key_pair;

                if ( $sp <= $key_sp && $key_sp <= $ep ) {
                    my $share = $ep - $key_sp + 1;
                    $share = $key_ep - $key_sp + 1 if ( $key_ep < $ep );
                    print "\t$share|$key_sp,$key_ep|$info";
                }
                elsif ( $key_sp <= $sp && $sp <= $key_ep ) {
                    my $share = $key_ep - $sp + 1;
                    $share = $ep - $sp + 1 if ( $ep < $key_ep );
                    print "\t$share|$key_sp,$key_ep|$info";
                }
                else {

                    # not in giving area. #
                }
            }

          LABLE_END:
            print "\n";
            select STDOUT;
        }
    }
    print STDERR "-- $dmr_file annotation done --\n";
    close FILE;

    # close FD #
    foreach my $db ( keys %hash_fd ) {
        close( $hash_fd{$db} );
    }
}

