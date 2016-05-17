package Module::hpcaller; 
use strict;
use warnings;
use Data::Dumper; 

# need to keep 
#getargs simple_histogram find_hist_max hp_converter
#getargs string_to_hist hist_to_string stats_basic hp_converter

use Exporter qw(import); 

our @EXPORT_OK = qw( getargs simple_histogram find_hist_max string_to_hist hist_to_string stats_basic hp_converter ); 

sub hp_converter {
    # input: T:9 | T9 => TTTTTTTTT
    # input: TTTTTTTTT => T9 
    my $hp = shift @_; 
    if ( $hp =~ /([ATCG])[:]{0,1}(\d+)/ ) {
        # T9 or T:9 => TTTTTTTTT
        my ( $nt, $l ) = ( $1, $2 ) ;
        my $hp_str; 
        until ( length ($hp_str) == $l ) {
            $hp_str .= $nt; 
        }
        return $hp_str; 
    } elsif ( $hp !~ /\d/ ) {
        my $l = length $hp; 
        my $nt = substr( $hp, 0, 1);
        return "$nt:$l"; 
    }
}

sub hist_to_string {
    # input: histogram as hashref
    # output: stringified histogram useful for printing
    # format: x1,y1;x2,y2;x3,y3
    my ( $h ) = @_ ;
    my @keys = sort { $a <=> $b } keys %$h; 
    my @pairs; 
    map { $h->{$_} ||= 0; push @pairs, "$_,$h->{$_}" } @keys; 
    return join(";", @pairs); 
}

sub string_to_hist {
    # takes a stringified histogram as input
    # returns a hash reference histogram
    # x1,y1;x2,y2;x3,y3 
    my ( $str ) = @_; 
    my @pairs = split(";", $str); 
    my $h; 
    map { my @e = split(",", $_); $h->{$e[0]} = $e[1] } @pairs; 
    return $h; 
}

sub getargs {
    my %opts; #should be an HoA
    my @args = @ARGV; 
    die "\n\nSome arguments required!\n\n" if (scalar(@args) == 0); 
    my $lastflag;
    for my $i (0..(scalar(@args)-1)) {
        chomp $args[$i];
        if ( $args[$i] =~ /^-/ ) { #it's a flag
            $lastflag = $args[$i]; 
            if ( defined($args[$i+1]) == 0 ) { # it's the last argument
                push @{$opts{$args[$i]}}, 1;
            } elsif ( $args[$i+1] =~ /^-/ ) { # it's a standalone flag 
                push @{$opts{$args[$i]}}, 1;
            } elsif ( $args[$i+1] !~ /^-/ ) { # it's a flag with an option
                next; 
            }
        } elsif ( $args[$i] !~ /^-/ ) { # it's not a flag 
            if ( $args[$i-1] =~ /^-/ ) { # it's an option to the previous flag
                push @{ $opts{ $args[$i-1] } }, $args[$i]; 
            } elsif ( $args[$i-1] !~ /^-/ ) { # append arg to last flag seen
                push @{ $opts{ $lastflag } }, $args[$i]; 
            }
        } elsif ( defined($args[$i]) == 0 ) {
            print "$i\t$args[$i]\n";
            die "\n\nGetargs subroutine has somehow reached an element of @ARGV that is undefined!\n\n"; 
        }
        #        print "$lastflag\t$args[$i]", "\n"; 
        

    }
    # convert arrays to scalars when there is only one value in the array
    for my $k (keys %opts) {
        if ( scalar @{ $opts{$k} } == 1 ) {
            my $value = shift @{ $opts{$k} } ; 
            $opts{$k} = $value; 
        } elsif ( scalar @{ $opts{$k} } > 1 ) {
            next; 
        }
    }
    
    return \%opts; 
}

sub simple_histogram {
    my $array = shift @_; 
    return 1 if not defined @$array; 
    my $hist = (); 
    for my $element (@{$array}) {
        $hist->{$element} ++; 
    }
    return $hist; 
}

sub find_hist_max {
    # finds the key for the maximum value in a hash
    # returns the key and the maximum value
    use List::Util qw(max sum); 
    use List::MoreUtils qw(minmax any none uniq); 
    my ( $hist ) = @_; 
    
    # new methodology
    my @k = sort keys %$hist; 
    my $max_value = $hist->{ $k[0] };
    map { $max_value = $hist->{ $_ } if $hist->{ $_ } > $max_value } @k; 
    my @modes;
    map { push @modes, $_ if $hist->{ $_ } == $max_value } @k; 
    return ( \@modes, $max_value ); 
    
}

sub stats_basic {
    use List::Util qw( max );
    # takes a histogram stored in hashref
    my ( $h ) = @_ ;
    my @error = (0, "na", "na", ["na"], "na", "na" ); 
    
    # keys
    my @k = sort keys %$h;
    if ( scalar @k == 0 ) { # hash is empty
        return @error;
    }
    
    # construct array
    my @raw;
    for my $k ( @k ) {
        next unless defined $h->{$k}; 
        my $i = $h->{$k}; 
        while ( $i > 0 ) {
            push @raw, $k;
            $i --; 
        }
    }
    @raw = sort @raw; 
    if ( scalar @raw == 0 ) { # has entries but values are 0 
        return @error; 
    }

    # n
    my $n = scalar @raw; 
    
    # mean 
    use List::Util qw (sum); 
    my $mean = ( sum @raw ) / $n ; 
        
    # modes 
    my $max_value = max( @{ $h }{ keys %$h } );
    my @modes;
    map { push @modes, $_ if $h->{ $_ } == $max_value } @k; 
    my ( $mode_ref, $mode_count ) = ( \@modes, scalar @modes ); 
    
    # mode fraction
    my $mode_frac;
    map { $mode_frac += $h->{$_} } @modes; 
    $mode_frac = $mode_frac/$n; 
    
    #median 
    my $median; 
    if ( $n % 2 == 1 ) { # odd number of elements
        my $i = sprintf "%u", $n/2 + 0.5; 
        $median = $raw[$i-1]; 
    } elsif ( $n % 2 == 0 ) { # even number of elements 
        my $i = sprintf "%u", $n/2; 
        $median = (( $raw[$i-1] + $raw[$i] )/2);
    } else {
        die "basic_stats failed in median\n";
    }
    
    # variance and stdev
    my $d;
    map { $d->{'dev'}->{$_} = ($_ - $mean) } @k; # dev->deviations
    map { $d->{'dev2'}->{$_} = ( $d->{'dev'}->{$_} * $d->{'dev'}->{$_} ) } @k; # square of deviations
    map { $d->{'sum'}->{$_} = ( $d->{'dev2'}->{$_} * $h->{$_} ) } @k; # multiply each squared deviation by the number of times it occurs
    #print STDERR Dumper($d);
    my @sum_dev2 = @{ $d->{'sum'} }{ @k }; # array of summed, squared deviations per key 
    @sum_dev2 = sort { $a <=> $b } @sum_dev2; # sort to minimize float errors in addition
    my $sum_of_dev2 = sum @sum_dev2;  # the sum of squared deviations
    my $variance = ( $sum_of_dev2 / $n ); 
    my $stdev = sqrt( $variance ); 
    
    
    
    return ( $n, $mean, $stdev, $mode_ref, $mode_frac, $median ); 
}

1; 