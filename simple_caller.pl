#! /usr/bin/perl

use strict;
use warnings FATAL => qw(all);

# takes the output of hp_aggregator.pl and does a simple depth and mode based calling

use Data::Dumper;
use My::Module qw(getargs simple_histogram hp_converter string_to_hist);
use List::Util qw(max min sum);
use List::MoreUtils qw(uniq none indexes firstidx any);

&main;
exit 1; 

sub main {
    my $opts = &getargs; 
    my $agg = $opts->{'-agg'}; 
    $opts->{'-loc_min_dp'} ||= 200;
    $opts->{'-loc_max_dp'} ||= 5000;
    $opts->{'-s_min_dp'} ||= 10;
    $opts->{'-s_max_dp'} ||= 100;
    
    
    my @ag_headers;
    my @ag_col_heads;
    
    open ( my $AG, $agg ) or die "could not open file: $agg\n"; 
    while (<$AG>) {
        chomp $_;
        if ( $_ =~ /^##/ 
            ) { push @ag_headers, $_;
                print $_, "\n";
        } 
        elsif ( $_ =~ /^#CHROM/ 
        ) { $_ =~ s/^#//;
            @ag_col_heads = split("\t", $_);
            print "#". $_, "\n";
        }
        else {
            my $lh = {};
            @{ $lh }{ @ag_col_heads } = split("\t", $_);
            &process1( \@ag_col_heads, $lh, $opts ); 
            
        }
    }
    
}

sub process1 {
    my ( $fields, $lh, $opts ) = @_; 
    
    $lh->{"REF"} =~ /([ATCG]):{0,1}(\d+)/; 
    my ( $ref_n, $ref_l ) = ( $1, $2 );
    
    
    my $allele_hash = {};
    $allele_hash->{ $ref_l } = "0/0";
    
    # fill in some defaults
    my $olh = {}; #output line hash
    # copy chrom, pos, id, ref, alt, qual, filter, info fields from input to output
    @{$olh}{ @$fields[0..7] } = @{$lh}{ @$fields[0..7] };
    $olh->{"REF"} = &hp_converter($olh->{"REF"});
    # fill the format field (defines format for sample fields)
    $olh->{"FORMAT"} = "GT:DP:HIST";
    # fill in N/N as the default genotype for all samples
    map { $olh->{$fields->[$_]} = "N/N:" . $lh->{$fields->[$_]} } ( 9 .. $#$fields );
    
    
    #check locus depth to see if within min and max
    # if no, print and proceed to next line
    $lh->{"INFO"} = /DP4=(\d+):/;
    my $loc_dp = $1;
    if ( $loc_dp < $opts->{'-loc_min_dp'} || $loc_dp > $opts->{'-loc_max_dp'} ) {
        print join("\t", @{$olh}{@$fields}), "\n";
        return 1; 
    } else {
        #iterate through samples
        for my $s_name ( @$fields[9..$#$fields] ) {
            my $s_str = $lh->{$s_name};
            my ($s_dp, $s_hist_str) = split(":", $s_str);
            my $s_hist = &string_to_hist($s_hist_str); #convert the histogram string to a hashref
            # check to see if sample meets depth requirements
            if ( $s_dp < $opts->{'-s_min_dp'} || $s_dp > $opts->{'-s_max_dp'} ) {
                # sample fails depth requirements; leave as N/N
            } else {
                # sample meets depth requirements; analyze more 
                #determine mode
                my $s_max_val = max ( @{$s_hist}{ keys %$s_hist } );
                my @s_modes;
                map { push @s_modes, $_ if $s_hist->{$_} >= $s_max_val; } ( keys %$s_hist );
                if ( scalar @s_modes != 1 ) {
                    # sample is not unimodal; call as N/N
                    next;
                } else {
                    # is mode different from the reference?
                    my $s_mode = $s_modes[0];
                    if ( $s_mode == $ref_l ) {
                        # sample has reference allele
                        $olh->{$s_name} = "$allele_hash->{$ref_l}:"  . "$lh->{$s_name}";
                        next; 
                    } else {
                        #print STDERR "sample has non-reference mode\n";
                        # sample mode is different from reference allele
                        my @observed_alleles = ( keys %$allele_hash );
                        if ( any { $_ eq $s_mode } @observed_alleles ) {
                            # sample has an already observed allele, assign that genotype
                            $olh->{$s_name} = "$allele_hash->{$s_mode}:" . "$lh->{$s_name}"; 
                        } else {
                            # sample has new genotype, add to allele hash and assign new genotype
                            my $current_allele_count = scalar ( keys %$allele_hash ); 
                            $allele_hash->{$s_mode} = join("/", $current_allele_count, $current_allele_count);
                            $olh->{$s_name} = "$allele_hash->{$s_mode}:" . "$lh->{$s_name}"; 
                        }
                    }
                }
            }
        }
        #finished iterating through samples
        # construct alternate allele list for ALT field
        my $allele_count = scalar ( keys %$allele_hash ); 
        my @alt_alleles = (); 
        if ( $allele_count > 1 ) {
            # list the alternate alleles in numerical order
            for my $i ( 1 .. $allele_count ) {
                # test each alternate allele to see if it is the right number
                for my $allele ( keys %$allele_hash ) {
                    if ( $allele_hash->{$allele} eq "$i/$i" ) {
                        # it's the right number, add it to the alt alleles list
                        push @alt_alleles, &hp_converter( "$ref_n" . "$allele" );
                    } else {
                        # do nothing
                    }
                }
            }
            $olh->{"ALT"} = join(",", @alt_alleles); 
        } else {
            # do nothing, the alt allele field is a .
        }
        
    }
    print join("\t", @{$olh}{@$fields}), "\n"; #print output line
    return 1;
}
