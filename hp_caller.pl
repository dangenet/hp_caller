#! /usr/bin/perl

# this is going to call samples relative to the locus without paying attention to the rest of the genome
# this is hp_caller3 renamed hp_caller for publication

use strict;
use warnings;
use Data::Dumper;
use My::Module qw( getargs string_to_hist hist_to_string hp_converter sample_ids cohort_ids); 
use List::Util qw( sum max min ); 
use List::MoreUtils qw( true any none uniq all firstidx ); 
use Scalar::Util qw(looks_like_number); 
use Math::Round qw( nearest_ceil nearest ); 
use Storable; 
use Math::BigInt; 
use Math::BigFloat;

my $opts = &getargs;

&main;
exit;

sub main {
    my $agg = $opts->{'-v'};
    die "input file required with -v\n" if not -e $agg;
    $opts->{'--min_sdp'} ||= 5;
    $opts->{'--max_sdp' } ||= 150;
    $opts->{'--nsr'} ||= 1;
    $opts->{'--uncallable'} ||= 3;
    $opts->{'--loc_min'} ||= 200;
    $opts->{'--loc_max'} ||= 2500;
    $opts->{'-mp'} ||= 5; #PHRED scaled prob of a sample being different from locus
    $opts->{'--asymmetry'} ||= 0.2; #maximum fractional asymmetry around the mode in a distribution
    $opts->{'--debug'} ||= 0;
    $opts->{'--variants_only'} ||= 0; 
    $opts->{'--minSGQ'} ||= 10; #PHRED scaled min sample distribution quality 
    $opts->{'--minQUAL'} ||= 10; #PHRED scaled quality of locus distribution quality
    # both QUAL scores are 100 * ( mode counts - 2nd most abundant counts) / total depth
    
    # print parameters
    
    open ( my $AGG_FILE, $agg ) or die "could not open $agg";
    while ( <$AGG_FILE> ) {
        chomp $_; 
        if ( $_ =~ /^##/
        ) {
            push @{$opts->{"HEADERS"}}, $_; 
        } elsif ( $_ =~ /^#CHROM/ 
        ) {
            print $_, "\n";
            @{$opts->{"FIELDS"}} = split("\t", $_);
            $opts->{"FIELDS"}->[0] =~ s/^#//; 
        } elsif ( $_ !~ /^#/ 
        ) { 
            debug("reading data line in main") if $opts->{'--debug'};
            my $lh = &locus_hashmaker( $opts, $_ );
            
            #really probably want to exclude 22582 from locus histogram
            # since it can actually be heterozygous
            # delete the 22582 histogram data 
            map { $lh->{"INFO"}->{"HIST"}->{$_} -= $lh->{22582}->{"HS"}->{$_} } (keys %{$lh->{'22582'}->{"HS"}});
            &locus_caller( $opts, $lh );
        }
    }
    close $AGG_FILE; 
    
    return 1; 
}

erR39pcwrMDR

sub locus_caller {

    my ($opts, $lh) = @_; 
    &debug("locus_caller") if $opts->{'--debug'} ;
    # pull some useful subsets of data out of the hashref
    my @fields = @{$opts->{"FIELDS"}};
    my @samples = @fields[9..$#fields];
    # make a hashref that just contains the histograms without the other stuff
    my $locus_hist->{"LOCUS"} = $lh->{"INFO"}->{"HIST"};
    map { $locus_hist->{$_} = $lh->{$_}->{"HS"} } @samples;
    # set default return values
    my $loc_dp = $lh->{"INFO"}->{"DP4"};
    $lh->{"QUAL"} = "na";
    $lh->{"INFO"}->{"MODE"} = "na"; 
    map { $lh->{$_}->{"MODE"} = "na"; $lh->{$_}->{"GQ"} = "na"; $lh->{$_}->{"GT"} = "N"; $lh->{$_}->{"MP"} = "na" } @samples;
    #start constructing the allele list
    my @loc_alleles;
    push @loc_alleles, $lh->{"REF"}->{"l"}; #ref genome allele gets index 0

    #determine if the locus is callable
    # the problem with this is that later I am only going to call samples as mutant if they
    # have a given GQ (probability of having come from the locus) greater than 5
    # however, I am going to call samples with GQ less than 5 as being the same as the locus
    # so effectively, those samples with GQ less than 5 are uncallable, because they can't be called as mutants
    # but I'm not taking that into account 
    my $locus_flag = &callable_locus( $locus_hist, \@samples, $opts );

#    if ( $locus_flag == 0 && $opts->{'--variants_only'} == 0 ) { #print and return early if locus is uncallable
#        &vcf_printer($opts, $lh);
#        return 1;
#    }
    
    #fill in the locus mode hashref
    my $loc_modes; 
    map { $loc_modes->{$_} = &just_modes($locus_hist->{$_}); }( "LOCUS", @samples);
    #fill in output values
    $lh->{"INFO"}->{"MODE"} = join(",", @{ $loc_modes->{"LOCUS"} } );
    $lh->{"QUAL"} = &locus_qual_score($opts, $locus_hist->{"LOCUS"} );
    #add the locus mode to the allele list if not the same as reference 
    # doing it now ensures that it will have allele number 1
    if ( $loc_modes->{"LOCUS"}->[0] != $loc_alleles[0] ) {
        push @loc_alleles, $loc_modes->{"LOCUS"}->[0]; 
    }
    
    # do a quick screen to see if any samples are callable and have modes different than the locus
    # don't do the per-sample genotype quality calculation unless there are some 

    my $pvf = 0; #potential variant flag

    for my $s ( @samples ) {
        next if $s == 22582;
        
        my $callable = &callable_sample( $opts, $lh->{$s}->{"DP"}, $locus_hist->{$s} );
        $lh->{$s}->{"MODE"} = join(",", @{$loc_modes->{$s}});
        
        if ( $callable == 1 
            && $loc_modes->{$s}->[0] != $loc_modes->{"LOCUS"}->[0] 
        ) { # sample is callable and has a mode different from the locus
            $pvf = 1;
        } elsif ( $callable == 1 && $loc_modes->{$s}->[0] == $loc_modes->{"LOCUS"}->[0] ) { 
            #sample is callable and has same mode as locus
            $lh->{$s}->{"GT"} = firstidx { $loc_modes->{$s}->[0] == $_ } @loc_alleles;
        } else { # sample is uncallable
            next;
        }
        $lh->{$s}->{"GQ"} = &locus_qual_score( $opts, $locus_hist->{$s} );

    }
    
#    if ( $pvf == 0 && $opts->{'--variants_only'} == 0 ) { #terminate early if no samples have mode different from locus
#        $lh->{"REF"} = $lh->{"REF"}->{"seq"};
#        &vcf_printer($opts, $lh);
#        return 1;
#    }
    
    # for lines containing potentially variant samples
    # calculate probabilities and mark variant samples for those 
    my $vsf; #variant sample flag hashref
    map { $vsf->{$_} = 0 } @samples; # set all to nonvariant initially
    my $sph; #sample probability hash
    map { $sph->{$_} = "na"; } @samples; #prefill with error values

    if ( $pvf == 1 && $locus_flag == 1 ) {
        for my $s ( @samples ) {
            next if $s == 22582;
            
            $sph->{$s} = &prob_sample_came_from_locus($loc_modes->{"LOCUS"}->[0], $locus_hist->{$s}, $locus_hist->{"LOCUS"} );
            if (&callable_sample( $opts, $lh->{$s}->{"DP"}, $locus_hist->{$s} ) 
            &&  $sph->{$s} <= 10 ** ( -1 * $opts->{'-mp'} )
            &&  $loc_modes->{$s}->[0] != $loc_modes->{"LOCUS"}->[0] ) {
                #sample is callable, has different mode and meets GQ threshold
                $vsf->{$s} = 1; 
                # add allele to the list if not already represented
                if ( none { $loc_modes->{$s}->[0] == $_ } @loc_alleles ) {
                    push @loc_alleles, $loc_modes->{$s}->[0];
                }
                #mark its genotype
                $lh->{$s}->{"GT"} = firstidx { $loc_modes->{$s}->[0] == $_ } @loc_alleles;
                #compute phred scaled probability that sample is !!!non!!!-locus
                # that is
                $lh->{$s}->{"MP"} = $sph->{$s}->copy()->bround(5)->blog(10)->bmul(-1);
            } elsif ( &callable_sample( $opts, $lh->{$s}->{"DP"}, $locus_hist->{$s} ) ) { 
                # sample is callable, but did not meet GQ threshold 
                # may or may not be mutant 
                # call as having locus mode
                $lh->{$s}->{"GT"} = firstidx { $loc_modes->{"LOCUS"}->[0] == $_ } @loc_alleles;
                # GQ is -log(10) of p(sample has indicated genotype)
                $lh->{$s}->{"MP"} = $sph->{$s}->copy()->bround(5)->blog(10)->bmul(-1);
            } else {
                # sample is uncallable 
            }
        }
    }
    
    #construct the alternate allele field
    my @loc_seqs;
    for my $i ( 0..$#loc_alleles ) {
        $loc_seqs[$i] = $lh->{"REF"}->{"nt"};
        while ( length $loc_seqs[$i] < $loc_alleles[$i] ) {
            $loc_seqs[$i] .= $lh->{"REF"}->{"nt"};
        }
    }
    
    $lh->{"ALT"} = join(",", @loc_seqs[1..$#loc_seqs]) if scalar @loc_seqs > 1;
    $lh->{"REF"} = $loc_seqs[0];
    
    if ( $opts->{'--variants_only'} == 1 ) {    
        &vcf_printer( $opts, $lh) if any { $vsf->{$_} == 1 } @samples;
    } elsif ( $opts->{'--variants_only'} == 0 ) {
        &vcf_printer( $opts, $lh);
    }
    
    #    if ( scalar @loc_alleles != 1 ) {
    #        print join("\t", @{$lh}{qw(CHROM POS)}), "\t";
    #        map { print join("\t", $_, $loc_seqs[$_] ), "\t" } ( 0 .. $#loc_alleles );
    #        print "\n";
    #    }

    return 1; 
}


sub vcf_printer {
    my ( $opts, $lh ) = @_;
    my @fields = @{$opts->{"FIELDS"}};
    my @samples = @fields[9..$#fields];

    
    if ( ref $lh->{"REF"} eq "HASH" ) {
        $lh->{"REF"} = $lh->{"REF"}->{"seq"};
    }
    
    my $info_str = "MODE=" . $lh->{"INFO"}->{"MODE"} . ":"; 
    $info_str .= "DP4=" . $lh->{"INFO"}->{"DP4"} . ":";
    $info_str .= "HIST=" . &hist_to_string($lh->{"INFO"}->{"HIST"}); 
    $lh->{"INFO"} = $info_str;

    my @format = qw( GT MODE DP GQ MP HS );
    $lh->{"FORMAT"} = join(":", @format);
    $lh->{"FORMAT"} =~ s/HS/HIST/;
    
    for my $s ( @samples ) {
        $lh->{$s}->{"HS"} = &hist_to_string( $lh->{$s}->{"HS"} );
        $lh->{$s}->{"GT"} .= "/$lh->{$s}->{'GT'}"; 
        $lh->{$s} = join(":", @{$lh->{$s}}{ @format } );
    }
    
    
    print join("\t", @{$lh}{@fields} ), "\n";

    
    return 1;
}

sub callable_locus {
    # applies criteria for determining locus is callable
    # returns boolean
    my ( $lh, $samples, $opts ) = @_;
    &debug("callable_locus") if ( $opts->{'--debug'} == 1 );
    my $locus_flag = 1; 
    my @locus_modes = @{ &just_modes( $lh->{"LOCUS"} ) };
    my $loc_dp = sum( @{$lh->{"LOCUS"}}{( keys %{$lh->{"LOCUS"}} ) } );
    
    if ( $locus_flag == 1 ) {
        $locus_flag = &mode_checker(\@locus_modes);
    }
    if ( $locus_flag == 1 ) {
        $locus_flag = &asymmetry($opts, $locus_modes[0], $lh->{"LOCUS"});
    }
    if ( $locus_flag == 1 ) {
        $locus_flag = &callable_sample_check( $opts, $lh, $samples );
    }
    if ( $locus_flag == 1 ) {
        $locus_flag = &locus_qual_check( $opts, $lh->{"LOCUS"} );
    }

    return $locus_flag;
}

sub locus_qual_check { 
    my ( $opts, $lh ) = @_;
    my $return = 0;
    my $qual = &locus_qual_score( $opts, $lh);
    if ( $qual >= $opts->{'--minQUAL'} ) {
        $return = 1;
    }
    return $return;
}

sub locus_qual_score {
    my ( $opts, $lh ) = @_; #opts, #locus hash histogram
    &debug("locus_qual_score") if $opts->{'--debug'};
    my $loc_modes = &just_modes($lh);
    my @values = sort {$b<=>$a} ( @{$lh}{(keys %$lh)} );
    push @values, 1; #if there is only one value in values, add one to the end
    my $dp = sum( @values );
    
    my $qual = 100 * ($values[0] - $values[1]) / $dp;
    # the mode fraction minus the next most abundant fraction
    # max score is almost 100; e.g., 100*(50-1)/50
    # a typical score might be ( 150 - 50 ) / 200 = 50
    # min score is 100 * ( 101-99)/200 = 1
    
    # for a sample, the score for the weakest possible distribution is 100 * (3-2)/5 = 20
    # 100 * ( 3 - 1 ) / 5 = 40
    # could have 100*(25-24)/49 = 2.04
    $qual = &nearest_ceil(0.1, $qual);
    return $qual;
}

sub callable_sample_check { 
    #determines whether the number of uncallable samples at a locus exceeds the threshold
    debug("callable_sample_check") if $opts->{'--debug'};
    my ( $opts, $lh, $samples ) = @_;
    my $return = 0; 
    my $uncallable_count = 0; 
    
    for my $s ( @$samples ) {
        next if $s == 22582;

        my $callable = &callable_sample( $opts, sum( @{ $lh->{$s} }{( keys %{$lh->{$s}} )} ), $lh->{$s} );
        $uncallable_count ++ unless $callable;
    }
    if ( $uncallable_count <= $opts->{'--uncallable'} ) {
        $return = 1;
    }
    return $return;
}

sub mode_checker { #returns true if unimodal
    &debug("mode_checker") if $opts->{'--debug'};
    my ($mr) = @_;
    my $return = 0;
    if ( scalar @$mr == 1 ) {
        $return = 1;
    }
    return $return;
}

sub locus_depth_checker {
    &debug("locus_depth_checker") if $opts->{'--debug'};
    my ( $dp ) = @_;
    my $return = 0;
    if ( $dp >= $opts->{'--loc_min'} || $dp <= $opts->{'--loc_max'} ) {
        $return = 1;
    }
    return $return;
}

sub asymmetry {
    # checks to see if the locus has a peak at the mode with decreasing counts going out from there
    my ( $opts, $mode, $lh ) = @_; 
    &debug("asymmetry") if $opts->{'--debug'};

    my $return = 0; # default to false
    my @keys = ( min(keys %$lh) - 1 .. max( keys %$lh ) + 1 );
    map { $lh->{$_} ||= 0 } @keys; #backfill the locus hash with zeros 
    my $loc_dp = sum( @{$lh}{@keys} );
    
    map { $lh->{$_}++ } @keys; #smooth
    my $ph;
    # locus mode frac is $ph->{$mode};
    map { $ph->{$_} = $lh->{$_}/($loc_dp + scalar @keys) } @keys; #faster than using construct_ph
    
    my $mean = sum ( map { $_ * $lh->{$_} } @keys ) / sum ( @{$lh}{@keys} );
    
    #check for pronounced asymmetry around the mode
    my $mf = $ph->{$mode};
    my $upper_frac = sum ( @{$ph}{( $mode + 1 .. $keys[-1])} );
    my $lower_frac = sum ( @{$ph}{( $keys[0] .. $mode - 1)} );
    my $diff = abs( $upper_frac - $lower_frac );
    
    if ( $upper_frac < 0.25 && $lower_frac < 0.25 ) {
        $return = 1;
    } else {
        if ( $diff > $opts->{'--asymmetry'} ) {
            $return = 0;
        } else { 
            $return = 1;
        }
        
    }
#    if ( $return == 0  ) {
#        print STDERR "mode\t$mode\tmf\t$ph->{$mode}\n";
#        map { print STDERR join("\t", $_, $lh->{$_}, $ph->{$_}), "\n"; } @keys;
#        print STDERR "diff\t$diff\tupper\t$upper_frac\tlower\t$lower_frac\n\n";
#        print $diff, "\n";
#
#    }
    return $return; 
}

sub prob_sample_came_from_locus { 
    my ($h, $sh, $tlh, $opts) = @_; #sample hashref, #total locus histogram
    &debug("prob_sample_came_from_locus") if $opts->{'--debug'};

    my $lh; #subtract the sample values from the total locus histogram
    map { $sh->{$_} ||= 0; $lh->{$_} = $tlh->{$_} - $sh->{$_} } ( keys %$tlh );
    
    #my $ph = &construct_ph( $lh );
    
    my $prob = &binomial_probability( $h, $sh, $lh );
    if ( $prob > 1 || $prob < 0 ) {
        print STDERR Dumper($h, $sh, $tlh);
        print STDERR $prob, "\n";
        die;
    }
    return $prob;
}

sub debug {
    print STDERR join("\n", @_), "\n";
    return;
}

sub binomial_probability {
    &debug("in binomial_probability") if $opts->{'--debug'};
    my ( $h, $sh, $lh ) = @_; #success key, sample hist, locus hist
    #construct the probability table for probability of h, prob of not h
    my $ph;

    $ph->{$h} = Math::BigFloat->new( $lh->{$h} );
#    print STDERR $h, "\t", $ph->{$h}, "\n";

    my $loc_dp = sum( @{$lh}{(keys %$lh)} );
    #    $loc_dp += 2; #smooth
    
    $ph->{$h}->badd(1)->bdiv($loc_dp + 1); 
#    print STDERR $ph->{$h}, "\n";
    $ph->{"other"} = Math::BigFloat->new( $loc_dp );
#    print STDERR $ph->{"other"}, "\n";
    $ph->{"other"}->bsub( $lh->{$h} )->badd(1);
#    print STDERR $ph->{"other"}, "\n";
    $ph->{"other"}->bdiv($loc_dp + 1);
#    print STDERR $ph->{"other"}, "\n";

    #print STDERR "$h\t$ph->{$h}\nother\t$ph->{'other'}\n";
    # count reads at h and not at h
    my $shc;
    $shc->{$h} = $sh->{$h};
    map { $shc->{"other"} += $sh->{$_} unless $_ == $h } ( keys %$sh );
    
    my $bnc = &binom_coeff( $h, $sh );
    
    my $prob = $bnc->copy();
    
    my $k_prob = Math::BigFloat->new( $ph->{$h} );
    $k_prob->bpow( $shc->{$h} );
    my $Nk_prob = Math::BigFloat->new( $ph->{"other"} );
    $Nk_prob->bpow( $shc->{"other"} );
    
    $prob->bmul($k_prob)->bmul($Nk_prob); 
    map { print STDERR join("\t",$_, $ph->{$_}), "\n" } ( keys %$ph ) if ( $prob < 0 || $prob > 1 );
    return $prob;
}

sub binom_coeff {
    my ( $h, $sh ) = @_;
    my $N = Math::BigFloat->new( sum( @{$sh}{(keys %$sh)} ) ); #get N
    my $k = Math::BigFloat->new( $sh->{ $h } ); #get k
    my $N_k = $N->copy(); 
    $N_k->bsub( $k ); #calculate N-k
    $k->bfac(); #calculate k! 
    $N_k->bfac(); #calculate (N-k)! 

    my $bnc = $N->copy(); 
    $bnc->bfac(); 
    $bnc->bdiv($N_k)->bdiv($k);
    
    return $bnc; 
}

sub construct_ph {
    print STDERR "in construct_ph\n";
    my ( $lh ) = @_;
    my $lh_dp = sum( @{$lh}{( keys %$lh )});
    
    my $ph; # make the locus probability value hash
    for my $k ( keys %$lh ) { 
        $ph->{$k} = Math::BigFloat->new($lh->{$k});
        $ph->{$k}->bdiv($lh_dp)->blog();
    }
    # values are stored as the ln of the probability

    return $ph;
}

sub callable_sample {
    my ($opts, $dp, $hist) = @_;
    my $callable = 1; 
    if ( $dp < $opts->{'--min_sdp'} ) {
        # samples must meet a minimum sample depth threshold
        $callable = 0;
    } 
    if ( $dp > $opts->{'--max_sdp'} ) {
        # sample must be below max threshold
        $callable = 0;
    }
    my @modes = @{ &just_modes( $hist ) };
    if ( scalar @modes > 1 ) {
        # multimodal samples are not considered callable
        $callable = 0;
    }
    my $nsr = &net_supporting_reads($opts,$hist);
    if ( $nsr < $opts->{'--nsr'} ) {
        # samples must meet net supporting reads threshold
        $callable = 0;
    }
    my $sgq = &locus_qual_score( $opts, $hist );
    if ( $sgq < $opts->{'--minSGQ'} ) {
        $callable = 0;
    }
    
    return $callable;
}

sub net_supporting_reads {
    my ( $opts, $hr ) = @_;
    my @keys = ( sort {$a<=>$b} ( keys %$hr ) );
    my @values = ( sort {$b<=>$a} @{$hr}{@keys} );
    $values[1] ||= 0;
    $values[0] ||= 0;
    my $nsr = $values[0] - $values[1];
    return $nsr;
}

sub just_modes {
    my ( $hr ) = @_;
    my @keys = ( sort {$a<=>$b} ( keys %$hr ));
    my $max_val = max( @{$hr}{@keys} );
    my @modes;
    map { push @modes, $_ if $hr->{$_} == $max_val } @keys;
    return \@modes;
}

sub locus_hashmaker {
    my ($opts, $line) = @_;
    &debug("locus_hashmaker") if $opts->{'--debug'};
    my @fields = @{$opts->{"FIELDS"}};
    my @samples = @fields[9..$#fields];
    my $lh;
    @{ $lh }{ @fields } = split("\t", $_);
    # split the ref field into its components
    my $refstr = $lh->{"REF"};
    delete $lh->{"REF"};
    @{ $lh->{"REF"} }{("nt", "l", "seq")} = &refseq_parser($refstr);
    # split the info field into its components
    $lh->{"INFO"} = &info_parser( $lh->{"INFO"} );
    # parse the samples into the line hash
    for my $s ( @samples ) {
        $lh->{$s} = &sample_lh_parser( @{$lh }{("FORMAT", $s)});
    }
    return $lh; # or print
}

sub sample_lh_parser {
    my ( $format, $str ) = @_;
    my $sh; 
    my @formats = split(":", $format);
    @{$sh}{ @formats } = split(":", $str) or die "$str";
    $sh->{"HS"} = &string_to_hist($sh->{"HS"}) unless $sh->{"HS"} eq "na";
    return $sh; 
}

sub refseq_parser {
    my ( $hp ) = @_; # in format A8
    my ($nt, $l) = ($1, $2) if $hp =~ /([ATCG]):{0,1}(\d+)/; 
    my $hpstr = &hp_converter($hp);
    return ($nt, $l, $hpstr);
}

sub info_parser {
    my ( $info ) = @_; 
    my @temp = split(":", $info);
    my $ih;
    for my $i ( @temp ) {
        my ($key, $val) = split("=", $i); 
        if ( $key eq "HIST" ) {
            $val = &string_to_hist( $val );
        }
        $ih->{$key} = $val;
    }
    return $ih; 
}

#sub fast_binomial {
#    # currently not used
#    # does everything in log space 
#    # does not use Math::BigFloat
#    my ( $h, $sh, $lh) = @_;
#    
#    # p(data) = n!/(xi! .. xk! ) * p1^x1 .. pk^xk
#    # log transform 
#    # log(p(data)) = log( multinomial coefficient ) + x1* log( p1 ) + .. + xk * log ( pk )
#    
#    #convert histograms to binary data
#    my $sbin->{$h} = $sh->{$h};
#    my $lbin->{$h} = $lh->{$h};
#    for my $i ( keys %$lh ) {
#        next if $i == $h;
#        $sbin->{"other"} += $sh->{$i};
#        $lbin->{"other"} += $lh->{$i};
#    }
#    map { $lbin->{$_} ++ } ( $h, "other" ); # smooth
#    my $loc_dp = sum( @{$lbin}{($h, "other")} ); 
#    
#    my $ph;
#    map { $ph->{$_} = log( $lbin->{$_} ) - log ( $loc_dp ) } ( $h, "other" );
#    
#    #calculate binomial coefficient (in log space)
#    #numerator
#    my @numerator = ( sort {$a<=>$b} ( 1 .. sum( @{$sbin}{(keys %$sbin)} ) ) );
#    while ( $numerator[0] == 1 ) {
#        shift @numerator;
#    }
#    map { $numerator[$_] = log( $numerator[$_] ) } ( 0 .. $#numerator );
#    #denominator
#    my @denominator; 
#    map { push @denominator, $_ } ( 1 .. $sbin->{$h});
#    map { push @denominator, $_ } ( 1 .. $sbin->{"other"});
#    @denominator = sort {$a<=>$b} @denominator; 
#    while ( $denominator[0] == 1 ) {
#        shift @denominator;
#    }
#    map { $denominator[$_] = log( $denominator[$_] ) } ( 0 .. $#denominator );
#    #binomial coefficient
#    my $bnc = sum( @numerator );
#    $bnc -= sum( @denominator );
#    
#    #calculate probability
#    my $prob = 0; 
#    for my $i ( keys %$sbin ) {
#        $prob += $sbin->{$i} * $ph->{$i} 
#    }
#    $prob += $bnc;
#    
#    return $prob; # is ln(prob)
#}

#sub cohort_variant_test {
#    #came out of locus caller subroutine
#    my $locus_flag = 0;
#    my $locus_hist;
#    
#    # compare cohorts to the locus
#    # come up with about 5-6 over the entire genome
#    # ignore for the moment
#    my $variant_cohort_flag = 0;
#    my $variant_cohort;
#    if ( $locus_flag == 1 && $variant_cohort_flag == 1 ) {
#        my $cph; #cohort probability hash
#        my $cohorts = {
#            "empty" => [ 22539, 22544 .. 22551 ],
#            "WT" => [ 22540, 22552 .. 22559 ],
#            "E125Q" => [ 22541, 22560 .. 22567 ],
#            "Y162A" => [ 22542, 22568 .. 22575 ],
#            "N169S" => [ 22543, 22576 .. 22583 ]
#            
#        };
#        my @cohort_names = qw(empty WT E125Q Y162A N169S);
#        for my $c ( @cohort_names ) {
#            my $ch; # cohort hash - fill 
#            for my $s ( @{$cohorts->{$c}} ) {
#                map { $ch->{$_} += $locus_hist->{$s}->{$_} } ( keys %{$locus_hist->{$s}} ); 
#            }
#            my $rol; #rest of locus hash
#            map { $ch->{$_} ||= 0; $rol->{$_} = $locus_hist->{"LOCUS"}->{$_} - $ch->{$_} } ( keys %{$locus_hist->{"LOCUS"}} );
#            my $rol_modes = &just_modes($rol); #recalculate the mode since it could change when subtracting a cohort
#            my $c_modes = &just_modes($ch); 
#            if ( $c_modes->[0] != $rol_modes->[0] ) {
#                $cph->{$c} = &binomial_probability( $rol_modes->[0], $ch, $rol );
#                #map { print STDERR join("\t", $c, $_, $ch->{$_}, $rol->{$_}, $locus_hist->{"LOCUS"}->{$_}) , "\n" } ( sort {$a<=>$b} ( keys %$rol ) );
#                #print STDERR $cph->{$c}->bround(4), "\n";
#                $cph->{$c}->bround(10)->blog(10)->bmul(-1) unless $cph->{$c} == 0;
#            } else {
#                $cph->{$c} = 0;
#            }
#            
#        }
#        print join("\t", @{$cph}{@cohort_names}), "\n"; # if any { $cph->{$_} > 0 } @cohort_names;
#    }
#
#}

#sub multinomial_probability {
#    # given a sample histogram and a probability table
#    # calculate the probability with the multinomial function
#    # p(data) = n!/(xi! .. xk! ) * p1^x1 .. pk^xk
#    # log transform 
#    # log(p(data)) = log( multinomial coefficient ) + x1* log( p1 ) + .. + xk * log ( pk )
#    
#    my ( $h, $sh, $ph ) = @_;
#    my $sdp = sum( @{$sh}{( keys %$sh )}); #N
#    my $multinom_coeff = Math::BigFloat->new($sdp);
#    $multinom_coeff->bfac()->blog(); #numerator
#    for my $k ( @{$sh}{(keys %$sh)} ) { 
#        #divide by each term in the denominator (subtract in log space)
#        my $d = Math::BigFloat->new($k);
#        $d->bfac()->blog();
#        $multinom_coeff->bsub($d);
#    }
#    
#    my $prob = Math::BigFloat->bzero();
#    for my $k ( keys %$sh ) {
#        my $int_prob = Math::BigFloat->new( $sh->{$k} );
#        $int_prob->bmul( $ph->{$k} );
#        $prob->badd( $int_prob );
#    }
#    $prob->badd( $multinom_coeff );
#    
#    return $prob;
#}

#sub locus_distribution {
#    my ( $opts, $mode, $lh ) = @_; 
#    &debug("distribution_shape") if $opts->{'--debug'};
#    
#    my $return = 0; # default to false
#    my @keys = ( min(keys %$lh) - 1 .. max( keys %$lh ) + 1 );
#    map { $lh->{$_} ||= 0 } @keys; #backfill the locus hash with zeros 
#    my $loc_dp = sum( @{$lh}{@keys} );
#    # check for decreasing as go up from mode
#    my @upper_array;
#    for my $k ( $mode .. $keys[-1] ) {
#        push @upper_array, $lh->{$k};
#    }
#    print STDERR join(", ", @upper_array), "\n";
#    #    my $upper_call = &monotonic($opts, \@upper_array);
#    # check for decreasing as go down from mode
#    my @lower_array;
#    for my $k ( sort {$b<=>$a} ( $keys[0] .. $mode ) ) { 
#        push @lower_array, $lh->{$k};
#    }
#    print STDERR join(", ", @lower_array), "\n";
#    #    my $lower_call = &monotonic($opts, \@lower_array);
#    
#    #    if ( $upper_call == 1 && $lower_call == 1 ) {
#    #        $return = 1;
#    #    } else { 
#    #        $return = 0;
#    #    }
#    return $return;
#    
#    
#}
#
#sub monotonic {
#    # returns true if the input values are monotonic decreasing within thresholds
#    
#    my ($opts, $array) = @_;
#    &debug("monotonic") if $opts->{'--debug'};
#    
#    print STDERR join(", ", @$array), "\n";
#    my $return = 1;
#    for my $i ( 0 .. $#$array - 1 ) {
#        if ( all { $array->[$_] < $opts->{'--noise_threshold'} } ( $i  .. $#$array ) ) {
#            last;  # all remaining values are below noise threshold
#        } elsif ( $array->[$i] > $array->[$i+1] ) {
#            next; # decreasing, ok
#        } elsif ( $array->[$i] < $array->[$i+1] ) {
#            # increased, check to see if met threshold
#            if ( $array->[$i+1] - $array->[$i] > $opts->{'--peak_height'} ) {
#                $return = 0; # determine as secondary peak, fail
#                last; 
#            } elsif ( $array->[$i+1] - $array->[$i] <= $opts->{'--peak_height'} ) {
#                next; # below threshold, ok
#            }
#        }
#    }
#    return $return;
#}
