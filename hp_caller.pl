#! /usr/bin/perl

# this is going to call samples relative to the locus without paying attention to the rest of the genome
# this is hp_caller3 renamed hp_caller for publication

use strict;
use warnings;
use Data::Dumper;
use My::Module qw( getargs string_to_hist hist_to_string hp_converter); 
use List::Util qw( sum max min ); 
use List::MoreUtils qw( true any none uniq all firstidx ); 
use Scalar::Util qw(looks_like_number);
use Math::Round qw( nearest_ceil nearest ); 
use Math::BigInt; 
use Math::BigFloat;

my $opts = &getargs;

&main;
exit;

sub main {
    &help if defined $opts->{'-h'};

    my $agg = $opts->{'-v'};

    if ( not -e $agg || not defined $opts->{'-v'} ) {
        &help && die "input file required with -v\n" ;
    }
    $opts->{'--min_sdp'} ||= 5;
    $opts->{'--max_sdp' } ||= 150;
    $opts->{'--nsr'} ||= 1;
    $opts->{'--uncallable'} ||= 3;
    $opts->{'--loc_min'} ||= 200;
    $opts->{'--loc_max'} ||= 2500;
    $opts->{'-mp'} ||= 5; #PHRED scaled prob of a sample being different from locus
    $opts->{'--asymmetry'} ||= 0.2; #maximum fractional asymmetry around the mode in a distribution
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
            # data line
            # extract the data and store as a hashref
            my $lh = &locus_hashmaker( $opts, $_ );
            
            #really probably want to exclude 22582 from locus histogram
            # since it can actually be heterozygous
            # delete the 22582 histogram data 
            map { $lh->{"INFO"}->{"HIST"}->{$_} -= $lh->{22582}->{"HS"}->{$_} } (keys %{$lh->{'22582'}->{"HS"}});
            
            # pass the hashref to the caller
            &locus_caller( $opts, $lh );
        }
    }
    close $AGG_FILE; 
    
    return 1; 
}


# these subroutines extract data from the input line
sub locus_hashmaker {
    my ($opts, $line) = @_;
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

sub sample_lh_parser {
    my ( $format, $str ) = @_;
    my $sh; 
    my @formats = split(":", $format);
    @{$sh}{ @formats } = split(":", $str) or die "$str";
    $sh->{"HS"} = &string_to_hist($sh->{"HS"}) unless $sh->{"HS"} eq "na";
    return $sh; 
}

# this subroutine is run on every input line (locus)
sub locus_caller {

    my ($opts, $lh) = @_; 
    
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
    my $locus_flag = &callable_locus( $locus_hist, \@samples, $opts );
    
    #make a hashref that contains only the mode of each sample distribution
    my $loc_modes; 
    map { $loc_modes->{$_} = &just_modes($locus_hist->{$_}); }( "LOCUS", @samples);
    
    #fill in output values
    $lh->{"INFO"}->{"MODE"} = join(",", @{ $loc_modes->{"LOCUS"} } ); #mode of the locus
    $lh->{"QUAL"} = &locus_qual_score($opts, $locus_hist->{"LOCUS"} ); #locus qual score
    
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
    
    # for lines containing potentially variant samples
    # calculate probabilities and mark variant samples for those 
    my $vsf; #variant sample flag hashref
    map { $vsf->{$_} = 0 } @samples; # set all to nonvariant initially
    my $sph; #sample probability hash
    map { $sph->{$_} = "na"; } @samples; #prefill with error values

    if ( $pvf == 1 && $locus_flag == 1 ) {
        for my $s ( @samples ) {
            next if $s == 22582; #skip the mutant sample
            
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
    
    return 1; 
}

# the following subroutines are called by locus_caller and are arranged in the order in which they are called
sub callable_locus {
    # applies criteria for determining locus is callable
    # returns boolean
    my ( $lh, $samples, $opts ) = @_;
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
    if ( $locus_flag == 1 ) {
        $locus_flag = &locus_depth_checker( $loc_dp );
    }
    return $locus_flag;
}

sub just_modes {
    my ( $hr ) = @_;
    my @keys = ( sort {$a<=>$b} ( keys %$hr ));
    my $max_val = max( @{$hr}{@keys} );
    my @modes;
    map { push @modes, $_ if $hr->{$_} == $max_val } @keys;
    return \@modes;
}

sub mode_checker { #returns true if unimodal
    my ($mr) = @_;
    my $return = 0;
    if ( scalar @$mr == 1 ) {
        $return = 1;
    }
    return $return;
}

sub asymmetry {
    # checks to see if the locus has a peak at the mode with decreasing counts going out from there
    my ( $opts, $mode, $lh ) = @_; 
    
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
    return $return; 
}

sub callable_sample_check { 
    #determines whether the number of uncallable samples at a locus exceeds the threshold
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

sub locus_qual_score {
    my ( $opts, $lh ) = @_; #opts, #locus hash histogram
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

sub locus_qual_check { 
    my ( $opts, $lh ) = @_;
    my $return = 0;
    my $qual = &locus_qual_score( $opts, $lh);
    if ( $qual >= $opts->{'--minQUAL'} ) {
        $return = 1;
    }
    return $return;
}

sub locus_depth_checker {
    my ( $dp ) = @_;
    my $return = 0;
    if ( $dp >= $opts->{'--loc_min'} && $dp <= $opts->{'--loc_max'} ) {
        $return = 1;
    }
    return $return;
}

#math is done here
sub prob_sample_came_from_locus { 
    my ($h, $sh, $tlh, $opts) = @_; #sample hashref, #total locus histogram
    
    my $lh; #subtract the sample values from the total locus histogram
    map { $sh->{$_} ||= 0; $lh->{$_} = $tlh->{$_} - $sh->{$_} } ( keys %$tlh );
    
    my $prob = &binomial_probability( $h, $sh, $lh );
    if ( $prob > 1 || $prob < 0 ) {
        print STDERR Dumper($h, $sh, $tlh);
        print STDERR $prob, "\n";
        die;
    }
    return $prob;
}

sub binomial_probability {
    my ( $h, $sh, $lh ) = @_; #success key, sample hist, locus hist
    #construct the probability table for probability of h, prob of not h
    my $ph;
    
    $ph->{$h} = Math::BigFloat->new( $lh->{$h} );
    
    my $loc_dp = sum( @{$lh}{(keys %$lh)} );
    
    $ph->{$h}->badd(1)->bdiv($loc_dp + 1); 
    $ph->{"other"} = Math::BigFloat->new( $loc_dp );
    $ph->{"other"}->bsub( $lh->{$h} )->badd(1);
    $ph->{"other"}->bdiv($loc_dp + 1);
    
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

#printing is done here
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

sub help  {
    
    my $help = <<'END_HELP';

USAGE: hp_caller.pl -v [data file from hp_aggregator]

    -v  data file from hp_aggregator (required)
    --min_sdp   
        minimum sample depth for calling, default 5
    --max_sdp   
        maximum sample depth for calling, default 150
    --nsr       
        net supporting reads or the minimum number of excess
        reads supporting a mutant call, default 1
    --uncallable
        maximum number of uncallable samples for a locus
        to be considered callable, default 3
    --loc_min   
        minimum reads for a locus to be considered callable,
        default 200
    --loc_max   
        maximum reads for a locus to be considered callable,
        default 2500
    -mp         
        PHRED-scaled threshold for calling a sample as mutant,
        default 5
    --asymmetry 
        maximum asymmetry of the locus distribution around its
        mode, as fraction of total reads, default 0.2
    --variants_only 
        print only variants, default 0 (print all)
    --minSGQ
        minimum sample genotype quality required for a sample
        to be considered callable, default 10
    --minQUAL
        minimum locus quality required for locus to be considered
        callable, default 10
            
    Both QUAL scores are:
       100 * ( mode counts - 2nd most abundant counts) / total depth

END_HELP
    
    print STDERR $help;
    exit 0; 
}







