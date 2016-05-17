#! /usr/bin/perl

# note that this was written for samtools v 0.19+ and will probably not work right with the htslib-based samtools

# this script analyzes homopolymer-spanning reads in a list of bamfiles 
# output is vcf-like format, containing per sample histograms for each locus provided in the bam file
# required inputs
# --bedfile    bed format file containing list of homopolymers
# --refseq     fasta format file containing reference sequences
# --bamfiles   list of bam files to use (will output one sample field per bam file), I use globs
# also outputs a pattern length report in text format, which tells you about the distribution
# of the lengths of flanking sequences used to determine where homopolymers end
# bed file should be CHROM START END HP SCORE STRAND where 
# HP is A:8 for the sequence AAAAAAAA

use warnings;
use strict;
use Data::Dumper; 
use Module::hpcaller qw(getargs simple_histogram find_hist_max hp_converter); 
use Time::HiRes qw(gettimeofday); 
use List::Util qw(max min sum); 
use List::MoreUtils qw(uniq none indexes firstidx any); 

&help;

sub help {
    my $helpmsg = <<'END_HELP'; 
    
Usage: hp_aggregator.pl --bedfile [required] --refseq [required] --bamfiles [required][whitespace separated list]
    
    note that this was written for samtools v 0.19+ and 
        may not work right with newer samtools
    
    this script analyzes homopolymer-spanning reads in a list of bamfiles 
    output is vcf-like format, containing per sample histograms for 
        each locus provided in the bam file
    required inputs
    --bedfile    bed format file containing list of homopolymers
    --refseq     fasta format file containing reference sequences
    --bamfiles   list of bam files to use (will output one sample 
        field per bam file), I use globs
    also outputs a pattern length report in text format, which tells you 
        about the distribution of the lengths of flanking sequences used
        to determine where homopolymers end
    bed file should be CHROM START END HP SCORE STRAND where 
    "HP" is A:8 for the sequence AAAAAAAA
    It doesnt matter what "SCORE" is. 
    
END_HELP
    print $helpmsg;
    exit if ( any { $_ =~ /-h/ } @ARGV ) ; 
    exit if scalar @ARGV == 0;
}

&main;
exit;

sub main {
    my $opts_ref = &getargs; 
    my %opts = %$opts_ref; 
    my $bedfile = $opts{'--bedfile'} or &help; # the bedfile to be used
    my $refseq = $opts{'--refseq'} or &help; 
    
    # progress monitoring stuff
    my $start_time = gettimeofday(); 
    my $total_bamfiles = scalar(@{$opts{'--bamfiles'}}); 

	# read the bed lines into memory
    my $all_bedlines; 
    open(my $BED, $bedfile);
    while (<$BED>) {
        chomp $_; 
        push @$all_bedlines, $_; 
    }
    close ($BED); 

    # collect the data in a big hash
    my $results = ();
    my $bedline_counter;
    my $oldfrac = -1; 
    my @samples; 
    print STDERR "processing bed lines... "; 
    for my $bedline (@$all_bedlines) { #go through the bed file line by line
        my $bamfile_counter = 0; 
        # print STDERR "\nprocessing bamfiles... ";
        for my $bamfile ( @{$opts{'--bamfiles'}} ) { # go through each bam file 
            my $sn = $& if ($bamfile =~ /\d+/); # extract the sample number from the bam file name
            push @samples, $sn;
            @samples = uniq @samples; 
            my ($c, $s, $e, $hp, $score, $strand) = split("\t", $bedline);
            my @bedline_arr = split("\t", $bedline);
            $s++; 
            # collect data from the bamfile
            %{$results->{$hp}->{$c}->{$s}->{$sn}} = &bam_processor($bamfile, \@bedline_arr, $refseq);
            # if there are no spanning reads, the only entry will be sn->reference hp length= 0
            $bamfile_counter ++; 
        }
        # progress monitoring stuff 
        $bedline_counter ++;
        my $frac = sprintf("%d", 100*$bedline_counter/(scalar @$all_bedlines));
        print STDERR "$frac%... "  if ( $frac > 1 + $oldfrac);
        $oldfrac = $frac if ($frac > 1 + $oldfrac); 
    }
    print STDERR "\n";

    my $sample_vcf = &by_sample( $all_bedlines, $results, \@samples );
    
    print "##BAM_FILES:", join(';', @{ $opts{'--bamfiles'} } ), "\n";
    print "##BED_FILE:", $bedfile, "\n";
    &vcf_printer ( $all_bedlines, $sample_vcf, \@samples, $opts_ref ); 
    
#    my ($strain_vcf, $reference_dist, $trouble_spots) = &by_locus($all_bedlines, $results);
#    print Dumper($results); 
    
    print STDERR "\n\n", gettimeofday()-$start_time, "seconds runtime\n";
    return; 
}

sub vcf_printer { 
    my ( $bedlines, $vcf, $samples, $opts ) = @_ ; 
#    print Dumper($vcf); 
    my @field_names = qw(CHROM POS ID REF ALT QUAL FILTER INFO FORMAT); 
    my @all_fields;
    map { push @all_fields, $_ } @field_names; 
    map { push @all_fields, $_ } @$samples; 
    
    print "#", join("\t", @all_fields), "\n"; 
    for my $bedline ( @$bedlines ) {
        my ($c, $s, $e, $hp, $score, $strand) = split("\t", $bedline);
        $s ++; 
        my $loc = "$c\t$s"; 
        print join("\t", @{ $vcf->{$loc} }{ @all_fields } ), "\n"; 
    }
    return 1; 
}

sub by_sample {
    my ( $bedlines, $data, $samples ) = @_; 
    @$samples = sort { $a <=> $b } @$samples; 
    my $probs; # container for counts
    for my $bedline ( @$bedlines ) {
        my ( $c, $s, $e, $hp, $score, $strand ) = split("\t", $bedline);
        $s++; 
        $e++; 
        my ( $nt, $ref_l ) = split(":", $hp); 
        my $hp_str = &hp_converter( $hp ); 
        $probs->{"$c\t$s"}->{"CHROM"} = $c; 
        $probs->{"$c\t$s"}->{'POS'} = $s;
        $probs->{"$c\t$s"}->{'ID'} = ".";
        $probs->{"$c\t$s"}->{'REF'} = $hp; 
        
        # look at the samples individually 
        my $locus_hash; 
        for my $sn ( @$samples ) {
            # add the sample data to the locus histogram
            my @lengths;
            map { push @lengths, $_ } ( keys %{ $data->{$hp}->{$c}->{$s}->{$sn} } );
            map { $locus_hash->{$_} += $data->{$hp}->{$c}->{$s}->{$sn}->{$_} } @lengths; 
            # make a copy
            my %slice; 
            @slice{@lengths} = @{ $data->{$hp}->{$c}->{$s}->{$sn} }{ @lengths }; 
            # check your inputs
            if ( ( max @slice{@lengths} ) == 0 ) { 
                # happens when no reads span the homopolymer
                $probs->{"$c\t$s"}->{$sn}->{'depth'} = 0 ;
                $probs->{"$c\t$s"}->{$sn}->{'hist'} = "$ref_l,0";
            } elsif ( ( max @slice{@lengths} ) > 0 ) {
                # calculate the depth
                $probs->{"$c\t$s"}->{$sn}->{'depth'} = sum( @slice{@lengths} ); 
                # construct the histogram string
                my @k = sort { $a <=> $b } ( keys %slice ); 
                my @v = ();
                map { push @v, "$_,$slice{$_}" } @k; 
                $probs->{"$c\t$s"}->{$sn}->{'hist'} = join(";", @v);
            }
            $probs->{"$c\t$s"}->{$sn} = $probs->{"$c\t$s"}->{$sn}->{'depth'} . ":" . $probs->{"$c\t$s"}->{$sn}->{'hist'};
        }
        
        # do something with @locus_lengths
        my @locus_hist = map { join(",", $_, $locus_hash->{$_} ) } ( sort {$a <=> $b} keys %$locus_hash ); 
        my $H = join(";", @locus_hist); 
        my $total_depth = sum( @{ $locus_hash }{ keys %$locus_hash } ); 

        @locus_hist = ();
        $locus_hash = {}; 
                
        $probs->{"$c\t$s"}->{'ALT'} ||= "."; 
        $probs->{"$c\t$s"}->{'QUAL'} =  "."; # not sure what to define this as 
        $probs->{"$c\t$s"}->{'FILTER'} = "."; 
        $probs->{"$c\t$s"}->{'INFO'} = "DP4=$total_depth:HIST=$H"; # total info for loci
        $probs->{"$c\t$s"}->{'FORMAT'} = "DP:HS"; 
        # genotype, mode(s), # of modes, mode fraction, sample depth
        $H = ''; 
    }
    return $probs; 
}

sub bam_processor {
    my ($bamfile, $bedline, $ref_file) = @_; 
#    my $bamfile = shift @_;
#    my $ref_file = pop @_; 
#    my $bedline = shift @_; 
    my $chr = $bedline->[0];
    my $start = $bedline->[1]+1;
    my $end = $bedline->[2]; 
    my ($min, $max) = ($start -100, $end +100); 
    $min = 1 if $min < 1; 
    my $us_length = $start - $min;
    
    # retrieve the reference sequence 100 nt before and after the homopolymer
    my $ref_cmd = "samtools faidx $ref_file $chr:$min-$max |";
    open(my $REF, $ref_cmd);
    my $ref_seq = '';
    while (<$REF>) {
        chomp $_;
        next if $_ =~ /^>/; 
        next if $_ eq ''; 
        $ref_seq .= $_;
    }
    close( $REF); 
    die if not defined $ref_seq; # shouldn't happen
    # define the sequences upstream and downstream of the homopolymer
    my $upstream = substr( $ref_seq, 0, $us_length ); 
    my $homopolymer = substr( $ref_seq, $us_length, $end-$start+1);
    my $downstream = substr( $ref_seq, $us_length + length( $homopolymer ), 100);
    # catch the error if the homopolymer happens to be at the end of the chromosome
    if ( $upstream eq '' || $downstream eq '' ) {
        print STDERR "\n", join("\t", @$bedline), "\n", "at end of chromosome, returning 0 depth\n";
        return ( length($homopolymer) => 0 ); 
    }
    # pull out the sequences that flank the homopolymer
    my $start_seq = substr($upstream, -4, 4);
    my $end_seq = substr($downstream, 0, 4); 
    my $n = substr($homopolymer, 0, 1); 
    # want the pattern
    #     /(upstream_sequence)(homopolymer nucleotide){1,}(downstream sequence)/
    # to be unique in the region being examined,
    # and to not be vulnerable to readthrough of the terminating nucleotide on either side
    # Hence, the next chunk is dedicated to finding good flanking sequences 
    # expand the end match sequences until they contain at least 3 of the 4 nucleotides
    my $start_3 = substr( $start_seq, -1, 1 ); 
    while ( $start_seq !~ /[^$n$start_3]/ ) {
        # while the start string does not contain any nucleotides other than the homopolymer nucleotide and the terminating nucleotide 
        last if length($start_seq) >= length( $upstream); 
        $start_seq = substr($upstream, -1*length($start_seq) -1, length($start_seq) + 1); 
    }
    my $end_5 = substr( $end_seq, 0, 1); 
    while ( $end_seq !~ /[^$n$end_5]/ ) {
        last if length($end_seq) >= length($downstream); 
        $end_seq = substr( $downstream, 0, length($end_seq) + 1 ); 
    }
    # extend the start and end patterns if they lead to off-target matches in the upstream and downstream regions
    my $start_count = () = $upstream =~ /(?=$start_seq)/g; 
    my $start_miscount = () = $upstream =~ /(?=$start_seq)$n+/g; # 0 when unique
    my $end_count = () = $downstream =~ /(?=$end_seq)/g; 
    my $end_miscount = () = $downstream =~ /$n+(?=$end_seq)/g; # 0 when unique
    my $match_count = () = $ref_seq =~ /($<=$start_seq)$n+(?=$end_seq)/g;
    if ($start_miscount > 1 || $end_miscount > 1 || $match_count > 1 ) {
        # increase the length of the start sequence until it is unambiguous in the upstream region
        my $i = length($start_seq); 
        while ( $start_miscount > 1 ) {
            $i ++; 
            $start_seq = substr($upstream, -1*$i, $i ); 
            $start_count = () = $upstream =~ /(?=$start_seq)/g;
            $start_miscount = () = $upstream =~ /(?=$start_seq)$n+/g;
        }
        $i = length($end_seq); 
        # increase the length of the end sequence until it is unambiguous in the downstream region
        while ( $end_miscount > 1 ) {
            $i++; 
            $end_seq = substr($downstream, 0, $i); 
            $end_count = () = $downstream =~ /(?=$end_seq)/g;
            # this will overcount potential readthrough events
            $end_miscount = () = $downstream =~ /$n+(?=$end_seq)/g; 
        }
        $match_count = () = $ref_seq =~ /($<=$start_seq)$n+(?=$end_seq)/g;
        $i = min( length( $start_seq), length($end_seq) );
        # increase the length of both flanking sequences until the whole pattern is unambiguous
        while ( $match_count > 1 ) {
            $i ++;
            $start_seq = substr($upstream, -1*$i, $i) if $i > length( $start_seq);
            $end_seq = substr($downstream, 0, $i) if $i > length( $end_seq ); 
            $match_count = () = $ref_seq =~ /($<=$start_seq)$n+(?=$end_seq)/g;
        }
    }
    # seems I need to re-determine the match count here otherwise it incorrectly comes up as zero, not sure why
    $start_count = () = $ref_seq =~ /(?<=$start_seq)$n+/g;
    $end_count = () = $ref_seq =~ /$n+(?=$end_seq)/g; 
    $match_count = () = $ref_seq =~ /(?<=$start_seq)$n+(?=$end_seq)/g; 
    my $n_count = () = $ref_seq =~ /$n+/g; 
    if ( $match_count == 0 ) {
        die "match count was zero\n";
    }
    # test explicitly to see if readthrough of the terminating nucleotide can lead to incorrect lengths
    # die if it can
    my $ref_sub = join('', substr($upstream, -30, 30), $homopolymer, substr($downstream, 0, 20 )); 
    my $us_copy = $upstream;
    substr( $us_copy, -1, 1 ) = $n; 
    my $rdthr = join('', $us_copy, $homopolymer, $downstream); 
    if ( $rdthr =~ /($<=$start_seq)$n+(?=$end_seq)/ ) {
        print STDERR "\n", join("\t", @$bedline), "\n";
        print STDERR $ref_seq, "\n";
        print STDERR $ref_sub, "\n";
        print STDERR $rdthr, "\n";
        print STDERR $start_seq, "\t", $end_seq, "\n";
        print STDERR "readthrough\t", $&, "\n"; 
        die;
    }
    my $ds_copy = $downstream;
    substr( $ds_copy, 0, 1) = $n; 
    $rdthr = join('', $upstream, $homopolymer, $ds_copy);
    if ( $rdthr =~ /($<=$start_seq)$n+(?=$end_seq)/ ) {
        print STDERR "\n", join("\t", @$bedline), "\n";
        print STDERR $ref_seq, "\n"; 
        print STDERR $ref_sub, "\n";
        print STDERR $rdthr, "\n";
        print STDERR $start_seq, "\t", $end_seq, "\n";
        print STDERR "readthrough\t", $&, "\n"; 
        die;
    }
    # write the lengths of the flanking sequences to the report text file
    open( my $LREPORT, ">>", "pattern_length_report.txt");
    print $LREPORT join("\n", length( $start_seq ), length( $end_seq ) ), "\n"; 
    close( $LREPORT ); 
    
    
    #moving on. retrieve reads from the bam files and count homopolymer lengths
    my $bam_cmd = "samtools view -f 2 -F 1028 $bamfile $chr:$start-$end |"; 
    # -f 2 requires reads be properly paired
    # -F 1284 excludes duplicates (1024), alignment is not primary (256), unmapped (4)
    my $lengths = (); 
    open (my $BAM, $bam_cmd); 
    while (<$BAM>) {
        chomp $_; 
        my @bamline = split ("\t", $_); 
        next unless ( $bamline[3] <= $start - 4  && $bamline[3] +100 >= $end + 4 ) ; # must span locus
        next unless ( $bamline[9] =~ /$start_seq$n+$end_seq/ ); # must match pattern
        
        my $read_seq = $bamline[9];

        if ( $read_seq =~ /$start_seq($n{3,})$end_seq/ ) {
            # add to list of observed lengths unless the length is 0 
            push @$lengths, length($1) unless ( length( $1 ) == 0 ); 
        }
    }
    close ($BAM); 
    # make a histogram
    my $hist = &simple_histogram($lengths);
    if ($hist == 1) {
        # returns a 0 count for the reference length
        $hist = { length($homopolymer) => 0 };  
    }
#    for my $l ($min_length..$max_length) {
#        $hist->{$l} ||=0; 
#    }
    return %$hist; 
}

