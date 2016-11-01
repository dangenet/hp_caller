#! /usr/bin/perl

use warnings;
use strict;

# feed this a refseq fasta file on stdin 
# returns a bed-format text file of locations of all homopolymers 
# with lengths btw 6 & 35
# field 4 is the "name" entry, in the format T:12 (a 12 T run)
# field 5 is the "score" entry, 28*length 
# field 6 is + for top strand


my %genome;

&main;
exit;

sub main {
    &chromosomes; #read the genome into memory
    my $bedlines_ref = &find_polyN;
    print join("\n", @$bedlines_ref), "\n"; 
}

sub chromosomes {
    my $chromosome_number = '';
    my $chromosome;
    my $stupid_flag = 0;
    while (<>){
    	chomp $_; 
        if ( $_ =~ /^>/ && $stupid_flag == 0) {
            $chromosome_number = $& if ( $_ =~ /^>chr\d+|chrMito|pYES2_AAGfl/);
            $chromosome_number =~ s/^>//g; 
            $stupid_flag = 1;
        } elsif ($_ !~ /^>/ ) {
            $chromosome .= $_;
        } elsif ( $_ =~ /^>/ && $stupid_flag == 1 ) {
            $genome{ $chromosome_number } = $chromosome;
            $chromosome_number = $& if ( $_ =~ /^>chr\d+|chrMito|pYES2_AAGfl/);
            $chromosome_number =~ s/^>//g; 
            $chromosome = '';
        }
    }
    close (ARGV);
    $genome{ $chromosome_number } = $chromosome;
}

sub find_polyN {
    my $bedlines_ref; 
    for my $k ('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',  'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16' ) { # removed pYES2 and chrMito
        my %chromosome_bedlines = ();
        for my $match_string ('A', 'T', 'G', 'C') {
            while ( $genome{$k} =~ /(($match_string){4,})/g ) {
                my $index = pos $genome{$k}; 
                my $length = length($1); 
                my $start = $index - $length; 
                # bed line format, 0 based counting, end not included in feature
                # 1chrom 2chromStart 3chromEnd 4name 5score 6strand  
                # score = 28*length, which would color longer hps darker in UCSC
                my $str = join("\t", $start, $index, "$match_string:$length", $length * 28, '+' ); 
                $chromosome_bedlines{ $start } = $str; 
            }
        }
        foreach ( sort { $a <=> $b } ( keys %chromosome_bedlines ) ) {
            my $bedline = "$k\t" . $chromosome_bedlines{$_}; 
            push @$bedlines_ref, $bedline; 
        }
    }
    return $bedlines_ref; 
}
