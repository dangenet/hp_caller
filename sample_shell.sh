#!/bin/bash

#run in a directory containing your reference genome sequence as a fasta file and your bam files
#developed using samtools 0.19+, untested with the htslib samtools 
#2016-11-01 D. Eyler

R="your_reference_genome.fasta"
B="put_the_name_you_want_your_homopolymer_bed_file_to_have_here.bed"

hp_finder_BED.pl < $R > $B 

hp_aggregator.pl --bedfile $B --refseq $R --bamfiles ./*.bam > "name_of_raw_data_file.agg"

hp_caller.pl -v "name_of_raw_data_file.agg" > "name_of_output_file.vcf"
