# hp_caller

Step 1 for calling: Run hp_aggregator.pl. Like this:

hp_aggregator.pl --bedfile [required] --refseq [required] --bamfiles [required][whitespace separated list]

Your refseq should be in fasta format and have been indexed with samtools faidx. It will 
look for the .fai file. 

Your bedfile should have the columns CHROM START END HP SCORE STRAND where 
    "HP" is A:8 for the sequence AAAAAAAA
It doesn't matter what you put for the score.

I used bowtie2 to make my bam files and picard to mark duplicates. 

Run it and direct STDOUT to the file of your choice. Note I wrote hp_aggregator pretty 
early and it keeps everything in memory and spits it out in the end. 

Step 2: Run hp_caller.pl -v [output file from step 1] and direct STDOUT to the file of 
your choice. Run hp_caller.pl -h for a list of options with brief descriptions. 

That's it. The help in each script gives a little more information, and you can run 
head -n50 on each of them to get their dependencies... Math::BigFloat, and
Math::BigInt are the ones that you might not already have. 

These scripts really are intended as a proof of concept for the way we call homopolymers
and not as tools ready for someone to drop into their pipeline. That said you can reach
me @dangenet on twitter or on github. 