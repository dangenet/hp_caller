# hp_caller

Requirements:

samtools v0.1.19+

Perl v5.12
- Math::Round
- Math::BigInt
- Math::BigFloat

1) Index your reference fasta sequence using samtools faidx. 

samtools faidx REFSEQ.fasta 

2) Find the homopolymers in your reference sequence. You can use the perl script I wrote (hp_finder_BED.pl), or you can use something else. The resulting bed file should have the tab-separated columns CHROM START END HP SCORE STRAND where "HP" is A:8 for the sequence AAAAAAAA It doesn't matter what you put for the score.

hp_finder_BED.pl < REFSEQ.fasta > HOMOPOLYMERS.bed 

3) Run the aggregator script, hp_aggregator.pl. 

hp_aggregator.pl --bedfile HOMOPOLYMERS.bed --refseq REFSEQ.fasta --bamfiles ./*.bam > DATA_FILE.agg

Note that hp_aggregator.pl keeps everything in memory before writing the output data. If you have a large genome or many bamfiles you may run out of memory. Break your bed file up into several smaller bed files if you have this problem. You can put all the variant calls into one file at the end.  

4) Run hp_caller to call variants. Minimal syntax is:

hp_caller.pl -v DATA_FILE.agg > VARIANTS.vcf 

An example command line with all option parameters specified: 

hp_caller.pl -v DATA_FILE.agg --min_sdp 5 --max_sdp 150 --nsr 1 --uncallable 3 --loc_min 200 --loc_max 2500 -mp 5 --asymmetry 0.2 --minSGQ 10 --minQUAL 10 --variants_only > VARIANTS.vcf

Meaning of the options for hp_caller:

--min_sdp [integer >=3]
	Describes the minimum depth requirement for a sample to be callable.

--max_sdp [integer] 
	Describes the maximum sample read depth for a sample to be callable. With many reads at a locus, a binomial probability calculation is not required; the sample can simply be called according to the mode of its distribution. The time required for the calculation scales with the number of reads in the sample. Additionally, floating point errors may result if sample depth is too high. 

--nsr [integer]
	Net supporting reads, the number of reads supporting one hypothetical length vs. the next most likely hypothetical length. Default is 1; is redundant with sample QUAL score. If set to 5, there must be 5 more reads supporting the hypothesis that the homopolymer is of length L1 vs. than 
the hypothesis that the homopolymer is of length L2. 

--uncallable [integer]
	The maximum number of uncallable samples at a locus before the entire locus is considered uncallable. 

--loc_min [integer]
	The minimum number of reads required for a locus to be considered callable. Recommended not to be set below 50. 

--loc_max [integer] 
	The maximum number of reads at which to attempt calling the locus; default 2500. This is used to avoid calling at repetitive loci with artificially elevated coverage with respect to the rest of the genome. 

-mp [integer]
	PHRED scaled probability that the sample has a different genotype than the locus. Default 5, or 10^-5. 

--asymmetry [fraction greater than zero and less than 1]
	Maximum fractional asymmetry around the mode in a distribution. Default 0.2. If a sample or locus distribution has greater asymmetry than this threshold, the sample or locus will not be called.

--variants_only
	If this option is specified, only loci with variants will be output. By default all loci are  output regardless of whether a variant is present. 

--minSGQ [integer]
	Default 10. PHRED-scaled minimum sample distribution quality score for callable samples. Samples with quality scores less than this threshold will not be called.

--minQUAL [integer]
	Default 10. PHRED-scaled minimum locus distribution quality score for callable loci. No samples will be called at loci with quality scores less than this threshold. 

************
These scripts are intended as proof of concept rather than as tools ready to drop into a pipeline. 


You can reach me at @dangenet on Twitter or github. 

