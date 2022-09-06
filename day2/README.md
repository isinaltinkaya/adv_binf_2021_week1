# Mapping Exercise
In this exercise you will get acquainted with the "bwa" mapper and the "samtools" program to interpret sam files. 

We will use a subset of the human genome, chromosome 21, as a reference genome, to keep the CPU load and waiting time lower. 
The sequencing files you will need to map are single-end NGS data of a human individual.

## The data are located on ricco.popgen.dk at /ricco/data/ex3/
~~~bash
NA19238.fastq.gz
chr21.fa.gz
~~~
With the first being the fastq file of the individual and the last one being the reference genome for the human chromosome 21 in fasta format.

All the files are compressed using gzip. You can view the first few lines using the
command zmore, which is like more, but decompressing the file first.

Before we start its important to note that you should try to run the command with no arguments - the help page will get an explanation of the way run the command.

You should all copy, using cp, the reference file (chr21.fa.gz) to your own home directory and build the index yourself.

## Indexing for BWA alignment
Now you should map the reads in the fastq files against the reference genome. First, you
have to make the index (the Burrows-Wheler transform) and then run the mapping. 

You use the bwa index command to create FM-index of the reference genome, such that we can align our paired end query sequences to this index. 
You just need one argument in this case.

~~~bash
bwa index reference_file
~~~

NOTE: During the index construction we see the "[bwa_index] Construct SA from BWT and Occ..." 

Once the index is made, the second step is to map the reads, the aligners are both capable of alinging paired-end and single-end sequencing.

# Mapping
We will try multiple different alignment methods all based on the Burrows-Wheeler tranform. Namely the different approaches of bwa

If you run the bwa alignment commands with no arguments to get info about how to use it, then the number
of options may be a bit overwhelming, but you can run it with no additional options.

To run with threads we need to consider the CPU load so I suggest you add "-t 2" or "-t 4" to have it run in 2 or 4 threads if your computer
has multiple cores. 

All of the commands reads the compressed fastq files directly, so you need not
decompress them. By default the result comes on stdout (in the terminal), so you have to redirect to a file using ">". 

## BWA ALN
First we will find the coordinates of our input reads files in the Suffix Array (SA) created in the index files

The help page for bwa can be found at http://bio-bwa.sourceforge.net/bwa.shtml and further information can be found by simply typing bwa in your terminal

Q1) Which command do we need in order to generate the Suffix Array given the input reads?

/~~~bash
/bwa aln -t 2 reference_file readfile_1 > readfile1_SA.sai
/~~~

Q2) Using the generated SA intervals, perform single end alignment using the reference genome and the provided sequencing reads
/Then we will generate alignments in the SAM format given the single-end read by searching in the SA intervals
/~~~bash
/bwa samse reference_file readfile1_SA.sai readfile_1 > output_aln.sam
/~~~

## BWA MEM
We can also perform alignment using bwa mem
the bwa mem mode, which is the most commonly used these days
~~~bash
bwa mem -t 2 reference_file readfile_1 > output_mem.sam
~~~

# SAM/BAM output
Now take a look at the sam file – the output – and see if it makes sense - with the examples below using the output from bwa aln. 

Several tools have been created to handle sam files, we will use the samtools program. 

## File conversion
We often need the compressed version of the sam format, which is called bam. To do this we use "samtools view" for converting between formats. Try
the command without arguments to see the options.

~~~bash
samtools view -Sb output_aln.sam > output_aln.bam
~~~
## QUESTION
1. Try to convert sam file to .cram file (a different .sam conversion)
2. What is the difference between the commands creating the .bam and .cram
~~~bash
samtools view -T chr21.fa.gz -SC output_aln.sam > output_aln.cram
~~~
Creating the cram input file need the "-C" parameter which requires the "-T" parameter.
The "-S" parameter is also described below.
~~~bash
-C       output CRAM (requires -T)
-T, --reference FILE
-S       ignored (input format is auto-detected)
~~~

## Sorting
Often we need a file sorted on genomic coordinates. You can use "samtools sort" for this
~~~bash
samtools sort -@ 2 output_aln.bam -o output_aln.sorted.bam
~~~
Try to use the "samtools view" to look at the bam files, with the 3rd and 4th column being the chromosomal position, so you can check if it is
sorted

## QUESTION
1. Based on column 3 and 4, have we succesfully sorted the reads?

At first glance the reads have been sorted succesfully according to the chromosome and the coordinate, which we can see if we look at the first couple of lines (3 reads shown)
~~~bash
samtools view output_aln.sorted.bam | head -3
SRR002996.10376720	16	chr21	9719768	25	36M	*	0	0	AATTCTGAGAAACTTCTTTGTGAGGGTTGGATTTTT	@@@@?@@@?@?@??@???>><=8>6849',%0.%5+	XT:A:U	NM:i:2	X0:i:1	X1:i:0	XM:i:2	XO:i:0	XG:i:0	MD:Z:33C0A1
SRR003003.12879416	0	chr21	9719769	37	36M	*	0	0	ATTCTGAGAAACTTCTTTGTGAGGGTTGGATTCATC	@99;<7:9=@3DB?ECEE@CBF@>;?@B@@?@A?@@	XT:A:U	NM:i:1	X0:i:1	X1:i:0	XM:i:1	XO:i:0	XG:i:0	MD:Z:35T0
SRR002994.3881533	0	chr21	9719772	37	36M	*	0	0	CTGAGAAACTTCTTTGTGAGGGTTGGATTCATTTCA	'&,/*2636B5;5@<?2E@?@=A?@@@@??>@A?@@	XT:A:U	NM:i:0	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:36
~~~
However if we look at the last couple of lines
~~~bash
samtools view output_aln.sorted.bam | tail -3
SRR002999.11626831	4	*	0	0	*	*	0	0	TAATGATGTAAAATGAGATTAATGATGTGTCATTTT	-?*?@8?@+@0;=B<;>)5>5<B6=3+808-2'-)3
SRR002999.11627025	4	*	0	0	*	*	0	0	CGATTGTAATGGAATGGAATAGAATGGAATGGATTG	,?C<@30?0?,512:3*(/5,-*..++.*2-')&,,
SRR002999.11627329	4	*	0	0	*	*	0	0	TAGCATAAAAAAACTAAAATTACTCTCATAGTGGTA	?<?7@:?4:;:5A(..9<'91(6:.+'+-*3))-&,
~~~
We can see these reads have no information regarding the chromosome and the coordinate (column 3 and 4). This is because the reads are unaligned. Which we need to filter out.

## Filtering
To perform further analysis on raw bam files, we often need to filter the reads.
Filtering can be done on mapping quality (column 5) or filtering based on the type of alignment which can be determined based on the 2nd column containing a FLAG. 
The flags can be further explained here (https://broadinstitute.github.io/picard/explain-flags.html)

~~~bash
samtools view -F 4 -Sb output_aln.bam > output_aln.filter.bam
~~~
Ideally we should perform filtering before sorting


NOTE: You can of course also perform similar filtering, sorting and file conversion on the sam file created by the bwa mem command

## QUESTION
1. Which reads are filtered out using the parameters "-F 4"

-F means which read we should not include, flag 4 (look at the explain flag link) means unmapped reads.

2. Which parameters do we need to include to filter reads with a mapping quality above 30

Type in "samtools view" to see all the possible parameters

 -q INT   only include reads with mapping quality >= INT [0]

## Statistics
We can get simple stats of the bam files using "flagstat"

~~~bash
samtools flagstat output_aln.bam
~~~

### QUESTION
1. Which of the alignment methods bwa aln and bwa mem generates the most aligned reads?
Just by looking at the raw sam files
~~~bash
samtools flagstat output_aln.sam 
5277758 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
4271033 + 0 mapped (80.93% : N/A)

samtools flagstat output_mem.sam 
5277758 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
3754281 + 0 mapped (71.13% : N/A)
~~~

So we see that the output from bwa aln (output_aln) has a higher number of aligned reads (4271033) compared to the output from bwa mem (3754281). This is because bwa aln is better for shorter sequence reads. The average read length can be found by using the command from day1
~~~bash
zcat NA19238.fastq.gz | awk 'NR%4==2{sum+=length($0)}END{print sum/(NR/4)}'
36
~~~


To generate more various statistics about the mapping we can use
 
~~~bash
samtools stats output_aln.bam > output_aln.stats
~~~
These statistics can be visualized using plot-bamstats


