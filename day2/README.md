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
We will try multiple different alignment methods all based on the Burrows-Wheeler tranform. Namely bwa and Bowtie2

If you run the bwa alignment commands with no arguments to get info about how to use it, then the number
of options may be a bit overwhelming, but you can run it with no additional options.

To run with threads we need to consider the CPU load so I suggest you add "-t 2" or "-t 4" to have it run in 2 or 4 threads if your computer
has multiple cores. 

All of the commands reads the compressed fastq files directly, so you need not
decompress them. By default the result comes on stdout (in the terminal), so you have to redirect to a file using ">". 

## BWA ALN
First we will find the coordinates of our input reads files in the Suffix Array (SA) created in the index files
~~~bash
bwa aln -t 2 reference_file readfile_1 > readfile1_SA.sai
~~~

Then we will generate alignments in the SAM format given the single-end read by searching in the SA intervals
~~~bash
bwa samse reference_file readfile1_SA.sai readfile_1 > output_aln.sam
~~~

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

## Sorting
Often we need a file sorted on genomic coordinates. You can use "samtools sort" for this
~~~bash
samtools sort -@ 2 output_aln.bam -o output_aln.sorted.bam
~~~
Try to use the "samtools view" to look at the bam files, with the 3rd and 4th column being the chromosomal position, so you can check if it is
sorted

## QUESTION
1. Based on column 3 and 4, have we succesfully sorted the reads?

~~~bash
SRR002996.10376720	16	chr21	9719768	25	36M	*	0	0	AATTCTGAGAAACTTCTTTGTGAGGGTTGGATTTTT	@@@@?@@@?@?@??@???>><=8>6849',%0.%5+	XT:A:U	NM:i:2	X0:i:1	X1:i:0	XM:i:2	XO:i:0	XG:i:0	MD:Z:33C0A1
SRR003003.12879416	0	chr21	9719769	37	36M	*	0	0	ATTCTGAGAAACTTCTTTGTGAGGGTTGGATTCATC	@99;<7:9=@3DB?ECEE@CBF@>;?@B@@?@A?@@	XT:A:U	NM:i:1	X0:i:1	X1:i:0	XM:i:1	XO:i:0	XG:i:0	MD:Z:35T0
SRR002994.3881533	0	chr21	9719772	37	36M	*	0	0	CTGAGAAACTTCTTTGTGAGGGTTGGATTCATTTCA	'&,/*2636B5;5@<?2E@?@=A?@@@@??>@A?@@	XT:A:U	NM:i:0	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:36
~~~
## Filtering
To perform further alignments on raw bam files, we often need to filter the reads.
Filtering can be done on mapping quality or filtering based on the type of alignment which can be determined based on the 2nd column containing a FLAG. 
The flags can be further explained here (https://broadinstitute.github.io/picard/explain-flags.html)

~~~bash
samtools view -F 4 -Sb output_aln.bam > output_aln.filter.bam
~~~
After filtering we can then perform sorting again

## QUESTION
1. Which reads are filtered out using the parameters "-F 4"
2. Which parameters do we need to include to filter reads with a mapping quality above 30

## Statistics
We can get simple stats of the bam files using "flagstat"

~~~bash
samtools flagstat output_aln.bam
~~~

### QUESTION
1. Which of the alignment methods bwa aln and bwa mem generates the most aligned reads?

To generate more various statistics about the mapping we can use
 
~~~bash
samtools stats output_aln.bam > output_aln.stats
~~~
These statistics can be visualized using plot-bamstats


