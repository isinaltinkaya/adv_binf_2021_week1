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

## BWA ALIGNMENT using suffix arrays
First we will find the coordinates of our input reads files in the Suffix Array (SA) created in the index files

The help page for bwa can be found at http://bio-bwa.sourceforge.net/bwa.shtml and further information can be found by simply typing bwa in your terminal

Q1) Which command do we need in order to generate the Suffix Array given the input reads?

Q2) Using the generated SA intervals, perform single end alignment using the reference genome and the provided sequencing reads

## BWA MEM
Q3) Use the bwa mem command to generate another alignment file

# SAM/BAM output
Now take a look at the sam file from the previous the output from Q2 and Q3.

## File conversion
We often need the compressed version of the sam format, which is called bam. To do this we use "samtools view" for converting between formats. Try
the command without arguments to see the options.

Q4) Convert the sam output files to both a bam format and cram format.

Q5) What is the major difference between the conversion from sam to bam and sam to cram 

## Sorting
Often we need a file sorted on genomic coordinates.

Q6) Which command do we need for sorting the aligned files and which columns can we use to identify that the sorting was succesful? 

Q7) Looking at the sorted reads according to their chromosome is there any problematic readÂs?

## Filtering
To perform further analysis on raw bam files, we often need to filter the reads.
Filtering can be done on mapping quality (column 5) or filtering based on the type of alignment which can be determined based on the 2nd column containing a FLAG. 
The flags can be further explained here (https://broadinstitute.github.io/picard/explain-flags.html)

Q8) Filter out the reads with a mapping quality below 35 and keep all of the aligned reads

Q9) What command would you use to filter out (based on the alignment flag) all of the reads aligning to the reverse strand

Q10) What command would you use to keep all of the reads aligning to the reverse strand.

NOTE: You can of course also perform similar filtering, sorting and file conversion on the sam file created by the bwa mem command

## Statistics

Q11) What command do we need to extract simple statistics of both the aligned files obtained using the two difference mapping approaches 

Q12) Which of the two alignment commands generates the most aligned reads?

Q13) What is the average read length of the read sequence file, and how can this be used to explain the results obtained in Q12? 
The average read length can be found by using the command from day1
