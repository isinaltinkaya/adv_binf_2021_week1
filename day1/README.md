# Day 1: NGS data - Workflow, formats and programs

**0. Getting started**

**1. Working with FASTA files**

**2. FASTQ files**

**3. Adapter trimming using fastp**

**4. Sequence Alignment/Map Formats: SAM, BAM and CRAM**

**5. Testing for Damage Patterns with mapDamage**

**6. Variant Call Formats: VCF and BCF**

## Getting started


### Connecting to the server via SSH



X11 forwarding method will allow you to start a graphical application on the remote system and forward this application's windows to your local system. We need to enable X11 forwarding to view the plots we will be generating for the exercises.

We use `-X` option to enable X11 forwarding over SSH:

```sh
$ ssh -X <your_username>@<server_name_or_ip>
```

Replace with your remote server username. For example:


```sh
$ ssh -X isin@ricco.popgen.dk
```




### Setting up the working environment


CHANGEME
We have 3 main directories.

```sh
├── data
│   ├── alignment
│   ├── fasta
│   ├── fastq
│   └── reference_fasta
└── exercises
    ├── alignment_formats
    ├── mapdamage
    └── trimming

```

```bash
/TEACHING/BIOINF21/
├── data
├── programs
└── github
```


```sh
$ git clone https://github.com/isinaltinkaya/adv_binf_2021_week1
$ cd adv_binf_2021_week1
```
___
___
___


## 1. Working with FASTA files

Working directory: `day1/data/fasta`


**File extensions:** `.fasta` or `.fa` for generic FASTA (See [this link](https://en.wikipedia.org/wiki/FASTA_formathttps://en.wikipedia.org/wiki/FASTA_format) for details on other FASTA extensions `.fna` `.ffn` `.faa` and `.frn`)


For each record, the FASTA file consists of two main fields:

Line | Contains | Description | Example
--- | --- | --- | ---
1 | Sequence ID | Starts with a `>` sign followed by a sequence identifier. | `>Chr21`
2+ | Sequence | One or more lines containing the sequence itself. | `CACCTCCCCTCAGGCCGCATTGCAGTGGGGGCTGAGAGGAGGAAGCACCATGGCCCACCTCTTCTCACCCCT`


```sh
$ cat example.fasta 
>sequence1
CACCTCCCCTCAGGCCGCATTGCAGTGGGGGCTGAGAGGAGGAAGCACCATGGCCCACCTCTTCTCACCCCT
>sequence2
CGCGCTGTCCGCGCTGAGCCACCTGCACGCGTGCCAGCTGCGAGTGGACCCGGCCAGCTTCCAGGTGAGCGGCTG
>sequence3
CCGTGCTGGGCCCCTGTCCCCGGGAGGGCCCCGGCGGGGTGGGTGCGGGGGGCGTGCGG
>sequence4
TGAGCCTTGAGCGCTCGCCGCAGCTCCTGGGCCACTGCCTGCTGGTAACCCTCGCCCGGCACTACCCCGGAGACT
CCGTGCTGGGCCCCTGTCCCCGGGAGGGCCCCGGCGGGGTGGGTGCGGGGGGCGTGCGGGGCGGGTGCAGGCGAG
ACAGA
```

**1.1. First, we need to index the FASTA file**

```sh
$ samtools faidx example.fasta
$ cat example.fasta.fai 
sequence1	72	11	72	73
sequence2	75	95	75	76
sequence3	59	182	59	60
sequence4	155	253	75	76
```

**1.2. Count the number of sequences in a FASTA file**

```sh
$ grep -c ">" example.fasta
4
```

**1.3. Fetching regions from FASTA files**

We use `samtools faidx` to fetch genomic sequences:

```sh
samtools faidx genome.fa <chromosome|sequence id>:<start>-<end>
```


- Coordinates are 1-based in FASTA files.
- Ranges in `samtools faidx` is inclusive, so that `sequence1:2-10` corresponds to `[2-10]` from `sequence1`.


```sh
$ samtools faidx example.fasta sequence1:2-10
>sequence1:2-10
ACCTCCCCT
```


**1.4. Find the kmer "CAGGTGAGC" in the FASTA file**

```sh
$ grep "CAGGTGAGC" example.fasta
CGCGCTGTCCGCGCTGAGCCACCTGCACGCGTGCCAGCTGCGAGTGGACCCGGCCAGCTTCCAGGTGAGCGGCTG
```


**1.5. Get sequence list from FASTA file and save it to a file called "sequences.txt"**
```sh
$ grep "^>" example.fasta| tr -d '>' >> sequences.txt
$ cat sequences.txt
sequence1
sequence2
sequence3
sequence4
```

**1.6. Get sequence names and their lengths**
```sh
$ cat example.fasta.fai  | cut -f1,2
sequence1	72
sequence2	75
sequence3	59
sequence4	155
```

**1.7. Download chromosome 21 of h19 reference genome**
```sh
wget 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr21.fa.gz' -O chr21.fa.gz
```

___
___
___


## 2. FASTQ files


Working directory: `day1/data/fastq`


**File extensions:** `.fastq` or `.fq`

Line | Contains | Description | Example
--- | --- | --- | ---
1 | Read ID | Begins with a `@` character and is followed by a sequence identifier and description (optional) |  `@NS500474:51:HK7FVAFXX:1:11101:5669:1056 1:N:0:GACACT+NTGACG`
2 | Base calls | Raw sequence letters | `CAGCTGGCGTCGGCCGACGTGATCACCTTCACGATCGGAAGAACACACGTCTGAACTCCAGTCACGACACTATCT`
3 | Seperator or additional information (optional) | Begins with a `+` character and is optionally followed by the sequence identifier and descriptions (if any) | `+`
4 | Sequencing base quality (1 per base call) |  encodes the quality values for the sequence in Line 2, and must contain the same number of symbols as letters in the sequence. Each letter in line 4 corresponds to the base quality score of the sequence with the same index in line 2. | `AAAAAEEEEEEEEEEEEEEEE6EEEEEEEEEEEEEEEEAEEEEE/EEEE/EEEEEAEEEEEEEEEEEEEEEE/<A`




**2.1. Get the first read from lane 1 read 1**
```sh
$ zcat DATA_L001_R1.fastq.gz | head -4
@NS500474:51:HK7FVAFXX:1:11101:5669:1056 1:N:0:GACACT+NTGACG
CAGCTGGCGTCGGCCGACGTGATCACCTTCACGATCGGAAGAACACACGTCTGAACTCCAGTCACGACACTATCT
+
AAAAAEEEEEEEEEEEEEEEE6EEEEEEEEEEEEEEEEAEEEEE/EEEE/EEEEEAEEEEEEEEEEEEEEEE/<A
```

**2.2. Get the first read in R1 and R2 from lane 1 and 2**

```sh
$ for R in DATA_L00{1,2}_R*.fastq.gz; do echo "FILE: $R"; zcat $R|head -4;done
FILE: DATA_L001_R1.fastq.gz
@NS500474:51:HK7FVAFXX:1:11101:5669:1056 1:N:0:GACACT+NTGACG
CAGCTGGCGTCGGCCGACGTGATCACCTTCACGATCGGAAGAACACACGTCTGAACTCCAGTCACGACACTATCT
+
AAAAAEEEEEEEEEEEEEEEE6EEEEEEEEEEEEEEEEAEEEEE/EEEE/EEEEEAEEEEEEEEEEEEEEEE/<A
FILE: DATA_L001_R2.fastq.gz
@NS500474:51:HK7FVAFXX:1:11101:5669:1056 2:N:0:GACACT+NTGACG
GGGGGGGGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGGGGGGGGGGGNNN
+
AAAAAEEE######################################################E<E/EEAE<AE###
FILE: DATA_L002_R1.fastq.gz
@NS500474:51:HK7FVAFXX:2:11101:7383:1057 1:N:0:GACACT+NTGACG
TTGTTCGACATTCTTGTCCGAGCCCAGAAGAAGCGATCGGAAGAGCACACGTCTGAACCAAGTCACGACACTAT
+
/AAAAEEAAEEAE/AAAAAEEEAEAEEEEEEEE<EAEEEEEEAEEEEEEEE/A6AEEEA/AA<EAA/EAEEE/E
FILE: DATA_L002_R2.fastq.gz
@NS500474:51:HK7FVAFXX:2:11101:7383:1057 2:N:0:GACACT+NTGACG
GGGGGGGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGGGGGGGGGGNNNNN
+
AAAAAEE######################################################///////<//#####
```

**2.3. Return only the first 5 reads from a FASTQ file**
```sh
$ zcat DATA_L001_R1.fastq.gz | sed -n '2~4p' | head -5
CAGCTGGCGTCGGCCGACGTGATCACCTTCACGATCGGAAGAACACACGTCTGAACTCCAGTCACGACACTATCT
ATGCGGGGCTGAGGCTGCTGGAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGACACTATCTCGTATGCCGT
GCCCGGCACGATCACACCGGGGCTAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGACTATCTCGTATGC
GCCGCGGAACTCGCAAAGGCGTTCGGCGTGGACGTCCCGACTGCCAAACGCAGAACGGAAGAGCACACGTCAGAA
GGCTGTTGCTGCCTCGGGAGCATCAATCTCGCGAGATCGGAAGAGCACACGTCTGAACTCGAGTCACGACACTAT
```
or

```sh
$ zcat DATA_L001_R1.fastq.gz | awk 'NR%4==2{print $0}'|head -5
CAGCTGGCGTCGGCCGACGTGATCACCTTCACGATCGGAAGAACACACGTCTGAACTCCAGTCACGACACTATCT
ATGCGGGGCTGAGGCTGCTGGAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGACACTATCTCGTATGCCGT
GCCCGGCACGATCACACCGGGGCTAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGACTATCTCGTATGC
GCCGCGGAACTCGCAAAGGCGTTCGGCGTGGACGTCCCGACTGCCAAACGCAGAACGGAAGAGCACACGTCAGAA
GGCTGTTGCTGCCTCGGGAGCATCAATCTCGCGAGATCGGAAGAGCACACGTCTGAACTCGAGTCACGACACTAT
```
**2.4. Calculate the mean length**

```sh
$ zcat DATA_L001_R1.fastq.gz | awk 'NR%4==2{sum+=length($0)}END{print sum/(NR/4)}'
75.37
```

**2.5. Which lane contains the highest and lowest number of raw sequence reads?**

```sh
$ zcat DATA_L001_R1.fastq.gz | awk 'NR%4==2{print $0}'|wc -l
1416689
```

or

```sh
$ zcat DATA_L001_R1.fastq.gz |  echo $((`wc -l` / 4))
1416689
```

or

```sh
$ zcat DATA_L001_R1.fastq.gz | grep -c "^@" 
1416689
```

**2.6. How many times do we see the motif "GATTACA"?**
```sh
$ zcat DATA_L001_R1.fastq.gz | grep -c "GATTACA" 
1147
```



**2.7. Get fragment size distribution and related statistics from a FASTQ file**

```sh
$ zcat DATA_L001_R1.fastq.gz | awk 'NR%4==2{print length($0)}' | datamash mean 1
75.369972520433
```



___
___
___


## 3. Adapter trimming using [fastp](https://github.com/OpenGene/fastp)

Working directory: `day1/exercises/trimming`

We will use `fastp` for adapter trimming. Some examples of other commonly used FASTQ preprocessing and adapter removal tools are [AdapterRemoval](https://adapterremoval.readthedocs.io/en/latest/) and [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic).

```sh
fastp --in1 ../../data/fastq/DATA_L001_R1.fastq.gz \ 
--in2 ../../data/fastq/DATA_L001_R2.fastq.gz \ 
--out1 DATA_L001_R1_trimmed.fastq.gz \ 
--out2 DATA_L001_R2_trimmed.fastq.gz \ 
--merge \ 
--merged_out DATA_L001_merged_trimmed.fastq.gz \ 
--unpaired1 DATA_L001_unpaired_R1.fastq.gz \
--unpaired2 DATA_L001_unpaired_R2.fastq.gz \ 
--length_required 30 \ 
--detect_adapter_for_pe
```
- `--out1` and `--out2` outputs the reads that cannot be merged successfully, but both pass all the filters.
- `--merged_out` specifies the file name to write the merged reads.
- `--unpaired1` outputs be the reads that cannot be merged, read1 passes filters but read2 doesn't.
- `--unpaired2` outputs be the reads that cannot be merged, read2 passes filters but read1 doesn't.
- By default, fastp uses overlap analysis to do adapter trimming. However, you can specify `--detect_adapter_for_pe` to enable adapter trimming by adapter sequence detection.


```sh
> fastp --in1 ../../data/fastq/DATA_L001_R1.fastq.gz --in2 ../../data/fastq/DATA_L001_R2.fastq.gz --out1 DATA_L001_R1_trimmed.fastq.gz --out2 DATA_L001_R2_trimmed.fastq.gz --merge --merged_out DATA_L001_merged_trimmed.fastq.gz --unpaired1 DATA_L001_unpaired_R1.fastq.gz --unpaired2 DATA_L001_unpaired_R2.fastq.gz --length_required 30 --detect_adapter_for_pe


Detecting adapter sequence for read1...
>Nextera_LMP_Read1_External_Adapter | >Illumina Multiplexing Index Sequencing Primer
GATCGGAAGAGCACACGTCTGAACTCCAGTCAC

Detecting adapter sequence for read2...
>Illumina TruSeq Adapter Read 2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

Read1 before filtering:
total reads: 1416689
total bases: 106775811
Q20 bases: 102917164(96.3862%)
Q30 bases: 100732807(94.3405%)

Read2 before filtering:
total reads: 1416689
total bases: 107129620
Q20 bases: 95391775(89.0433%)
Q30 bases: 91822908(85.712%)

Merged and filtered:
total reads: 505004
total bases: 23375425
Q20 bases: 22800929(97.5423%)
Q30 bases: 22534869(96.4041%)

Filtering result:
reads passed filter: 1092634
reads failed due to low quality: 206032
reads failed due to too many N: 776
reads failed due to too short: 1533936
reads with adapter trimmed: 2428136
bases trimmed due to adapters: 108248691
reads corrected by overlap analysis: 33214
bases corrected by overlap analysis: 54428

Duplication rate: 0.629749%

Insert size peak (evaluated by paired-end reads): 0

Read pairs merged: 505004
% of original read pairs: 35.6468%
% in reads after filtering: 100%


JSON report: fastp.json
HTML report: fastp.html
```

Following files will be generated:

```sh
├── DATA_L001_merged_trimmed.fastq.gz
├── DATA_L001_R1_trimmed.fastq.gz
├── DATA_L001_R2_trimmed.fastq.gz
├── DATA_L001_unpaired_R1.fastq.gz
├── DATA_L001_unpaired_R2.fastq.gz
├── fastp.html
└── fastp.json
```

You can view the output html report by using
```sh
firefox --no-remote fastp.html
```

Alternatively, you can use `scp` to download the output html file to your local machine
```sh
scp isin@ricco.popgen.dk:/PATH_TO/fastp.html .
```
**3.1. What is the proportion of reads with residual adapters?**
```sh
2428136/(1416689*2)
[1] 0.8569757
```

**3.2. Compare the average fragment length of merged reads, read 1 and read 2 after trimming**2
```sh
$ for FILE in *trimmed.fastq.gz;do echo ${FILE}; zcat ${FILE}| awk 'NR%4==2{sum+=length($0)}END{print sum/(NR/4)}';done
DATA_L001_merged_trimmed.fastq.gz
46.2876
DATA_L001_R1_trimmed.fastq.gz
68.5516
DATA_L001_R2_trimmed.fastq.gz
72.4403
```
___
___
___


## 4. Sequence Alignment/Map Formats: SAM, BAM and CRAM

Working directory: `day1/exercises/alignment_formats`


The SAM format consists of two main sections: the header and the alignment itself. 


**1. Header:** The header consists of the lines beginning with a `@` symbol, and includes details regarding the reference sequence used in alignment.
```sh
$ samtools view -H DATA.bam
@HD	VN:1.3	SO:coordinate
@SQ	SN:chr21	LN:48129895
@PG	ID:bwa	PN:bwa	VN:0.7.15-r1142-dirty	CL:bwa mem ../reference_fasta/chr21.fa.gz ../../exercises/trimming/DATA_L001_merged_trimmed.fastq.gz
@PG	ID:bwa-7E1842B3	PN:bwa	VN:0.7.15-r1142-dirty	CL:bwa mem ../reference_fasta/chr21.fa.gz ../../exercises/trimming/DATA_L001_R1_trimmed.fastq.gz ../../exercises/trimming/DATA_L001_R2_trimmed.fastq.gz
```

The line beginning with `@HD` can contain header information including the version (`VN:`) and the sorting order (`SO:`) of the alignments.

There can be one or more lines that start `@SQ` to describe the reference sequences in the file. In our case, as our sample was aligned to only the chromosome 21, we see only one line. Each `@SQ` line can have multiple entries:
- `SN:` the reference sequence name
- `LN:` the reference sequence length
- `AN:` alternate names for this sequence
- `AS:` how the reference was assembled
- `M5:` the [MD5sum](https://en.wikipedia.org/wiki/Md5sum) of the sequence
- `SP:` the species, and
- `UR:` the URL for the sequence.

**2. Alignment:** Each line contains 11 fields (all required) and the fields are separated by tab symbols. The fields are:

Column index | Field Name | Type | Description
--- | --- | --- | --- 
0 | QNAME | String | Query sequence name
1 | FLAG | Integer | Bitwise flag
2 | RNAME | String | Reference sequence name
3 | POS | Integer | Left most mapping position (1-based)
4 | MAPQ | Integer | Mapping quality score
5 | CIGAR | String | The CIGAR string for the alignment
6 | RNEXT | String | The reference name of the mate (or next read)
7 | PNEXT | Integer | The position of the mate (or next read)
8 | TLEN | Integer | The observed template length
9 | SEQ | String | The sequence of the segment
10 | QUAL | String | The ASCII representation of the base quality score


**Bitwise flags:** The FLAGs in the second column (column 1) are comprised of a bitwise combination of numbers. Details can be found at [Decoding SAM Flags webpage](https://broadinstitute.github.io/picard/explain-flags.html).


**4.1. How many reads do we have for each different FLAG in our BAM file?**
```sh
$ samtools view DATA.bam |cut -f2|sort -n|uniq -c
   1852 0
 501278 4
   1874 16
     13 69
     12 73
  41104 77
      1 81
     80 83
      2 97
     73 99
      4 113
     13 117
     11 121
     12 133
     13 137
  41104 141
      2 145
     73 147
      1 161
     80 163
      4 177
     11 181
     13 185
      1 2048
      1 2161
      1 2233
```

**4.2. What properties does a read with FLAG value of 82 have?**

If we have 82 as the FLAG, since `2+16+64=82` we can understand that that the sequence is properly aligned, is reverse complemented, and is the first in pair.



**CIGAR string:** The CIGAR string is a representation of the sequence alignment in an abbreviated form. The letters correspond to:

Code | Description
--- | ---
M | Alignment match (but could be a sequence match or mismatch!)
I | Insertion relative to the reference
D | Deletion from the reference (i.e. insertion relative to the query!)
N | Reference skipped
S | soft clipping
H | Hard clipping
P | Padding
= | Sequence match
X | Sequence mismatch

**4.3. What properties does a read with CIGAR string of 3M1D2M have?**

Thus the CIGAR string 2M1D3M means there are three matches, 1 insertion in the query (1 deletion in the reference), and two more matches.



**4.4. Sort SAM files by coordinates**
```sh
$ samtools sort DATA_L001_R1_R2.sam >> DATA_L001_R1_R2_sorted.bam
$ samtools sort DATA_L001_merged.sam >> DATA_L001_merged_sorted.bam
```

**4.5. Merge BAM files into one file called `DATA.bam`**
```sh
$ samtools merge DATA.bam DATA_L001_merged_sorted.bam DATA_L001_R1_R2_sorted.bam
```

**4.6 Index the bam file**
```sh
$ samtools index DATA.bam
```


**4.7. Get mapped reads and write to a BAM file named `DATA_mapped.bam`, and a CRAM file named `DATA_mapped.cram`. What are the file sizes of the BAM and CRAM files? Are they different, if so, why?** 

We need to use `-T` to provide the reference sequence used in alignment with `-C` to convert a BAM file into a CRAM file. 

`-F` corresponds to "exclude reads with flag \<INT\>"
```sh
$ samtools view -F 4 -h -b DATA.bam > DATA_mapped.bam
$ samtools view -F 4 -h -T ../../data/reference_fasta/chr21.fa.gz -C DATA.bam > DATA_mapped.cram
$ du -sh DATA_mapped.{bam,cram}
912K	DATA_mapped.bam
100K	DATA_mapped.cram
```

**4.8. Count the number of unmapped reads** 

`-f` corresponds to "only include reads with flag \<INT\>"
```sh
$ samtools view -f 4 -c DATA.bam
```


**4.9. Get the HLCS gene (chr21:38,120,926-38,362,511) and write to a CRAM file**

```sh
$ samtools view -h -T ../../data/reference_fasta/chr21.fa.gz -C DATA.bam chr21:38120926-38362511  > DATA_HLCS.cram
```



**4.10. How many reads overlap with the HLCS gene?**
```sh
$ samtools view -c DATA_HLCS.cram
40
```

**4.11. Inspect the alignment in region `chr21:14338386-14338400` visually using `samtools tview`**

You can save the output as an html file and view it using `firefox --no-remote`
```sh
$ samtools tview DATA.bam -p chr21:14338386-14338400 -d html >> DATA_chr21_14338386-14338400_tview.html
$ firefox --no-remote DATA_chr21_14338386-14338400_tview.html 
```
or you can directly view it in command-line:
```sh
$ samtools tview DATA.bam -p chr21:14338386-14338400
```


**4.12. What is the mean depth of the file `DATA.bam` given a minimum mapping quality of 30 and a minimum base quality of 20?**

There are multiple ways to calculate the mean depth (average number of bases mapped to positions on the reference).
`samtools depth` returns three columns: reference, position, and coverage at position.
```sh
$ samtools depth -Q 30 -q 20 DATA.bam | head
chr21	9416492	1
chr21	9416493	1
chr21	9416494	1
chr21	9416495	1
chr21	9416496	1
chr21	9416497	1
chr21	9416498	1
chr21	9416499	1
chr21	9416500	1
chr21	9416501	1
```

We can use `datamash` to calculate the mean depth

```sh
$ samtools depth -Q 30 -q 20 DATA.bam | cut -f3 | datamash mean 1
1.0677556818182
```

- Can this number be lower than 1? Why not? (Hint: See the option `-a` in the manual)



**4.13. Get FLAG statistics using `samtools flagstat`**
```sh
$ samtools flagstat DATA.cram
587633 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
3 + 0 supplementary
0 + 0 duplicates
4098 + 0 mapped (0.70% : N/A)
82626 + 0 paired in sequencing
41313 + 0 read1
41313 + 0 read2
306 + 0 properly paired (0.37% : N/A)
320 + 0 with itself and mate mapped
49 + 0 singletons (0.06% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```

**4.14. Get alignment statistics using `samtools stats` and write to a file called `DATA.stats`, then plot it using `plot-bamstats`**
```sh
$ samtools stats DATA.bam >> DATA_stats.txt
```
You can use `grep` to fetch statistics of interest:
```sh
$ cat DATA_stats.txt | grep "^SN"
SN	raw total sequences:	587630
SN	filtered sequences:	0
SN	sequences:	587630
SN	is sorted:	1
SN	1st fragments:	546317
SN	last fragments:	41313
SN	reads mapped:	4095
SN	reads mapped and paired:	320	# paired-end technology bit set + both mates mapped
SN	reads unmapped:	583535
SN	reads properly paired:	306	# proper-pair bit set
SN	reads paired:	82626	# paired-end technology bit set
SN	reads duplicated:	0	# PCR or optical duplicate bit set
SN	reads MQ0:	2669	# mapped and MQ=0
SN	reads QC failed:	0
SN	non-primary alignments:	0
SN	total length:	29200223	# ignores clipping
SN	bases mapped:	174032	# ignores clipping
SN	bases mapped (cigar):	150514	# more accurate
SN	bases trimmed:	0
SN	bases duplicated:	0
SN	mismatches:	3315	# from NM fields
SN	error rate:	2.202453e-02	# mismatches / bases mapped (cigar)
SN	average length:	49
SN	maximum length:	121
SN	average quality:	34.4
SN	insert size average:	609.0
SN	insert size standard deviation:	2072.8
SN	inward oriented pairs:	33
SN	outward oriented pairs:	59
SN	pairs with other orientation:	4
SN	pairs on different chromosomes:	0

# for coverage distribution
$ cat DATA_stats.txt | grep "^COV"
COV	[1-1]	1	141896
COV	[2-2]	2	12155
COV	[3-3]	3	857
COV	[4-4]	4	443
COV	[5-5]	5	205
COV	[6-6]	6	93
COV	[7-7]	7	55
COV	[8-8]	8	15
COV	[9-9]	9	22
COV	[10-10]	10	10
COV	[14-14]	14	20
COV	[15-15]	15	6
COV	[16-16]	16	2
COV	[17-17]	17	21
COV	[18-18]	18	24
```
```sh
$ plot-bamstats DATA_stats.txt -p DATA_stats
```

**4.15. Check if the bam file is corrupted or not**


```sh
$ samtools quickcheck DATA.bam && echo "Data looks OK" || echo "Data is corrupted"
Data looks OK

$ samtools quickcheck -vvv DATA.bam 
verbosity set to 3
checking DATA.bam
opened DATA.bam
DATA.bam is sequence data
DATA.bam has 1 targets in header.
DATA.bam has good EOF block.

$ samtools quickcheck DATA_bad.bam && echo "Data looks OK" || echo "Data is corrupted"
Data is corrupted

$ samtools quickcheck -vvv DATA_bad.bam 
verbosity set to 3
checking DATA_bad.bam
opened DATA_bad.bam
DATA_bad.bam is sequence data
DATA_bad.bam has 1 targets in header.
DATA_bad.bam was missing EOF block when one should be present.
DATA_bad.bam
```
**4.16. Summarizing the base calls: pileup format**

Pileup format is a text-based format and is useful for summarizing the base calls of aligned reads to a reference sequence. This format enables the visual investigation of SNP/indel calling and alignment.



Using `samtools mpileup` we can obtain a summary of the coverage of mapped reads on a reference sequence at a single base pair resolution.


```sh
$ samtools mpileup DATA.bam >> DATA.mpileup
```

Each line in a pileup file represents a single genomic position and consists of following columns:

Column index | Contains | Details | Example
--- | --- | --- | ---
0 | Sequence ID | Name of the sequence | chr21
1 | Coordinate | 1-based coordinate of each site | 9416492
2 | Reference base  | Corresponding base in the reference genome (generated by `-f` option)| g 
3 | Read count |  Number of reads covering the position | 1
4 | Read bases | Read bases column contains information on whether if a read base was matched, mismatched, was inserted or deleted with respect to the reference. It also contains information about whether the read base was on the positive or negative strand with respect to the reference, whether a read base was at the start or at the end of a read. `.` corresponds to match to the reference base on the positive strand, `,` corresponds to match on the negative strand.`^]`  indicates that the base was at the beginning of a read, and the `$` indicates that the base was at the end of a read.| ^],
5 | Base qualities | Base quality score at position | E
6 | Mapping quality | Alignment mapping quality (generated by `-s` option) | ]



- We can also retain the actual nucleotide at each site by using `-f` option to specify the reference file used in alignment.
- `-O (or --output-BP)` outputs the base positions on reads.
- `-s (or --output-MQ)` outputs the mapping quality, `-a` outputs all positions (including zero depth), and  
- `-a -a (or -aa)` outputs all positions including unused reference sequences.
```sh
$ samtools mpileup -f ../../data/reference_fasta/chr21.fa.gz -O -s DATA.bam >> DATA_basePos_MQ_withRef.mpileup
$ samtools mpileup -f ../../data/reference_fasta/chr21.fa.gz -O -s -a DATA.bam >> DATA_all_basePos_MQ_withRef.mpileup
```
```sh
$ head *mpileup
==> DATA_all_basePos_MQ_withRef.mpileup <==
chr21	1	N	0	*	*	*	*
chr21	2	N	0	*	*	*	*
chr21	3	N	0	*	*	*	*
chr21	4	N	0	*	*	*	*
chr21	5	N	0	*	*	*	*
chr21	6	N	0	*	*	*	*
chr21	7	N	0	*	*	*	*
chr21	8	N	0	*	*	*	*
chr21	9	N	0	*	*	*	*
chr21	10	N	0	*	*	*	*

==> DATA_basePos_MQ_withRef.mpileup <==
chr21	9416492	g	1	^],	E	]	1
chr21	9416493	t	1	,	E	]	2
chr21	9416494	t	1	,	E	]	3
chr21	9416495	a	1	,	E	]	4
chr21	9416496	t	1	,	E	]	5
chr21	9416497	a	1	,	E	]	6
chr21	9416498	t	1	,	E	]	7
chr21	9416499	a	1	,	E	]	8
chr21	9416500	t	1	,	E	]	9
chr21	9416501	t	1	,	E	]	10

==> DATA.mpileup <==
chr21	9416492	N	1	^]g	E
chr21	9416493	N	1	t	E
chr21	9416494	N	1	t	E
chr21	9416495	N	1	a	E
chr21	9416496	N	1	t	E
chr21	9416497	N	1	a	E
chr21	9416498	N	1	t	E
chr21	9416499	N	1	a	E
chr21	9416500	N	1	t	E
chr21	9416501	N	1	t	E
```

Let's see the distribution of depth
```sh
$ cat DATA.mpileup | cut -f4 | sort -n | uniq -c
      3 0
 129418 1
   4718 2
    539 3
    296 4
    208 5
     94 6
     48 7
     18 8
     14 9
     10 10
$ cat DATA.mpileup | cut -f4 | sort -n | uniq -c  > mpileup_depth.txt
$ R
> jpeg("depth_dist.jpeg")
> barplot(read.table("mpileup_depth.txt")[1:10,1],names=1:10,xlab="Depth",ylab="Number of sites")
> dev.off()
```

You can view the plot using `feh`
```sh
$ feh depth_dist.jpeg
```



___
___
___


## 5. Testing for Damage Patterns with [mapDamage](https://ginolhac.github.io/mapDamage/)

Working directory: `day1/exercises/mapdamage`

```bash
$ mapDamage -i ../../data/alignment/DATA.bam -r ../../data/reference_fasta/chr21.fa.gz --no-stats
```

Generates the following files:
```sh
DATA.mapDamage/
├── dnacomp.txt
├── Fragmisincorporation_plot.pdf
├── Length_plot.pdf
├── lgdistribution.txt
├── misincorporation.txt
└── Runtime_log.txt
````

Inspect the damage patterns. 
- What would you expect to see if the sample was modern?
- What would you expect to see if the sample was an ancient sample that has a library preparation with [USER treatment](https://international.neb.com/products/m5505-user-enzyme#Product%20Information)?


___
___
___


## 6. Variant Call Formats: VCF and BCF



BCF format is the binary equivalent of the VCF file format.

```sh
$ bcftools mpileup --fasta-ref ../../data/reference_fasta_hs37d5/chr1.fa.gz NA06985.bam -r 1 > NA06985.mpileup

$ bcftools call -c NA06985.mpileup > NA06985.vcf

# Call only for variant sites
$ bcftools call -v -c NA06985.mpileup > NA06985_variant.vcf

$ head NA06985.vcf 
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##bcftoolsVersion=1.3.1-98-ga6a7829+htslib-1.3.1-64-g74bcfd7
##bcftoolsCommand=mpileup --fasta-ref ../../data/reference_fasta_hs37d5/chr1.fa.gz -r 1 NA06985.bam
##reference=file://../../data/reference_fasta_hs37d5/chr1.fa.gz
##contig=<ID=1,length=249250621>
... skipping the contig lines ...
##ALT=<ID=*,Description="Represents allele(s) other than observed.">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
##INFO=<ID=IDV,Number=1,Type=Integer,Description="Maximum number of reads supporting an indel">
##INFO=<ID=IMF,Number=1,Type=Float,Description="Maximum fraction of reads supporting an indel">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)",Version="3">
##INFO=<ID=RPB,Number=1,Type=Float,Description="Mann-Whitney U test of Read Position Bias (bigger is better)">
##INFO=<ID=MQB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality Bias (bigger is better)">
##INFO=<ID=BQB,Number=1,Type=Float,Description="Mann-Whitney U test of Base Quality Bias (bigger is better)">
##INFO=<ID=MQSB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better)">
##INFO=<ID=SGB,Number=1,Type=Float,Description="Segregation based metric.">
##INFO=<ID=MQ0F,Number=1,Type=Float,Description="Fraction of MQ0 reads (smaller is better)">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=AF1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele frequency (assuming HWE)">
##INFO=<ID=AF2,Number=1,Type=Float,Description="Max-likelihood estimate of the first and second group ALT allele frequency (assuming HWE)">
##INFO=<ID=AC1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele count (no HWE assumption)">
##INFO=<ID=MQ,Number=1,Type=Integer,Description="Root-mean-square mapping quality of covering reads">
##INFO=<ID=FQ,Number=1,Type=Float,Description="Phred probability of all samples being the same">
##INFO=<ID=PV4,Number=4,Type=Float,Description="P-values for strand bias, baseQ bias, mapQ bias and tail distance bias">
##INFO=<ID=G3,Number=3,Type=Float,Description="ML estimate of genotype frequencies">
##INFO=<ID=HWE,Number=1,Type=Float,Description="Chi^2 based HWE test P-value based on G3">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
##bcftools_callVersion=1.3.1-98-ga6a7829+htslib-1.3.1-64-g74bcfd7
##bcftools_callCommand=call -c NA06985.mpileup; Date=Sun Sep  5 21:54:10 2021
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA06985
1       13999984        .       T       .       28.2394 .       DP=1;MQ0F=0;AC1=2;DP4=0,0,0,0;MQ=0;FQ=-30.0004  GT:PL   0/0:0
```
```sh
$ bcftools call -m -v -Ob -o NA06985_multiallelic_variant.bcf NA06985.mpileup
```
Let's fill the allele frequency tags

```sh
$ bcftools +fill-tags NA06985.vcf -Ob -o NA06985_with_allelefreqs.bcf
```
```sh
$ bcftools stats NA06985.vcf > NA06985_vcfstats.txt
$ bcftools stats NA06985_variant.vcf > NA06985_variant_vcfstats.txt
$ plot-vcfstats NA06985_variant_vcfstats.txt -p NA06985_variant
$ plot-vcfstats NA06985_vcfstats.txt -p NA06985
```

```sh
$ bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' NA06985_with_allelefreqs.bcf > query.txt
$ head query.txt 
1	13999984	T	.	.
1	13999985	C	.	.
1	13999986	C	.	.
1	13999987	G	.	.
1	13999988	C	.	.
1	13999989	A	.	.
1	13999990	G	.	.
1	13999991	T	.	.
1	13999992	C	.	.
1	13999993	C	.	.
```


___
___
___



### Extras

**bgzip (Blocked GNU Zip Format):** Sequence files are usually zipped by bgzip (e.g. `vcf.gz`). bgzip files can be uncompressed with using gunzip

If you still have some time left, try to answer these questions:
- How many As, Cs, Gs, Ts and Ns are there in each FASTQ file?
- How can we calculate the GC content using a FASTQ file?
- How long is the reference genome? 
- How many reads has their mates were unmapped?
- How many reads are there that was aligned to a region included in 1000G sites with a mapping quality 30 as minimum?
- How can we sort a BAM file by coordinates?


___

Class working directory including the results of analyses above: `/TEACHING/BIOINF21/`

