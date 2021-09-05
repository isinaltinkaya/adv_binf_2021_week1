# Day 1: NGS data - Workflow, formats and programs


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
data
├── fastq
│   ├── DATA_L001_R1.fastq.gz
│   └── DATA_L001_R2.fastq.gz
└── reference_fasta
    ├── chr21.fa.gz
    ├── chr21.fa.gz.amb
    ├── chr21.fa.gz.ann
    ├── chr21.fa.gz.bwt
    ├── chr21.fa.gz.pac
    └── chr21.fa.gz.sa
```

```bash
/TEACHING/BIOINF21/
├── data
├── programs
└── 
```


```sh
$ git clone https://github.com/isinaltinkaya/adv_binf_2021_week1
$ cd adv_binf_2021_week1
```


## 1. Working with FASTA files


A file in FASTA format consists of:
- One line starting with a ">" sign followed by a sequence identifier.
- One or more lines containing the sequence itself. 


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
$ cat sequences
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




## 2. FASTQ


A FASTQ file consists of four lines per sequence:

- Line 1 begins with a `@` character and is followed by a sequence identifier and an optional description.
- Line 2 is the raw sequence letters.
- Line 3 begins with a `+` character and is optionally followed by the sequence identifier and descriptions (if any).
- Line 4 encodes the quality values for the sequence in Line 2, and must contain the same number of symbols as letters in the sequence. Each letter in line 4 corresponds to the base quality score of the sequence with the same index in line 2. 

**2.1 Get the first read from lane 1 read 1**
```sh
$ zcat DATA_L001_R1.fastq.gz | head -4
@NS500474:51:HK7FVAFXX:1:11101:5669:1056 1:N:0:GACACT+NTGACG
CAGCTGGCGTCGGCCGACGTGATCACCTTCACGATCGGAAGAACACACGTCTGAACTCCAGTCACGACACTATCT
+
AAAAAEEEEEEEEEEEEEEEE6EEEEEEEEEEEEEEEEAEEEEE/EEEE/EEEEEAEEEEEEEEEEEEEEEE/<A
```

**2.2 Get the first read in R1 and R2 from lane 1 and 2**

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

**2.3 Return only the first 5 reads from a FASTQ file**
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
**2.4 Calculate the mean length**

```sh
$ zcat DATA_L001_R1.fastq.gz | awk 'NR%4==2{sum+=length($0)}END{print sum/(NR/4)}'
75.37
```

**2.5 Which lane contains the highest and lowest number of raw sequence reads?**

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

**2.6 How many times do we see the motif "GATTACA"?**
```sh
$ zcat DATA_L001_R1.fastq.gz | grep -c "GATTACA" 
1147
```



**2.7 Get fragment size statistics and distribution from a FASTQ file**

```sh
$ zcat DATA_L001_R1.fastq.gz | awk 'NR%4==2{print length($0)}' >> fragment_lengths_L001_R1
```


```sh
$ R
> d<- scan("fragment_lengths_L001_R1")
Read 1416689 items
> hist(d)
> summary(d)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  35.00   75.00   76.00   75.37   76.00   76.00 
```

## Adapter trimming using fastp


https://github.com/OpenGene/fastp



- `--out1` and `--out2` outputs the reads that cannot be merged successfully, but both pass all the filters.
- `--merged_out` specifies the file name to write the merged reads.
- `--unpaired1` outputs be the reads that cannot be merged, read1 passes filters but read2 doesn't.
- `--unpaired2` outputs be the reads that cannot be merged, read2 passes filters but read1 doesn't.
- By default, fastp uses overlap analysis to do adapter trimming. However, you can specify `--detect_adapter_for_pe` to enable adapter trimming by adapter sequence detection.



```sh
cd day1/exercises/trimming

fastp --in1 ../../data/fastq/DATA_L001_R1.fastq.gz --in2 ../../data/fastq/DATA_L001_R2.fastq.gz --out1 DATA_L001_R1_trimmed.fastq.gz --out2 DATA_L001_R2_trimmed.fastq.gz --merge --merged_out DATA_L001_merged_trimmed.fastq.gz --unpaired1 DATA_L001_unpaired_R1.fastq.gz --unpaired2 DATA_L001_unpaired_R2.fastq.gz --length_required 30 --detect_adapter_for_pe
```

It will output the following:

```sh
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

And the following files are generated:

```sh
.
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










Sort
```sh
samtools sort .bam -o .bam
```

Index a bam file
```sh
samtools index SRR1234567/SRR1234567.sorted.bam
```

Get HLCS gene (chr21:38,120,926-38,362,511) and write to a new bam file
```sh
samtools view sample.bam chr21:38120926-38362511 -b >s HLCS.bam
```







### SAM, BAM and CRAM files


samtools tveiw
cram
index

count the number of reads in bam files
```sh
samtools view -c sample.bam
```

count the number of mapped reads
```sh
samtools view -F 4 -c SRR1234567/SRR1234567.bam
```

count the number of unmapped reads in sam/bam files
```sh
samtools view -f 4 -c SRR1234567/SRR1234567.bam
```

show statistics of sam/bam files
```sh
samtools stats SRR1234567/SRR1234567.bam
```

check whether my alignment file is sorted lexicographically (list all the chromosomes)
```sh
samtools view SRR1234567/SRR1234567.bam | cut -f3 | uniq
```



### alignment/mapping


### bcf/vcf


### Mpileup
```bash
samtools mpileup ${FILE}
```
#generates vcf
./bcftools/bcftools mpileup -b bams.list --no-reference |less
./bcftools/bcftools mpileup -b bams.list --fasta-ref hs37d5.fa.gz |./bcftools/bcftools call -m -v|less

./bcftools/bcftools mpileup -b bams.list --fasta-ref hs37d5.fa.gz -r 1 |./bcftools/bcftools call -m -v -Ob -o res.bcf
bcftools +fill-tags res.bcf -Ob -o res.with.af.bcf
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' res.with.af.bcf
./angsd/angsd -b bams.list -domajorminor 1 -gl 1 -domaf 1 -r 1 -snp_pval 1e-6
gunzip -c angsdput.mafs.gz




### Coverage/Depth
```bash
samtools depth ${FILE}
```

### Stats
```bash
samtools stats ${FILE}
```

```bash
samtools quickcheck -vvv ${FILE}
```

### Using Flags
https://broadinstitute.github.io/picard/explain-flags.html


### Damage
```bash
 mapDamage -i ${FILE} -r reference.fasta
```

### Extract regions
```bash
samtools view ${FILE}
```


### Using 1000G bed
```bash
samtools view ${FILE} in.bed
```



### vcf/bcf
### mpileup


# generates vcf
./bcftools/bcftools mpileup -b bams.list --no-reference |less
./bcftools/bcftools mpileup -b bams.list --fasta-ref hs37d5.fa.gz |./bcftools/bcftools call -m -v|less

./bcftools/bcftools mpileup -b bams.list --fasta-ref hs37d5.fa.gz -r 1 |./bcftools/bcftools call -m -v -Ob -o res.bcf
bcftools +fill-tags res.bcf -Ob -o res.with.af.bcf
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' res.with.af.bcf
./angsd/angsd -b bams.list -domajorminor 1 -gl 1 -domaf 1 -r 1 -snp_pval 1e-6
gunzip -c angsdput.mafs.gz




```sh
echo hi
```





### Extract regions
```bash
samtools view ${FILE}
```



bgzip (Blocked GNU Zip Format)
BAM, BCF and VCF file formats are typically bgzip compressed.
bgzip files can be uncompressed with gzip

### Other

- Base counts: How many As, Cs, Gs, Ts and Ns are there in the fastq file?
- GC content
- What is the coverage when we use 1000G
- Which chromosome has the most reads aligned to it? `samtools view -c`
- Mean read length, read length distribution
- How long is the reference genome? 

- How many reads has their mates were unmapped?
- How many insertions and deletions are there in the alignment?
- How many reads are there that was aligned to a region included in 1000G sites with a mapping quality 30 as minimum?
- Sort by name, sort by coordinates
- Cigar strings
- Write a code for counting the number of lines in each fastq file and write the results into a new file.

Class working directory: `/TEACHING/BIOINF21/`

