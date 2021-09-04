# Day 1: NGS data - Workflow, formats and programs




# Getting started


## Connecting to the server via SSH





X11 forwarding method will allow you to start a graphical application on the remote system and forward this application's windows to your local system. We need to enable X11 forwarding to view the plots we will be generating for the exercises.

We use `-X` option to enable X11 forwarding over SSH:

```sh
ssh -X <your_username>@<server_name_or_ip>
```

Replace with your remote server username. For example:


```sh
ssh -X isin@ricco.popgen.dk
```




## Working environment setup


Class working directory: `/TEACHING/BIOINF21/`
Student working directory: `/TEACHING/BIOINF21/students/<student_id>`


First, we go to the class directory:
```sh
cd /TEACHING/BIOINF21/
```

We have 3 main directories.

```bash
/TEACHING/BIOINF21/
├── data
├── programs
└── students
```




```sh
cd students
mkdir <student_id>
cd <student_id>
git clone https://github.com/{REPLACE_ME}
```





## Exercises



- Sort
```sh
samtools sort .bam -o .bam
```

- Index a bam file
```sh
samtools index SRR1234567/SRR1234567.sorted.bam
```

- Get HLCS gene (chr21:38,120,926-38,362,511) and write to a new bam file
```sh
samtools view sample.bam chr21:38120926-38362511 -b >s HLCS.bam
```



- Count the total number of bases in a FASTA file
```sh
grep -v ">" imputfile.fasta | wc | awk '{print $3 - $1}'
```


- Calculate the fragment lengths in fastq files 
```sh
zcat DATA_L001_R1.fastq.gz | awk '{if(NR%4==2) print length($1)}' > DATA_L001_R1_length.txt
```



### fasta

FASTA format consists of:

- One line starting with a ">" sign followed by a sequence identifier.
- One or more lines containing the sequence itself. 
- Calculate the fragment lengths in fasta files 


```sh
cat ref.fasta | awk '{if(NR%4==0) print length($1)}' > fasta_length.txt
```


- Extract ids from fasta file
```sh
grep -o -E "^>\w+" file.fasta | tr -d ">"
```



- Get a histogram of sequence lengths from FASTA files 
cat ref.fasta | awk '{if(NR%4==0) print length($1)}' | sort -n | uniq -c


Download chr21 
```sh
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr21.fa.gz' -O chr21.fa.gz
```



- Extract region from fasta

```sh
samtools faidx genome.fa chr:X-Y
```





### fastq


- Extract only the sequences from a fastq file
```sh
sed -n '2~4p' file.fq
```

- Mean length
```sh
awk 'NR%4==2{sum+=length($0)}END{print sum/(NR/4)}' input.fastq
```

- Convert SAM file to BAM file
```sh
samtools view -h -b -S aln-pe.sam > aln-pe.bam
```


- Which lane contains the highest and lowest number of raw sequence reads?


- How many reads in each fastq file contain the motif "GATTACA"?

- 
- Get a histogram of sequence lengths from FASTQ files 
zcat DATA_L001_R1.fastq.gz | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c



- Count the number of reads in fastq files

```sh
cat sample.fastq | echo $((`wc -l` / 4))
```

or

```sh
cat DATA_L001_R1.fastq.gz | grep -c @
```

- Count the number of reads in gzipped fastq files
- 
- 
```shell
zcat sample.fastq.gz | echo $((`wc -l` / 4))
```


### Trimming



                    
                    
Using fastp we can
1. Assess the quality of our sequence fastq data (in the output HTML report file) 
2. Trim off low quality bases
3. Detect and remove the adapter sequencenes

```sh
fastp -i sample1_r1.fastq -I sample1_r2.fastq -o sample1_out.R1.fq.gz -O sample1_out.R2.fq.gz --html \
 sample1_results.html --json sample1_results.json --report_title sample1_results
``` 
 
 


```sh
fastp --in1 DATA_L001_R1.fastq.gz --in2 DATA_L001_R2.fastq.gz --out1 "DATA_L001_R1_trimmed.fastq.gz" --out2 "DATA_L001_R1_trimmed.fastq.gz" -A -g -Q -L -w ${task.cpus} --json "DATA_L001.json" 
```


```sh
fastp --merge --in1 R1.fq --in2 R2.fq --out1 R1.trimmed.fq --out2 R2.trimmed.fq --merged_out merged.fq --html fastp.html --json fastp.json --report_title "SAMPLE" --thread 4
```
                    
                    



### sam/bam/cram



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

 get the most common bitwise flags to identify alignment issues
```sh
cat $FILE | grep -v "^@" | cut -f 2 | sort | uniq -c | sort -nr | head -20
```






- Check if BAM file is sorted
```bash
samtools view -H test.bam | grep @HD
```
Will show `SO:coordinate` if sorted by coordinates.




### bcf/vcf
```sh
./bcftools/bcftools mpileup -b bams.list --no-reference |less
./bcftools/bcftools mpileup -b bams.list --fasta-ref hs37d5.fa.gz |./bcftools/bcftools call -m -v|less

./bcftools/bcftools mpileup -b bams.list --fasta-ref hs37d5.fa.gz -r 1 |./bcftools/bcftools call -m -v -Ob -o res.bcf
bcftools +fill-tags res.bcf -Ob -o res.with.af.bcf
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' res.with.af.bcf
./angsd/angsd -b bams.list -domajorminor 1 -gl 1 -domaf 1 -r 1 -snp_pval 1e-6
gunzip -c angsdput.mafs.gz
```



### Mpileup
```bash
samtools mpileup ${FILE}
```


### Coverage/Depth
```bash
samtools depth ${FILE}
```

### Stats
```bash
samtools stats ${FILE}
plot-bamstats
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

- Extract regions
```bash
samtools view ${FILE}
```


- Using 1000G bed
```bash
samtools view ${FILE} in.bed
```



### vcf/bcf
### mpileup




### Other-ideas

- Base counts: How many As, Cs, Gs, Ts and Ns are there in the fastq file?
- GC content
- What is the coverage when we use 1000G
- Which chromosome has the most reads aligned to it? `samtools view -c`
- Mean read length, read length distribution
- How long is the reference genome? 
```bash
samtools faidx ${FILE}
```
- How many reads has their mates were unmapped?
- How many insertions and deletions are there in the alignment?
- How many reads are there that was aligned to a region included in 1000G sites with a mapping quality 30 as minimum?
- Sort by name, sort by coordinates
- Cigar strings
- Write a code for counting the number of lines in each fastq file and write the results into a new file.





```sh
bwa aln -t 4 reference.fasta R1.fq -f R1.sai
bwa aln -t 4 reference.fasta R2.fq -f R2.sai
bwa sampe -r "@RG\\tID:ILLUMINA-${libraryid}\\tSM:${libraryid}\\tPL:illumina\\tPU:ILLUMINA-${libraryid}-${seqtype}" reference.fasta R1.sai R2.sai R1.fq R2.fq | samtools sort -@ 4 -O bam - > SAMPLE.mapped.bam
samtools index SAMPLE.mapped.bam
```

```sh
# PE collapsed, or SE data 
bwa aln -t 2 reference.fasta R1.fq  -f collapsed.sai
bwa samse -r "@RG\\tID:ILLUMINA-${libraryid}\\tSM:${libraryid}\\tPL:illumina\\tPU:ILLUMINA-${libraryid}-${seqtype}" reference.fasta R1.sai R1.fq  | samtools sort -@ 4 -O bam - > SAMPLE.mapped.bam
samtools index SAMPLE.mapped.bam
```
  
  
  
 
 

```

bwa mem reference/chrM.fa sample1_out.R1.fq.gz sample1_out.R2.fq.gz -R '@RG\tID:sample1\tSM:sample1\tLB:sample1\tPL:ILLUMINA' -o sample1.sam 



```



- Sort a BAM file
```bash
samtools sort -o test_sorted.bam -T tmp test.bam
```
- Extract run ID, flow cell ID and Lane number
`@NS500474:51:HK7FVAFXX:1:11101:5669:1056 1:N:0:GACACT+NTGACG` 




- Extract sample name
```bash
samtools view -H f.bam | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq
```



- Change sample name 
For all read groups in a BAM file
```bash
samtools view -H f.bam  | sed "s/SM:[^\t]*/SM:NEW_NAME/g" | samtools reheader - test.bam > test_SM.bam
```



- Calculate coverage for each position of a bed
```sh
bedtools coverage -d -a my_bed.bed -b my_bam.bam
```

- Find where reads map in a BAM file

only report non-zero coverage
and create contiguous regions with similar coverage.
only keep regions with at least 20 reads:
```sh
bedtools genomecov -bg -ibam test.bam | awk '$4>=20'

```

- Converting a SAM file to a BAM file

``` sh
samtools view -h i.bam > i.sam
```


check size of SAM file
compare sam bam fsizes

```sh
ls -lh i.sam
```


``` sh
ls -lh i.bam
```

``` sh
samtools coverage i.bam
```


