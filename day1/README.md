# Day 1 - Draft

## Files
fastq
bam/cram
fasta
vcf/bcf
mpileup

## Exercises

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



### Extract regions
```bash
samtools view ${FILE}
```


### Extract region from fasta

```bash
samtools faidx genome.fa chr:X-Y
```




### Other

- Base counts: How many As, Cs, Gs, Ts and Ns are there in the fastq file?
- GC content
- What is the coverage when we use 1000G
- Which chromosome has the most reads aligned to it? `samtools view -c`
- Mean read length, read length distribution
- How long is the reference genome? 
- ```bash
samtools faidx ${FILE}
```
- How many reads has their mates were unmapped?
- How many insertions and deletions are there in the alignment?
- How many reads are there that was aligned to a region included in 1000G sites with a mapping quality 30 as minimum?
- Sort by name, sort by coordinates
- Cigar strings

## Todo
- Download 1000G bed
- trim,map->bam/cram
