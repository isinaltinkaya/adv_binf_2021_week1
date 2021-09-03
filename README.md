# {ClassID}: {FullClassName}


Schedule, date, time, who is teaching when and what




Friday(10/9/2021, Lei Zhao): 
Given that the BWT algorithm is introduced, I will introduce the MapQ score, a measurement to assess alignment. Based on the On the accuracy of short read mapping.
Once the short reads fragments have been mapped, different read bases will be piled up at the each genome positions. But even the read bases at the same genome position can vary a lot due to uncertainty introduced by different error schemes. I will introduce an important concept, genotype likelihood, which helps to integrate the uncertainty into the further genomic inferences. For reference, please have a read at The genome analysis toolkit: a map reduce framework for analyzing next-generation dna sequencing data, Patterns of damage in genomic DNA sequences from a Neandertal and ReQON: a Bioconductor package for recalibrating quality scores from next-generation sequencing data.

#### Draft repository structure:

```bash
{REPOSITORY_NAME}
├── day1
│   ├── exercises
│   │   ├── exercise_1
│   │   └── exercise_2
│   └── class_notes
├── day2
│   ├── exercises
│   └── class_notes
├── data
└── README.md

```




#### Draft class directory structure:


```bash
/PATH/TO/{CLASSID}
├── student-id-1
│   └── {REPOSITORY}
├── student-id-2
│   └── {REPOSITORY}
└── student-id-3
    └── {REPOSITORY}
```

```bash
ssh {SERVER}
```

Class working directory: `PATH/TO/CLASSID`


```bash
cd /PATH/TO/{CLASSID}
mkdir {STUDENT-ID}
cd {STUDENT-ID}
```

Git is a version-control system software for tracking changes in a set of files, and is useful for coordinating work among multiple people working in collaboration. GitHub is one of the most popular git repository hosting services. Git is not covered within the scope of this class, and will only be used for downloading the exercise materials. To learn more about git and GitHub, see https://lab.github.com/ (not required).

```bash
git clone {REPO_LINK}
```


