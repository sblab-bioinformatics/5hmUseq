```
-----------A-----A------A-A-A-A
...<<------T-----T------T-T-T-T
           *     *             
positions  35    43

* = 12% level of incorporation hmU/T

in the reverse strand - check T->C conversion:

- hmU --(hmU-seq)--> fU --(read as)--> C
- T --(hmU-seq)--> T --(read as)--> T

in the forwad strand - check A->G conversion:

- hmU --(hmU-seq)--> fU --(read as)--> G
- T --(hmU-seq)--> T --(read as)--> A
```

- 25% incorporation libraries:
  - F8-6-M7v2-25pct-ox-1   ox
  - F8-6-M7v2-25pct-ox-2   ox
  - F8-6-M7v2-25pct-ox-3   ox
  - F8-6-M7v2-25pct-nooxctrl-1   nooxctrl
  - F8-6-M7v2-25pct-nooxctrl-2   nooxctrl
  - F8-6-M7v2-25pct-nooxctrl-3   nooxctrl



## Software requirements
Essential:
- Standard unix tools: mv, mkdir, cd, echo ...
- [FastQC v0.11.3](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [cutadapt v1.12](http://cutadapt.readthedocs.io/en/stable/guide.html)
- 

Optional:
- [Slurm](https://slurm.schedmd.com/overview.html) cluster job scheduling system



## Processing sequencing reads

### Renaming fastq files

Adding the sequencing run id (180109_M00886_0188_000000000_BGP8K) to the file name.

```bash
cd fastq # fastq = directory containing the six *.fastq.gz files shown above

for fastq in *.fastq.gz
do
  mv $fastq 180109_M00886_0188_000000000_BGP8K_$fastq
done
```


### Quality check

```bash
mkdir ../fastqc

for fq in *.fastq.gz
do
  bname=${fq%_L001_R1_001.fastq.gz}
  sbatch -J $bname -o ../fastqc/$bname.log --mem 4096 --wrap "fastqc --noextract --nogroup -q -o ../fastqc $fq"
done
```


### Trim Illumina adaptors

```bash
mkdir ../fastq_trimmed

for fq in *.fastq.gz
do
  bname=${fq%_L001_R1_001.fastq.gz}
  sbatch -J $bname -o ../fastq_trimmed/$bname.log --mem 4096 --wrap "cutadapt -a AGATCGGAAGAGC -m 15 -q 20 -o ../fastq_trimmed/$fq $fq > ../fastq_trimmed/$bname.txt"
done
```



## Alignment

### Prepare and index reference template

```bash
mkdir ../reference
cd ../reference
echo -e ">template\nGCTCGCTTTGTTGGTTTCCTTGTTCTCTGTGCCCACTGCCTGACGGGCGGAAAGCAGCGCGAGCAAGCGAGACAGGACAC" > template.fa
bwa index template.fa
```






