
## Software requirements

Essential:
- Standard unix tools
- [FastQC v0.11.3](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [cutadapt v1.12](http://cutadapt.readthedocs.io/en/stable/guide.html)
- [python v2.7.12](https://www.python.org/). Libraries:
  - [os](https://docs.python.org/2/library/os.html)
  - [gzip](https://docs.python.org/2/library/gzip.html)
  - [biopython](http://biopython.org/)
  - [collections](https://docs.python.org/2/library/collections.html)
  - [re](https://docs.python.org/2/library/re.html)

Optional:
- [slurm](https://slurm.schedmd.com/overview.html) cluster job scheduling system



## Tables

*Trypanosoma brucei* chromosome 2:

- Single-base resolution 5hmU [sites](Tryp_chr2_5hmUsites.bed)
- Chemical-enrichment 5hmU [regions](Tryp_chr2_5hmUregions.bed)



## Single-base resolution 5hmU sites

### Processing sequencing reads

#### Quality check

```bash
mkdir ../fastqc

for fq in *.fastq.gz
do
  bname=${fq%_L001_R1_001.fastq.gz}
  sbatch -J $bname -o ../fastqc/$bname.log --mem 4096 --wrap "fastqc --noextract --nogroup -q -o ../fastqc $fq"
done
```

#### Trim Illumina adaptors

```bash
mkdir ../fastq_trimmed

for fq in *.fastq.gz
do
  bname=${fq%_L001_R1_001.fastq.gz}
  sbatch -J $bname -o ../fastq_trimmed/$bname.log --mem 4096 --wrap "cutadapt -a AGATCGGAAGAGC -m 15 -q 20 -o ../fastq_trimmed/$fq $fq > ../fastq_trimmed/$bname.txt"
done
```


### Alignment

#### Full genome conversion and indexing

Conversion:

```python
ifile = open("reference/TriTrypDB-9.0_TbruceiTREU927_Genome.fasta", "r")
ilines = ifile.readlines()
ifile.close()

lines_t2c = []
lines_a2g = []

for line in ilines:
  if line.startswith(">"):
    line_t2c = line.split(" | ")[0] + "_TC\n"
    line_a2g = line.split(" | ")[0] + "_AG\n"
  else:
    line_t2c = line.replace("T", "C").replace("t", "C")
    line_a2g = line.replace("A", "G").replace("a", "G")
  lines_t2c.append(line_t2c)
  lines_a2g.append(line_a2g)

len(ilines) == len(lines_t2c) == len(lines_a2g) # True

ofile = open("reference/TriTrypDB-9.0_TbruceiTREU927_Genome.t2c.python.fasta", "w")
ofile.write("%s" % "".join(lines_t2c) + "".join(lines_a2g))
ofile.close()
```

Indexing:

```bash
bwa index reference/TriTrypDB-9.0_TbruceiTREU927_Genome.t2c.python.fasta
```


## Chemical-enrichment 5hmU regions

### Processing sequencing reads

Under construction ..

### Alignment

Under construction ..

