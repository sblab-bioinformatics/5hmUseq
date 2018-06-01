
## Software requirements

Essential:
- Standard unix tools
- [FastQC v0.11.3](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [cutadapt v1.12](http://cutadapt.readthedocs.io/en/stable/guide.html)
- [python v2.7.12](https://www.python.org/). Libraries:
  - [os](https://docs.python.org/2/library/os.html)
  - [re](https://docs.python.org/2/library/re.html)
  - [sys](https://docs.python.org/2/library/sys.html)
  - [string](https://docs.python.org/2/library/string.html)  
- [bwa v0.7.15-r1140](http://bio-bwa.sourceforge.net/)
- [samtools v1.3.1](http://samtools.sourceforge.net/)
- [pysamstats v1.0.1](https://github.com/alimanfoo/pysamstats) (pysam 0.11.2.2)
- [bioawk v20110810](https://github.com/lh3/bioawk)
- [R v3.3.2](https://www.r-project.org/). Libraries:
  - [data.table v1.10.4](https://cran.r-project.org/web/packages/data.table/index.html)
  - [limma v3.30.11](http://bioconductor.org/packages/release/bioc/html/limma.html)
  - [ggplot2 v2.2.1](http://ggplot2.org/)
- [tableCat.py](https://github.com/dariober/bioinformatics-cafe/blob/master/tableCat/tableCat.py)

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

#### Read conversion

```bash
cd fastq_trimmed
mkdir ../fastq_trimmed_converted

for fq in *.fastq.gz
do
bname=`basename $fq`
zcat $fq \
| paste - - - - \
| python -c "import sys
n= 0
for line in sys.stdin:
    n+=1
    line= line.strip().split('\t')
    line[0]= '@' + str(n) + ':' + line[1]
    line[1]= line[1].replace('T', 'C').replace('t', 'c')
    print '\n'.join(line)
sys.stderr.write('N:' + str(n) + '\n')
" | gzip > ../fastq_trimmed_converted/${bname%%.*}.t2c.fastq.gz
done
```

#### Align

```bash
cd fastq_trimmed_converted
mkdir ../sam

for fq in *.fastq.gz
do
  bname=${fq%_R1_001.t2c.fastq.gz}
  sbatch -J $bname -o ../sam/$bname.log --mem 16384 --wrap "bwa mem -M -t 20 ../reference/TriTrypDB-9.0_TbruceiTREU927_Genome.t2c.python.fasta $fq > ../sam/$bname.sam"
done
```

#### Clean

```bash
cd sam

for sam in *.sam
do
bname=${sam%.sam}
nohup samtools view -h -F 2820 -q 10 $sam > $bname.clean.sam &
done
```

#### Rename chroms and fix sequence

```python
import sys, string, re, os

all_files = os.listdir("sam/")

for f in sorted(all_files):
  if "clean.sam" in f:
    print f
    ifile = open("sam/%s" % f, "r")
    ilines = ifile.readlines()
    ifile.close()
    lines_sam = []
    for line in ilines:
        if line.startswith("@"):
            if line.startswith("@SQ\t") and "_AG\t" in line:
                continue
            if line.startswith("@SQ\t") and "_TC\t" in line:
                line= line.replace("_TC\t", "\t")
            lines_sam.append(line.strip() + '\n')
            continue
        line= line.strip().split("\t")
        line[2]= re.sub("_TC$|_AG$", "", line[2])
        rname, rseq= line[0].split(":")
        if len(rseq) != len(line[9]):
            sys.stderr.write(line)
            sys.exit(1)
        line[0]= rname
        line[9]= rseq
        if int(line[1]) & 16 == 16:
            line[9]= line[9].upper().translate(string.maketrans("ACTGN", "TGACN"))[::-1]
        lines_sam.append("\t".join(line) + '\n')
    ofile = open("sam/%s" % (f.replace(".sam", "") + ".renamed.fixed.sam"), "w")
    ofile.write("%s" % "".join(lines_sam))
    ofile.close()
```

#### Convert sam to bam

```bash
cd sam

for sam in *.clean.renamed.fixed.sam
do
bname=${sam%.sam}
nohup samtools view -b $sam -o $bname.bam &
done
```

#### Merge and sort

```bash
cd sam

for id in `ls *.clean.renamed.fixed.bam | cut -d "_" -f 1-2 | sort | uniq`
do
sbatch -J $id --wrap "samtools merge $id.bam $id*.clean.renamed.fixed.bam && \
samtools sort -T ../tmp/$id -o $id.sorted.bam $id.bam"
done
```


### Extract 5hmU

#### Index and run pysamstats

```bash
cd bam
mkdir ../pysamstats

for bam in *.sorted.bam
do
  bname=${bam%.bam}
  sbatch -J $bname -o ../pysamstats/$bname.log --mem 16384 --wrap "samtools index $bam &&
  pysamstats -D 1000000 -f ../reference/TriTrypDB-9.0_TbruceiTREU927_Genome.fasta --type variation_strand $bam > ../pysamstats/$bname.pss.txt"
done
```

#### Count T (A) and C (G) in reference T fwd and A rev

```bash
cd pysamstats

for txt in *.pss.txt
do
  bname=`basename $txt .pss.txt`
  nohup bioawk -tc hdr '$ref == "T" || $ref == "A" { \
      if($ref == "T"){
          T=$T_fwd; C=$C_fwd
      } else if($ref == "A"){
          T=$A_rev; C=$G_rev
      }
      print \
        $chrom, \
        $pos, \
        $ref, \
        T,
        C
  }' $txt > ${bname}.depth.txt &
done
```

#### Defining significant sites

```r
library(data.table)
library(limma)
library(ggplot2)

# Set width
options(width = 300)

# Load input files and set column names
txt <- fread("tableCat.py -i *.sorted.depth.txt -r .sorted.depth.txt")
setnames(txt, names(txt), c('chrom', 'pos', 'ref', 'cnt_T', 'cnt_C', 'library_id'))

# Redefine library_id column fields
txt[library_id == "nooxctrl-1", library_id := "noox_1"]
txt[library_id == "nooxctrl-2", library_id := "noox_2"]
txt[library_id == "nooxctrl-3", library_id := "noox_3"]
txt[library_id == "5hmUseq-1", library_id := "ox_1"]
txt[library_id == "5hmUseq-2", library_id := "ox_2"]
txt[library_id == "5hmUseq-3", library_id := "ox_3"]

# Select chromosomes of interest
txt <- txt[chrom == "Tb927_01_v5.1" | chrom == "Tb927_02_v5.1" | chrom == "Tb927_03_v5.1" | chrom == "Tb927_04_v5.1" | chrom == "Tb927_05_v5.1" | chrom == "Tb927_06_v5.1" | chrom == "Tb927_07_v5.1" | chrom == "Tb927_08_v5.1" | chrom == "Tb927_09_v5.1" | chrom == "Tb927_10_v5.1" | chrom == "Tb927_11_v5.1" | chrom == "Tb927_11_RH_fork_v5.1" | chrom == "Tb927_11_Homologues_1_v5.1" | chrom == "Tb927_11_Homologues_2_v5.1" | chrom == "Tb927_11_Homologues_3_v5.1" | chrom == "Tb927_11_bin_v5.1"]

# Define pct conversion
txt[, pct_t2c := cnt_C / (cnt_C + cnt_T)]

# Define logit function
logit <- function(p){
    stopifnot(p >= 0); stopifnot(p <= 1);
    p<- ifelse(p == 0, min(p[p > 0]), p)
    p<- ifelse(p == 1, max(p[p < 1]), p)
    return( -log(1/p -1) )
}

# Cast table and select all positions
txt2_hmUseq_nooxctrl <- dcast.data.table(data = txt, chrom + pos ~ library_id, value.var= 'pct_t2c')

# Select positions with counts for all
txt2_hmUseq_nooxctrl <- txt2_hmUseq_nooxctrl[complete.cases(txt2_hmUseq_nooxctrl)]
dim(txt2_hmUseq_nooxctrl) #  10938119        8

pos_hmUseq_nooxctrl <- txt2_hmUseq_nooxctrl[, list(chrom, pos)]
txt2_hmUseq_nooxctrl[, chrom := NULL]
txt2_hmUseq_nooxctrl[, pos := NULL]

# Logit, design, fitting and top sites
mat_hmUseq_nooxctrl <- data.table(apply(data.frame(txt2_hmUseq_nooxctrl), 2, logit))
design_hmUseq_nooxctrl <- cbind(intercept = 1, hmUseq_nooxctrl = ifelse(grepl('^ox', names(mat_hmUseq_nooxctrl)), 1, 0)) # 1 = ox, 0 = noox
fit_hmUseq_nooxctrl <- lmFit(mat_hmUseq_nooxctrl, design_hmUseq_nooxctrl)
fit_hmUseq_nooxctrl <- eBayes(fit_hmUseq_nooxctrl)
tt_hmUseq_nooxctrl <- data.table(topTable(fit_hmUseq_nooxctrl, coef = 2, number = Inf, sort.by = 'none'))
tte_hmUseq_nooxctrl <- cbind(pos_hmUseq_nooxctrl, tt_hmUseq_nooxctrl, txt2_hmUseq_nooxctrl * 100)

tte_hmUseq_nooxctrl[, ox_avg := (ox_1 + ox_2 + ox_3)/3]
tte_hmUseq_nooxctrl[, noox_avg := (nooxl_1 + noox_2 + noox_3)/3]

tte_hmUseq_nooxctrl[order(P.Value)][logFC > 0 & adj.P.Val < 0.1 & hmUseq_avg > 5 & nooxctrl_avg < 5] # 43434
tte_hmUseq_nooxctrl[order(P.Value)][logFC > 0 & adj.P.Val < 0.1 & hmUseq_avg > 6 & nooxctrl_avg < 5] # 14600
tte_hmUseq_nooxctrl[order(P.Value)][logFC > 0 & adj.P.Val < 0.1 & hmUseq_avg > 7 & nooxctrl_avg < 5] # 5867
tte_hmUseq_nooxctrl[order(P.Value)][logFC > 0 & adj.P.Val < 0.1 & hmUseq_avg > 8 & nooxctrl_avg < 5] # 2898
tte_hmUseq_nooxctrl[order(P.Value)][logFC > 0 & adj.P.Val < 0.1 & hmUseq_avg > 9 & nooxctrl_avg < 5] # 1717
tte_hmUseq_nooxctrl[order(P.Value)][logFC > 0 & adj.P.Val < 0.1 & hmUseq_avg > 10 & nooxctrl_avg < 5] # 1189

tte_hmUseq_nooxctrl[order(P.Value)][logFC > 0 & adj.P.Val < 0.1 & hmUseq_avg > 8 & nooxctrl_avg < 4] # 2898
tte_hmUseq_nooxctrl[order(P.Value)][logFC > 0 & adj.P.Val < 0.1 & hmUseq_avg > 8 & nooxctrl_avg < 3] # 2898
tte_hmUseq_nooxctrl[order(P.Value)][logFC > 0 & adj.P.Val < 0.1 & hmUseq_avg > 8 & nooxctrl_avg < 2] # 2893
tte_hmUseq_nooxctrl[order(P.Value)][logFC > 0 & adj.P.Val < 0.1 & hmUseq_avg > 8 & nooxctrl_avg < 1] # 2884
tte_hmUseq_nooxctrl[order(P.Value)][logFC > 0 & adj.P.Val < 0.1 & hmUseq_avg > 8 & nooxctrl_avg < 6] # 2898
tte_hmUseq_nooxctrl[order(P.Value)][logFC > 0 & adj.P.Val < 0.1 & hmUseq_avg > 8] # 2975
tte_hmUseq_nooxctrl[order(P.Value)][logFC > 0 & adj.P.Val < 0.1 & hmUseq_avg > 8 & nooxctrl_avg > 90] # 41

tte_hmUseq_nooxctrl_adjpval0.1_s_5 <- tte_hmUseq_nooxctrl[order(P.Value)][logFC > 0 & adj.P.Val < 0.1 & hmUseq_avg > 5 & nooxctrl_avg < 5]
tte_hmUseq_nooxctrl_adjpval0.1_s_6 <- tte_hmUseq_nooxctrl[order(P.Value)][logFC > 0 & adj.P.Val < 0.1 & hmUseq_avg > 6 & nooxctrl_avg < 5]
tte_hmUseq_nooxctrl_adjpval0.1_s_7 <- tte_hmUseq_nooxctrl[order(P.Value)][logFC > 0 & adj.P.Val < 0.1 & hmUseq_avg > 7 & nooxctrl_avg < 5]
tte_hmUseq_nooxctrl_adjpval0.1_s_8 <- tte_hmUseq_nooxctrl[order(P.Value)][logFC > 0 & adj.P.Val < 0.1 & hmUseq_avg > 8 & nooxctrl_avg < 5]
tte_hmUseq_nooxctrl_adjpval0.1_s_9 <- tte_hmUseq_nooxctrl[order(P.Value)][logFC > 0 & adj.P.Val < 0.1 & hmUseq_avg > 9 & nooxctrl_avg < 5]
tte_hmUseq_nooxctrl_adjpval0.1_s_10 <- tte_hmUseq_nooxctrl[order(P.Value)][logFC > 0 & adj.P.Val < 0.1 & hmUseq_avg > 10 & nooxctrl_avg < 5]

# Save significant sites in bed files
tte_hmUseq_nooxctrl_adjpval0.05_s_8[, start := as.integer(pos-1)]
tte_hmUseq_nooxctrl_adjpval0.05_s_8[, end := as.integer(start+1)]

write.table(tte_hmUseq_nooxctrl_adjpval0.05_s_8[, c("chrom", "start", "end", "logFC", "P.Value", "adj.P.Val", "ox_avg", "noox_avg")], "../bed/hmUseq_nooxctrl_nocountfilter_adjpval0.05_hmUseqavg8_nooxctrlavg5.bed", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
system('bedtools sort -i ../bed/hmUseq_nooxctrl_nocountfilter_adjpval0.05_hmUseqavg8_nooxctrlavg5.bed > ../bed/hmUseq_nooxctrl_nocountfilter_adjpval0.05_hmUseqavg8_nooxctrlavg5.sorted.bed')

```



## Chemical-enrichment 5hmU regions

### Processing sequencing reads

Under construction ..

### Alignment

Under construction ..

