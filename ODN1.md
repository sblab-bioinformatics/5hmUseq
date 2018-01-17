## Software requirements
Essential:
- Standard unix tools: zcat, paste, awk, tr ...
- [FastQC v0.11.3](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [cutadapt v1.12](http://cutadapt.readthedocs.io/en/stable/guide.html)
- [bwa v0.7.15-r1140](http://bio-bwa.sourceforge.net/)
- [samtools v1.3.1](http://samtools.sourceforge.net/)
- [bedtools v2.26.0](http://bedtools.readthedocs.io/en/latest/)
- [pysamstats v1.0.1](https://github.com/alimanfoo/pysamstats) (pysam 0.11.2.2)
- [python v2.7.12](https://www.python.org/). Libraries:
  - [os](https://docs.python.org/2/library/os.html)
- [R v3.3.2](https://www.r-project.org/). Libraries:
  - [data.table](https://cran.r-project.org/web/packages/data.table/index.html)
  - [edgeR](http://bioconductor.org/packages/release/bioc/html/edgeR.html)
  - [ggplot2](http://ggplot2.org/)
  - [ggrepel](https://cran.r-project.org/web/packages/ggrepel/index.html)

Optional:
- [slurm](https://slurm.schedmd.com/overview.html) cluster job scheduling system



## Processing sequencing reads

### Quality check

Just like in [ODN2.md](https://github.com/sblab-bioinformatics/5hmUseq/blob/master/ODN2.md#quality-check)


### Trim Illumina adaptors

```bash
cd fastq # directory containing the *.fastq.gz files
mkdir ../fastq_trimmed

for fq in *.fq.gz
do
    bsub -oo $fq.log "cutadapt -a AGATCGGAAGAGC --trimmed-only -M 90 $fq | gzip > ../fastq_trimmed/$fq"
done
```

Extract reads with both tags in forward in one file and both tags in reverse in another file:

```bash
cd ../fastq_trimmed

for fq in *.fq.gz
do
    bname=${fq/.*/}
    zcat $fq | paste - - - - | awk -v OFS='\t' -v FS='\t' '$2 ~ "ATCGAG" && $2 ~ "CGTGTC"' | tr '\t' '\n' > ${bname}.fw.fq
    zcat $fq | paste - - - - | awk -v FS='\t' '$2 ~ "CTCGAT" && $2 ~ "GACACG"' | tr '\t' '\n' > ${bname}.rc.fq
done
```

Now we need to get the random barcode and put it in the header. We also deduplicate based on barcode.

```bash
for fq in fk*.fw.fq
do
    bname=`basename $fq .fq`
    python ~/Tritume/getInsertAndBarcode.py $fq \
    | paste  - - - - \
    | sort -k1,1 \
    | groupBy -g 1 -c 2,3,4 -o first,first,first \
    | tr '\t' '\n' > $bname.ins.fq
done

## For reversed reads we revcomp first. Then extract insert and revecomp back.
for fq in fk*.rc.fq
do
    bname=`basename $fq .fq`
    seqtk seq -r $fq \
    | python ~/Tritume/getInsertAndBarcode.py - \
    | paste  - - - - \
    | sort -k1,1 \
    | groupBy -g 1 -c 2,3,4 -o first,first,first \
    | tr '\t' '\n' \
    | seqtk seq -r - > $bname.ins.fq
done
```

