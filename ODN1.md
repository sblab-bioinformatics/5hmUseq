- Normal PCR, dATP dilution 1   rep1   fk457
- Normal PCR, dATP dilution 1   rep2   fk458
- Normal PCR, dATP dilution 1   rep3   fk187

- Bst, 37 deg Tpol, 2 mM MgSO4, dATP dilution 500   rep1   fk382
- Bst, 37 deg Tpol, 2 mM MgSO4, dATP dilution 500   rep2   fk391
- Bst, 37 deg Tpol, 2 mM MgSO4, dATP dilution 500   rep3   fk173

- Bst, 37 deg, Tpol, 10 mM MgSO4, dATP dilution 500   rep1   fk386
- Bst, 37 deg, Tpol, 10 mM MgSO4, dATP dilution 500   rep2   fk449
- Bst, 37 deg, Tpol, 10 mM MgSO4, dATP dilution 500   rep3   fk450

- Bst, 37 deg, Tpol, 10 mM MgSO4, dATP dilution 50   rep1   fk380
- Bst, 37 deg, Tpol, 10 mM MgSO4, dATP dilution 50   rep2   fk381

- Klenow exo-, 37 deg, Tpol, dATP dilution 500   rep1   fk419
- Klenow exo-, 37 deg, Tpol, dATP dilution 500   rep2   fk420

- VentR exo-, 37deg, Tpol, 10 mM MgSO4, dATP dilution 500   rep1   fk423
- VentR exo-, 37deg, Tpol, 10 mM MgSO4, dATP dilution 500   rep2   fk424

- Sulfolobus pol, 37deg, Tpol, 10 mM MgSO4, dATP dilution 500   rep1   fk459
- Sulfolobus pol, 37deg, Tpol, 10 mM MgSO4, dATP dilution 500   rep3   fk461

- Bst, 37 deg, Tpol, 10 mM MgSO4, noox ctrl, dATP dilution 500   rep1   fk421
- Bst, 37 deg, Tpol, 10 mM MgSO4, noox ctrl, dATP dilution 500   rep2   fk422
- Bst, 37 deg, Tpol, 10 mM MgSO4, noox ctrl, dATP dilution 500   rep3   fk277



## Software requirements
Essential:
- Standard unix tools: zcat, paste, awk, tr, sort, sed ...
- [FastQC v0.11.3](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [cutadapt v1.12](http://cutadapt.readthedocs.io/en/stable/guide.html)
- [getInsertAndBarcode.py](getInsertAndBarcode.py)
- [bedtools v2.26.0](http://bedtools.readthedocs.io/en/latest/)
- [seqtk 1.0-r45](https://github.com/lh3/seqtk)
- [bwa v0.7.15-r1140](http://bio-bwa.sourceforge.net/)
- [samtools v1.3.1](http://samtools.sourceforge.net/)
- [pysamstats v1.0.1](https://github.com/alimanfoo/pysamstats) (pysam 0.11.2.2)
- [tableCat.py](https://github.com/dariober/bioinformatics-cafe/blob/master/tableCat/tableCat.py)
- [python v2.7.12](https://www.python.org/). Libraries:

Optional:
- [lsf](https://www.ibm.com/us-en/marketplace/hpc-workload-management) cluster job scheduling system



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
    python getInsertAndBarcode.py $fq \
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
    | python getInsertAndBarcode.py - \
    | paste  - - - - \
    | sort -k1,1 \
    | groupBy -g 1 -c 2,3,4 -o first,first,first \
    | tr '\t' '\n' \
    | seqtk seq -r - > $bname.ins.fq
done
```



## Alignment

### Prepare and index reference template

```bash
echo '>mod10' > mod10.fa
echo 'ATCGAGAATCCCGGTGCCGATACCHACTCTTGHAGAA' >> mod10.fa
bwa index mod10.fa
```


### Align and count 

```bash
for fq in fk*.ins.fq
do
    bname=`basename $fq .fq`
    bsub "
    bwa mem -k 5 -T 10 -L 200 mod10.fa $fq \
    | samtools sort -o - -O bam -T aln.tmp.bam > ${bname}.bam &&
    samtools index ${bname}.bam &&
    pysamstats -D 1000000 -d -f mod10.fa --type variation_strand ${bname}.bam > ${bname}.var.txt"
done
```


### Prepare data files

```bash
tableCat.py -H -i *.ins.var.txt -r '.ins.var.txt' | sed 's/\r//g' > mod10.20160725.txt
```
