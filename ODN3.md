
- F8-6-M7v2-25pct-ox-1   ox   id:171214_M00886_0185_000000000-BGP8V   (Fig. 2c)
- F8-6-M7v2-25pct-ox-2   ox   id:171214_M00886_0185_000000000-BGP8V   (Fig. 2c)
- F8-6-M7v2-25pct-ox-3   ox   id:171214_M00886_0185_000000000-BGP8V   (Fig. 2c)
- F8-6-M7v2-25pct-nooxctrl-1   nooxctrl   id:171214_M00886_0185_000000000-BGP8V   (Fig. 2c)
- F8-6-M7v2-25pct-nooxctrl-2   nooxctrl   id:171214_M00886_0185_000000000-BGP8V   (Fig. 2c)
- F8-6-M7v2-25pct-nooxctrl-3   nooxctrl   id:171214_M00886_0185_000000000-BGP8V   (Fig. 2c)



## Software requirements
Essential:
- Standard unix tools: mv, mkdir, cd, cp, echo ...
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



## Processing sequencing reads

### Renaming fastq files

Adding the sequencing run id (171214_M00886_0185_000000000-BGP8V) to the file name.

```bash
cd fastq # fastq = directory containing the *.fastq.gz files shown above

for fastq in *.fastq.gz
do
  mv $fastq 171214_M00886_0185_000000000-BGP8V_$fastq
done
```


### Quality check

In `clust1-headnode`,

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



## Check dinucleotide frequency around 5hmU site

```python
import gzip
import re
from Bio.Seq import Seq
from collections import Counter
import os

files = os.listdir("../fastq_trimmed")

for f in sorted(files):
  if "fastq.gz" in f:
    print f
    ifile = gzip.open('../fastq_trimmed/%s' % f, 'rb')
    ireads = ifile.read().split('\n')[1:-1:4]
    ifile.close()
    print "Total number of reads:", len(ireads)
    fwd = []
    fwd_d = {}
    rev = []
    rev_d = {}
    for read in ireads:
      if "N" not in read:
        #print read
        fwd_m = re.search(r"CTGTGGCTCTGCGTCCTTGTCCT[ACGT]{6}ACACAGCGCA", read)
        rev_m = re.search(r"TGCGCTGTGT[ACGT]{6}AGGACAAGGACGCAGAGCCACAG", read)
        if len(read) >= 42:
          if fwd_m and not rev_m and fwd_m.end() < 146:
            fwd_dinucleotide = "%s%s" % (read[fwd_m.end()], read[fwd_m.end() + 2])
            fwd_5hmU = read[fwd_m.end() + 1]
            fwd_barcode = read[(fwd_m.end() - 16):(fwd_m.end() - 10)]
            fwd.append((fwd_dinucleotide, fwd_5hmU, fwd_barcode))
            #print fwd_dinucleotide, fwd_5hmU, fwd_barcode
            if (fwd_dinucleotide, fwd_barcode) in fwd_d:
              fwd_d[(fwd_dinucleotide, fwd_barcode)].append(fwd_5hmU)
            else:
              fwd_d[(fwd_dinucleotide, fwd_barcode)] = [fwd_5hmU]
          if rev_m and not fwd_m and rev_m.start() > 3:
            rev_dinucleotide = "%s%s" % (str(Seq(read[rev_m.start() - 1]).reverse_complement()), str(Seq(read[rev_m.start() - 3]).reverse_complement()))
            rev_5hmU = str(Seq(read[rev_m.start() - 2]).reverse_complement())
            rev_barcode = str(Seq(read[(rev_m.start() + 10):(rev_m.start() + 16)]).reverse_complement())
            rev.append((rev_dinucleotide, rev_5hmU, rev_barcode))
            #print rev_dinucleotide, rev_5hmU, rev_barcode
            if (rev_dinucleotide, rev_barcode) in rev_d:
              rev_d[(rev_dinucleotide, rev_barcode)].append(rev_5hmU)
            else:
              rev_d[(rev_dinucleotide, rev_barcode)] = [rev_5hmU]
    print "Total number of matching fwd reads:", len(fwd)
    print "Total number of matching fwd dinucleotide and barcode combinations:", len(fwd_d)
    fwd_dr = {}
    for p in fwd_d:
      random.seed(1)
      fwd_dr[p] = [random.choice(fwd_d[p])]
    print "Total number of matching fwd dinucleotide and barcode combinations (random):", len(fwd_dr)
    print "dn\tlen\tT\tC\tTtoC\tlen_r\tT_r\tC_r\tTtoC_r"
    for dn in sorted(set([p[0] for p in fwd_d.keys()])):
      l = []
      l2 = []
      for p in fwd_d:
        if p[0] == dn:
          l += fwd_d[p]
          l2 += fwd_dr[p]
      print "%s\t%i\t%i\t%i\t%f\t%i\t%i\t%i\t%f" % (dn, len(l), l.count("T"), l.count("C"), round(100*float(l.count("C"))/(l.count("C")+l.count("T")), 2), len(l2), l2.count("T"), l2.count("C"), round(100*float(l2.count("C"))/(l2.count("C")+l2.count("T")), 2))
    print "Total number of matching rev reads:", len(rev)
    print "Total number of matching rev dinucleotide and barcode combinations:", len(rev_d)
    rev_dr = {}
    for p in rev_d:
      random.seed(1)
      rev_dr[p] = [random.choice(rev_d[p])]
    print "Total number of matching rev dinucleotide and barcode combinations (random):", len(rev_dr)
    print "dn\tlen\tT\tC\tTtoC\tlen_r\tT_r\tC_r\tTtoC_r"
    for dn in sorted(set([p[0] for p in rev_d.keys()])):
      l = []
      l2 = []
      for p in rev_d:
        if p[0] == dn:
          l += rev_d[p]
          l2 += rev_dr[p]
      print "%s\t%i\t%i\t%i\t%f\t%i\t%i\t%i\t%f" % (dn, len(l), l.count("T"), l.count("C"), round(100*float(l.count("C"))/(l.count("C")+l.count("T")), 2), len(l2), l2.count("T"), l2.count("C"), round(100*float(l2.count("C"))/(l2.count("C")+l2.count("T")), 2))
    print "\n"

#171214_M00886_0185_000000000-BGP8V_F7-192-M20-nooxctrl-1_S4_L001_R1_001.fastq.gz
#Total number of reads: 605962
#Total number of matching fwd reads: 180719
#Total number of matching fwd dinucleotide and barcode combinations: 51615
#Total number of matching fwd dinucleotide and barcode combinations (random): 51615
#dn	len	T	C	TtoC	len_r	T_r	C_r	TtoC_r
#AA	19054	18959	63	0.330000	3668	3647	13	0.360000
#AC	11708	11673	18	0.150000	3357	3343	6	0.180000
#AG	19480	19351	88	0.450000	3701	3674	19	0.510000
#AT	11293	11226	29	0.260000	3342	3325	10	0.300000
#CA	9430	9344	26	0.280000	3093	3058	10	0.330000
#CC	5543	5525	10	0.180000	2588	2574	7	0.270000
#CG	9807	9746	23	0.240000	3206	3183	11	0.340000
#CT	5525	5509	10	0.180000	2631	2624	4	0.150000
#GA	16014	15902	34	0.210000	3538	3506	7	0.200000
#GC	9588	9561	12	0.130000	3149	3142	3	0.100000
#GG	15847	15741	57	0.360000	3546	3522	12	0.340000
#GT	9379	9338	29	0.310000	3149	3140	5	0.160000
#TA	11353	11213	37	0.330000	3344	3288	12	0.360000
#TC	6477	6460	8	0.120000	2819	2812	3	0.110000
#TG	12966	11702	311	2.590000	3506	3129	101	3.130000
#TT	7255	7236	15	0.210000	2978	2972	6	0.200000
#Total number of matching rev reads: 185741
#Total number of matching rev dinucleotide and barcode combinations: 52325
#Total number of matching rev dinucleotide and barcode combinations (random): 52325
#dn	len	T	C	TtoC	len_r	T_r	C_r	TtoC_r
#AA	20073	20033	6	0.030000	3700	3696	0	0.000000
#AC	12005	11959	7	0.060000	3411	3399	0	0.000000
#AG	19282	19177	50	0.260000	3705	3681	14	0.380000
#AT	11723	11657	18	0.150000	3402	3386	5	0.150000
#CA	9750	9708	3	0.030000	3166	3151	1	0.030000
#CC	5571	5543	8	0.140000	2642	2635	2	0.080000
#CG	9724	9674	13	0.130000	3173	3162	3	0.090000
#CT	5781	5768	3	0.050000	2682	2678	1	0.040000
#GA	16804	16732	6	0.040000	3637	3628	0	0.000000
#GC	9851	9831	3	0.030000	3199	3195	1	0.030000
#GG	16787	16728	23	0.140000	3618	3610	3	0.080000
#GT	9617	9603	4	0.040000	3199	3194	0	0.000000
#TA	11766	11639	8	0.070000	3417	3369	2	0.060000
#TC	6785	6766	4	0.060000	2873	2865	2	0.070000
#TG	12889	11761	256	2.130000	3535	3202	90	2.730000
#TT	7333	7299	9	0.120000	2966	2950	3	0.100000
#
#
#171214_M00886_0185_000000000-BGP8V_F7-192-M20-nooxctrl-2_S5_L001_R1_001.fastq.gz
#Total number of reads: 318502
#Total number of matching fwd reads: 91346
#Total number of matching fwd dinucleotide and barcode combinations: 40247
#Total number of matching fwd dinucleotide and barcode combinations (random): 40247
#dn	len	T	C	TtoC	len_r	T_r	C_r	TtoC_r
#AA	9822	9779	31	0.320000	3171	3158	10	0.320000
#AC	5881	5857	11	0.190000	2595	2586	2	0.080000
#AG	9783	9722	38	0.390000	3138	3123	10	0.320000
#AT	5753	5707	17	0.300000	2634	2615	9	0.340000
#CA	4663	4612	17	0.370000	2361	2335	8	0.340000
#CC	2851	2834	7	0.250000	1802	1793	2	0.110000
#CG	4888	4862	14	0.290000	2432	2417	10	0.410000
#CT	2866	2861	3	0.100000	1827	1823	2	0.110000
#GA	8034	7980	24	0.300000	2958	2933	9	0.310000
#GC	4941	4928	2	0.040000	2416	2409	0	0.000000
#GG	8050	8013	17	0.210000	3001	2985	7	0.230000
#GT	4795	4774	13	0.270000	2381	2374	6	0.250000
#TA	5614	5548	17	0.310000	2605	2573	8	0.310000
#TC	3207	3196	7	0.220000	1969	1962	3	0.150000
#TG	6614	5998	129	2.110000	2831	2546	60	2.300000
#TT	3584	3570	11	0.310000	2126	2117	8	0.380000
#Total number of matching rev reads: 105036
#Total number of matching rev dinucleotide and barcode combinations: 42649
#Total number of matching rev dinucleotide and barcode combinations (random): 42649
#dn	len	T	C	TtoC	len_r	T_r	C_r	TtoC_r
#AA	11179	11149	8	0.070000	3248	3236	3	0.090000
#AC	6665	6649	1	0.020000	2751	2744	0	0.000000
#AG	10981	10927	25	0.230000	3290	3276	8	0.240000
#AT	6741	6695	3	0.040000	2828	2808	1	0.040000
#CA	5441	5417	1	0.020000	2565	2550	1	0.040000
#CC	3162	3145	2	0.060000	1927	1916	1	0.050000
#CG	5568	5540	8	0.140000	2575	2564	5	0.190000
#CT	3359	3349	2	0.060000	1981	1975	2	0.100000
#GA	9615	9576	3	0.030000	3128	3115	2	0.060000
#GC	5498	5484	1	0.020000	2576	2570	0	0.000000
#GG	9397	9357	8	0.090000	3107	3096	3	0.100000
#GT	5718	5707	2	0.040000	2595	2592	0	0.000000
#TA	6513	6436	8	0.120000	2741	2706	2	0.070000
#TC	3753	3746	1	0.030000	2133	2131	0	0.000000
#TG	7304	6670	149	2.190000	2937	2655	69	2.530000
#TT	4142	4123	3	0.070000	2267	2254	1	0.040000
#
#
#171214_M00886_0185_000000000-BGP8V_F7-192-M20-nooxctrl-3_S6_L001_R1_001.fastq.gz
#Total number of reads: 843895
#Total number of matching fwd reads: 245208
#Total number of matching fwd dinucleotide and barcode combinations: 55822
#Total number of matching fwd dinucleotide and barcode combinations (random): 55822
#dn	len	T	C	TtoC	len_r	T_r	C_r	TtoC_r
#AA	25634	25526	75	0.290000	3833	3809	18	0.470000
#AC	15603	15551	23	0.150000	3579	3567	5	0.140000
#AG	26297	26133	112	0.430000	3858	3835	13	0.340000
#AT	15139	15026	44	0.290000	3558	3533	7	0.200000
#CA	12878	12770	36	0.280000	3469	3443	9	0.260000
#CC	7703	7668	17	0.220000	2981	2969	6	0.200000
#CG	13037	12942	37	0.290000	3471	3445	10	0.290000
#CT	7525	7498	14	0.190000	2948	2939	3	0.100000
#GA	22224	22075	45	0.200000	3782	3762	4	0.110000
#GC	13142	13104	22	0.170000	3445	3437	2	0.060000
#GG	22088	21925	83	0.380000	3774	3752	13	0.350000
#GT	12968	12918	32	0.250000	3477	3458	9	0.260000
#TA	15260	15089	37	0.240000	3588	3528	17	0.480000
#TC	8587	8556	17	0.200000	3095	3084	7	0.230000
#TG	17659	16029	374	2.280000	3745	3337	108	3.130000
#TT	9464	9424	23	0.240000	3219	3206	8	0.250000
#Total number of matching rev reads: 274202
#Total number of matching rev dinucleotide and barcode combinations: 57308
#Total number of matching rev dinucleotide and barcode combinations (random): 57308
#dn	len	T	C	TtoC	len_r	T_r	C_r	TtoC_r
#AA	29273	29208	8	0.030000	3884	3875	4	0.100000
#AC	17593	17546	8	0.050000	3694	3685	1	0.030000
#AG	28874	28728	62	0.220000	3912	3894	8	0.210000
#AT	17492	17363	26	0.150000	3681	3656	4	0.110000
#CA	14282	14227	4	0.030000	3511	3495	1	0.030000
#CC	8393	8365	4	0.050000	3056	3046	0	0.000000
#CG	14383	14317	13	0.090000	3527	3506	1	0.030000
#CT	8711	8692	5	0.060000	3146	3137	2	0.060000
#GA	24776	24673	8	0.030000	3877	3856	2	0.050000
#GC	14216	14169	9	0.060000	3538	3525	0	0.000000
#GG	24540	24449	27	0.110000	3830	3813	3	0.080000
#GT	14500	14475	4	0.030000	3572	3565	0	0.000000
#TA	17420	17238	18	0.100000	3697	3649	4	0.110000
#TC	9592	9571	4	0.040000	3230	3221	2	0.060000
#TG	19174	17496	353	1.980000	3791	3392	74	2.140000
#TT	10983	10931	10	0.090000	3362	3346	6	0.180000
#
#
#171214_M00886_0185_000000000-BGP8V_F7-192-M20-ox-1_S1_L001_R1_001.fastq.gz
#Total number of reads: 750489
#Total number of matching fwd reads: 194742
#Total number of matching fwd dinucleotide and barcode combinations: 52454
#Total number of matching fwd dinucleotide and barcode combinations (random): 52454
#dn	len	T	C	TtoC	len_r	T_r	C_r	TtoC_r
#AA	20799	13526	5957	30.580000	3712	2388	1072	30.980000
#AC	11559	8825	2062	18.940000	3315	2558	563	18.040000
#AG	22223	15485	5047	24.580000	3781	2646	835	23.990000
#AT	11566	7332	3548	32.610000	3327	2112	999	32.110000
#CA	10959	6410	3464	35.080000	3275	1926	1005	34.290000
#CC	6103	4533	1065	19.020000	2684	1969	479	19.570000
#CG	11309	7388	2631	26.260000	3360	2186	801	26.820000
#CT	6212	3571	1950	35.320000	2724	1564	845	35.080000
#GA	18229	11481	4611	28.650000	3702	2364	901	27.600000
#GC	10388	8033	1704	17.500000	3235	2514	515	17.000000
#GG	17432	12627	3738	22.840000	3633	2635	768	22.570000
#GT	10842	6779	3218	32.190000	3350	2086	1008	32.580000
#TA	11294	7914	2725	25.610000	3283	2288	781	25.450000
#TC	6172	5369	643	10.700000	2708	2356	288	10.890000
#TG	12739	8978	2575	22.290000	3517	2493	677	21.360000
#TT	6916	4703	1912	28.900000	2848	1940	782	28.730000
#Total number of matching rev reads: 256559
#Total number of matching rev dinucleotide and barcode combinations: 56599
#Total number of matching rev dinucleotide and barcode combinations (random): 56599
#dn	len	T	C	TtoC	len_r	T_r	C_r	TtoC_r
#AA	26061	25999	8	0.030000	3865	3856	0	0.000000
#AC	16024	15972	6	0.040000	3582	3576	1	0.030000
#AG	26541	26419	48	0.180000	3867	3848	3	0.080000
#AT	16495	16399	14	0.090000	3624	3602	1	0.030000
#CA	13057	12987	5	0.040000	3480	3463	0	0.000000
#CC	7693	7670	5	0.070000	3008	2995	5	0.170000
#CG	13318	13264	15	0.110000	3492	3473	3	0.090000
#CT	8045	8025	6	0.070000	3077	3071	2	0.070000
#GA	23668	23583	12	0.050000	3802	3790	0	0.000000
#GC	13872	13819	9	0.070000	3505	3489	4	0.110000
#GG	23103	23029	19	0.080000	3813	3805	2	0.050000
#GT	14620	14599	6	0.040000	3543	3535	2	0.060000
#TA	16275	16082	17	0.110000	3671	3617	5	0.140000
#TC	9209	9175	8	0.090000	3167	3154	4	0.130000
#TG	18073	16413	410	2.440000	3749	3362	94	2.720000
#TT	10505	10468	7	0.070000	3354	3343	2	0.060000
#
#
#171214_M00886_0185_000000000-BGP8V_F7-192-M20-ox-2_S2_L001_R1_001.fastq.gz
#Total number of reads: 843667
#Total number of matching fwd reads: 222737
#Total number of matching fwd dinucleotide and barcode combinations: 54126
#Total number of matching fwd dinucleotide and barcode combinations (random): 54126
#dn	len	T	C	TtoC	len_r	T_r	C_r	TtoC_r
#AA	23430	15260	6598	30.190000	3778	2456	1069	30.330000
#AC	13075	10040	2291	18.580000	3397	2554	621	19.560000
#AG	25380	17693	5794	24.670000	3842	2657	899	25.280000
#AT	13005	8267	3923	32.180000	3415	2187	1008	31.550000
#CA	12636	7517	3815	33.670000	3419	2013	1045	34.170000
#CC	6912	5113	1232	19.420000	2824	2085	517	19.870000
#CG	12738	8327	2966	26.260000	3432	2248	775	25.640000
#CT	7066	4218	2110	33.340000	2883	1739	840	32.570000
#GA	21327	13549	5296	28.100000	3770	2417	946	28.130000
#GC	11879	9132	1942	17.540000	3379	2598	581	18.280000
#GG	20332	14651	4382	23.020000	3726	2687	787	22.650000
#GT	12681	7978	3765	32.060000	3447	2193	1025	31.850000
#TA	12958	9087	3128	25.610000	3412	2361	840	26.240000
#TC	6924	6057	721	10.640000	2849	2497	293	10.500000
#TG	14477	10365	2781	21.150000	3551	2511	700	21.800000
#TT	7917	5386	2152	28.550000	3002	2016	857	29.830000
#Total number of matching rev reads: 280034
#Total number of matching rev dinucleotide and barcode combinations: 57518
#Total number of matching rev dinucleotide and barcode combinations (random): 57518
#dn	len	T	C	TtoC	len_r	T_r	C_r	TtoC_r
#AA	28647	28587	10	0.030000	3901	3891	2	0.050000
#AC	16609	16558	3	0.020000	3631	3616	0	0.000000
#AG	29037	28888	59	0.200000	3904	3884	9	0.230000
#AT	18238	18117	25	0.140000	3711	3698	4	0.110000
#CA	14317	14250	7	0.050000	3518	3498	3	0.090000
#CC	8395	8363	3	0.040000	3036	3021	2	0.070000
#CG	14744	14687	12	0.080000	3612	3598	3	0.080000
#CT	8737	8717	5	0.060000	3121	3114	1	0.030000
#GA	25821	25700	13	0.050000	3839	3822	2	0.050000
#GC	15020	14966	4	0.030000	3581	3564	0	0.000000
#GG	25815	25721	19	0.070000	3861	3841	4	0.100000
#GT	16111	16085	10	0.060000	3635	3628	2	0.060000
#TA	17704	17480	21	0.120000	3721	3673	1	0.030000
#TC	9267	9231	4	0.040000	3225	3216	0	0.000000
#TG	20005	18305	371	1.990000	3802	3442	82	2.330000
#TT	11567	11525	8	0.070000	3420	3402	1	0.030000
#
#
#171214_M00886_0185_000000000-BGP8V_F7-192-M20-ox-3_S3_L001_R1_001.fastq.gz
#Total number of reads: 168499
#Total number of matching fwd reads: 36996
#Total number of matching fwd dinucleotide and barcode combinations: 23957
#Total number of matching fwd dinucleotide and barcode combinations (random): 23957
#dn	len	T	C	TtoC	len_r	T_r	C_r	TtoC_r
#AA	3956	1992	1634	45.060000	2103	1036	891	46.240000
#AC	2205	1485	532	26.380000	1475	1002	355	26.160000
#AG	4367	2519	1382	35.430000	2254	1301	705	35.140000
#AT	2197	1090	934	46.150000	1502	747	634	45.910000
#CA	2208	964	949	49.610000	1497	662	628	48.680000
#CC	1238	762	338	30.730000	977	595	269	31.130000
#CG	2233	1176	740	38.620000	1522	786	507	39.210000
#CT	1275	534	554	50.920000	1020	424	454	51.710000
#GA	3403	1621	1239	43.320000	1957	941	685	42.130000
#GC	2019	1335	497	27.130000	1425	951	352	27.010000
#GG	3211	1932	1013	34.400000	1894	1127	606	34.970000
#GT	2089	994	866	46.560000	1473	720	588	44.950000
#TA	2104	1200	743	38.240000	1483	846	529	38.470000
#TC	1054	864	151	14.880000	889	728	129	15.050000
#TG	2266	1412	606	30.030000	1530	938	422	31.030000
#TT	1171	634	467	42.420000	956	521	375	41.850000
#Total number of matching rev reads: 50187
#Total number of matching rev dinucleotide and barcode combinations: 29057
#Total number of matching rev dinucleotide and barcode combinations (random): 29057
#dn	len	T	C	TtoC	len_r	T_r	C_r	TtoC_r
#AA	5127	5116	0	0.000000	2446	2439	0	0.000000
#AC	3136	3126	0	0.000000	1878	1871	0	0.000000
#AG	5167	5130	17	0.330000	2485	2468	7	0.280000
#AT	3324	3305	6	0.180000	1916	1906	3	0.160000
#CA	2454	2439	2	0.080000	1611	1601	2	0.120000
#CC	1496	1487	0	0.000000	1150	1145	0	0.000000
#CG	2558	2544	5	0.200000	1650	1639	4	0.240000
#CT	1530	1523	3	0.200000	1177	1171	3	0.260000
#GA	4821	4797	0	0.000000	2334	2320	0	0.000000
#GC	2670	2662	0	0.000000	1709	1703	0	0.000000
#GG	4596	4576	5	0.110000	2242	2230	3	0.130000
#GT	2820	2818	1	0.040000	1770	1769	1	0.060000
#TA	3094	3055	5	0.160000	1857	1834	3	0.160000
#TC	1825	1822	0	0.000000	1336	1333	0	0.000000
#TG	3572	3241	68	2.060000	2083	1879	42	2.190000
#TT	1997	1984	0	0.000000	1413	1403	0	0.000000
```

The hashed data above is the output from the dinucleotide context analysis:

- Each block from the python printing statement above corresponds to one of the six files.
- Each block contains two tables, first fwd and second rev (where the rev dinucleotides have been reverse complemented).
  - Each table has two set of columns:
    - first four (len, T, C and TtoC) counting all Ts and Cs in 5hmU sites for all reads without selecting from reads that contain the same barcode and dinucleotide context
    - second four (len_r, T_r, C_r and TtoC_r) counting at 5hmU sites after selecting one read at random from the pool of reads containing the same barcode and dinucleotide context



