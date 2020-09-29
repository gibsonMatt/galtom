# Galapagos ddRAD sequence pre-processing

Matt Gibson
Indiana University
2020


## Data pre-processing
* Data was generated using ddRAD at the IU Center for Genomics and Bioinformatics. The two enzymes used were EcoRI and PstI. Individuals were split over two runs of the NextSeq. Sequenced paired-end 300 (2x150 reads). Mid output. 

## Pre-process data

### Filter adapters, fix bases in read overlap
* Seems easiest to filter adapters prior to running process_radtags. `fastp` is run prior ro Stacks, but without trimming for quality.


```
#PBS -k oe #Keep the output and errors from the run
#PBS -l nodes=1:ppn=6,vmem=50gb,walltime=20:00:00
#PBS -M gibsomat@iu.edu
#PBS -m abe # Get notifications for the job
#PBS -N trim_adapted_pimp_poolA

cd /N/dc2/projects/gibsonTomato/galtom/data/rawdata/poolA


#Group 1
./N/dc2/projects/gibsonTomato/fastp -i GSF1973-PoolA-group1_S1_R1_001.fastq -I GSF1973-PoolA-group1_S1_R2_001.fastq
 -o ../../adapter_clean_rawdata/GSF1973-PoolA-group1_S1_R1_001_noadapters.fastq -O ../../adapter_clean_rawdata/GSF1973-PoolA-group1_S1_R2_001_noadapters.fastq-Q -L -w 6

#Group 2
./N/dc2/projects/gibsonTomato/fastp -i GSF1973-PoolA-group2_S2_R1_001.fastq -I GSF1973-PoolA-group2_S2_R2_001.fastq
 -o ../../adapter_clean_rawdata/GSF1973-PoolA-group2_S2_R1_001_noadapters.fastq -O ../../adapter_clean_rawdata/GSF1973-PoolA-group2_S2_R2_001_noadapters.fastq-Q -L -w 6

#Group 3
./N/dc2/projects/gibsonTomato/fastp -i GSF1973-PoolA-group3_S3_R1_001.fastq -I GSF1973-PoolA-group3_S3_R2_001.fastq
 -o ../../adapter_clean_rawdata/GSF1973-PoolA-group3_S3_R1_001_noadapters.fastq -O ../../adapter_clean_rawdata/GSF1973-PoolA-group3_S3_R2_001_noadapters.fastq-Q -L -w 6

#Group 4
./N/dc2/projects/gibsonTomato/fastp -i GSF1973-PoolA-group4_S4_R1_001.fastq -I GSF1973-PoolA-group4_S4_R2_001.fastq
 -o ../../adapter_clean_rawdata/GSF1973-PoolA-group4_S4_R1_001_noadapters.fastq -O ../../adapter_clean_rawdata/GSF1973-PoolA-group4_S4_R2_001_noadapters.fastq-Q -L -w 6

#Group 5
./N/dc2/projects/gibsonTomato/fastp -i GSF1973-PoolA-group5_S5_R1_001.fastq -I GSF1973-PoolA-group5_S5_R2_001.fastq
 -o ../../adapter_clean_rawdata/GSF1973-PoolA-group5_S5_R1_001_noadapters.fastq -O ../../adapter_clean_rawdata/GSF1973-PoolA-group5_S5_R2_001_noadapters.fastq-Q -L -w 6


#Group6
./N/dc2/projects/gibsonTomato/fastp -i GSF1973-PoolsAB-group6_S6_R1_001.fastq -I GSF1973-PoolslAB-group6_S6_R2_001.fastq
 -o ../../adapter_clean_rawdata/GSF1973-PoolsAB-group6_S6_R1_001_noadapters.fastq -O ../../adapter_clean_rawdata/GSF1973-PoolsAB-group6_S6_R2_001_noadapters.fastq-Q -L -w 6

#Same done for poolB
```



### Fix RE site with `recoverRE.py`

Base quality at restriction sites is typically terrible due to low sequence diversity. We can correc these to improve mapability of the reads.  `recoverRE.py` corrects the restriction site sequences, allowing only one mismatch between the true/expected overhangs (TGCAG and AATTC).

```
python recoverRE.py --truesite TGCAG --startpos 6 --endpos 10 --out /N/dc2/projects/gibsonTomato/galtom/data/recoverREsite/GSF1973-PoolA-group1_S1_R1_001_noadapters_refix.fastq /N/dc2/projects/gibsonTomato/galtom/data/adapter_clean_rawdata/GSF1973-PoolA-group1_S1_R1_001_noadapters.fastq &

python recoverRE.py --truesite TGCAG --startpos 6 --endpos 10 --out /N/dc2/projects/gibsonTomato/galtom/data/recoverREsite/GSF1973-PoolA-group2_S2_R1_001_noadapters_refix.fastq /N/dc2/projects/gibsonTomato/galtom/data/adapter_clean_rawdata/GSF1973-PoolA-group2_S2_R1_001_noadapters.fastq &

python recoverRE.py --truesite TGCAG --startpos 6 --endpos 10 --out /N/dc2/projects/gibsonTomato/galtom/data/recoverREsite/GSF1973-PoolA-group3_S3_R1_001_noadapters_refix.fastq /N/dc2/projects/gibsonTomato/galtom/data/adapter_clean_rawdata/GSF1973-PoolA-group3_S3_R1_001_noadapters.fastq &

python recoverRE.py --truesite TGCAG --startpos 6 --endpos 10 --out /N/dc2/projects/gibsonTomato/galtom/data/recoverREsite/GSF1973-PoolA-group4_S4_R1_001_noadapters_refix.fastq /N/dc2/projects/gibsonTomato/galtom/data/adapter_clean_rawdata/GSF1973-PoolA-group4_S4_R1_001_noadapters.fastq &

python recoverRE.py --truesite TGCAG --startpos 6 --endpos 10 --out /N/dc2/projects/gibsonTomato/galtom/data/recoverREsite/GSF1973-PoolA-group5_S5_R1_001_noadapters_refix.fastq /N/dc2/projects/gibsonTomato/galtom/data/adapter_clean_rawdata/GSF1973-PoolA-group5_S5_R1_001_noadapters.fastq &

python recoverRE.py --truesite TGCAG --startpos 6 --endpos 10 --out /N/dc2/projects/gibsonTomato/galtom/data/recoverREsite/GSF1973-PoolAB-group6_S6_R1_001_noadapters_refix.fastq /N/dc2/projects/gibsonTomato/galtom/data/adapter_clean_rawdata/GSF1973-PoolAB-group6_S6_R1_001_noadapters.fastq &



python recoverRE.py --truesite AATTC --startpos 1 --endpos 5 --out /N/dc2/projects/gibsonTomato/galtom/data/recoverREsite/PoolA/GSF1973-PoolA-group1_S1_R2_001_noadapters_refix.fastq /N/dc2/projects/gibsonTomato/galtom/data/adapter_clean_rawdata/poolB/GSF1973-PoolA-group1_S1_R2_001_noadapters.fastq &

python recoverRE.py --truesite AATTC --startpos 1 --endpos 5 --out /N/dc2/projects/gibsonTomato/galtom/data/recoverREsite/PoolA/GSF1973-PoolA-group2_S2_R2_001_noadapters_refix.fastq /N/dc2/projects/gibsonTomato/galtom/data/adapter_clean_rawdata/poolB/GSF1973-PoolA-group2_S2_R2_001_noadapters.fastq &

python recoverRE.py --truesite AATTC --startpos 1 --endpos 5 --out /N/dc2/projects/gibsonTomato/galtom/data/recoverREsite/PoolA/GSF1973-PoolA-group3_S3_R2_001_noadapters_refix.fastq /N/dc2/projects/gibsonTomato/galtom/data/adapter_clean_rawdata/poolB/GSF1973-PoolA-group3_S3_R2_001_noadapters.fastq &

python recoverRE.py --truesite AATTC --startpos 1 --endpos 5 --out /N/dc2/projects/gibsonTomato/galtom/data/recoverREsite/PoolA/GSF1973-PoolA-group4_S4_R2_001_noadapters_refix.fastq /N/dc2/projects/gibsonTomato/galtom/data/adapter_clean_rawdata/poolB/GSF1973-PoolA-group4_S4_R2_001_noadapters.fastq &

python recoverRE.py --truesite AATTC --startpos 1 --endpos 5 --out /N/dc2/projects/gibsonTomato/galtom/data/recoverREsite/PoolA/GSF1973-PoolA-group5_S5_R2_001_noadapters_refix.fastq /N/dc2/projects/gibsonTomato/galtom/data/adapter_clean_rawdata/poolB/GSF1973-PoolA-group5_S5_R2_001_noadapters.fastq &

python recoverRE.py --truesite AATTC --startpos 1 --endpos 5 --out /N/dc2/projects/gibsonTomato/galtom/data/recoverREsite/PoolA/GSF1973-PoolAB-group6_S6_R2_001_noadapters_refix.fastq /N/dc2/projects/gibsonTomato/galtom/data/adapter_clean_rawdata/poolB/GSF1973-PoolAB-group6_S6_R2_001_noadapters.fastq &
wait
```
- Same for poolB



## Demultiplex with `process_radtags`

### PoolA

```
process_radtags -P -p /N/dc2/projects/gibsonTomato/galtom/data/recoverREsite/poolA/ -b /N/dc2/projects/gibsonTomato/galtom/data/barcodes/poolA/poolA_stacks_barcodes.txt -o ./demultiplexed/poolA/ -r --inline_index --renz_1 pstI --renz_2 ecoRI -y gzfastq -s 5 --len_limit 30
```


### PoolB
```
process_radtags -P -p /N/dc2/projects/gibsonTomato/galtom/data/recoverREsite/poolB/ -b /N/dc2/projects/gibsonTomato/galtom/data/barcodes/poolB/poolB_stacks_barcodes.txt -o ./demultiplexed/poolB/ -r --inline_index --renz_1 pstI --renz_2 ecoRI -y gzfastq -s 5 --len_limit 30
```


## Merge individuals split over the two pools

```
rsync -a -P *.gz data/demultiplexed/
```



## Map with `bwa`


```
for r1 in *.1.fq.gz
do
	echo "Mapping..."
	echo "$r1"
	toremove='.1.'
	newsubstr='.2.'
	r2="${r1/$toremove/$newsubstr}"
	
	echo "$r2"

	outname="/N/dc2/projects/gibsonTomato/galtom/stacks/map/mapped_reads/${r1::-8}.sam"
	echo "$outname"

	bwa mem -t 12 /N/dc2/projects/gibsonTomato/pimpGEA/data/genome/S_lycopersicum_chromosomes.3.00.fa $r1 $r2 > $outname &

done
wait

```


## Filter Bams

We filter bams based on mapping quality (MQ > 3), no supplementary reads, and proper paired-end sequencing
```
#!/usr/bin/env bash

cd $1

for f in *.sam
do
	echo "$f" 1>&2
	outname="/N/dc2/projects/gibsonTomato/galtom/stacks/filter_mapped/filtered_bams/${f::-4}.filtered.bam"
	
	samtools view -f 2 -q 4 -h -b -o $outname $f &
done
wait

cd /N/dc2/projects/gibsonTomato/galtom/stacks/filter_mapped/filtered_bams/

for f in *.filtered.bam
do

	outname="/N/dc2/projects/gibsonTomato/galtom/stacks/filter_mapped/filtered_bams/${f::-4}.filtered.bam"
	samtools flagstat $f > "$outname.stats" &
done
wait

```


## Sort bams

```
for f in *.bam
do
	echo "$f"
	outname="/N/dc2/projects/gibsonTomato/galtom/stacks/sort_bams/sorted_bams/${f::-4}.sorted.bam"
	
	samtools sort $f > $outname
done
```


## Stack Reference Mapping Pipeline

### `ref_map.pl`

```
/N/dc2/projects/gibsonTomato/bin/ref_map.pl -T 6 --samples /N/dc2/projects/gibsonTomato/pimpGEA/STACKS/sorted_bams/ --popmap /N/dc2/projects/gibsonTomato/pimpGEA/STACKS/gstacks/popmap.txt -o /N/dc2/projects/gibsonTomato/pimpGEA/STACKS/gstacks/outfiles/
```

# VCF filtering

## Remove .filtered.sorted from indv names
```
sed -i '' 's/.filtered.sorted//g' populations.snps.vcf
```

## Initial stats

Total # records: 150,128



## Filter 1: Remove chromosome 0 (unassembled regions)
* See local `./associations/data/vcfs/stacks/populations.snps.filter1.vcf.stats`

Total # snps: 146,287

```
grep -Ev 'SL3.0ch00' populations.snps.sorted.vcf > populations.snps.filter1.vcf
```

## Filter 2: For genotype allele depth greater than 8

Total # of snps: 

Doesn't remove sites, but changes genotypes to missing if DP < 8

```
vcftools --vcf populations.snps.filter1.vcf --recode --recode-INFO-all --out populations.snps.filter2 --minDP 8

vcftools --vcf populations.snps.filter2.recode.vcf --recode --recode-INFO-all --out populations.snps.filter3 --max-missing-count 88
```

* Max-missing count: 35 (80%) ->  sites
* Max-missing count: 88 (50%) -> 6282 sites

## Filter 4: Remove MG120-35 because I dont know what it is

```
vcftools --vcf populations.snps.filter3.recode.vcf --remove-indv MG120-35 --recode --recode-INFO-all --out populations.snps.filter4.recode
```

## Filter 5: Remove SNPs in full LD

```
bcftools +prune -l 0.9 -w 1000 populations.snps.filter4.recode.recode.vcf -Ov -o populations.snps.filter5.vcf
```