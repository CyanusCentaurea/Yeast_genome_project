# Analysis of the genetic changes effect on the manifestation of the prion-like factor [MCS +] in *Saccharomyces cerevisiae* strains
**Student:**   
Sidorenko Oksana, Bioinformatics Institute   

**Supervisors:**   
Lavrentii Danilov,    
Department of Genetics and Biotechnology, SPbU

Alexandr Rubel,    
Laboratory of Amyloid Biology, Department of Genetics and Biotechnology,
SPbU

# Table of contents
1. [Introduction](#introduction)
2. [Aim](#aim)
3. [Tasks](#tasks)
4. [Raw data](#raw_data)
5. [Workflow](#workflow)
6. [Versions of the programs used](#versions)
7. [Preparing the raw reds](#preparing_reads)
8. [Part 1: reference selection, alignment and variant calling](#part_1)
   * [Step 1: alignment – map to reference](#step_1)
   * [Step 2: mark duplicates + sort](#step_2)
   * [Step 3: collect alignment & insert size metrics](#step_3)
   * [Step 4: variant calling](#step_4)
   * [Step 5: extract SNPs & indels](#step_5)
   * [Step 6: filter SNPs](#step_6)
   * [Step 7: filter indels](#step_7)
   * [Step 8: exclude filtered variants](#step_8)
   * [Step 9: base Quality Score Recalibration (BQSR)](#step_9)
   * [Step 10: apply BQSR](#step_10)
   * [Step 11: variant calling](#step_11)
   * [Step 12: extract SNPs & indels](#step_12)
   * [Step 13: filter SNPs](#step_13)
   * [Step 14: filter indels](#step_14)
   * [Step 15: annotate SNPs and predict effects](#step_15)
   * [Overlaps between VCF files](#overlaps)
9. [Part 2: genome assembly](#part_2)
   * [Genome size estimation: raw reads](#genome_size_raw)
   * [Genome size estimation: corrected reads](#genome_size_cor)
   * [Annotation](#annotation)
   * [D-genies plotting](#d_genies)
   * [Characteristics of genes in the inverted region](#inversion_genes)
10. [Conclusion and further plans](#conclusion)
11. [Literature](#literature)

## Introduction <div id='introduction'/>
The increased interest in the study of amyloids is due to their
association with the development of a number of diseases in humans and animals, for example,
Alzheimer's and Parkinson's disease, type II diabetes, some types of cancer
diseases. Infectious amyloids (prions) capable of being transmitted both between cells of the same organism, and between
organisms are isolated into a separate group. In amyloid and prion diseases, extracellular or
intracellular accumulation of protein aggregates occurs. In addition to pathological
amyloids in various organisms, functional amyloids have been identified,
performing important biological functions. There are about 10 known
prions in the yeast *Saccharomyces cerevisiae*. In research led
Professor Yu. O. Chernoff in the yeast S. cerevisiae revealed a previously unknown
cytoplasmically inherited factor that exhibits characteristic properties
for prions. This prion-like factor was named [MCS + ]. Presence of this
factor in cells leads to a decrease in the accuracy of translation termination and
reading premature stop codons as significant. Structural protein of
the [MCS+] factor is currently unknown. At the same time, it is shown that its
manifestation in cells is associated with the presence in cells of unknown genetic
changes. That is, in strains that do not carry genetic changes, the factor
[MCS + ] is also contained, but has no manifestation. Thus, the identification
differences in the genomes of strains that differ in the manifestation of the factor [MCS + ],
will shed light on the nature of genetic factors influencing the manifestation of
[MCS + ], and may also contribute to the discovery of the nature of the [MCS + ] factor itself.   

## Aim: <div id='aim'/>
Find changes in the *S. cerevisiae* genome (including genomic mutations) that affect
manifestation of yeast prion-like factor [MCS + ].  

## Tasks: <div id='tasks'/>
1. Conduct data quality control;
2. Assemble the genomes of two strains;
3. Find SNPs unique to each strain;
4. Assess the presence of possible chromosomal rearrangements in the genomes of 2 strains.

## Raw data <div id='raw_data'/>
The available data at the start of the project were: four samples of Illumina paired end reads (two samples for each *S. cerevisiae* strain (with Rub_115 and Rub_117 prefixes, respectively)).

## Workflow <div id='workflow'/>

The project combined two major tasks:   
1. Reference selection, alignment and variant calling;
2. Genome assembly.

![img_1.png](Figures/img_1.png)

## Versions of the programs used <div id='versions'/>
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) v0.11.9
* [fastp](https://github.com/OpenGene/fastp?ysclid=lhgdfykpzu319786469) v. 0.20.1
* [GATK4](https://gatk.broadinstitute.org/hc/en-us) v4.1.3.0
* [BWA-MEM](https://bio-bwa.sourceforge.net/) v. 0.7.17-r1188
* [Picard Tools](https://broadinstitute.github.io/picard/) v. 2.20.5
* [R](https://www.r-project.org/) version 3.6.0 (2019-04-26) (variant-calling-pipeline-gatk4 [Docker container](https://hub.docker.com/r/gencorefacility/variant-calling-pipeline-gatk4))
* R version 4.2.1
* [Samtools](https://www.htslib.org/) v. 1.9 (variant-calling-pipeline-gatk4 Docker container)
* Samtools v. 1.16.1
* [SnpEff](https://pcingola.github.io/SnpEff/) v. 4.3i
* [SPAdes genome assembler](https://cab.spbu.ru/software/spades/) v3.15.4
* [Jellyfish mer counter](https://genome.umd.edu/jellyfish.html) v. 2.3.0
* [QUAST](https://cab.spbu.ru/software/quast/) v5.2.0
* [RepeatMasker](https://www.repeatmasker.org/) v. 4.1.5
* [ProteinOrtho](https://www.bioinf.uni-leipzig.de/Software/proteinortho/) v. 6.1.7 
* [BUSCO](https://busco.ezlab.org/) v. 5.4.4


## Preparing the raw reads <div id='preparing_reads'/>

The quality of raw paired reads was assessed in FastQC. About 2% of reads in each sample were trimmed by fastp (standard parameters) due to low quality and too many N.   

![img_2.png](Figures/img_2.png)

![img_3.png](Figures/img_3.png)

![img_4.png](Figures/img_4.png)

![img_5.png](Figures/img_5.png)

## Part 1: reference selection, alignment and variant calling <div id='part_1'/>

We used [variant calling pipeline](https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/) with GATK4 published by Mohammed Khalfan on 2020-03-25.    
Here we give an example of commands for one sample. All commands were executed similarly for all four samples.
![img_6.png](Figures/img_6.png)

### Step 1: alignment – map to reference <div id='step_1'/>

Tool: BWA-MEM   

Input: .fastq files, reference genome

Output: aligned_reads.sam 

Notes:   
-Y tells BWA to use soft clipping for supplementary alignments;    
-K tells BWA to process INT input bases in each batch regardless of nThreads (for reproducibility);     
-R Readgroup header line such as '@RG\tID:foo\tSM:bar'. This information is key for downstream GATK functionality. GATK will not work without a read group tag.    

Command example:     
`$ bwa mem -t 8 -K 100000000 -Y -R '@RG\tID:Rub115_ATTCAGAA-CCTATCCT_L001\tLB:Rub115_ATTCAGAA-CCTATCCT_L001\tPL:ILLUMINA\tPM:HISEQ\tSM:Rub115_ATTCAGAA-CCTATCCT_L001' GCA_014898935.1_ASM1489893v1_genomic.fna.gz fastp_Rub115_ATTCAGAA-CCTATCCT_L001_R1_001.fastq fastp_Rub115_ATTCAGAA-CCTATCCT_L001_R2_001.fastq > Rub115_ATTCAGAA-CCTATCCT_L001_alignment.sam`


### Step 2: mark duplicates + sort <div id='step_2'/>

Tool: GATK4 MarkDuplicatesSpark   

Input: aligned_reads.sam  

Output:   
sorted_dedup_reads.bam   
sorted_dedup_reads.bam.bai   
dedup_metrics.txt   

Notes:   
In GATK4, the Mark Duplicates and Sort Sam steps have been combined into one step using the MarkDuplicatesSpark tool. In addition, the BAM index file (.bai) is created as well by default. The tool is optimized to run on queryname-grouped alignments (that is, all reads with the same queryname are together in the input file). The output of BWA is query-grouped, however if provided coordinate-sorted alignments, the tool will spend additional time first queryname sorting the reads internally. Due to MarkDuplicatesSpark queryname-sorting coordinate-sorted inputs internally at the start, the tool produces identical results regardless of the input sort-order. That is, it will flag duplicates sets that include secondary, and supplementary and unmapped mate records no matter the sort-order of the input. This differs from how Picard MarkDuplicates behaves given the differently sorted inputs. (i.e. coordinate sorted vs queryname sorted). 

Command example:   
`$ gatk MarkDuplicatesSpark -I Rub115_ATTCAGAA-CCTATCCT_L001_alignment.sam -M Rub115_ATTCAGAA-CCTATCCT_L001_dedup_metrics.txt -O Rub115_ATTCAGAA-CCTATCCT_L001_sorted_dedup_reads.bam`

### Step 3: collect alignment & insert size metrics <div id='step_3'/>

Tool: Picard Tools, R, Samtools 

Input:   
sorted_dedup_reads.bam   
reference genome  

Output:   
alignment_metrics.txt     
insert_metrics.txt    
insert_size_histogram.pdf    
depth_out.txt   

Commands example:       
`$ java -jar apps/picard/2.17.11/picard-2.17.11.jar CollectAlignmentSummaryMetrics R=GCA_014898935.1_ASM1489893v1_genomic.fna I=Rub115_ATTCAGAA-CCTATCCT_L001_sorted_dedup_reads.bam O=Rub115_ATTCAGAA-CCTATCCT_L001_alignment_metrics.txt`

`$ java -jar apps/picard/2.17.11/picard-2.17.11.jar CollectInsertSizeMetrics INPUT=Rub115_ATTCAGAA-CCTATCCT_L001_sorted_dedup_reads.bam OUTPUT=Rub115_ATTCAGAA-CCTATCCT_L001_insert_metrics.txt HISTOGRAM_FILE=Rub115_ATTCAGAA-CCTATCCT_L001_insert_size_histogram.pdf`

`$ samtools depth -a Rub115_ATTCAGAA-CCTATCCT_L001_sorted_dedup_reads.bam > Rub115_ATTCAGAA-CCTATCCT_L001_depth_out.txt`

![img_7.png](Figures/img_7.png)

### Step 4: variant calling <div id='step_4'/>

Tool: GATK4   

Input:   
sorted_dedup_reads.bam   
reference genome   

Output:   
raw_variants.vcf   

Notes:   
First round of variant calling. The variants identified in this step will be filtered and provided as input for Base Quality Score Recalibration (BQSR)

Commands example:    
`$ gatk CreateSequenceDictionary -R GCA_014898935.1_ASM1489893v1_genomic.fna # creating the FASTA sequence dictionary file`

`$ gatk HaplotypeCaller -R GCA_014898935.1_ASM1489893v1_genomic.fna -I Rub115_ATTCAGAA-CCTATCCT_L001_sorted_dedup_reads.bam -O Rub115_ATTCAGAA-CCTATCCT_L001_raw_variants.vcf`

### Step 5: extract SNPs & indels <div id='step_5'/>

Tool: GATK4   

Input:    
raw_variants.vcf   
reference genome   

Output:    
raw_indels.vcf    
raw_snps.vcf   

Notes:   
This step separates SNPs and Indels so they can be processed and used independently

Commands example:   
`$ gatk SelectVariants -R GCA_014898935.1_ASM1489893v1_genomic.fna -V Rub115_ATTCAGAA-CCTATCCT_L001_raw_variants.vcf --select-type-to-include SNP -O Rub115_ATTCAGAA-CCTATCCT_L001_raw_snps.vcf` # separating SNPs

`$ gatk SelectVariants -R GCA_014898935.1_ASM1489893v1_genomic.fna -V Rub115_ATTCAGAA-CCTATCCT_L001_raw_variants.vcf --select-type-to-include INDEL -O Rub115_ATTCAGAA-CCTATCCT_L001_raw_indels.vcf` # separating indels

### Step 6: filter SNPs <div id='step_6'/>

Tool: GATK4   

Input:   
raw_snps.vcf   
reference genome   

Output:   
filtered_snps.vcf   
filtered_snps.vcf.idx   

Notes:
QD < 2.0: This is the variant confidence (from the QUAL field) divided by the unfiltered depth of non-hom-ref samples. This annotation is intended to normalize the variant quality in order to avoid inflation caused when there is deep coverage. For filtering purposes it is better to use QD than either QUAL or DP directly.

SNPs which are ‘filtered out’ at this step will remain in the filtered_snps.vcf file, however they will be marked as ‘_filter’, while SNPs which passed the filter will be marked as ‘PASS’. We need to extract and provide only the passing SNPs to the BQSR tool, we do this in step 9. 

Command example:   
`$ gatk VariantFiltration -R GCA_014898935.1_ASM1489893v1_genomic.fna -V Rub115_ATTCAGAA-CCTATCCT_L001_raw_snps.vcf -O Rub115_ATTCAGAA-CCTATCCT_L001_filtered_snps.vcf -filter-name "QD_filter" -filter "QD < 2.0"`

### Step 7: filter indels <div id='step_7'/>

Tool: GATK4   

Input:   
raw_indels.vcf   
reference genome  

Output:   
filtered_indels.vcf   
filtered_indels.vcf.idx   

Notes:   
QD < 2.0: This is the variant confidence (from the QUAL field) divided by the unfiltered depth of non-hom-ref samples. This annotation is intended to normalize the variant quality in order to avoid inflation caused when there is deep coverage. For filtering purposes it is better to use QD than either QUAL or DP directly.

Indels which are ‘filtered out’ at this step will remain in the filtered_indel.vcf file, however they will be marked as ‘_filter’, while Indels which passed the filter will be marked as ‘PASS’. We need to extract and provide only the passing Indels to the BQSR tool, we do this in step 9. 

Command example:   
`$ gatk VariantFiltration -R GCA_014898935.1_ASM1489893v1_genomic.fna -V Rub115_ATTCAGAA-CCTATCCT_L001_raw_indels.vcf -O Rub115_ATTCAGAA-CCTATCCT_L001_filtered_indels.vcf -filter-name "QD_filter" -filter "QD < 2.0"`

### Step 8: exclude filtered variants <div id='step_8'/>  

Tool: GATK4   

Input:   
filtered_snps.vcf   
filtered_indels.vcf   

Output:   
bqsr_snps.vcf   
bqsr_indels.vcf   

Notes:   
We need to extract only the passing variants and provide this as input to BQSR (next step). 

Commands example:   
`$ gatk SelectVariants --exclude-filtered -V Rub115_ATTCAGAA-CCTATCCT_L001_filtered_snps.vcf -O Rub115_ATTCAGAA-CCTATCCT_L001_bqsr_snps.vcf`

`$ gatk SelectVariants --exclude-filtered -V Rub115_ATTCAGAA-CCTATCCT_L001_filtered_indels.vcf -O Rub115_ATTCAGAA-CCTATCCT_L001_bqsr_indels.vcf`

### Step 9: base Quality Score Recalibration (BQSR) <div id='step_9'/>

Tool: GATK4   

Input:   
sorted_dedup_reads.bam (from step 2)   
bqsr_snps.vcf   
bqsr_indels.vcf   
reference genome   

Output:   
recal_data.table

Notes:   
BQSR is performed twice. The second pass is optional, only required to produce a recalibration report. We preferred to skip the second pass.

Command example:   
`$ gatk BaseRecalibrator -R GCA_014898935.1_ASM1489893v1_genomic.fna -I Rub115_ATTCAGAA-CCTATCCT_L001_sorted_dedup_reads.bam --known-sites Rub115_ATTCAGAA-CCTATCCT_L001_bqsr_snps.vcf --known-sites Rub115_ATTCAGAA-CCTATCCT_L001_bqsr_indels.vcf -O Rub115_ATTCAGAA-CCTATCCT_L001_recal_data.table`

### Step 10: apply BQSR <div id='step_10'/>

Tool: GATK4   

Input:   
recal_data.table   
sorted_dedup_reads.bam   
reference genome   

Output:   
recal_reads.bam   

Notes:    
This step applies the recalibration computed in the first BQSR step to the bam file. This recalibrated bam file is now analysis-ready. 

Command example:   
`$ gatk ApplyBQSR -R GCA_014898935.1_ASM1489893v1_genomic.fna -I Rub115_ATTCAGAA-CCTATCCT_L001_sorted_dedup_reads.bam -bqsr Rub115_ATTCAGAA-CCTATCCT_L001_recal_data.table -O Rub115_ATTCAGAA-CCTATCCT_L001_recal_reads.bam`

### Step 11: variant calling <div id='step_11'/>

Tool: GATK4   

Input:    
recal_reads.bam   
reference genome   

Output:    
raw_variants_recal.vcf    

Notes:    
Second round of variant calling performed using recalibrated (analysis-ready) bam

Command example:    
`$ gatk HaplotypeCaller -R GCA_014898935.1_ASM1489893v1_genomic.fna -I Rub115_ATTCAGAA-CCTATCCT_L001_recal_reads.bam -O Rub115_ATTCAGAA-CCTATCCT_L001_raw_variants_recal.vcf`

### Step 12: extract SNPs & indels <div id='step_12'/>

Tool: GATK4   

Input:   
raw_variants_recal.vcf   
reference genome   

Output:   
raw_indels_recal.vcf   
raw_snps_recal.vcf   

Notes:    
This step separates SNPs and Indels so they can be processed and analyzed independently

Commands example:    
`$ gatk SelectVariants -R GCA_014898935.1_ASM1489893v1_genomic.fna -V Rub115_ATTCAGAA-CCTATCCT_L001_raw_variants_recal.vcf -select-type-to-include SNP -O Rub115_ATTCAGAA-CCTATCCT_L001_raw_snps_recal.vcf` # separating SNPs

`$ gatk SelectVariants -R GCA_014898935.1_ASM1489893v1_genomic.fna -V Rub115_ATTCAGAA-CCTATCCT_L001_raw_variants_recal.vcf -select-type-to-include INDEL -O Rub115_ATTCAGAA-CCTATCCT_L001_raw_indels_recal.vcf` # separating indels

### Step 13: filter SNPs <div id='step_13'/>

Tool: GATK4   

Input:   
raw_snps_recal.vcf    
reference genome   

Output:  
filtered_snps_final.vcf  
filtered_snps_final.vcf.idx  

Notes:  
QD < 2.0: This is the variant confidence (from the QUAL field) divided by the unfiltered depth of non-hom-ref samples. This annotation is intended to normalize the variant quality in order to avoid inflation caused when there is deep coverage. For filtering purposes it is better to use QD than either QUAL or DP directly.

SNPs which are ‘filtered out’ at this step will remain in the filtered_snps.vcf file, however they will be marked as ‘_filter’, while SNPs which passed the filter will be marked as ‘PASS’.

Command example:   
`$ gatk VariantFiltration -R GCA_014898935.1_ASM1489893v1_genomic.fna -V Rub115_ATTCAGAA-CCTATCCT_L001_raw_snps_recal.vcf -O Rub115_ATTCAGAA-CCTATCCT_L001_filtered_snps_final.vcf -filter-name "QD_filter" -filter "QD < 2.0"`

### Step 14: filter indels <div id='step_14'/>

Tool: GATK4   

Input:  
raw_indels_recal.vcf  
reference genome  

Output:   
filtered_indels_final.vcf   
filtered_indels_final.vcf.idx   

Notes:  
QD < 2.0: This is the variant confidence (from the QUAL field) divided by the unfiltered depth of non-hom-ref samples. This annotation is intended to normalize the variant quality in order to avoid inflation caused when there is deep coverage. For filtering purposes it is better to use QD than either QUAL or DP directly.

Indels which are ‘filtered out’ at this step will remain in the filtered_indels.vcf file, however they will be marked as ‘_filter’, while Indels which passed the filter will be marked as ‘PASS’.

Command example:   
`$ gatk VariantFiltration -R GCA_014898935.1_ASM1489893v1_genomic.fna -V Rub115_ATTCAGAA-CCTATCCT_L001_raw_indels_recal.vcf -O Rub115_ATTCAGAA-CCTATCCT_L001_filtered_indels_final.vcf -filter-name "QD_filter" -filter "QD < 2.0"`

### Step 15: annotate SNPs and predict effects <div id='step_15'/>

Tool: SnpEff

Input:   
filtered_snps_final.vcf   

Output:   
filtered_snps_final.ann.vcf   
snpeff_summary.html   
snpEff_genes.txt   

Commands example:   
`$ echo "scerevisiae.genome : scerevisiae" >> snpEff.config` # config editing

`$ mkdir -p data/scerevisiae` # creating folder for the database

`$ cp GCF_000146045.2_R64_genomic.gtf ../new_gatk/data/scerevisiae/genes.gtf` # putting there the .gbk reference file (unzipped and renamed to genes.gtf)

`$ cp ../alignment/GCA_014898935.1_ASM1489893v1_genomic.fna ../../apps/snpeff/4_3i/data/scerevisiae/sequences.fa` # putting there the reference .fasta file (renamed to sequences.fa)

`$ java -jar ../apps/snpeff/4_3i/snpEff.jar build -gtf22 -v scerevisiae` # creating the database

`$ ll ../../apps/snpeff/4_3i/data/scerevisiae/`                          
genes.gtf   
sequences.fa   
snpEffectPredictor.bin   

`$ java -jar ../../apps/snpeff/4_3i/snpEff.jar ann scerevisiae Rub115_ATTCAGAA-CCTATCCT_L001_filtered_snps_final.vcf -v > Rub115_ATTCAGAA-CCTATCCT_L001_filtered_snps_final.ann.vcf` # running SnpEff

This step was also repeated using gatk CNNScoreVariants:   
`$ gatk CNNScoreVariants -V Rub115_ATTCAGAA-CCTATCCT_L001_raw_variants_recal.vcf -R GCA_014898935.1_ASM1489893v1_genomic.fna -O Rub115_ATTCAGAA-CCTATCCT_L001_CNN_annotated.vcf`

The results matched for these two tools: there were 12315, 12352, 12239 and 12279 SNPs detected in Rub_115_L001, Rub_115_L002, Rub_117_L001, Rub_117_L002 samples, respectively.

We went through all the above steps with the reference *S. cerevisiae* strain S288C, but found about 25000 SNPs for each sample, so we changed the reference to *S. cerevisiae* strain 74-D694.

![img_8.png](Figures/img_8.png)

### Overlaps between VCF files <div id='overlaps'/>

We were interested in identifying SNPs unique to the two strains. In other words, by what variants do the strains differ from each other?   
We used bcftools v. 1.13 `index` and `isec` commands to intersect VCF files with SNPs received after pipeline execution. There were 96 SNPs unique for the strain 115 and 95 SNPs unique for the strain 117 detected. Of these, 5 and 11 SNPs were located in genes, passed filters and were not synonymous for strain 115 and strain 117, respectively.

115 unique SNPs (in genes, passed filters and not synonymous):

| Chromosome            | Position | Annotated gene    | Effect            |
|-----------------------|----------|-------------------|-------------------|
| JADBMI010000006.1 (V) | 95498    | *YEL030W (ECM10)* | GTC>CTC (Val>Leu) |
| CM026509.1 (VI)       | 16639    | *YNL331C (AAD14)* | CTG>GTG (Leu>Val) |
| CM026510.1 (VII)      | 1131344  | *YGR295C (COS6)*  | CTG>ATG (Leu>Met) |
| CM026510.1 (VII)      | 1131352  | *YGR295C (COS6)*  | CTT>CCT (Leu>Pro) |
| CM026515.1 (XIII)     | 597989   | *YMR173W-A* (-)   | AAT>AAG (Asn>Lys) |

117 unique SNPs (in genes, passed filters and not synonymous):

| Chromosome                    | Position | Annotated gene    | Effect            |
|-------------------------------|----------|-------------------|-------------------|
| CM026506.1 (I)                | 18862    | *YAL063C (FLO9)*  | GAC>GTC (Asp>Val) |
| CM026506.1 (I)                | 18949    | *YAL063C (FLO9)*  | TTT>TCT (Phe>Ser) |
| CM026506.1 (I)                | 205450   | *YHR211W (FLO5)*  | ACC>ATC (Thr>Ile) |
| CM026507.1 (II)               | 803305   | *YBR301W (PAU24)* | GTC>GCC (Val>Ala) |
| CM026508.1 (IV)               | 2055     | *YDL248W (COS7)*  | GTC>ATC (Val>Ile) |
| CM026508.1 (IV)               | 2060     | *YDL248W (COS7)*  | TGG>TGA (Trp>*)   |   
| CM026516.1 (XIV)              | 8540     | *YNL332W (THI12)* | ATG>ATT (Met>Ile) |
| CM026516.1 (XIV)              | 93745    | *YDR210W-B*       | GGA>GAA (Gly>Glu) |
| CM026518.1 (XVI)              | 881755   | *YPR181C (SEC23)* | GAG>GAT (Glu>Asp) |
| CM026519.1 (Plasmid 2-micron) | 1244     | not annotated     | AGT>AAT (Ser>Asn) |  
| JADBMI010000028.1 (Plasmid?)  | 14425    | *YHL040C (ARN1)*  | AAT>GAT (Asn>Asp) |

Of the detected SNPs, the following ones could be distinguished:

* *YEL030W (ECM10)* (strain 115):    
**[*Saccharomyces* genome database (SGD)](https://www.yeastgenome.org/)** [[1]](#1) [[2]](#2):   
Heat shock protein of the Hsp70 family; localized in mitochondrial nucleoids, plays a role in protein translocation, interacts with Mge1p in an ATP-dependent manner; overexpression induces extensive mitochondrial DNA aggregations; ECM10 has a paralog, SSC1, that arose from the whole genome duplication.  
**[The National Center for Biotechnology Information (NCBI)](https://www.ncbi.nlm.nih.gov/)** [[3]](#3) [[4]](#4):   
Predicted to enable several functions, including ATP binding activity; ATP hydrolysis activity; and misfolded protein binding activity. Involved in protein targeting to mitochondrion. Located in mitochondrial nucleoid. Human ortholog(s) of this gene implicated in Parkinson's disease and autosomal dominant sideroblastic anemia 4. Orthologous to human HSPA9 (heat shock protein family A (Hsp70) member 9).

* *YMR173W-A* (strain 115):    
**SGD**:   
Dubious open reading frame; unlikely to encode a functional protein, based on available experimental and comparative sequence data; overlaps the verified gene DDR48/YML173W.   
*DDR48/YML173W*, **SGD**:   
DNA damage-responsive protein; expression is increased in response to heat-shock stress or treatments that produce DNA lesions; contains multiple repeats of the amino acid sequence NNNDSYGS; protein abundance increases in response to DNA replication stress.  
*DDR48/YML173W*, **NCBI**:  
Enables ATP hydrolysis activity and GTPase activity. Involved in DNA repair. Located in cytosol.

## Part 2: genome assembly <div id='part_2'/>

*De novo* genome assembly was performed with SPAdes.   

Command example:   
`$ SPAdes-3.15.4-Linux/bin/spades.py -1 Rub115_ATTCAGAA-CCTATCCT_L001_R1_001.fastq -2 Rub115_ATTCAGAA-CCTATCCT_L001_R2_001.fastq -o ./assembly/Rub115_ATTCAGAA-CCTATCCT_L001`

### Genome size estimation: raw reads <div id='genome_size_raw'/>

We used three approaches to estimate the genome size:  
1. With Jellyfish and R;  
2. With [Genomescope web tool](http://www.genomescope.org/);  
3. With formula.

The first approach involves creating a .histo file and k-mer profile plotting in accordance with [this tutorial](https://koke.asrc.kanazawa-u.ac.jp/HOWTO/kmer-genomesize.html).

First, we run jellyfish with parameters:   
`-m` or “mer” specifies the length (23)   
`-C` tells it to ignore directionality (it treats each read the same as its reverse complement)   
`-t` number of threads   
`-s` is an initial estimate for the size of the hash table jellyfish uses   
`-o` specifies the name of the output file, we chose a name with the k-mer length (23) in it.   

Command example:  
`$ jellyfish count -C -m 23 -t 10 -s 2G -o Rub115_ATTCAGAA-CCTATCCT_L001_23.jf Rub115_ATTCAGAA-CCTATCCT_L001_R1_001.fastq Rub115_ATTCAGAA-CCTATCCT_L001_R2_001.fastq`

Then we made tables for histograms.  

Command example:  
`$ jellyfish histo Rub115_ATTCAGAA-CCTATCCT_L001_23.jf > Rub115_ATTCAGAA-CCTATCCT_L001_23.histo`

After that we run the R script (can be found in Scripts/k_mer_profile_raw_reads.Rmd) on the .histo files obtained as a result of executing the previous commands (can be found in Data_for_scripts/Raw_reads_genome_size_plotting).

The second approach involves using the Genomescope web tool on these .histo files.

The third approach involves calculations by the formula:  
![img_9.png](Figures/img_9.png)  
where:  
*M* = k-mer peak;   
*K* = K-mer size;   
*L* = average read length;   
*T* = total bases  

M was obtained from the R script above, L and T were obtained as a result of FastQC analysis.  

### Genome size estimation: corrected reads <div id='genome_size_cor'/>
SPAdes works in two-step mode: error correction and assembly. It contains corrected reads in the “corrected” folder. We repeated the calculations and compared the results with the ones for uncorrected reads. 

R script can be found in Scripts/k_mer_profile_corrected_reads.Rmd and the .histo files can be found in Data_for_scripts/Corrected_reads_genome_size_plotting.

Here we present a comparison of the results obtained.

![img_10.png](Figures/img_10.png)

![img_11.png](Figures/img_11.png)

We see that, as a rule (but not in all cases), the calculations for corrected reads are slightly closer to the actual genome size. Also, for corrected reads, we give a QUAST estimate of the genome size.  

### Annotation <div id='annotation'/>

Repeat masking was performed with RepeatMasker.  
Command example:   

`$ RepeatMasker -species 'Saccharomyces cerevisiae' Rub115_ATTCAGAA-CCTATCCT_L001_cor_contigs.fasta`

Annotation was performed with [Augustus web interface](http://bioinf.uni-greifswald.de/webaugustus) (default parameters).

We obtained 5374, 5377, 5371, 5369 protein sequences for Rub_115_L001, Rub_115_L002, Rub_117_L001, Rub_117_L002 masked contigs predicted, respectively.

These sequences were mapped to the known S288C genes with ProteinOrtho and BUSCO. 

[Protein fasta](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_protein.faa.gz) was downloaded for mapping with ProteinOrtho clustering. Contains 6017 protein sequences.   
ProteinOrtho command example:   
`$ perl ./usr/local/bin/proteinortho6.pl Rub115_L001_cor_contigs_masked_predictions.aa GCA_000146045.2_R64_protein.faa`

We have found 5138 (about 85.4% of reference genome), 5143 (about 85.5% of reference genome), 5141 (about 85.4% of reference genome), 5145 (about 85.5% of reference genome) genes for Rub_115_L001, Rub_115_L002, Rub_117_L001, Rub_117_L002 masked contigs predicted, respectively.

BUSCO command example:   
`$ docker run -u 1000 -v busco:v5.4.4_cv1 busco -i ./busco_input -l saccharomycetes_odb10 -m proteins -o busco_output`  
parameters:  
-u for user id (obtained with `$ id -u`)  
-l for [lineages](https://busco.ezlab.org/list_of_lineages.html)  

![img_12.png](Figures/img_12.png)

where BUSCOs - Benchmarking Universal Single-Copy Orthologs.

Genome assembly: summary

| Sample                                                                                                        | Rub115_L001                  | Rub115_L001                  | Rub117_L001                  | Rub117_L002                  |   
|---------------------------------------------------------------------------------------------------------------|------------------------------|------------------------------|------------------------------|------------------------------|
| Number of contigs*                                                                                            | 263                          | 260                          | 257                          | 248                          |
| Largest contig, bases*                                                                                        | 514 399                      | 551 284                      | 551 369                      | 520 144                      |
| Total length, bases*                                                                                          | 11 509 482                   | 11 511 662                   | 11 503 821                   | 11 505 033                   |
| N50*                                                                                                          | 130 154                      | 133 828                      | 126 323                      | 117 401                      |
| Number of S288C genes found (%)**                                                                             | 5138 (85.4%)                 | 5143 (85.5%)                 | 5141 (85.4%)                 | 5145 (85.5%)                 |
| Number of BUSCOs (Benchmarking Universal Single-Copy Orthologs): complete BUSCOs: single-copy + duplicated*** | 2094 <br/>(2049 + 45), 98.0% | 2093 <br/>(2048 + 45), 97.9% | 2095 <br/>(2050 + 45), 98.0% | 2091 <br/>(2046 + 45), 97.8% |
\* Calculated with Quast.   
** Calculated as one-to-one orthologs with ProteinOrtho clustering.   
*** Calculated with BUSCO.   


### D-Genies plotting <div id='d_genies'/>
[D-Genies](https://dgenies.toulouse.inra.fr/) - an online tool designed to compare two genomes - was used for comparing reference strain 74-D694 and assemblies with dot plots.   

![img_15.png](Figures/img_15.png)

![img_16.png](Figures/img_16.png)

A long inversion was detected (chromosome CM026512.1 (IX), 148411 - 304998) in sample Rub_117_L001. A coverage plot was built to make sure there was no duplication or deletion.  
Coverage file was obtained with command:  
`samtools depth -a Rub117_GAATTCGT-CCTATCCT_L001.bam  > Rub117_GAATTCGT-CCTATCCT_L001.coverage`  

This file can be found in 117_001_coverage/Rub117_GAATTCGT-CCTATCCT_L001_recal_reads.zip. The plot was obtained by processing this file with the R script which can be found in Scripts/117_1_coverage_plot.Rmd

![img_14.png](Figures/img_14.png) 

There are no signs of changes in the copyness of chromosomes or their sections on the plot.

### Characteristics of genes in the inverted region <div id='inversion_genes'/>

Characteristics of genes in the inverted region:

| Locus tag   | Gene      | Coordinates (Chromosome: IX; CM026512.1) | Gene description                                                                                        |
|-------------|-----------|------------------------------------------|---------------------------------------------------------------------------------------------------------|
| *YIL107C*   | *PFK26*   | 148453-150936                            | 6-phosphofructo-2-kinase                                                                                |
| *YIL106W*   | *MOB1*    | 151590-152616                            | Mob1p                                                                                                   |
| *YIL105C*   | *SLM1*    | 152753-154819                            | phosphatidylinositol 4,5-bisphosphate-binding protein                                                   |
| *YIL104C*   | *SHQ1*    | 155162-156685                            | Hsp90 cochaperone SHQ1                                                                                  |
| *YIL103W*   | *DPH1*    | 156930-158207                            | 2-(3-amino-3-carboxypropyl)histidine synthase                                                           |
| *YIL101C*   | *XBP1*    | 160519-162465                            | Xbp1p                                                                                                   |
| *YIL099W*   | *SGA1*    | 163215-164864                            | glucan 1,4-alpha-glucosidase                                                                            |
| *YIL098C*   | *FMC1*    | 164981-165448                            | Fmc1p                                                                                                   |   
| *YIL097W*   | *FYV10*   | 165636-167186                            | glucose-induced degradation complex subunit FYV10                                                       |
| *YIL096C*   | *BMT5*    | 167326-168336                            | 25S rRNA (uracil2634-N3)-methyltransferase                                                              |
| *YIL095W*   | *PRK1*    | 169144-17157                             | serine/threonine protein kinase PRK1                                                                    |
| *YIL094C*   | *LYS12*   | 171724-172839                            | homoisocitrate dehydrogenase                                                                            |  
| *YIL093C*   | *RSM25*   | 173207-174001                            | mitochondrial 37S ribosomal protein RSM25                                                               |
| *YIL091C*   | *UTP25*   | 176247-178412                            | rRNA-binding ribosome biosynthesis protein UTP25                                                        |
| *YIL090W*   | *ICE2*    | 178810-180285                            | Ice2p                                                                                                   |
| *YIL088C*   | *AVT7*    | 182773-184245                            | Avt7p                                                                                                   |
| *YIL087C*   | *AIM19*   | 184489-184962                            | Aim19p                                                                                                  |
| *YIL085C*   | *KTR7*    | 185333-186886                            | putative mannosyltransferase                                                                            |
| *YIL084C*   | *SDS3*    | 187119-188102                            | Sds3p                                                                                                   |
| *YIL083C*   | *CAB2*    | 188398-189495                            | phosphopantothenate--cysteine ligase CAB2                                                               |
| *YIL079C*   | *AIR1*    | 190330-191412                            | TRAMP complex RNA-binding subunit                                                                       |
| *YIL078W*   | *THS1*    | 191905-194109                            | threonine--tRNA ligase THS1                                                                             |
| *YIL077C*   | *RCI37*   | 194396-195358                            | Rci37p                                                                                                  |
| *YIL076W*   | *SEC28*   | 196063-196953                            | coatomer subunit epsilon                                                                                |
| *YIL075C*   | *RPN2*    | 197267-200104                            | proteasome regulatory particle base subunit RPN2                                                        |
| *YIL074C*   | *SER33*   | 1200485-201894                           | phosphoglycerate dehydrogenase SER33                                                                    |
| *YIL073C*   | *SPO22*   | 202341-205358                            | Spo22p                                                                                                  |
| *YIL072W*   | *HOP1*    | 206006-207823                            | Hop1p                                                                                                   |
| *YIL071C*   | *PCI8*    | 208064-209098                            | Pci8p                                                                                                   |
| *YIL070C*   | *MAM33*   | 209676-210476                            | Mam33p                                                                                                  |
| *YER074W*   | *RPS24A*  | 210955-211772                            | ribosomal 40S subunit protein S24A                                                                      |
| *YIL068C*   | *SEC6*    | 212460-214877                            | SNARE-binding exocyst                                                                                   |
| *YIL067C*   | -         | 215127-217163                            | -                                                                                                       |
| *YIL066C*   | *RNR3*    | 217502-220159                            | ribonucleotide-diphosphate reductase subunit RNR3                                                       |
| *YIL065C*   | *FIS1*    | 220706-221173                            | Fis1p                                                                                                   |
| *YIL064W*   | *EFM4*    | 221365-222114                            | Efm4p                                                                                                   |
| *YIL063C*   | *YRB2*    | 222159-223142                            | Yrb2p                                                                                                   |
| *YIL062C*   | *ARC15*   | 223396-223860                            | Arc15p                                                                                                  |
| *YIL061C*   | *SNP1*    | 224055-224957                            | U1 snRNP complex subunit SNP1                                                                           |
| *YIL057C*   | *RGI2*    | 226965-227459                            | Rgi2p                                                                                                   |
| *YIL056W*   | *VHR1*    | 229057-230979                            | Vhr1p                                                                                                   |
| *YIL055C*   | -         | 231108-232991                            | -                                                                                                       |
| *YIL053W*   | *GPP1*    | 234181-234933                            | glycerol-1-phosphatase RHR2                                                                             |
| *YIL052C*   | *RPL34B*  | 235293-236131                            | ribosomal 60S subunit protein L34B                                                                      |
| *YIL051C*   | *MMF1*    | 236912-237349                            | isoleucine biosynthesis protein MMF1                                                                    |
| *YIL050W*   | *PCL7*    | 237982-238839                            | Pcl7p                                                                                                   |
| *YIL048W*   | *NEO1*    | 240502-243957                            | putative aminophospholipid-translocating P4-type ATPase NEO1                                            |
| *YIL047C*   | *YIL047C* | 244180-246888                            | Syg1p                                                                                                   |
| *YIL046W*   | *YIL046W* | 247717-249639                            | ubiquitin-binding SDF ubiquitin ligase complex subunit MET30                                            |
| *YIL045W*   | *PIG2*    | 250227-251843                            | putative protein phosphatase regulator PIG2                                                             |
| *YIL044C*   | *AGE2*    | 252015-252911                            | GTPase-activating protein AGE2                                                                          |
| *YIL043C*   | *CBR1*    | 253134-253988                            | cytochrome-b5 reductase                                                                                 |
| *YIL042C*   | *PKP1*    | 254170-255354                            | protein kinase PKP1                                                                                     |
| *YIL041W*   | *GVP36*   | 255587-256567                            | Gvp36p                                                                                                  |
| *YIL039W*   | *TED1*    | 257490-258911                            | Ted1p                                                                                                   |
| *YIL038C*   | *NOT3*    | 259207-264064                            | CCR4-NOT core subunit NOT3                                                                              |
| *YIL036W*   | *CST6*    | 264730-266496                            | Cst6p                                                                                                   |
| *YIL035C*   | *CKA1*    | 266857-267975                            | casein kinase 2 catalytic subunit CKA1                                                                  |
| *YIL034C*   | *CAP2*    | 268294-269157                            | F-actin-capping protein subunit beta                                                                    |
| *YIL033C*   | *BCY1*    | 269486-270736                            | cAMP-dependent protein kinase regulatory subunit BCY1                                                   |
| *YIL031W*   | *ULP2*    | 271700-274804                            | SUMO protease ULP2                                                                                      |
| *YIL030C*   | *SSM4*    | 275117-279076                            | E3 ubiquitin-protein ligase SSM4                                                                        |
| *YIL027C*   | *EMC5*    | 284366-284791                            | Emc5p                                                                                                   |
| *YIL026C*   | *IRR1*    | 285164-288616                            | Irr1p                                                                                                   |
| *YIL024C*   | -         | 289291-289860                            | uncharacterized protein                                                                                 |
| *YIL023C*   | *YKE4*    | 290075-291115                            | Zn(2+) transporter YKE4                                                                                 |
| *YIL022W*   | *TIM44*   | 291857-293152                            | protein translocase subunit TIM44                                                                       |
| *YIL021W*   | *RPB3*    | 293597-294553                            | DNA-directed RNA polymerase                                                                             |
| *YIL020C*   | *HIS6*    | 294726-295511                            | 1-(5-phosphoribosyl)-5- ((5-phosphoribosylamino)methylideneamino)imidazole-4-carboxamide isomerase HIS6 |
| *YIL019W*   | *FAF1*    | 295817-296824                            | Faf1p                                                                                                   |
| *YFR031C-A* | *RPL2A*   | 297265-298624                            | ribosomal 60S subunit protein L2A                                                                       |
| *YIL017C*   | *VID28*   | 298892-301657                            | glucose-induced degradation complex subunit VID28                                                       |
| *YIL016W*   | *SNL1*    | 302161-302625                            | Snl1p                                                                                                   |
| *YIL015W*   | *BAR1*    | 303034-304797                            | aspartyl protease BAR1                                                                                  |

Of these genes, the following ones could be distinguished:   

*YIL101C (XBP1)*:   
**NCBI**:   
Enables DNA-binding transcription factor activity, RNA polymerase II-specific. Involved in cellular response to oxidative stress; negative regulation of transcription by RNA polymerase II; and positive regulation of histone deacetylation. Located in nucleus.   
**SGD**:   
Transcriptional repressor; binds promoter sequences of cyclin genes, CYS3, and SMF2; not expressed during log phase of growth, but induced by stress or starvation during mitosis, and late in meiosis; represses 15% of all yeast genes as cells transition to quiescence; important for maintaining G1 arrest and for longevity of quiescent cells; member of Swi4p/Mbp1p family; phosphorylated by Cdc28p; relative distribution to nucleus increases upon DNA replication stress.

It was shown that prion [SWI] affects expression of FLO genes, variants of which we have identified in the strain 117 [[5]](#5).

*RPS24A*:   
**NCBI**:  
Predicted to be a structural constituent of ribosome. Involved in maturation of SSU-rRNA from tricistronic rRNA transcript (SSU-rRNA, 5.8S rRNA, LSU-rRNA). Located in mitochondrion. Human ortholog(s) of this gene implicated in Diamond-blackfan anemia 3. Orthologous to human RPS24 (ribosomal protein S24).   
**SGD**:  
Protein component of the small (40S) ribosomal subunit; homologous to mammalian ribosomal protein S24, no bacterial homolog; RPS24A has a paralog, RPS24B, that arose from the whole genome duplication.

*YIL023C (YKE4)*:    
**NCBI**:   
Enables zinc ion transmembrane transporter activity. Involved in zinc ion transport. Located in endoplasmic reticulum. Human ortholog(s) of this gene implicated in Ehlers-Danlos syndrome and Ehlers-Danlos syndrome spondylodysplastic type 3. Orthologous to several human genes including SLC39A7 (solute carrier family 39 member 7).  
**SGD**:   
Zinc transporter; localizes to the ER; null mutant is sensitive to calcofluor white, leads to zinc accumulation in cytosol; ortholog of the mouse KE4 and member of the ZIP (ZRT, IRT-like Protein) family.
It was shown that that members of the prion gene family represent an ancient branch of a larger family of ZIP metal ion transporters [[6]](#6)[[7]](#7).

*YFR031C-A (RPL2A)*:  
**NCBI**:   
Enables ATP binding activity; DNA binding activity; and chromatin binding activity. Contributes to ATP hydrolysis activity. Involved in chromosome organization and negative regulation of meiotic DNA double-strand break formation. Located in mitochondrion. Part of condensin complex. Orthologous to human SMC2 (structural maintenance of chromosomes 2).   
**SGD**:   
Ribosomal 60S subunit protein L2A; homologous to mammalian ribosomal protein L2 and bacterial L2; RPL2A has a paralog, RPL2B, that arose from the whole genome duplication.

*YIL016W (SNL1)*:   
**NCBI**:    
Enables ribosome binding activity. Involved in nuclear pore organization and protein folding. Located in endoplasmic reticulum membrane and nuclear envelope.  
**SGD**:   
Ribosome-associated protein; proposed to act in protein synthesis, nuclear pore complex biogenesis and maintenance as well as protein folding and prion maintenance; has similarity to the mammalian BAG-1 protein.  

Also, in connection with the revealed effect of the prion-like factor [MSC+] (presence of this factor in cells leads to a decrease in the accuracy of translation termination and reading premature stop codons as significant), we propose for additional studies all the genes associated with translation (*YER074W, YFR031C-A, YIL043C, YIL052C, YIL103W, YIL038C, YIL051C, YIL093C, YIL078W*).

### Conclusion and further plans <div id='conclusion'/>

Thus, in our work, we have identified a number of candidate genes for further testing. We can offer to study protein - protein interaction by the co-immunoprecipitation
method, DNA-protein interaction by the Chip-seq method or gene knockdown to verify our findings and for
further experiments. It is also necessary to study the identified SNPs, apply statistical approaches and obtain gene deletions to confirm or refute the involvement of these findings in the manifestation if [MCS + ].

If one of these candidates is highlighted as a result of additional research, it may make sense to analyze the differential expression and its structure for enrichment with certain amino acids and isolate a hypothetical prion domain, as well as check possible prion-prion interactions, since it is well-documented that a single yeast cell can harbor more than one prion element and that when two prions co-exist, they influence not only each other's de novo appearance but also propagation [[5]](#5).

## Literature <div id='literature'/>
1. Cherry, J. M., Hong, E. L., Amundsen, C., Balakrishnan, R., Binkley, G., Chan, E. T., Christie, K. R., Costanzo, M. C., Dwight, S. S., Engel, S. R., Fisk, D. G., Hirschman, J. E., Hitz, B. C., Karra, K., Krieger, C. J., Miyasato, S. R., Nash, R. S., Park, J., Skrzypek, M. S., … Wong, E. D. (2012). Saccharomyces Genome Database: the genomics resource of budding yeast. <i>Nucleic Acids Research</i>, <i>40</i>(Database issue). https://doi.org/10.1093/NAR/GKR1029 <div id='1'/> 
2. Cherry, J. M., Ball, C., Weng, S., Juvik, G., Schmidt, R., Adler, C., Dunn, B., Dwight, S., Riles, L., Mortimer, R. K., &#38; Botstein, D. (1997). Genetic and physical maps of Saccharomyces cerevisiae. <i>Nature</i>, <i>387</i>(6632 Suppl), 67. https://doi.org/10.1038/387s067 <div id='2'/> 
3. Sayers, E. W., Bolton, E. E., Brister, J. R., Canese, K., Chan, J., Comeau, D. C., Connor, R., Funk, K., Kelly, C., Kim, S., Madej, T., Marchler-Bauer, A., Lanczycki, C., Lathrop, S., Lu, Z., Thibaud-Nissen, F., Murphy, T., Phan, L., Skripchenko, Y., … Sherry, S. T. (2022). Database resources of the national center for biotechnology information. <i>Nucleic Acids Research</i>, <i>50</i>(D1), D20–D26. https://doi.org/10.1093/NAR/GKAB1112 <div id='3'/> 
4. Gene [Internet]. Bethesda (MD): National Library of Medicine (US), National Center for Biotechnology Information; 2004 – [cited 2023 12 05]. Available from: https://www.ncbi.nlm.nih.gov/gene/ <div id='4'/>
5. Goncharoff, D. K., Du, Z., &#38; Li, L. (2018). A brief overview of the Swi1 prion—[SWI+]. <i>FEMS Yeast Research</i>, <i>18</i>(6), 61. https://doi.org/10.1093/FEMSYR/FOY061 <div id='5'/>  
6. Ehsani, S., Mehrabian, M., Pocanschi, C. L., &#38; Schmitt-Ulms, G. (2012). The ZIP-prion connection. <i>Prion</i>, <i>6</i>(4), 317. https://doi.org/10.4161/PRI.20196 <div id='6'/>  
7. Schmitt-Ulms, G., Ehsani, S., Watts, J. C., Westaway, D., &#38; Wille, H. (2009). Evolutionary descent of prion genes from the ZIP family of metal ion transporters. <i>PloS One</i>, <i>4</i>(9). https://doi.org/10.1371/JOURNAL.PONE.0007208 <div id='7'/> 