# Analysis of the genetic changes effect on the manifestation of the prion-like factor [MCS +] in *Saccharomyces cerevisiae* strains
**Student:**   
Sidorenko Oksana, Bioinformatics Institute   

**Supervisors:**   
Lavrentii Danilov,    
Department of Genetics and Biotechnology, SPbSU

Alexandr Rubel,    
Laboratory of Amyloid Biology, Department of Genetics and Biotechnology,
SPbSU
## Introduction 
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
will shed light on the nature of genetic factors influencing the manifestation
[MCS + ], and may also contribute to the discovery of the nature of the [MCS + ] factor itself.   

## Aim:
Find changes in the *S. cerevisiae* genome (including genomic mutations) that affect
manifestation of yeast prion-like factor [MCS + ].  

## Tasks:
1. Conduct data quality control;
2. Assemble the genomes of two strains;
3. Find SNPs unique to each strain;
4. Assess the presence of possible chromosomal rearrangements in the genomes of 2 strains.

## Raw data
The available data at the start of the project were: four samples of Illumina paired end reads (two samples for each *S. cerevisiae* strain (with Rub_115 and Rub_117 prefixes, respectively)).

## Workflow

The project combined two major tasks:   
1. Reference selection, alignment and variant calling;
2. Genome assembly.

![img_1.png](img_1.png)

## Preparing the raw reads

The quality of raw paired reads was assessed in [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). About 2% of reads in each sample were trimmed by [fastp](https://github.com/OpenGene/fastp?ysclid=lhgdfykpzu319786469) (v. 0.20.1, standard parameters) due to low quality and too many N.
![img_2.png](img_2.png)

![img_3.png](img_3.png)

![img_4.png](img_4.png)

![img_5.png](img_5.png)

## Part 1: reference selection, alignment and variant calling

We used [variant calling pipeline](https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/) with [GATK4](https://gatk.broadinstitute.org/hc/en-us) published by Mohammed Khalfan on 2020-03-25.    
Here we give an example of commands for one sample. All commands were executed similarly for all four samples.
![img_6.png](img_6.png)

### Step 1: alignment – map to reference

Tool: [BWA-MEM](https://bio-bwa.sourceforge.net/)    

Input: .fastq files, reference genome

Output: aligned_reads.sam 

Notes:   
-Y tells BWA to use soft clipping for supplementary alignments;    
-K tells BWA to process INT input bases in each batch regardless of nThreads (for reproducibility);     
-R Readgroup header line such as '@RG\tID:foo\tSM:bar'. This information is key for downstream GATK functionality. GATK will not work without a read group tag.    

Command example:     
`$ bwa mem -t 8 -K 100000000 -Y -R '@RG\tID:Rub115_ATTCAGAA-CCTATCCT_L001\tLB:Rub115_ATTCAGAA-CCTATCCT_L001\tPL:ILLUMINA\tPM:HISEQ\tSM:Rub115_ATTCAGAA-CCTATCCT_L001' GCA_014898935.1_ASM1489893v1_genomic.fna.gz fastp_Rub115_ATTCAGAA-CCTATCCT_L001_R1_001.fastq fastp_Rub115_ATTCAGAA-CCTATCCT_L001_R2_001.fastq > Rub115_ATTCAGAA-CCTATCCT_L001_alignment.sam`


### Step 2: mark duplicates + sort

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

### Step 3: collect alignment & insert size metrics

Tool: [Picard Tools](https://broadinstitute.github.io/picard/), [R](https://www.r-project.org/), [Samtools](https://www.htslib.org/)  

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

![img_7.png](img_7.png)

### Step 4: variant calling

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

### Step 5: extract SNPs & indels

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

### Step 6: filter SNPs

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

### Step 7: filter indels

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

### Step 8: exclude filtered variants   

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

### Step 9: base Quality Score Recalibration (BQSR)

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

### Step 10: apply BQSR

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

### Step 11: variant calling

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

### Step 12: extract SNPs & indels

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

### Step 13: filter SNPs

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

### Step 14: filter indels

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

### Step 15: annotate SNPs and predict effects

Tool: [SnpEff](https://pcingola.github.io/SnpEff/)

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

![img_8.png](img_8.png)

## Overlaps between VCF files

We were interested in identifying SNPs unique to the two strains. In other words, by what variants do the strains differ from each other?   
We used bcftools v. 1.13 `index` and `isec` commands to intersect VCF files with SNPs received after pipeline execution. There were 96 SNPs unique for the strain 115 and 95 SNPs unique for the strain 117 detected. Of these, 5 and 11 SNPs were located in genes, passed filters and were not synonymous for strain 115 and strain 117, respectively.