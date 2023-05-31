# Scripts information
## * **117_1_coverage_plot.Rmd**   
R script to obtain a coverage plot for the sample Rub_117_1 to make sure there was no duplication or deletion.
The coverage file Rub117_GAATTCGT-CCTATCCT_L001.coverage from the directory [Data_for_scripts/117_001_coverage](https://github.com/CyanusCentaurea/Yeast_genome_project/tree/main/Data_for_scripts/117_001_coverage) is required.    
To run the script, the file Rub117_GAATTCGT-CCTATCCT_L001.coverage must be pre-unzipped to the same directory where it is located (Data_for_scripts
/117_001_coverage/). This can be done using [unzip](https://www.linux.org/docs/man1/unzip.html) command which can be installed with:   
`$ sudo apt install unzip`

...and runned with   
`$ unzip Data_for_scripts/117_001_coverage/Rub117_GAATTCGT-CCTATCCT_L001_recal_reads.zip -d Data_for_scripts/117_001_coverage/`

## * **k_mer_profile_corrected_reads.Rmd**   
R script to estimate the genome size after SPAdes correction step (corrected reads).   
.histo files from the directory [Data_for_scripts/Corrected_reads_genome_size_plotting](https://github.com/CyanusCentaurea/Yeast_genome_project/tree/main/Data_for_scripts/Corrected_reads_genome_size_plotting) are required.

## * **k_mer_profile_raw_reads.Rmd**   
R script to estimate the genome size before SPAdes correction step (raw reads).   
.histo files from the directory [Data_for_scripts/Raw_reads_genome_size_plotting](https://github.com/CyanusCentaurea/Yeast_genome_project/tree/main/Data_for_scripts/Raw_reads_genome_size_plotting) are required.