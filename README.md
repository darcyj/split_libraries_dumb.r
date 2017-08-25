# split_libraries_dumb.r
Script to demultiplex illumina raw fastq files WITHOUT filtering

 by Jack Darcy
 24 AUG 2017
 Script to demultiplex illumina raw fastq files WITHOUT filtering

 I wrote this program because I want to demultiplex, THEN join my paired-end fastq files. 
 That's impossible just using qiime - you have to join first, since split_libraries_fastq.py
 filters R1 and R2 indipendently. Then, you get two fastqs of different lengths, and now
 you can't use vsearch anymore if you want to join them together. You COULD use the fastx 
 toolkit to do this, but it's SLOW and hard to understand. This is FAST. WAY FAST.
 
 Required R packages: Rscript, data.table
 
 -Jack 

 usage: split_libraries_dumb.r --r1 r1.fastq --r2 r2.fastq -i index.fastq -m mappintgfile.txt
