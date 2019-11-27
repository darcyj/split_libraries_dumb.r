# split_libraries_dumb.r
Script to demultiplex illumina raw fastq files WITHOUT filtering

 by John L. Darcy
 26 NOV 2019
 
I wrote this program because I want to demultiplex, THEN join my paired-end fastq files. 
 That's impossible just using qiime 1 - you have to join first, since split_libraries_fastq.py
 filters R1 and R2 indipendently. Then, you get two fastqs of different lengths, and now
 you can't use vsearch anymore if you want to join them together. You COULD use the fastx 
 toolkit to do this, but it's SLOW and hard to understand. Qiime 2 has a solution for this
 problem, but if you're like me you don't like how Qiime 2 keeps all files as a baby-friendly
 "qiime zipped archive" object. 
 
 This is a no-nonsense demultiplexer that works only with EXACT BARCODE-INDEX matches, meaning
 there's no wiggle room if your index read is off by a NT. You don't want that anyway because 
 with illumina, the index read is the most high-quality part of a read, so if it's bad then your
 R1 and certainly R2 reads will be garbage for sure. 
 
 Required R packages: Rscript, data.table
 
 -JLD
 
## Note for Windows users:
 This definitely won't read gzipped files on windows. Should work every time on OS X and Linux,
 as long as you have gunzip installed. Actually, do they even have rscript for Windows???

## Usage: 
Easy mode:
```split_libraries_dumb.r --r1 r1.fastq.gz --r2 r2.fastq.gz -i index.fastq.gz -m mappintgfile.txt```

options (run with ```--help``` to see these):
* ```--r1```: R1 reads filepath in fastq format. If gzipped, must end with .gz. Optional.
* ```--r2```: see above. 
* ```--r3```: see above. 
* ```--r4```: see above. 
* ```-i / --index```: Index reads filepath in fastq format. If gzipped, must end in .gz. Required.
* ```-m / --map```: Metadata filepath in tab separated format. Required.
* ```--skip```: Used to skip "first" or "last" character of index reads. Default is "none".
* ```rc_barcodes```: Flag. If used, reverse-complements your barcodes before anything else.
* ```--add_Cas1.8_data```: Flag. If used, adds Cas1.8 data to output files. 
* ```split```: Flag. If used, separate output files will be written per sample per read. Files will be named sample_Rx.fastq.gz.
* ```prefix```: Prefix for output fastq files. Will be appended to \_Rx.fastq. Overridden by --split. Files will be named prefix_Rx.fastq.gz. Default is "demuxed".
* ```--bc_col```: Column number of metadata file containing barcodes. If 0 (default), finds "BarcodeSequence" in column labels.
* ```--samp_col```: Col number of metadata file containing sample IDs. Default=1.
* ```-h / --help```: Show the above, but explained more tersely.


