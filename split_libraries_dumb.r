#!/usr/bin/Rscript

# Jack Darcy
# 24 AUG 2017
# Script to demultiplex illumina raw fastq files WITHOUT filtering
# I wrote this program because I want to demultiplex, THEN join my paired-end files. 
# That's impossible just using qiime - you have to join first, since split_libraries_fastq.py
# filters R1 and R2 indipendently. Then, you get two fastqs of different lengths, and now
# you can't use vsearch anymore if you want to join them together. You COULD use the fastx 
# toolkit to do this, but it's SLOW and hard to understand. This is FAST. WAY FAST. -Jack 

# usage: split_libraries_dumb.r --r1 r1.fastq --r2 r2.fastq -i index.fastq -m mappintgile.txt

suppressPackageStartupMessages(require(optparse))

option_list <- list(
	make_option(c("a", "--r1"), action="store", default=NA, type='character',
			help="R1 reads in fastq format"),
	make_option(c("b", "--r2"), action="store", default=NA, type='character',
			help="R2 reads in fastq format"),
	make_option(c("-i", "--index"), action="store", default=NA, type='character',
			help="Index reads in fastq format"),
	make_option(c("-m", "--map"), action="store", default=NA, type='character',
			help="QIIME format mapping file"),
	make_option(c("-o", "--output_directory"), action="store", default="", type='character',
			help="Directory where you want output fastq files to be. Default=pwd."),
	make_option("--skip", action="store", default="none", type='character',
			help="Used to skip 'first' or 'last' character of index read."),
	make_option("--prefix", action="store", default="filtered", type='character',
		help="Prefix for output fastq files. Will be appended to _r[12].fastq")
# help option -h/--help is included by optparse by default
)
opt = parse_args(OptionParser(option_list=option_list))


skip <- opt$skip
r1_fp <- opt$r1
r2_fp <- opt$r2
in_fp <- opt$index
mf_fp <- opt$map
outprefix <- opt$prefix
outdir <- opt$output_directory

# strip terminal "/" characters from outdir
# don't actually need to do this, but it's nice when things look right
while(substr(outdir, start=nchar(outdir), stop=nchar(outdir)) == "/"){
	outdir <- substr(outdir, start=1, stop=nchar(outdir) - 1)
}

# test that all input files have the right number of lines BEFORE
# reading them into memory. This is only tested on linux. 
print("Verifyinng file lengths.")
linecount <- function(fp){
	command <- paste("wc -l", fp, "| cut -f 1 -d ' '")
	return(as.numeric(system(command, intern=TRUE)))
}
counts <- c(linecount(r1_fp), linecount(r2_fp), linecount(in_fp))
if( ! all(counts == counts[1])){stop("CRITICAL ERROR: Input fastq files are of unequal length.")}
if( ! all(counts %% 4 == 0)){stop("CRITICAL ERROR: Invalid fastq files (length not multiple of 4)")}	


library(data.table)
library(parallel)

# read in index reads and mapping file
print("Reading index file.")
index <- fread(in_fp, header=F, stringsAsFactors=FALSE, sep='\t')

# get only barcode lines from index
print("Simplifying index file")
index <- index[(1:nrow(index) + 2) %% 4 == 0]
if(skip == "last"){
	startbp <- 1
	stopbp  <- nchar(index[[1]][1]) - 1
	index[[1]] <- substr(index[[1]], start=startbp, stop=stopbp)
}else if (skip == "first"){
	startbp <- 2
	stopbp <- nchar(index[[1]][1])
	index[[1]] <- substr(index[[1]], start=startbp, stop=stopbp)
}

# read in mapping file
print("Reading mapping file.")
map <- fread(mf_fp, header=T, stringsAsFactors=FALSE, sep='\t')

# compare index reads to barcodes
print("Comparing index reads to barcodes.")
good_seqs <- index[[1]] %in% map$BarcodeSequence

# function to translate barcode into sampleid
# good_index_df is just data.frame(good_seqs, index)[i]
# trans_df is just data.frame(map[[1]], map$BarcodeSequence)[i]
bc2sampid <- function(good_index, trans_df){
	if(good_index[1] == TRUE){
		samp <- trans_df[[2]] [ trans_df[[1]] == good_index[2] ]	
	}else{
		samp <- "NA"
	}
	return(samp)
}
sampids <- character(length(good_seqs))
for(i in 1:nrow(map)){
	bc_i <- map$BarcodeSequence[i]
	samp_i <- map[[1]] [i]
	sampids[which(index[[1]] == bc_i)] <- samp_i
}

rm(index)
sampids_good <- sampids[good_seqs]
good_lines <- rep(good_seqs, each=4)

# status message
print(paste("Found", sum(good_seqs), "hits out of", length(good_seqs), "total."))

# get ready for filtering
r1_outname <- paste(outdir, "/", outprefix, "_r1.fastq", sep="")
r2_outname <- paste(outdir, "/", outprefix, "_r2.fastq", sep="")

# filter r1
print("Reading in R1 fastq file.")
r1_fastq <- fread(r1_fp, header=F, stringsAsFactors=FALSE, sep='\t')
print("Filtering R1 fastq file.")
r1_fastq <- r1_fastq[[1]][good_lines]
r1_newnames <- paste(">", sampids_good, "_R1-", (1:length(sampids_good)), sep="")
r1_fastq[(1:length(r1_fastq) + 3) %% 4 == 0] <- r1_newnames
print("Writing out R1 fastq file.")
fwrite(list(r1_fastq), file=r1_outname, quote=F)
rm(r1_fastq)

# filter r2
print("Reading in R2 fastq file.")
r2_fastq <- fread(r2_fp, header=F, stringsAsFactors=FALSE, sep='\t')
print("Filtering r2 fastq file.")
r2_fastq <- r2_fastq[[1]][good_lines]
r2_newnames <- paste(">", sampids_good, "_r2-", (1:length(sampids_good)), sep="")
r2_fastq[(1:length(r2_fastq) + 3) %% 4 == 0] <- r2_newnames
print("Writing out r2 fastq file.")
fwrite(list(r2_fastq), file=r2_outname, quote=F)
rm(r2.fastq)

print("All done. Files written:")
print(r1_outname)
print(r2_outname)
