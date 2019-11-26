#!/usr/bin/Rscript

# Jack Darcy
# 25 JAN 2019
# script to filter QIIME format OTU maps

# usage: split_libraries_dumb.r --r1 r1.fastq --r2 r2.fastq -i index.fastq -m mappintgile.txt

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(parallel))

option_list <- list(
	make_option(c("a", "--r1"), action="store", default=NA, type='character',
			help="R1 reads in fastq format"),
	make_option(c("b", "--r2"), action="store", default=NA, type='character',
			help="R2 reads in fastq format"),
	make_option(c("-i", "--index"), action="store", default=NA, type='character',
			help="Index reads in fastq format"),
	make_option(c("-m", "--map"), action="store", default=NA, type='character',
			help="QIIME format mapping file"),
	make_option("--skip", action="store", default="none", type='character',
			help="Used to skip 'first' or 'last' character of index read."),
	make_option("--rc_barcodes", action="store_true", default=FALSE, type='logical',
		help="Reverse-complements your barcodes before anything else."),
	make_option("--add_Cas1.8_data", action="store_true", default=FALSE, type='logical',
		help="Adds Casava 1.8 tags to sequence names (1:N:0:ATGATATGATGA)"),
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
rc_barcodes <- opt$rc_barcodes
cas18names <- opt$add_Cas1.8_data

# functions to reverse complement a string

## reverse-complement function
# x is a string of IUPAC nucleotide characters, all upper case
rc <- function(x){
	# map for complementary nucleotides
	comp_map <- setNames(
		c("A", "T", "G", "C", "N", "R", "Y", "S", "W", "K", "M", "B", "V", "D", "H", "-"), 
		c("T", "A", "C", "G", "N", "Y", "R", "S", "W", "M", "K", "V" ,"B" ,"H", "D", "-"))
	# turn string x into character vector
	x <- unlist(strsplit(x, split=""))
	# reverse-complement x
	x_rc <- rev(comp_map[unlist(x)])
	# turn "NA"s into Ns
	x_rc[is.na(x_rc)] <- "N"
	# return x_rc as string
	return(paste(x_rc, collapse=""))
}

library(data.table)
library(parallel)

fread_fq_gz <- function(fp){
	if(endsWith(fp, "gz")){
	        in_cmd <- paste("gunzip -c", fp)
	        return(fread(cmd=in_cmd, header=FALSE, sep=NULL)[[1]])
	}else{
		return(fread(fp, header=FALSE, sep=NULL)[[1]])
	}
}

# read in index reads and mapping file
print("Reading index file.")
index <- fread_fq_gz(in_fp)
nlines_index <- length(index)

# get only barcode lines from index
print("Simplifying index")
index <- index[(1:nlines_index) %% 4 == 2 ]

if(skip == "last"){
	startbp <- 1; stopbp  <- nchar(index[1]) - 1
	index <- substr(index, start=startbp, stop=stopbp)
}else if (skip == "first"){
	startbp <- 2; stopbp <- nchar(index[1])
	index <- substr(index, start=startbp, stop=stopbp)
}

# read in mapping file
print("Reading mapping file.")
map <- fread(mf_fp, header=T, stringsAsFactors=FALSE, sep='\t')

if(rc_barcodes){
	print("Reverse-complementing barcodes within mapping file.")
	for(i in 1:nrow(map)){
		map$BarcodeSequence[i] <- rc(map$BarcodeSequence[i])	
	}
}

# compare index reads to barcodes
print("Comparing index reads to barcodes.")
good_seqs <- index %in% map$BarcodeSequence

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

# for each read, figure out its sampleID and its index sequence
sampids <- indices <- character(length(good_seqs))
for(i in 1:nrow(map)){
	bc_i <- map$BarcodeSequence[i]
	samp_i <- map[[1]] [i]
	sampids[which(index == bc_i)] <- samp_i
	indices[which(index == bc_i)] <- bc_i
}
rm(index)
sampids_good <- sampids[good_seqs]
indices_good <- indices[good_seqs]
good_lines <- rep(good_seqs, each=4)

# status message
print(paste("Found", sum(good_seqs), "hits out of", length(good_seqs), "total."))

# get ready for filtering
r1_outname <- paste(outprefix, "_r1.fastq", sep="")
r2_outname <- paste(outprefix, "_r2.fastq", sep="")

# filter r1
print("Reading in R1 fastq file.")
r1_fastq <- fread_fq_gz(r1_fp)
if(length(r1_fastq) != nlines_index){ stop("ERROR: r1 and index files have different numbers of lines.") }
print("Filtering R1 fastq file.")
r1_fastq <- r1_fastq[good_lines]
r1_newnames <- paste("@", sampids_good, "_R1-", (1:length(sampids_good)), sep="")
if(cas18names){r1_newnames <- paste(r1_newnames, " 1:N:2:", indices_good, sep="")}
r1_fastq[(1:length(r1_fastq) + 3) %% 4 == 0] <- r1_newnames
print("Writing out R1 fastq file.")
fwrite(list(r1_fastq), file=r1_outname, quote=F)
rm(r1_fastq)

# filter r2
print("Reading in R2 fastq file.")
r2_fastq <- fread_fq_gz(r2_fp)
if(length(r2_fastq) != nlines_index){ stop("ERROR: r2 and index files have different numbers of lines.") }
print("Filtering r2 fastq file.")
r2_fastq <- r2_fastq[good_lines]
r2_newnames <- paste("@", sampids_good, "_R2-", (1:length(sampids_good)), sep="")
if(cas18names){r2_newnames <- paste(r2_newnames, " 2:N:2:", indices_good, sep="")}
r2_fastq[(1:length(r2_fastq) + 3) %% 4 == 0] <- r2_newnames
print("Writing out r2 fastq file.")
fwrite(list(r2_fastq), file=r2_outname, quote=F)
rm(r2_fastq)
