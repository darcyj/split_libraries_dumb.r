#!/usr/bin/Rscript

# John L. Darcy
# 26 NOV 2019
# script to split barcoded illumina libraries by sample

# simple usage: split_libraries_dumb.r --r1 r1.fastq --r2 r2.fastq -i index.fastq -m mappingfile.txt

# thwow errors if packages are missing:
packages = c("optparse", "R.utils", "data.table")
for(p in packages){
	p_ok <- suppressWarnings(suppressPackageStartupMessages(require(p, character.only=T)))
	if(!p_ok){
		stop(paste0("Package \"", p, "\" is required but isn't installed."))
	}
}

# check if installed version of data.table is new enough for gzip support:
dt_vers <- as.character(packageVersion("data.table"))
if(compareVersion("1.12.6", dt_vers) == 1){
	stop("You are using an old version of data.table. Please update it to 1.12.6 or higher.")
}

option_list <- list(
	make_option("--r1", action="store", default=NA, type="character",
		help="R1 reads in fastq format (optional)"),
	make_option("--r2", action="store", default=NA, type="character",
		help="R2 reads in fastq format (optional)"),
	make_option("--r3", action="store", default=NA, type="character",
		help="R3 reads in fastq format (optional)"),
	make_option("--r4", action="store", default=NA, type="character",
		help="R4 reads in fastq format (optional)"),
	make_option(c("-i", "--index"), action="store", default=NA, type="character",
		help="Index reads in fastq format"),
	make_option(c("-m", "--map"), action="store", default=NA, type="character",
		help="metadata map file (must include sampleids and barcodes"),
	make_option("--allowed_mismatch", action="store", default=0, type="integer",
		help="Num allowed char mismatches between read and barcode. Default=0."),
	make_option("--skip", action="store", default="none", type="character",
		help="Used to skip 'first' or 'last' character of index reads. Default=none."),
	make_option("--rc_barcodes", action="store_true", default=FALSE, type="logical",
		help="Reverse-complements your barcodes before anything else."),
	make_option("--add_Cas1.8_data", action="store_true", default=FALSE, type="logical",
		help="Adds Casava 1.8 tags to sequence names (1:N:0:ATGATATGATGA)"),
	make_option("--split", action="store_true", default=FALSE, type="logical",
		help="Splits outputs by sample ID, creating separate files with sample ID prefixes."),
	make_option("--prefix", action="store", default="demuxed", type="character",
		help="Prefix for output fastq files. Will be appended to _r[12].fastq. Not used with --split. Default=demuxed."),
	make_option("--samp_col", action="store", default=1, type="integer",
		help="Col # of map containing sample IDs. Default=1."),
	make_option("--bc_col", action="store", default=0, type="integer",
		help="Col # of map containing barcodes. If 0, finds \"BarcodeSequence\" in column labels. Default=0."),
	make_option(c("-q", "--quiet"), action="store_true", default=FALSE, type="logical",
		help="Quiets output."),
	make_option(c("-s", "--summary"), action="store", default="none", type="character",
		help="Summary counts table filepath, not written if 'none'. Default=none.")
	# help option -h/--help is included by optparse by default
)
opt = parse_args(OptionParser(option_list=option_list))

# some stuff for debugging; make your own opt!
if(FALSE){
	opt <- list()
	opt$r1 = "../../rawfastq/miseq2014/Undetermined_S0_L001_R1_001.fastq.gz"\
	opt$r2 = "../../rawfastq/miseq2014/Undetermined_S0_L001_R2_001.fastq.gz"
	opt$r3 = NA
	opt$r4 = NA
	opt$index = "../../rawfastq/miseq2014/Undetermined_S0_L001_I1_001.fastq.gz"
	opt$map = "mapping_2014_small.txt"
	opt$allowed_mismatch = 0
	opt$skip = "last"
	opt$rc_barcodes <- TRUE
	opt$add_Cas1.8_data <- TRUE
	opt$split <- TRUE
	opt$prefix <- "demuxed"
	opt$samp_col <- 1
	opt$bc_col <- 0
	opt$quiet <- FALSE
	opt$summary <- "counts.txt"

}


## reverse-complement function
# simple function avoids dependency balogna
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

# message function wrapper for compatibility with quiet mode
msg <- function(x){ if( ! opt$quiet ){ message(x) } }

library(data.table)

fread_fq_gz <- function(fp){
	return(fread(fp, header=FALSE, sep=NULL)[[1]])
}

# read in index reads and mapping file
msg("Reading index file.")
index <- fread_fq_gz(opt$index)
nlines_index <- length(index)

# simplify index to include only nucleotide lines
msg("Simplifying index.")
index <- index[(1:nlines_index) %% 4 == 2 ]

# prune extra NTs from index reads if requested
if(opt$skip == "last"){
	startbp <- 1; stopbp  <- nchar(index[1]) - 1
	index <- substr(index, start=startbp, stop=stopbp)
}else if (opt$skip == "first"){
	startbp <- 2; stopbp <- nchar(index[1])
	index <- substr(index, start=startbp, stop=stopbp)
}

# read in and simplify mapping file
msg("Reading mapping file.")
map <- read.delim(opt$map, header=T, stringsAsFactors=FALSE, sep='\t', comment.char="")
if(opt$bc_col == 0 && "BarcodeSequence" %in% colnames(map)){
	bcs <- map$BarcodeSequence
}else if(opt$bc_col == 0){
	stop("Error: bc_col set to 0 but \"BarcodeSequence\" not in column names of map.")
}else{
	bcs <- map[,opt$bc_col]
}
map <- data.frame(
	sampleid=as.character(map[[opt$samp_col]]), #using [[col]] because it's a data.table
	barcode=bcs, stringsAsFactors=FALSE
)

# reverse-compliment barcodes if requested
if(opt$rc_barcodes){
	msg("Reverse-complementing barcodes within mapping file.")
	map$barcode <- sapply(X=map$barcode, FUN=rc)
}

# Figure out which reads belong to each barcode in map
msg("Comparing index reads to barcodes.")
msg(paste0("(Run with ",opt$allowed_mismatch, " allowed mismatches)"))
if(opt$allowed_mismatch <= 0){
	sampids <- character(length(index))
	for(i in 1:nrow(map)){
		sampids[which(index == map$barcode[i])] <- map$sampleid[i]
		# only output if i is divisible by 3, or if it's nrow(map)
		if(3 %% i == 0 | i == nrow(map)){
			msg(paste0("  ", i, " / ", nrow(map), " samples done"))
		}
	}
}else{
	# alternate comparison with wiggle room
	# function that compares string of the same length
	strdiff <- function(s1, s2){
		if(nchar(s1) != nchar(s2)){
			return(nchar(s1))
		}else{
			return( sum(strsplit(s1, split="")[[1]] != strsplit(s2, split="")[[1]]) )
		}
	}
	sampids <- character(length(index)) # stores sample destinations for each read
	hitsums <- integer(length(index))   # stores how many samples each read "hits"
	for(i in 1:nrow(map)){
		diffs_i <- sapply(X=index, FUN=strdiff, s2=map$barcode[i] )
		sampids[diffs_i <= opt$allowed_mismatch] <- map$sampleid[i]
		hitsums[diffs_i <= opt$allowed_mismatch] <- hitsums[diffs_i <= opt$allowed_mismatch] + 1
		rm(diffs_i)
		if(3 %% i == 0 | i == nrow(map)){
			msg(paste0("  ", i, " / ", nrow(map), " samples done"))
		}
	}
	# get rid of seqs that hit multiple barcodes
	n_mult_hit <- sum(hitsums > 1)
	sampids[hitsums > 1] <- ""
	msg(paste0(n_mult_hit, " seqs matched multiple barcodes"))
	msg("(Those hits were dicarded)")
	rm(hitsums)
}

# status msg
msg(paste("Found", sum(sampids != ""), "good hits out of", length(sampids), "total."))
# write summary file
if(opt$summary != "none"){
	msg("Tabulating summary counts table.")
	sum_out <- table(factor(sampids, levels=map$sampleid))
	# remove empty category (seqs that weren't in any given sample)
	sum_out <- sum_out[names(sum_out) != ""]
	sum_out <- data.frame(
		Sample = names(sum_out),
		n_seqs = as.numeric(sum_out)
	)
	msg("Writing counts table to file.")
	write.table(sum_out, file=opt$summary, quote=FALSE, sep="\t", row.names=FALSE)
}

# function to read in fastq, add sampleids, make names, and return a table of fq entries (rows)
	# each entry has 4 cols: 1=sample, 2=name, 3=seq, 4=plus, 5=qual
	# fqfp : fastq file path
	# samps : vector of sample ids, with "" meaning unassigned.
	# r : read number (1 or 2.... or 3?)
	# addcas : add casava data? T/F
	# ind : index reads (only needed if addcas=TRUE)
process_reads <- function(fqfp, samps, r, ind=NULL, addcas=F){
	# read in data
	fastq <- fread_fq_gz(fqfp)
	# check that length matches
	if(length(fastq) != length(samps) * 4){ stop("ERROR: r1 and index files have different numbers of lines.") }
	# fix names
	newnames <- paste0("@", samps, "_R", r, "-", (1:length(samps)))
	if(addcas){newnames <- paste0(newnames, " 1:N:2:", index)}
	# turn into table and return
	output <- cbind(
		samp=samps,
		name=newnames,
		seq=fastq[(1:length(fastq) + 2) %% 4 == 0],
		plus=rep("+", length(samps)),
		qual=fastq[(1:length(fastq) + 0) %% 4 == 0] 
	)
	# remove unassigned seqs
	output <- output[samps != "", ]
	return(output)
}

# writes reads processed with process_reads
	# readstable : output of process_reads
	# r : read number (1 or 2.... or 3?)
	# split : if true, write 1 file per sample
	# prefix : prefix for output file (if split is F)
write_reads <- function(readstable, r, split, prefix){
	# function that turns a readstable (from process_reads) into a character array of fq lines, in order.
	rtab2char <- function(rtab){ as.vector(t(rtab[,2:5])) } # this is hacky but it works.
	if(split == TRUE){
		samps <- unique(readstable[,1])
		for(s in samps){
			out_fp_s <- paste0(s, "_r", r, ".fastq.gz")
			rtab_s <- readstable[readstable[,1] == s, ,drop=F]
			fwrite(list(rtab2char(rtab_s)), file=out_fp_s, quote=F, compress="gzip")
			rm(out_fp_s, rtab_s)
		}
	}else{
		out_fp <- paste0(prefix, "_r", r, ".fastq.gz")
		fwrite(list(rtab2char(readstable)), file=out_fp, quote=F, compress="gzip")
	}
}

# process R1
if(! is.na(opt$r1)){
	msg(paste0("Processing R1 (", opt$r1, ")"))
	rt <- process_reads(fqfp=opt$r1, sampids, r=1, ind=index, addcas=opt$add_Cas1.8_data)
	write_reads(rt, r=1, split=opt$split, prefix=opt$prefix)
	rm(rt)
}

# process R2
if(! is.na(opt$r2)){
	msg(paste0("Processing R2 (", opt$r2, ")"))
	rt <- process_reads(fqfp=opt$r2, sampids, r=2, ind=index, addcas=opt$add_Cas1.8_data)
	write_reads(rt, r=2, split=opt$split, prefix=opt$prefix)
	rm(rt)
}

# process R3
if(! is.na(opt$r3)){
	msg(paste0("Processing R3 (", opt$r3, ")"))
	rt <- process_reads(fqfp=opt$r3, sampids, r=3, ind=index, addcas=opt$add_Cas1.8_data)
	write_reads(rt, r=3, split=opt$split, prefix=opt$prefix)
	rm(rt)
}

# process R4
if(! is.na(opt$r4)){
	msg(paste0("Processing R4 (", opt$r4, ")"))
	rt <- process_reads(fqfp=opt$r4, sampids, r=4, ind=index, addcas=opt$add_Cas1.8_data)
	write_reads(rt, r=4, split=opt$split, prefix=opt$prefix)
	rm(rt)
}

msg("All done.")
