#!/usr/bin/env Rscript

# load library
suppressPackageStartupMessages(library("optparse", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))
suppressPackageStartupMessages(library("data.table", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))

##Specify desired options in a list

option_list <- list(
    make_option(c("-m","--mut-file"), help="the mutation file"),
    make_option(c("-b","--bed-file"), help="the bed region"),
    make_option(c("-c","--cds-file"), help="the cds region"),
    make_option(c("-t","--muttype-file"), help="mutation type"),
    make_option(c("-s","--sample-name"), help="sample name"),
    make_option(c("-o","--output-dir"), help="the output dir, use to save the result")

)

# Get command line options
arguments <- parse_args(OptionParser(usage = "%prog [options]", option_list = option_list), positional_arguments = 0)
opt <- arguments$options

MutFile <- opt$`mut-file`
BedFile <- opt$`bed-file`
CdsFile <- opt$`cds-file`
MutTypeFile <- opt$`muttype-file`
SampleName <- opt$`sample-name`
OutputDir <- opt$`output-dir`
# def
TMB_calculate_fun <- function(mut_file, bed_file, cds_file, muttype_file, SampleName, output_dir){
	# read value
	mut_value <- fread(mut_file)
	bed_value <- fread(bed_file)
	cds_value <- fread(cds_file)
	muttype_value <- fread(muttype_file, header = F, sep = "\t")
	# bed expand
	bed_value$Id <- 1:nrow(bed_value)
	bed_expand_value <- bed_value[,.(V1 = V1, Position = V2:V3), by = Id]
	bed_expand_value$Chr_pos <- paste(bed_expand_value$V1, bed_expand_value$Position, sep = "_")
	# cds expand
	cds_value$Id <- 1:nrow(cds_value)
	cds_expand_value <- cds_value[,.(V1 = V1, Position = V2:V3), by = Id]
	cds_expand_value$Chr_pos <- paste(cds_expand_value$V1, cds_expand_value$Position, sep = "_")
	# bed select cds
	bed_cds_value <- bed_expand_value[bed_expand_value$Chr_pos %in% cds_expand_value$Chr_pos,]
	bed_cds_length <-  nrow(bed_cds_value)
	# mut select bed_cds
	mut_sel_value <- mut_value[paste(mut_value$Chromosome, mut_value$Start, sep = "_") %in% bed_cds_value$Chr_pos, ]
	# mut select muttype 
	for(i in 1:nrow(muttype_value)){
                muttype <- unlist(strsplit(as.character(muttype_value[i,]), ","))
		mut_sel_value1 <- mut_sel_value[mut_sel_value$Mutation_type %in% muttype,]
		mut_sel_value2 <- mut_sel_value1[mut_sel_value1$Function_region %in% muttype,]
		mut_sel_all_value <- rbind(mut_sel_value1, mut_sel_value2)
		SampleID_TMB <- mut_sel_all_value[,.(TMB=length(Chromosome) / bed_cds_length * 1000000), by = SampleID]
		write.table(SampleID_TMB, paste(output_dir, "/", SampleName, ".muttype", i, ".txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)
	}
}
##
#mut_file <- "/disk/gvcgroup/fushun/tmux_survive/cell_data/Data_mutlist.value.sort.header.ann.txt"
#bed_file <- "/disk/gvcgroup/fushun/tmux_survive/data/CDS.txt"
#cds_file <- "/disk/gvcgroup/fushun/tmux_survive/data/CDS.txt"
#muttype_file <- "/disk/gvcgroup/fushun/tmux_survive/cell_data/Data_mut_type.txt"
#output_dir <- "/disk/gvcgroup/fushun/tmux_survive/cell_data/TMB/test"
#TMB_calculate_fun(mut_file, bed_file, cds_file, muttype_file, output_dir)
TMB_calculate_fun(MutFile, BedFile, CdsFile, MutTypeFile, SampleName, OutputDir)

