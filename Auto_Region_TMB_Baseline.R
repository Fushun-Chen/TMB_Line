#!/usr/bin/env Rscript

# load library
suppressPackageStartupMessages(library("optparse", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))
suppressPackageStartupMessages(library("data.table", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))
##Specify desired options in a list

option_list <- list(
    make_option(c("-t","--tmb-file"), help="the sample tmb file"),
    make_option(c("-r","--response-file"), help="the response file"),
    make_option(c("-b","--start-end"), help="the region for start end,eg:50,55"),
    make_option(c("-l","--step-length"), help="increment of the sequence,eg:1"),
    make_option(c("-s","--sample-name"), help="sample name"),
    make_option(c("-o","--output-dir"), help="the output dir, use to save the result")

)

# Get command line options
arguments <- parse_args(OptionParser(usage = "%prog [options]", option_list = option_list), positional_arguments = 0)
opt <- arguments$options

TmbFile <- opt$`tmb-file`
ResponseFile <- opt$`response-file`
Region <- opt$`start-end`
Step <- as.numeric(opt$`step-length`)
SampleName <- opt$`sample-name`
OutputDir <- opt$`output-dir`
#
TMB_value <- fread(TmbFile)
colnames(TMB_value) <- c("SampleID", "TMB")
Response_value <- fread(ResponseFile)

TMB_Response_value <- merge(TMB_value, Response_value, by = "SampleID")

Start_End <- as.numeric(unlist(strsplit(Region, ",")))
Positions <- seq(Start_End[1], Start_End[2], Step)

pvalue_list <- list()
for(Position in Positions){
	TMB_Response_value_cut1 <- TMB_Response_value[TMB_Response_value$TMB < Position, ]
	TMB_Response_value_cut2 <- TMB_Response_value[TMB_Response_value$TMB >= Position, ]

	df_mat <- matrix(c(sum(TMB_Response_value_cut2$Response == 1, na.rm = T), sum(TMB_Response_value_cut1$Response == 1, na.rm = T), sum(TMB_Response_value_cut2$Response == 0, na.rm = T), sum(TMB_Response_value_cut1$Response == 0, na.rm = T)), nrow = 2) 
	pvalue <- fisher.test(df_mat)$p.value
	pvalue_list[[Position]] <- pvalue
}

df_pvalue <- data.frame(position = Positions, pvalue = unlist(pvalue_list))
cut_num <- df_pvalue$position[which.min(df_pvalue$pvalue)]

TMB_value$TMB_Scope <- paste("<=", cut_num, sep = "")
TMB_value$TMB_Scope[TMB_value$TMB > cut_num] <- paste(">", cut_num, sep = "") 
TMB_value$TMB_State <- "low"
TMB_value$TMB_State[TMB_value$TMB > cut_num] <- "high"

write.table(TMB_value, paste(OutputDir, "/", SampleName, ".regionCut.txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)

