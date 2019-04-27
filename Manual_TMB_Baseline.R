#!/usr/bin/env Rscript

# load library
suppressPackageStartupMessages(library("optparse", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))
suppressPackageStartupMessages(library("data.table", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))

##Specify desired options in a list

option_list <- list(
    make_option(c("-t","--tmb-file"), help="the sample tmb file"),
    make_option(c("-n","--cut-num"), help="the cut num"),
    make_option(c("-s","--sample-name"), help="sample name"),
    make_option(c("-o","--output-dir"), help="the output dir, use to save the result")

)

# Get command line options
arguments <- parse_args(OptionParser(usage = "%prog [options]", option_list = option_list), positional_arguments = 0)
opt <- arguments$options

TmbFile <- opt$`tmb-file`
CutNum <- as.numeric(opt$`cut-num`)
SampleName <- opt$`sample-name`
OutputDir <- opt$`output-dir`
#
TMB_value <- fread(TmbFile)
colnames(TMB_value) <- c("SampleID", "TMB")

if(CutNum == 2){
	tmb_median <- round(quantile(TMB_value$TMB, 0.5), 2)
	TMB_value$TMB_Scope <- 0
	TMB_value$TMB_Scope[TMB_value$TMB > tmb_median] <- paste(">", tmb_median, sep = "")
	TMB_value$TMB_Scope[TMB_value$TMB <= tmb_median] <- paste("<=", tmb_median, sep = "")
	TMB_value$TMB_State <- "low"
	TMB_value$TMB_State[TMB_value$TMB > tmb_median] <- "high"

}

if(CutNum == 3){
	tmb_low <- round(quantile(TMB_value$TMB, 1/3), 2)
	tmb_high <- round(quantile(TMB_value$TMB, 2/3), 2)
	
	TMB_value$TMB_Scope <- paste(tmb_low, "<=x<=", tmb_high, sep = "")
	TMB_value$TMB_Scope[TMB_value$TMB < tmb_low] <- paste("<", tmb_low, sep = "")
	TMB_value$TMB_Scope[TMB_value$TMB > tmb_high] <- paste(">", tmb_high, sep = "")	

	TMB_value$TMB_State <- "median"
	TMB_value$TMB_State[TMB_value$TMB > tmb_high] <- "high"
	TMB_value$TMB_State[TMB_value$TMB < tmb_low] <- "low"

}

if(CutNum == 4){
	q1 <- round(quantile(TMB_value$TMB, 1/4), 2)
	q2 <- round(quantile(TMB_value$TMB, 2/4), 2)
	q3 <- round(quantile(TMB_value$TMB, 3/4), 2)
	
	TMB_value$TMB_Scope <- paste("<", q1, sep = "")
	n2_logical <- TMB_value$TMB >= q1 & TMB_value$TMB < q2
	TMB_value$TMB_Scope[n2_logical] <- paste(q1, "<=x<", q2, sep = "")
	n3_logical <- TMB_value$TMB >= q2 & TMB_value$TMB < q3
	TMB_value$TMB_Scope[n3_logical] <- paste(q2, "<=x<=", q3, sep = "")
	TMB_value$TMB_Scope[TMB_value$TMB > q3] <- paste(">", q3, sep = "")

	TMB_value$TMB_State <- "q1"
	TMB_value$TMB_State[n2_logical] <- "q2"
	TMB_value$TMB_State[n3_logical] <- "q3"
	TMB_value$TMB_State[TMB_value$TMB > q3] <- "q4"

}

write.table(TMB_value, paste(OutputDir, "/", SampleName, ".CutNum.", CutNum, ".txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)


