#!/usr/bin/env Rscript

# load library
suppressPackageStartupMessages(library("optparse", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))
suppressPackageStartupMessages(library("data.table", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))
suppressPackageStartupMessages(library("pROC", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))
##Specify desired options in a list

option_list <- list(
    make_option(c("-t","--tmb-file"), help="the sample tmb file"),
    make_option(c("-r","--response-file"), help="the response file"),
    #make_option(c("-c","--choice-method"), help="choice method of dividing base line "),
    make_option(c("-s","--sample-name"), help="sample name"),
    make_option(c("-o","--output-dir"), help="the output dir, use to save the result")

)

# Get command line options
arguments <- parse_args(OptionParser(usage = "%prog [options]", option_list = option_list), positional_arguments = 0)
opt <- arguments$options

TmbFile <- opt$`tmb-file`
ResponseFile <- opt$`response-file`
SampleName <- opt$`sample-name`
OutputDir <- opt$`output-dir`
#
TMB_value <- fread(TmbFile)
colnames(TMB_value) <- c("SampleID", "TMB")
Response_value <- fread(ResponseFile)

TMB_Response_value <- merge(TMB_value, Response_value, by = "SampleID")
TMB_Response_Roc <- roc(TMB_Response_value$Response, TMB_Response_value$TMB)

png(paste(OutputDir, "/", SampleName, ".ROC.png", sep = ""))
p <- plot(TMB_Response_Roc, print.auc=TRUE,  print.thres=TRUE, legacy.axes = FALSE)
print(p)
dev.off()


TMB_threshold <- round(coords(TMB_Response_Roc, "best", ret=c("threshold")), 2)
TMB_value$TMB_Scope <- paste("<=", TMB_threshold, sep = "")
TMB_value$TMB_Scope[TMB_value$TMB > TMB_threshold] <- paste(">", TMB_threshold, sep = "") 

TMB_value$TMB_State <- "low"
TMB_value$TMB_State[TMB_value$TMB > TMB_threshold] <- "high"

write.table(TMB_value, paste(OutputDir, "/", SampleName, ".ROC.txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)

