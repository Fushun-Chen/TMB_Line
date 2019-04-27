#!/usr/bin/env Rscript

# load library
suppressPackageStartupMessages(library("optparse", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))
suppressPackageStartupMessages(library("data.table", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))
suppressPackageStartupMessages(library("mgcv", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))
##Specify desired options in a list

option_list <- list(
    make_option(c("-t","--tmb-file"), help="the sample tmb file"),
    make_option(c("-r","--response-file"), help="the response file"),
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

OR_list <- list()
i <- 1
for(tmb in TMB_Response_value$TMB){

	TMB_Response_value_low <- TMB_Response_value[TMB_Response_value$TMB <= tmb,]
	TMB_Response_value_high <- TMB_Response_value[TMB_Response_value$TMB > tmb,]
	tmb_high_response <- sum(TMB_Response_value_high$Response == 1, na.rm = T)
	tmb_high_noresponse <- sum(TMB_Response_value_high$Response == 0, na.rm = T)
	
	tmb_low_response <- sum(TMB_Response_value_low$Response == 1, na.rm = T)
	tmb_low_noresponse <- sum(TMB_Response_value_low$Response == 0, na.rm = T)

	OR <- c(tmb_high_response * tmb_low_noresponse) / c(tmb_low_response * tmb_high_noresponse)
	OR_list[[i]] <- OR
	i <- i + 1
}
OR_value <- do.call(rbind, OR_list)

TMB_Response_value$OR <- OR_value
TMB_Response_value$logOR <- log(OR_value)
TMB_Response_value_sel <- TMB_Response_value[!is.na(TMB_Response_value$logOR),]
TMB_Response_value_sel2 <- TMB_Response_value_sel[!is.infinite(TMB_Response_value_sel$logOR),]
write.table(TMB_Response_value_sel2, paste(OutputDir, "/", SampleName, ".GAM.TMB_ResponseNum.txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)

# run GAM 
TMB_logOR_GAM <- gam(logOR ~ s(TMB), data = TMB_Response_value_sel2)
png(paste(OutputDir, "/", SampleName, ".GAM.png", sep = ""))
p <- plot(TMB_logOR_GAM, ylab = "logOR")
print(p)
dev.off()

