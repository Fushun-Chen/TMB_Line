#!/usr/bin/env Rscript

# load library
suppressPackageStartupMessages(library("optparse", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))
suppressPackageStartupMessages(library("data.table", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))
suppressPackageStartupMessages(library("ggpubr", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))
suppressPackageStartupMessages(library("survival", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))
suppressPackageStartupMessages(library("survminer", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))

##Specify desired options in a list

option_list <- list(
    make_option(c("-t","--tmbState-file"), help="the sample tmb file, the file need have 4 columns(SampleID TMB TMB_Scope TMB_State)"),
    make_option(c("-r","--response-file"), help="the response file"),
    make_option(c("-s","--sample-name"), help="sample name"),
    make_option(c("-o","--output-dir"), help="the output dir, use to save the result")

)

# Get command line options
arguments <- parse_args(OptionParser(usage = "%prog [options]", option_list = option_list), positional_arguments = 0)
opt <- arguments$options

TmbFile <- opt$`tmbState-file`
ResponseFile <- opt$`response-file`
SampleName <- opt$`sample-name`
OutputDir <- opt$`output-dir`
#
TMB_value <- fread(TmbFile)
Response_value <- fread(ResponseFile)

TMB_Response_value <- merge(TMB_value, Response_value, by = "SampleID")

Low_TMB_Response_value <- TMB_Response_value[TMB_Response_value$TMB_State %in% c("low", "median", "q1", "q2", "q3"), ]
High_TMB_Response_value <- TMB_Response_value[TMB_Response_value$TMB_State %in% c("high", "q4"), ]


df_mat <- matrix(c(sum(High_TMB_Response_value$Response == 1, na.rm = T), sum(Low_TMB_Response_value$Response == 1, na.rm = T), sum(High_TMB_Response_value$Response == 0, na.rm = T), sum(Low_TMB_Response_value$Response == 0, na.rm = T)), nrow = 2) 
colnames(df_mat) <- c("Response", "NoResponse")
rownames(df_mat) <- c("High", "Low")
pvalue <- fisher.test(df_mat)$p.value
df_mat_value <- data.frame(df_mat)
df_mat_value$pvalue <- pvalue
write.table(df_mat_value, paste(OutputDir, "/", SampleName, ".ResponseNum.pvalue.txt", sep = ""), col.names = T, row.names = T, sep = "\t", quote = F)

Low_TMB_Response_value$TMB_State <- "TMB_Low"
High_TMB_Response_value$TMB_State <- "TMB_High"

All_TMB_Response_value <- rbind(High_TMB_Response_value, Low_TMB_Response_value)
All_TMB_Response_value$State[All_TMB_Response_value$State %in% c(2, "Dead")] <- 2
All_TMB_Response_value$State[All_TMB_Response_value$State %in% c(1, "Alive")] <- 1
All_TMB_Response_value$State <- as.numeric(All_TMB_Response_value$State)
#boxplot
png(paste(OutputDir, "/", SampleName, ".OS.boxplot.png", sep = ""))
p <- ggboxplot(All_TMB_Response_value, x = "TMB_State", y = "OS", color = "TMB_State", palette = "jco", add = "jitter") + stat_compare_means(label.x = 1.5)
print(p)
dev.off()
# survive
fit <- survfit(Surv(OS, State) ~ TMB_State, data = All_TMB_Response_value)
png(paste(OutputDir, "/", SampleName, ".TMB.survive.png", sep = ""))
p <- ggsurvplot(fit,conf.int=F,pval = TRUE, pval.coord = c(quantile(All_TMB_Response_value$OS, 0.95), 0.9),ggtheme=theme_minimal())
print(p)
dev.off()

