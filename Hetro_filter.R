#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

Het_data = read.table(args[1],header = TRUE)
means = mean(Het_data$F)
SD = sd (Het_data$F) 
keep_samples = subset(Het_data, F <= means +3* SD & F >= means-3 * SD)

write.table(keep_samples[,c(1,2)], args[2], quote=F, row.names=F)
