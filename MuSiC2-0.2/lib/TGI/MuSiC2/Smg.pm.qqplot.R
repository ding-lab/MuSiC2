####
# Generate qq plot for SMG results 
#
###
# for example 
#R --slave --args < qqplot.R music2_smg_test_detailed out.pdf
#
# Fetch command line arguments
args = commandArgs();
input = as.character(args[4])
output = as.character(args[5])

pdf( output )
# plot 1: 
#read.table( input, header = TRUE )[,11]->p
read.table( input, header = F )[,9]->p
p = p[p>0]
p = p[p<1]
p = p[!is.na(p)]
OBS = sort(-log10(p))
EXP = sort(-log10(1:length(p)/length(p)))
plot(EXP, OBS, col="red", pch=20); abline(a=0,b=1, col="lightgray", lty=1, lwd=2)
title("SMG test qq-plot")

dev.off()

