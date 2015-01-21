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

#input="smgs_detailed"
#output="qqplot.pdf"

read.table(input,header=T )->z

gc=function(p)
{
p=1-p
#pm=median(p)
pm=median(p[p>0 & p<1])
lambda=qchisq(pm,1)/qchisq(0.5,1)
x2=qchisq(p,1)/lambda
p=pchisq(x2,1)
p=1-p
p
}

z$P_CT_corrected=gc(z[,9])
z$FDR_CT_corrected=p.adjust(z$P_CT_corrected,method="fdr")

pdf(output,10,7 )

par(mfrow=c(1,2))

# plot 1: uncorrected
p=z[,9]
p = p[p>0]
p = p[p<1]
p = p[!is.na(p)]
OBS = sort(-log10(p))
EXP = sort(-log10(1:length(p)/length(p)))
plot(EXP, OBS, col="red", pch=20); abline(a=0,b=1, col="lightgray", lty=1, lwd=2)
title("SMG test qq-plot")


# plot 2: GC-corrected
p=z$P_CT_corrected
p = p[p>0]
p = p[p<1]
p = p[!is.na(p)]
OBS = sort(-log10(p))
EXP = sort(-log10(1:length(p)/length(p)))
plot(EXP, OBS, col="red", pch=20); abline(a=0,b=1, col="lightgray", lty=1, lwd=2)
title("SMG test qq-plot(GC corrected)")

dev.off()

#write.csv(z,paste(input,"corrected",sep="."),row.names=F,quote=F)
write.csv(z, file="smgs_detailed.corrected", row.names=F, quote=F)

