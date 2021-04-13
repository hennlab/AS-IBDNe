# Plot as-IBDne
# Source: Browning et al https://doi.org/10.1371/journal.pgen.1007385
# Usage: Rscript plot_ibdne.R {params.prefix} {params.suffix} {output} {comma separate list of ref pops} {number of ref ancestries}
# Example: Rscript plot_ibdne.R results/IBDne/americans.anc _2cM.ibdne.ne results/plots/americans.ibdne.pdf CEU,PEL,YRI 3 blue,green,orange,purple
# where the each ancestry file is: results/IBDne/americans.anc1_2cM.ibdne.ne results/IBDne/americans.anc2_2cM.ibdne.ne results/IBDne/americans.anc3_2cM.ibdne.ne


args = commandArgs(trailingOnly=TRUE)

pops = unlist(strsplit(args[4], split=","))
nanc <- args[5]
colors = unlist(strsplit(args[6], split=","))


pdf(args[3],height=3,width=10)
par(mfrow=c(1,3),mar=c(4.5,4.5,3.5,1.5)+0.1)
for(anc in 1:nanc){
  filename=paste(args[1],anc,args[2],sep="")
  print(anc)
  print(filename)
  x=read.table(filename,header=T)
  col2="#00000080"
  plot(x[,1],x[,2],type="l",log="y",xlim=c(0,100),ylim=c(1e3,1e7),xlab="g (generations before present)",ylab="Ne (effective population size)",xaxs="i",yaxt="n",las=1,lwd=2.5,lty=32,col="lightgray",main=paste(pops[anc],"ancestry"))
  abline(h=c(1e2,1e3,1e4,1e5,1e6,1e7),v=seq(20,80,20),col="lightgray",lty="dotted")
  axis(2,at=c(1e3,1e4,1e5,1e6,1e7),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6),expression(10^7)),las=1)
  polygon(c(x[,1],x[nrow(x):1,1]),c(x[,3],x[nrow(x):1,4]),col=colors[anc],density=-1,border=col2)
  lines(x[,1],x[,2],lwd=1)
}
dev.off()
