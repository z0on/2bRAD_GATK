getwd()
setwd("/Users/c-monster/Documents/MethEcoGeno_fall2014/rad_gatk/genomeScans")

# install.packages("gridExtra")
# install.packages("ggplot2")

library(ggplot2)
library(gridExtra)

###########################

nn=load("LD_quantized.RData")
head(ld.q)
library(plyr)
ldp=ldply(si)
ldp$lpv=-log(pvals,10)
names(ldp)=c("CHR","POS","pLD")
ldp=ldp[!is.na(ldp$pLD),]
ldp$POS=as.numeric(ldp$POS)
head(ldp,30)
save(ldp,ld.q,ldq,pvals,si,file="LD_quantized_done.RData")
############################
load("LD_quantized_done.RData")

head(ldp)

# putting bayescan table together
bs=read.table("OK.baye_fst.txt",header=T)
system('grep -E "^#" OK.recode.vcf | wc -l') # put the number after skip= in the next line
vcf=read.table("OK.recode.vcf",skip=49)
bscan=data.frame(cbind("POS"=vcf[,2],"Q"=bs[,3],"Fst"=bs[,5]))
bscan=cbind("CHROM"=vcf[,1],bscan)
bscan$Q=-log(bscan$Q+1e-5,10)

head(bscan)

fst=read.table("OK.fst",header=T)
names(fst)[3]="Fst"
kpi=read.table("K.pi",header=T)  
opi=read.table("O.pi",header=T)
ktd=read.table("K.td",header=T)
otd=read.table("O.td",header=T)
ktd=ktd[ktd$N_SNPS>1,]
otd=otd[otd$N_SNPS>1,]
ktd$POS=ktd$BIN_START+37500
otd$POS=otd$BIN_START+37500
ktd$pop="K"
otd$pop="O"
table(ktd$N_SNPS>3)

# finding interesting contigs
longs=unique(ldp$CHR[ldp$POS>100000])
hild=unique(ldp$CHR[ldp$pLD>2])
hifst=unique(fst$CHROM[fst$Fst>0.3])
hibs=unique(bscan$CHROM[bscan$Q>1.33])
hitd=unique(ktd$CHROM[ktd$TajimaD>1.8])
i1=intersect(hifst,longs)
i2=intersect(i1,longs)
i1a=intersect(hibs,longs)
i2a=intersect(i1a,hild)

# function to force simple formatting of axis labels
fmt <- function(){
    function(x) format(x,nsmall = 1,scientific = FALSE)
}

plot(density(ktd$TajimaD),col="cyan3",xlim=c(-3.5,3.5),ylim=c(0,0.5),lwd=1.5,bty="n",main=NA,mgp=c(2.3,1,0),xlab="Tajima's D in 75 kB bins")
lines(density(otd$TajimaD),col="coral",lwd=1.5)
legend("topleft",bty="n",lty=1,col=c("coral","cyan3"),legend=c("Orph","Kepp"),cex=0.9,lwd=1.5)

c="c2013605"

#c="c1855353"
#c="c2379886"


c1=subset(kpi,CHROM==c)
c2=subset(opi,CHROM==c)
f=subset(fst,CHROM==c)
l=subset(ldp,CHR==c)
kt=subset(ktd,CHROM==c)
ot=subset(otd,CHROM==c)
bbs=subset(bscan,CHROM==c)

plot1 <- ggplot(f, aes(POS,Fst))+geom_point(shape=1,cex=2,colour="grey50")+geom_smooth(colour="blueviolet",se=F,span=0.35,lwd=1)+theme_bw()+theme(legend.position="none")+xlim(0,max(f$POS))+theme(axis.title.x = element_blank())

plot2 <- ggplot(bbs, aes(POS,Q))+geom_point(shape=1,cex=2,colour="grey50")+geom_smooth(colour="steelblue",se=F,span=0.35,lwd=1)+theme_bw()+theme(legend.position="none")+xlim(0,max(f$POS))+theme(axis.title.x = element_blank())+ylab("bayescan")+scale_y_continuous(labels = fmt())

plot3 <- ggplot(c1, aes(POS,PI))+geom_point(shape=1,cex=2,colour="grey50")+geom_smooth(colour="cyan3",span=0.35,se=F,lwd=1)+theme_bw()+theme(legend.position="none")+xlim(0,max(c1$POS))+theme(axis.title.x = element_blank())+ylab("PI:Keppels")

plot4 <- ggplot(c2, aes(POS,PI))+geom_point(shape=1,cex=2,colour="grey50")+geom_smooth(colour="coral",span=0.35,se=F,lwd=1)+theme_bw()+theme(legend.position="none")+xlim(0,max(c2$POS))+theme(axis.title.x = element_blank())+ylab("PI:Orpheus")

plot5=ggplot(kt, aes(POS,TajimaD))+geom_bar(stat="identity",fill="cyan3")+theme_bw()+theme(legend.position="none")+theme(axis.title.x = element_blank())+ylab("TajD:Keppels")+geom_abline(intercept=0,slope=0,lwd=0.3,lty=3)+scale_y_continuous(labels = fmt())

plot6=ggplot(ot, aes(POS,TajimaD))+geom_bar(stat="identity",fill="coral")+theme_bw()+theme(legend.position="none")+theme(axis.title.x = element_blank())+ylab("TajD:Orpheus")+geom_abline(intercept=0,slope=0,lwd=0.3,lty=3)+scale_y_continuous(labels = fmt())

#plot5 <- ggplot(kt, aes(POS,TajimaD))+geom_point(shape=1,cex=2,colour="grey50")+geom_smooth(colour="cyan3",se=F,span=0.5,lwd=1,lty=2)+theme_bw()+theme(legend.position="none")+xlim(0,max(f$POS))+geom_abline(intercept=0,slope=0,lwd=0.3,lty=3)+theme(axis.title.x = element_blank())+ylab("TajD:Keppels")

#plot6 <- ggplot(ot, aes(POS,TajimaD))+geom_point(shape=1,cex=2,colour="grey50")+geom_smooth(colour="coral",se=F,span=0.5,lwd=1,lty=2)+theme_bw()+theme(legend.position="none")+xlim(0,max(f$POS))+geom_abline(intercept=0,slope=0,lwd=0.3,lty=3)+theme(axis.title.x = element_blank())+ylab("TajD:Orpheus")

plot7 <- ggplot(l, aes(POS,pLD))+geom_point(shape=1,cex=2,colour="grey50")+geom_smooth(colour="chartreuse",se=F,span=0.3,lwd=1)+scale_y_continuous(labels = fmt())+theme_bw()+theme(legend.position="none")+xlim(0,max(f$POS))+ylab("-log10(P) LD")

quartz()
grid.arrange(plot1, plot2, plot3, plot4,plot5,plot6, plot7,nrow=7)


