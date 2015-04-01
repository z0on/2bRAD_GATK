ld=read.table("out.geno.ld",header=T)  
ld$dist=(ld$POS2-ld$POS1)
ld$coord=paste(ld$CHR,ld$POS1,sep=".")
ld$coord2=paste(ld$CHR,ld$POS2,sep=".")
ld=ld[!is.na(ld$R.2),]
str(ld)
ld.q=c()
w=500
for (d in 1:nrow(ld)) {
	dis=ld$dist[d]
	r2=ld$R.2[d]
	if (dis>15000 ) {next}
	bin=ceiling(dis/w)
	ld.q=data.frame(rbind(ld.q,c("r2"=r2,"bin"=bin)))
}

ldq=c()
for (b in 1:30) {
	s=subset(ld.q,bin==b)
	qu=	quantile(1-s$r2,probs = c(0.0001,0.001,seq(0.01,1,0.01)))
	ldq=data.frame(cbind(ldq,1-qu))
}

ldtest=function(r2,dist,w=500,quantiles=ldq) {
	bin=ceiling(dist/w)
	if (bin>30) {
		return(NA)
	} else {
		q=sum(r2>=quantiles[,bin])-1
		if (q==101) { 
			return(0.0001) 
		}
		if (q==100) { 
			return(0.001) 
		} else {
			return(1-q/100)
		}
	}
}

pvals=c()
coords=unique(c(ld$coord,ld$coord2))
for (co in coords) {
	s=rbind(subset(ld,coord==co),subset(ld,coord2==co))
	ps=c()
	for(i in 1:nrow(s)) {
		if (s$dist[i]>15000	| is.na(s$R.2[i]) ) { next }
		ps=append(ps,ldtest(s$R.2[i],s$dist[i]))
	}
	if(length(ps)>0 && !is.na(ps)) {
		pv=NA
		if (length(ps)>1) {
			ftest=-2*sum(log(ps))
			df=2*length(ps)
			pv=1-pchisq(q=ftest,df=df)
		} else {
			pv=ps
		}
		pvals=append(pvals,pv)
	} else {
		pvals=append(pvals,NA)
	}
}
length(coords)

si=strsplit(as.character(coords),split=".",fixed=T)
save(si,pvals,ldq,ld.q,file="LD_quantized.RData")

# run these commands after the script is done:

# load("LD_quantized.RData")
# library(plyr)
# ldp=ldply(si)
# ldp$lpv=-log(pvals,10)
# names(ldp)=c("CHR","POS","pLD")
# ldp=ldp[!is.na(ldp$pLD),]
# ldp$POS=as.numeric(ldp$POS)
# ldp$coord=paste(ldp$CHR,ldp$POS,sep=".")
# save(ldp,ldq,ld.q,file="LD_Quantized.RData")
