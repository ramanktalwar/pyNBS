
library(survival)

args <- commandArgs(trailingOnly = TRUE)

dir=args[1]

sink(paste0(dir,"/coxph.txt"))

d=read.delim(paste0(dir,"/pat2surv2labels.txt"),header=T,row.names=1)

d=d[which(d$K2!="NaN"),]
d=d[which(d$DFS_MONTH!="NA"),]
d=d[which(d$DFS_MONTH>0),]

# output file format: K, OS_simple, DFS_simple
out=as.data.frame(matrix(ncol=3))
v=colnames(d)[c(-1,-2,-3,-4,-5,-6)]
for (i in 1:length(v)) {
    d2=d[,c(1,2,3,4,5,6,i+6)]
    d2$cluster=factor(d2[[v[i]]])
    d2$cluster_relevel <- relevel(d2$cluster, ref = names(which.max(summary(d2$cluster))))
    cr_OS_simple=tryCatch(coxph(Surv(OS_MONTHS,OS_STATUS)~cluster_relevel,data=d2,control = coxph.control(iter.max = 500, toler.chol = 1e-210)), error = function(e) e)
    cr_DFS_simple=tryCatch(coxph(Surv(DFS_MONTHS,DFS_STATUS)~cluster_relevel,data=d2,control = coxph.control(iter.max = 500, toler.chol = 1e-210)), error = function(e) e)
    out[i,1]=v[i]
    out[i,2]=summary(cr_OS_simple)$logtest[3]
    out[i,3]=summary(cr_DFS_simple)$logtest[3]
    pdf(paste0(dir,"/",v[i],"_OS.pdf"))
    fit <- survfit(Surv(OS_MONTHS,OS_STATUS)~cluster, data=d2)
    print(fit)
    print("#############################################################################################################################################################################################")
    plot(fit,col=rainbow(i+1),xlab="Overall survival (months)",ylab="Survival probability",xlim=c(0,120),mark.time=TRUE)
    legend("topright",levels(factor(d2[,7])),col=rainbow(i+1),lty=rep(1,i+1))
    dev.off()
    pdf(paste0(dir,"/",v[i],"_DFS.pdf"))
    fit <- survfit(Surv(DFS_MONTHS,DFS_STATUS)~cluster, data=d2)
    print(fit)
    print("#############################################################################################################################################################################################")
    plot(fit,col=rainbow(i+1),xlab="Disease free survival (months)",ylab="Survival probability",xlim=c(0,120),mark.time=TRUE)
    legend("topright",levels(factor(d2[,7])),col=rainbow(i+1),lty=rep(1,i+1))
    dev.off()
}
write.table(out,paste0(dir,"/coxph_FTest.txt"),sep="\t",row.name=F, quote =F,col.name=F)

sink()

