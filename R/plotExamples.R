
plotExampleSingleGeneEffects <- function(gene,features,mainEffects,ylim) {

  # For each target gene, features and replicate, the mean over all siRNA
  # designs is computed  
  MT1 = apply(mainEffects$target,c(1,4,5),mean,na.rm=TRUE)

  # The main effects are devided by the median deviation per feature to make
  # main effects comparable accross phenotypes.
  for (i in seq_len(dim(MT1)[2])) {
    for (j in seq_len(dim(MT1)[3])) {
      MT1[,i,j] = MT1[,i,j] / mad(MT1[,i,j],center=0.0)
    }
  }
  MT = apply(MT1,c(1,2),mean,na.rm=TRUE)

  # The bars show the mean over the replicate measurements
  MT2 = MT[,features]
  colnames(MT2) = HD2013SGI:::humanReadableNames[features]

  if (missing(ylim)) {
    ylim=range(MT2[gene,])
    if (ylim[1] > 0.0) { ylim[1] = 0.0 }
    if (ylim[2] < 0.0) { ylim[2] = 0.0 }
  }

  par(mar=c(12,4,4,1)+0.1,xpd=NA)
  bp = barplot(MT2[gene,],ylab="z-score",ylim=ylim,las=2,main=gene)
  abline(h=0.0,xpd=FALSE)

  # Circle symbols represent individual replicate measurements
  points(bp,MT1[gene,features,1])
  points(bp,MT1[gene,features,2])

}

plotExampleInteractions <- function(feature,target,query,Interactions,mainEffects) {
  PI = Interactions$piscore
  PI = PI[target,,query,,feature,]
  dim(PI) = c(prod(dim(PI)[1:2]),dim(PI)[3])
  
  # mean and standard dev. over replicates
  X = apply(PI,1,mean)
  names(X) = c("pi11","pi21","pi12","pi22")  
  SD = apply(PI,1,sd)

  par(xpd=NA,mar=c(8,5,4,1))
  bp=barplot(X,las=2, #col=col,
             main=sprintf("target=%s\nquery=%s\n%s",target,query,
                            HD2013SGI:::humanReadableNames[feature]),
             names.arg = rep("",4))
  abline(h=0.0,xpd=FALSE)
  
  # Circle symbols represent individual replicate measurements
  points(bp,PI[,1])
  points(bp,PI[,2])
  
  
  # Text (+,-) indicating which siRNA designs are used
  r = range(X)
  if (r[1] > 0) { r[1] = 0.0 }
  if (r[2] < 0) { r[2] = 0.0 }
  d = diff(r)
  A = c("+")
  B = c("+")
  C = c("-")
  text(x=bp,y=rep(r[1]-2*d/7,4),c(A,C,A,C))
  text(x=bp,y=rep(r[1]-3*d/7,4),c(C,A,C,A))
  text(x=bp,y=rep(r[1]-4*d/7,4),c(B,B,C,C))
  text(x=bp,y=rep(r[1]-5*d/7,4),c(C,C,B,B))
  text(rep(0.0,4),r[1] - (2:5)*d/7,adj=c(1,0.5),
       sprintf("%s #%d",c(target,target,query,query),c(1,2,1,2)))
}

