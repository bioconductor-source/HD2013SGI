
HD2013SGIselectByStability <- function(subsample,
                               preselect = "count",
                               Rdim = 40, 
                               verbose = TRUE) {

  n = 1
  I = phenotype = 5
  m = 3

  D = subsample$D
  phenotype = subsample$phenotype
  
  Sel = match(preselect, phenotype)
  d = dim(D)
  dimnames(D) = list(NULL,NULL,phenotype)
  correlation = rep(NA_real_, Rdim)
  correlationAll = list()
  ratioPositive = rep(NA_real_, Rdim)
  selected = rep("", Rdim)
  DSel = D[,1,]
  DSel[] = NA_real_
  Dall = D
  for (k in 1:Rdim) {
    if (k > 1) {
      for (i in 2:dim(D)[3]) {
        model = lm(as.vector(D[,1,i]) ~ DSel[,1:(k-1),drop=FALSE]+0)
        D[,1,i] = model$residuals
        model = lm(as.vector(D[,2,i]) ~ DSel[,1:(k-1),drop=FALSE]+0)
        D[,2,i] = model$residuals
      }
    }
    C = apply(D, 3, function(x) {
      a = x[,1]
      b = x[,2]
      I = which(is.finite(a) & is.finite(b))
      cor(a[I],b[I]) 
      } )
    if (k > length(preselect)) {
      I = names(C)[which.max(C)]
    } else {
      I = preselect[k]
    }
    correlation[k] = C[I]
    ratioPositive[k] = sum(C > 0, na.rm=TRUE) / length(C)
    selected[k] = I
    if ( verbose ) {
      message("k=",k," selected = ",selected[k], " cor = ", correlation[k], " r = ", ratioPositive[k])
    }
    correlationAll[[k]] = C
    DSel[,k] = apply(D[,,I,drop=FALSE],1,mean,na.rm=TRUE)
    D = D[,,dimnames(D)[[3]] != I,drop=FALSE]
  }

  res = list(selected = selected, correlation = correlation, ratioPositive = ratioPositive, correlationAll = correlationAll)
  res
}

