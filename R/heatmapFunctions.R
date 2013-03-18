
trsf = function(x) {
  exp(-x)-exp(-1)
}

toRaster = function(x, cuts, col) {
  cux =cut(x,cuts,include.lowest = TRUE, labels=FALSE)
  rv = x
  rv[] = col[cux]
  return(rv)
}

toMatrix = function(x) {
  dim(x) = c(dim(x)[1], prod(dim(x)[-1]))
  return(x)
}

HD2013SGIorderDim = function (x, i) {
  px = toMatrix(aperm(x, c(i, setdiff(seq(along = dim(x)), 
                                      i))))
  row.names(px) = dimnames(x)[[i]]
  theCorr = cor(t(px), use = "pairwise.complete.obs")
  hc = hclust(as.dist(trsf(theCorr)))
  return(hc)
}

HD2013SGIHeatmapHuman <- function (x, cuts, col, colnames = TRUE, rownames = FALSE, 
                            mrow=10, mcol=10, cexrow=1, cexcol=1,border=0.1,space=0.05) {
  stopifnot(is.array(x), length(dim(x)) == 3)
  mrow = unit(cexrow*mrow,"lines")
  mcol = unit(cexcol*mcol,"lines")
  rx = toRaster(x, cuts = cuts, col = col)
  width = rep(1,2*dim(x)[3]-1)
  width[seq(2,length(width),by=2)] = space
  u = rep("null",2*dim(x)[3]-1)
  u[seq(2,length(width),by=2)] = "lines"
  width=unit(width,u)
  pushViewport(viewport(x = unit(1,"npc")-unit(border,"lines"), y = unit(border,"lines"),just=c("right","bottom"),
                        width = unit(1, "npc")-unit(2*border,"lines")-mrow,
                        height = unit(1, "npc")-unit(2*border,"lines")-mcol,
                        layout=grid.layout(nrow=1,ncol=2*dim(x)[3]-1,
                                           widths=width,heights=unit(1,"npc"))))
  grid.rect()
  for (i3 in 1:dim(x)[3]) {
    pushViewport(viewport(layout.pos.row=1,layout.pos.col=i3*2-1))
    grid.raster(rx[, , i3], x = unit(0.5, "npc"), width = unit(1, "npc"),
                y = unit(0.5, "npc"), height = unit(1, "npc"), interpolate = FALSE)
    popViewport()
  }
  
  if (rownames) {
    pushViewport(viewport(layout.pos.row=1,layout.pos.col=1,clip="off",yscale=c(dim(rx)[1]+0.5,0.5)))
    grid.text(dimnames(rx)[[1]],x=unit(rep(-0.25,dim(rx)[1]),"lines"),y=unit(1:dim(rx)[1],"native"),
              gp=gpar(cex=cexrow),just=c("right","center"))
    popViewport()
  }
  
  if (colnames) {
    for (i3 in 1:dim(x)[3]) {
      pushViewport(viewport(layout.pos.row=1,layout.pos.col=i3*2-1))
      grid.text(dimnames(rx)[[3]][i3],x=unit(0.5,"npc"),y=unit(1,"npc")+unit(0.25,"lines"),gp=gpar(cex=cexcol),just=c("left","center"),rot=90)
      popViewport()
    }
  }
  
  popViewport()
}

