
parsePlateBarcodes <- function(plates) {
  platenr = substr(plates,1,3)

  # Each plate contains the same set of taerget siRNAs. In all wells of one
  # plate the same query siRNA is pipetted. The query siRNA on the
  # \Sexpr{length(platenr)} plates are grouped in sample (double gene knock
  # down) and negative control to measure the main effects (single knock down
  # effects) of the target genes.

  queryGroup = rep("sample",length(plates))
  queryGroup[grep("N",plates)] = "negControl"

  # The remainder of the plate barcodes contain the targetDesign (CI or CII).

  r = substr(plates,4,10000)
  #  print(head(r))

  S = which(queryGroup == "sample")
#  N = 161:168
  targetDesign = sapply(strsplit(r,split="[QN]"),function(x) { x[1] } )
  targetDesign[targetDesign == "CI"] = 1
  targetDesign[targetDesign == "CII"] = 2
  targetDesign = as.integer(targetDesign)

  # The remainder of the plate barcodes contain the query gene.

  r = sapply(strsplit(r,split="[QN]"),function(x) { x[2] } )
  #  print(head(r))
  
  queryGene = rep("NegControl",length(plates))
  queryGene[S] = substr(r[S],1,2)

  # The remainder of the plate barcodes contain the query design.

  r[S] = substr(r[S],3,100)
  #  print(head(r))
  
  queryDesign = sapply(strsplit(r,split="[R]"),function(x) { x[1] } )
  queryDesign[queryDesign == "I"] = 1
  queryDesign[queryDesign == "II"] = 2
  queryDesign = as.integer(queryDesign)

  # The remainder of the plate barcodes contain the biological replicate.
  replicate = sapply(strsplit(r,split="[R]"),function(x) { x[2] } )
  replicate[replicate == "I"] = 1
  replicate[replicate == "II"] = 2
  replicate = as.integer(replicate)

  # The plate annotation is summarized in a table.
  PlateAnnotation = data.frame(plate = plates,
                               targetDesign = targetDesign,
                               queryGroup = queryGroup,
                               queryGene = queryGene,
                               queryDesign = queryDesign, 
                               replicate = replicate,
                               stringsAsFactors=FALSE)

  PlateAnnotation
}

