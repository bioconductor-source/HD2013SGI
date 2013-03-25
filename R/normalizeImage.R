
normalizeImage <- function(Img, minInt, maxInt) {
  Img[Img > maxInt] = maxInt
  Img[] = Img[] - minInt
  Img[Img < 0] = 0
  Img[] = Img[] / (maxInt - minInt)
  Img[Img > 1] = 1
  Img
}

permuteColor <- function(x, normalize = TRUE) {
  M <- max(x)
  R1 <- sample(M)
  R2 <- sample(M)
  R3 <- sample(M)
  ch1 = x
  ch2 = x
  ch3 = x
  ch1[ch1 > 0] <- R1[ch1[ch1>0]]
  ch2[ch2 > 0] <- R2[ch2[ch2>0]]
  ch3[ch3 > 0] <- R3[ch3[ch3>0]]
  Img = Image(data=combine(ch1,ch2,ch3),colormode="Color")
  if (normalize) {
    Img <- normalize(Img)
  }
  Img
}
