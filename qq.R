find_conf_intervals = function(row){
  i = row[1]
  len = row[2]
  if (i < 10000 | i %% 100 == 0){
    return(c(-log10(qbeta(0.95,i,len-i+1)), -log10(qbeta(0.05,i,len-i+1))))
  } else { # Speed up
    return(c(NA,NA))
  }
}

confidence.intervals <- function(e){
  xspace = 0.078
  print("1")
  ci = apply(cbind( 1:length(e), rep(length(e),length(e))), MARGIN=1, FUN=find_conf_intervals)
  print("2")
  bks = append(seq(10000,length(e),100),length(e)+1)
  print("3")
  for (i in 1:(length(bks)-1)){
    ci[1, bks[i]:(bks[i+1]-1)] = ci[1, bks[i]]
    ci[2, bks[i]:(bks[i+1]-1)] = ci[2, bks[i]]
  }
  colnames(ci) = names(e)
  ## Extrapolate to make plotting prettier (doesn't affect intepretation at data points)
  slopes = c((ci[1,1] - ci[1,2]) / (e[1] - e[2]), (ci[2,1] - ci[2,2]) / (e[1] - e[2]))
  print("4")
  extrap_x = append(e[1]+xspace,e) ## extrapolate slightly for plotting purposes only
  extrap_y = cbind( c(ci[1,1] + slopes[1]*xspace, ci[2,1] + slopes[2]*xspace), ci)
  print("5")
  polygon(c(extrap_x, rev(extrap_x)), c(extrap_y[1,], rev(extrap_y[2,])),
          col = "grey81", border = "grey81")
}

quant.subsample <- function(y, m=100, e=1) {
  ## m: size of a systematic sample
  ## e: number of extreme values at either end to use
  x <- sort(y)
  n <- length(x)
  quants <- (1 + sin(1:m / (m+1) * pi - pi/2))/2
  sort(c(x[1:e], quantile(x, probs=quants), x[(n+1-e):n]))
  ## Returns m + 2*e sorted values from the EDF of y
}

get.points <- function(pv) {
  suppressWarnings(as.numeric(pv))
  names(d) = names(pv)
  d = d[!is.na(d)]
  d = d[d>0 & d<1]
  d = d[order(d,decreasing=F)]
  y = -log10(d)
  x = -log10( ppoints(length(d) ))
  m <- 0.001 * length(x)
  e <- floor(0.0005 * length(x))
  return(list(x=quant.subsample(x, m, e), y=quant.subsample(y, m, e)))
}

fqq <- function(x, y, ...) {
  plot(0,
       col=FALSE,
       xlim=range(x),
       ylim=range(y),
       xlab=expression(Expected~~-log[10](italic(p))),
       ylab=expression(Observed~~-log[10](italic(p))),
       ...)
  abline(0,1,col=2)
  points(x,y, ...)
}

args <- commandArgs(trailingOnly = TRUE)
pv.f = args[1]
qq.f = args[2]
nrows = as.numeric(args[3])
message(Sys.time())
message("READING")
d <- read.table(pv.f, header=TRUE, sep=" ", nrows=nrows, colClasses=c("numeric"))
message(Sys.time())
message("LAMBDA")
chisq <- qchisq(1-d$P_VAL,1)
lambda = median(chisq)/qchisq(0.5,1)
message(Sys.time())
message("PLOTING")
p <- get.points(d$P_VAL)
png(file=qq.f)
fqq(p$x, p$y, main=paste(pv.f, lambda, sep="\n"), cex.axis=1.5, cex.lab=1.5)
dev.off()
message(Sys.time())
