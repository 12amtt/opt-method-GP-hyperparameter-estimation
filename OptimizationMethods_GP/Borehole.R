#########   borehole function  #########

###  function that scales functions to unit input
unit.scale=function(fun,ranges){
  scaled.fun=function(x){
    xx = x*(ranges[,2]-ranges[,1])+ranges[,1]
    fun(xx)
  }
  return(scaled.fun)
}

borehole <- function(xx){
  # OUTPUT AND INPUT:
  # y  = water flow rate
  # xx = c(rw, r, Tu, Hu, Tl, Hl, L, Kw)
  
  rw <- xx[1]
  r  <- xx[2]
  Tu <- xx[3]
  Hu <- xx[4]
  Tl <- xx[5]
  Hl <- xx[6]
  L  <- xx[7]
  Kw <- xx[8]
  
  frac1 <- 2 * pi * Tu * (Hu-Hl)
  
  frac2a <- 2*L*Tu / (log(r/rw)*rw^2*Kw)
  frac2b <- Tu / Tl
  frac2 <- log(r/rw) * (1+frac2a+frac2b)
  
  y <- frac1 / frac2
  return(y)
}


bore.ranges=matrix(c(.05,.15,
                     100,50000,
                     63070,115600,
                     990,1110,
                     63.1,116,
                     700,820,
                     1120,1680,
                     9855,12045),ncol=2,byrow=TRUE)

#### generating function
bore=unit.scale(borehole,bore.ranges)

