####### piston function #######

###  function that scales functions to unit input
unit.scale=function(fun,ranges){
  scaled.fun=function(x){
    xx = x*(ranges[,2]-ranges[,1])+ranges[,1]
    fun(xx)
  }
  return(scaled.fun)
}


pistonfun <- function(xx)
{
  # OUTPUT AND INPUT:
  # C = cycle time
  # xx = c(M, S, V0, k, P0, Ta, T0)
  
  M  <- xx[1]
  S  <- xx[2]
  V0 <- xx[3]
  k  <- xx[4]
  P0 <- xx[5]
  Ta <- xx[6]
  T0 <- xx[7]
  
  Aterm1 <- P0 * S
  Aterm2 <- 19.62 * M
  Aterm3 <- -k*V0 / S
  A <- Aterm1 + Aterm2 + Aterm3
  
  Vfact1 <- S / (2*k)
  Vfact2 <- sqrt(A^2 + 4*k*(P0*V0/T0)*Ta)
  V <- Vfact1 * (Vfact2 - A)
  
  fact1 <- M
  fact2 <- k + (S^2)*(P0*V0/T0)*(Ta/(V^2))
  
  C <- 2 * pi * sqrt(fact1/fact2)
  return(C)
}

piston.ranges=matrix(c(30, 60,
                       0.005, 0.020,
                       0.002, 0.010,
                       1000, 5000,
                       90000, 110000,
                       290, 296,
                       340, 360),ncol=2,byrow=TRUE)

piston=unit.scale(pistonfun,piston.ranges)
