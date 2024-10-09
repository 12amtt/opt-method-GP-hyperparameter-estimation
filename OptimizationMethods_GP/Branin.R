
### generate from original Branin function
branin <- function(xx, a=1, b=5.1/(4*pi^2), c=5/pi, r=6, s=10, t=1/(8*pi), nugget)
{
  x1 <- xx[,1]
  x2 <- xx[,2]
  
  term1 <- a * (x2 - b*x1^2 + c*x1 - r)^2
  term2 <- s*(1-t)*cos(x1)
  
  y <- term1 + term2 + s
  if (nugget > 0) {
    y <- y + rnorm(n, sd = sqrt(nugget)) #nugget
  }
  return(y)
}


### generate data from rescaled Branin function
braninsc <- function(xx, nugget)
{
  x1 <- xx[,1]
  x2 <- xx[,2]
  
  x1bar <- 15*x1 - 5
  x2bar <- 15 * x2
  
  term1 <- x2bar - 5.1*x1bar^2/(4*pi^2) + 5*x1bar/pi - 6
  term2 <- (10 - 10/(8*pi)) * cos(x1bar)
  
  y <- (term1^2 + term2 - 44.81) / 51.95
  if (nugget > 0) {
    y <- y + rnorm(n, sd = sqrt(nugget)) #nugget
  }
  return(y)
}





