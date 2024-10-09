
check_converge <- function(array_nll, algs){
  result <- array(NA, dim = dim(array_nll))
  
  for(i in 1:dim(array_nll)[3]){
    for(j in 1:dim(array_nll)[2]){
      list_nll <- array_nll[,j,i]
      med_nll <- median(list_nll, na.rm = TRUE)
      if(med_nll >= 0){
        result[,j,i] <- array_nll[,j,i] <= 1.001*med_nll
      }else{ 
        #if final NLL is below zero
        result[,j,i] <- array_nll[,j,i] <= 0.999*med_nll
      }
    }
  }
  dimnames(result)[[1]] <- algs
  return(result)
}

check_converge_10000 <- function(array_nll, algs){
  result <- array(NA, dim = dim(array_nll))
  
  for(i in 1:dim(array_nll)[3]){
    for(j in 1:dim(array_nll)[2]){
      list_nll <- array_nll[,j,i]
      med_nll <- median(list_nll, na.rm = TRUE)
      if(med_nll >= 0){
        result[,j,i] <- array_nll[,j,i] <= 1.0001*med_nll
      }else{ 
        #if final NLL is below zero
        result[,j,i] <- array_nll[,j,i] <= 0.9999*med_nll
      }
    }
  }
  dimnames(result)[[1]] <- algs
  return(result)
}

check_converge_norep <- function(array_nll, algs){
  result <- array(NA, dim = dim(array_nll))
  
  for(j in 1:dim(array_nll)[2]){
    list_nll <- array_nll[,j]
    med_nll <- median(list_nll, na.rm = TRUE)
    if(med_nll >= 0){
      result[,j] <- array_nll[,j] <= 1.001*med_nll
    }else{ 
      #if final NLL is below zero
      result[,j] <- array_nll[,j] <= 0.999*med_nll
    }
  }
  
  dimnames(result)[[1]] <- algs
  return(result)
}