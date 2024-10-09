##### Generate data and fit GP with Matern () #####

### necessary packages
library(RandomFields)
# install.packages("remotes")
# library(remotes)
# install_version("RandomFieldsUtils", "1.2.5")
# install_version("RandomFields", "3.3.14")
library(gpboost)

##### generating function #####
#### simulate data from GP with Matern kernel
gen_data <- function(coords, n, rho, sigma2, sigma2_error, nu){
  if (rho == 0) {
    iid_no_GP <- TRUE
  } else {
    iid_no_GP <- FALSE
  }
  
  
  if (iid_no_GP) {
    eps <- rnorm(n, sd = sqrt(sigma2))
  } else {
    if (nu == 0.5) {
      RFmodel <- RMexp(var=sigma2, scale=rho)
    } else if (nu > 1e3) {
      RFmodel <- RMgauss(var=sigma2, scale=rho)
    } else {
      RFmodel <- RMmatern(var=sigma2, scale=rho, nu=nu)
    }
    sim <- RFsimulate(RFmodel, x=coords) # ignore warning
    # Note: check whether this works for d > 2 and if not, find work-around
    eps <- sim$variable1
    # eps <- eps - mean(eps)
  }
  y <- eps + rnorm(n, sd=sqrt(sigma2_error))
  return(y)
}

##### simulating y function #####
simulate_response_variable <- function (lp, rand_eff, likelihood) {
  ## Function that simulates response variable for various likelihoods
  n <- length(rand_eff)
  if (likelihood == "gaussian") {
    xi <- sqrt(0.1) * rnorm(n) # error term, variance = 0.1
    y <- lp + rand_eff + xi
  } else if (likelihood == "bernoulli_probit") {
    probs <- pnorm(lp + rand_eff)
    y <- as.numeric(runif(n) < probs)
  } else if (likelihood == "bernoulli_logit") {
    probs <- 1/(1+exp(-(lp + rand_eff)))
    y <- as.numeric(runif(n) < probs)
  } else if (likelihood == "poisson") {
    mu <- exp(lp + rand_eff)
    y <- qpois(runif(n), lambda = mu)
  } else if (likelihood == "gamma") {
    mu <- exp(lp + rand_eff)
    y <- qgamma(runif(n), scale = mu, shape = 1)
  } else if (likelihood == "negative_binomial") {
    mu <- exp(lp + rand_eff)
    y <- qnbinom(runif(n), mu = mu, size = 1.5)
  }
  return(y)
}

##### generate initial values ######
gen_init <- function(coords, y, nu){
  if(nu == Inf){
    denom = sqrt(3) 
  } else if(nu == 0.5){
    denom = 3
  } else if(nu == 1.5){
    denom = 4.7/sqrt(3)
  } else if(nu == 2.5){
    denom = 5.9/sqrt(5)
  } 
  if (length(y)>1000){
    idx = sample(1:n, 1000, replace = FALSE)
    init_cov_pars = c(rep(var(y[idx])/2, 2), mean(dist(coords[idx,]))/denom)
  } else {
    init_cov_pars = c(rep(var(y)/2, 2), mean(dist(coords))/denom)
  }
  return(init_cov_pars)
}

##### fitting function with GPBoost optimization methods #####
fit_matern <- function(y, inputs, optimizer_cov, X = NULL, nu, init_cov_pars=NULL,
                       use_nesterov_acc = FALSE, trace = FALSE, gp_approx = NULL,
                       delta_rel_conv = 1e-8, likelihood='gaussian', maxit=1000){
  ## dimensions
  n = nrow(inputs)
  d = ncol(inputs)
  
  if (nu == Inf){
    cov_function = 'gaussian'
  } else {
    cov_function = 'matern'
  }

  ## fit and record time
  start.time <- Sys.time()
  if (n <= 2000){
    gp_model = fitGPModel(gp_coords = inputs, cov_function = cov_function, 
                          X=X, likelihood=likelihood, y=y, cov_fct_shape = nu, 
                          params = list(optimizer_cov = optimizer_cov,
                                        use_nesterov_acc = use_nesterov_acc, 
                                        trace = trace, init_cov_pars = init_cov_pars,
                                        delta_rel_conv = delta_rel_conv, maxit = maxit))
  }
  else {
    # if sample size > 2000, use "gp_approx" option to change approximation method
    gp_model = fitGPModel(gp_coords = inputs, cov_function = cov_function, 
                          X=X, likelihood=likelihood, y=y, cov_fct_shape = nu,
                          gp_approx = gp_approx, ind_points_selection = "random",
                          params = list(optimizer_cov = optimizer_cov,
                                        use_nesterov_acc = use_nesterov_acc, 
                                        trace = trace, init_cov_pars = init_cov_pars,
                                        delta_rel_conv = delta_rel_conv, maxit = maxit))
  }
  end.time <- Sys.time()
  time.taken <- round(as.numeric(difftime(time1 = end.time, time2 = start.time, units = "secs")), 3)
  return(c(gp_model$get_cov_pars(), 
           nll = gp_model$get_current_neg_log_likelihood(),
           time_converge = time.taken, 
           convergence = gp_model$get_num_optim_iter(), 
           coef = gp_model$get_coef()))

}


##### fitting function with optim() built-in optimization methods #####

fit_matern_optim = function(y, inputs, optimizer, X=NULL, nu, init_cov_pars=NULL, reltol = 1e-8,
                            trace = 0, likelihood = 'gaussian', gp_approx = NULL){
  ## dimensions
  n = nrow(inputs)
  d = ncol(inputs)

  if (nu == Inf){
    cov_function = 'gaussian'
  } else {
    cov_function = 'matern'
  }

  if (n <= 2000) {
    gp_model <- GPModel(gp_coords=inputs, cov_function=cov_function,
                        likelihood=likelihood, cov_fct_shape=nu)
  }
  else {
    gp_model <- GPModel(gp_coords=inputs, cov_function=cov_function,
                        likelihood=likelihood, cov_fct_shape=nu, gp_approx=gp_approx,
                        ind_points_selection = "random")
  }


  #### initial covariate coefficient
  if(!missing(X)){
    init_coef <- rep(0,ncol(X))
  } else{
    init_coef <- NULL
  }


  eval_nll <- function(pars, gp_model, y, X, likelihood) {

    if (likelihood == "gaussian") {
      cov_pars <- exp(pars[1:3])
    } else {
      cov_pars <- exp(pars[1:2])
    }
    if (likelihood == "gamma") {
      aux_pars <- exp(pars[3])
    } else {
      aux_pars <- NULL
    }
    if(!is.null(X)){
      fixed_effects <- as.numeric(X %*% pars[-c(1:3)])
      return(neg_log_likelihood(gp_model, cov_pars=cov_pars, y=y, aux_pars=aux_pars,
                                fixed_effects = fixed_effects))
    }
    else{
      return(neg_log_likelihood(gp_model, cov_pars=cov_pars, y=y, aux_pars=aux_pars))
    }
  }

  eval_nll_lbfgsb <- function(pars, gp_model, y, X, likelihood) {
    if (likelihood == "gaussian") {
        cov_pars <- pars[1:3]
      } else {
        cov_pars <- pars[1:2]
      }
      if (likelihood == "gamma") {
        aux_pars <- pars[3]
      } else {
        aux_pars <- NULL
      }

    if(!is.null(X)){
      fixed_effects <- as.numeric(X %*% pars[-c(1:3)])
      return(neg_log_likelihood(gp_model, cov_pars=cov_pars, y=y, aux_pars=aux_pars,
                         fixed_effects = fixed_effects))
    }
    else{
      return(neg_log_likelihood(gp_model, cov_pars=cov_pars, y=y, aux_pars=aux_pars))
    }
  }

  start.time <- Sys.time()
  if(optimizer == 'L-BFGS-B'){
    if(is.null(X)){
      opt <- optim(par = init_cov_pars, fn = eval_nll_lbfgsb, y=y, X=X, gp_model=gp_model,
                   likelihood=likelihood, method = optimizer, lower = rep(1e-10, 3),
                   control = list(trace = trace, factr=reltol))
    } else {
      opt <- optim(par = c(init_cov_pars, init_coef), fn = eval_nll_lbfgsb, y=y, X=X, gp_model=gp_model,
                   likelihood=likelihood, method = optimizer, lower = c(rep(1e-10, 3), -Inf),
                   control = list(trace = trace, factr=reltol))
    }
  }
  else if(optimizer == 'L-BFGS'){
    opt <- optim(par = c(log(init_cov_pars), init_coef), fn = eval_nll, y=y, X=X, gp_model=gp_model,
                 likelihood=likelihood, method = "L-BFGS-B",
                 control = list(trace = trace, factr=reltol))
  }
  else {
    opt <- optim(par = c(log(init_cov_pars), init_coef), fn = eval_nll, y=y, X=X, gp_model=gp_model,
                 likelihood=likelihood, method = optimizer,
                 control = list(trace = trace, reltol = reltol))
  }
  end.time <- Sys.time()
  time.taken <- round(as.numeric(difftime(time1 = end.time, time2 = start.time, units = "secs")), 3)
  if(optimizer == 'L-BFGS-B'){
    return(c(Error_term = opt$par[1],
             GP_var = opt$par[2],
             GP_range = opt$par[3],
             nll = opt$value,
             time_converge = time.taken,
             convergence = opt$convergence,
             coef = opt$par[-c(1:3)]))
  } else {
    return(c(Error_term = exp(opt$par)[1],
             GP_var = exp(opt$par)[2],
             GP_range = exp(opt$par)[3],
             nll = opt$value,
             time_converge = time.taken,
             convergence = opt$convergence,
             coef = opt$par[-c(1:3)]))
  }
  # estimated parameters, time and convergence
}

###### fitting module for binary data ######
fit_matern_binary <- function(y, inputs, optimizer_cov, X = NULL, nu, init_cov_pars=NULL,
                              use_nesterov_acc = FALSE, trace = FALSE, gp_approx = NULL,
                              delta_rel_conv = 1e-8, likelihood='binary', maxit=1000){
  ## dimensions
  n = nrow(inputs)
  d = ncol(inputs)
  
  if (nu == Inf){
    cov_function = 'gaussian'
  } else {
    cov_function = 'matern'
  }
  ## denominator of initial range parameter under different smoothness
  if(nu == Inf){
    denom = sqrt(3) 
  } else if(nu == 0.5){
    denom = 3
  } else if(nu == 1.5){
    denom = 4.7
  } else if(nu == 2.5){
    denom = 5.9
  } else {
    denom = 2
  }
  ## default initial parameters for matern covariance
  #...
  
  
  ## fit and record time
  start.time <- Sys.time()
  
  if (n <= 2000){
    gp_model = fitGPModel(gp_coords = inputs, cov_function = cov_function, 
                          X=X, likelihood=likelihood, y=y, cov_fct_shape = nu, 
                          params = list(optimizer_cov = optimizer_cov, 
                                        use_nesterov_acc = use_nesterov_acc, 
                                        trace = trace, init_cov_pars = init_cov_pars,
                                        delta_rel_conv = delta_rel_conv, maxit = maxit))
  }
  else if (gp_approx == "vecchia"){
    gp_model = fitGPModel(gp_coords = inputs, cov_function = cov_function, 
                          X=X, likelihood=likelihood, y=y, cov_fct_shape = nu,
                          gp_approx = gp_approx, matrix_inversion_method = "iterative",
                          params = list(optimizer_cov = optimizer_cov, 
                                        use_nesterov_acc = use_nesterov_acc, 
                                        trace = trace, init_cov_pars = init_cov_pars,
                                        delta_rel_conv = delta_rel_conv, maxit = maxit))
  }
  else {
    gp_model = fitGPModel(gp_coords = inputs, cov_function = cov_function, 
                          X=X, likelihood=likelihood, y=y, cov_fct_shape = nu,
                          gp_approx = gp_approx,
                          params = list(optimizer_cov = optimizer_cov, 
                                        use_nesterov_acc = use_nesterov_acc, 
                                        trace = trace, init_cov_pars = init_cov_pars,
                                        delta_rel_conv = delta_rel_conv, maxit = maxit))
  }
  
  end.time <- Sys.time()
  time.taken <- round(as.numeric(difftime(time1 = end.time, time2 = start.time, units = "secs")), 3)
  return(c(gp_model$get_cov_pars(), 
           nll = gp_model$get_current_neg_log_likelihood(),
           time_converge = time.taken, 
           convergence = gp_model$get_num_optim_iter(), 
           coef = gp_model$get_coef()))
  
}