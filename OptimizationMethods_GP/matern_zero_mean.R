##### comparison on simulated GP with matern covariance #####

source('generate_fit_matern.R')

### load packages for experiments
library(dplyr)

### dimensions
d <- 2
sigma2 <- 1 # marginal variance

#### settings of true parameters
n_s <- c(200, 500, 1000, 2000)  #sample size
rho_s <- c(0, 0.01, 0.05, 0.1, 0.5, sqrt(d))  #range 
stn_s <- c(0.1, 0.5, 1, 2, 5, 10)  #signal-to-noise ratio
nu_s <- c(0.5, 1.5, 2.5, Inf)  #smoothness

alg_names <- c("lbfgs_not_profile_out_nugget", "lbfgs",
               "L-BFGS(optim)", "L-BFGS-B(optim)", "BFGS(optim)")  #five bfgs algorithms
# alg_names <- c("Gradient Descent", "GD with Nesterov", "Fisher-scoring",
#                "Nelder-Mead", "Newton", "LBFGS")  #six algorithms

comb <- data.frame() #combinations of parameters
comb <- rbind(comb, data.frame(sample_size = n_s, range = rep(0.1, length(n_s)), 
                               signal_to_noise = rep(2, length(n_s)), smoothness = rep(1.5, length(n_s))))
comb <- rbind(comb, data.frame(sample_size = rep(500, length(rho_s)), range = rho_s, 
                               signal_to_noise = rep(2, length(rho_s)), smoothness = rep(1.5, length(rho_s))))
comb <- rbind(comb, data.frame(sample_size = rep(500, length(stn_s)), range = rep(0.1, length(stn_s)),
                               signal_to_noise = stn_s, smoothness = rep(1.5, length(stn_s))))
comb <- rbind(comb, data.frame(sample_size = rep(500, length(nu_s)), range = rep(0.1, length(nu_s)), 
                               signal_to_noise = rep(2, length(nu_s)), smoothness = nu_s))
comb <- distinct(comb)



##### simulation for all parameter settings #####

reps <- 10

## Initial a result table 
results <- array(dim = c(6, length(alg_names), dim(comb)[1], reps))
dimnames(results)[[1]] <- c('error_var', 'GP_var', 'range', 'NLL', 'time', 'convergence')
dimnames(results)[[2]] <- alg_names

set.seed(2024)
for(rep in 1:reps){
  for(i in 1:dim(comb)[1]){
    n <- as.numeric(comb[i,][1])
    rho <- as.numeric(comb[i,][2])
    stn <- as.numeric(comb[i,][3])
    sigma2_error <- sigma2/stn
    nu <- as.numeric(comb[i,][4])
    
    coords <- matrix(runif(n*d), ncol=d)
    y = gen_data(coords, n, rho, sigma2, sigma2_error, nu)
    
    ###### 1. default initial values 
    init_pars <- gen_init(coords, y, nu)
    
    ###### 2. true value as starting points
    
    # ## to make sure the data sample is the same in all initialization strategies
    # default_init <- gen_init(coords, y, nu)
    # 
    # ## set initial values as true parameters
    # init_pars <- c(sigma2_error, sigma2, rho)
    # if(rho == 0){
    #   init_pars <- c(sigma2_error, sigma2, 0.01)
    # }
    
    ###### 3. distant points as starting points
    
    # ## to make sure the data sample is the same in all initialization strategies
    # default_init <- gen_init(coords, y, nu)
    #
    # ## set as unreasonable (large) values
    # init_pars <- 10*default_init
    
    ####### fitting models
    
    # results[,1,i,rep] = fit_matern(y, coords, optimizer_cov='gradient_descent',
    #                                init_cov_pars = init_pars, nu=nu, trace = FALSE)[1:6]
    # print(paste0('1/',length(alg_names),' algorithms finished!'))
    # 
    # results[,2,i,rep] = fit_matern(y, coords, optimizer_cov='gradient_descent',
    #                                use_nesterov_acc=TRUE,
    #                                init_cov_pars = init_pars, nu=nu, trace = FALSE)[1:6]
    # print(paste0('2/',length(alg_names),' algorithms fini shed!'))
    # 
    # results[,3,i,rep] = fit_matern(y, coords, optimizer_cov='fisher_scoring',
    #                                init_cov_pars = init_pars, nu=nu, trace = FALSE)[1:6]
    # print(paste0('3/',length(alg_names),' algorithms finished!'))
    # 
    # results[,4,i,rep] = fit_matern(y, coords, optimizer_cov='nelder_mead',
    #                                init_cov_pars = init_pars, nu=nu, trace = FALSE)[1:6]
    # print(paste0('4/',length(alg_names),' algorithms finished!'))
    # 
    # # if(i!=6){
    # results[,5,i,rep] = fit_matern(y, coords, optimizer_cov='newton',
    #                                init_cov_pars = init_pars, nu=nu, trace = FALSE)[1:6]
    # print(paste0('5/',length(alg_names),' algorithms finished!'))
    # # }
    # 
    # results[,6,i,rep] = fit_matern(y, coords, optimizer_cov='lbfgs',
    #                                init_cov_pars = init_pars, nu=nu, trace = FALSE)[1:6]
    # print(paste0('6/',length(alg_names),' algorithms finished!'))
    
    results[,1,i,rep] = fit_matern(y, coords, optimizer_cov='lbfgs_not_profile_out_nugget',
                                   init_cov_pars = init_pars, nu=nu, trace = FALSE)[1:6]
    print(paste0('1/',length(alg_names),' algorithms finished!'))
    
    results[,2,i,rep] = fit_matern(y, coords, optimizer_cov='lbfgs',
                                   init_cov_pars = init_pars, nu=nu, trace = FALSE)[1:6]
    print(paste0('2/',length(alg_names),' algorithms finished!'))
    
    results[,3,i,rep] = fit_matern_optim(y, coords, optimizer='L-BFGS',
                                         init_cov_pars = init_pars, nu=nu, trace = 0)
    print(paste0('3/',length(alg_names),' algorithms finished!'))
    
    results[,4,i,rep] = fit_matern_optim(y, coords, optimizer='L-BFGS-B',
                                         init_cov_pars = init_pars, nu=nu, trace = 0)
    print(paste0('4/',length(alg_names),' algorithms finished!'))
    
    results[,5,i,rep] = fit_matern_optim(y, coords, optimizer='BFGS',
                                         init_cov_pars = init_pars, nu=nu, trace = 0)
    print(paste0('5/',length(alg_names),' algorithms finished!'))
    
    print(paste0("Round ",rep," with comb ",i," finished!"))
  }
  print(paste0("Round ",rep," finished!"))
}

# save(d, comb, sigma2, alg_names, results, file='results/matern_bfgs_0606.RData')


####### check convergence #######
source('check_convergence.R')

converge <- check_converge(results["NLL",,,], alg_names)
for(i in 1:dim(converge)[2]){
  for(j in 1:dim(converge)[1]){
    if(any(!converge[j,i,], na.rm = TRUE)){
      print(paste0(alg_names[j]," algorithm fails to converge in repetition ",toString(which(!converge[j,i,]))," of the combination ",i))
    }
  }
}



####### large sample size (>2000) experiment #######
gp_approx_s = c("vecchia", "fitc_stable")
n_l = c(5000, 20000)

# alg_names <- c("lbfgs_not_profile_out_nugget", "lbfgs",
#                "L-BFGS(optim)", "L-BFGS-B(optim)", "BFGS(optim)")  #five bfgs algorithms
alg_names <- c("Gradient Descent", "GD with Nesterov", "Fisher-scoring",
               "Nelder-Mead", "Newton", "LBFGS")  #six algorithms

d <- 2
rho <- 0.1
stn <- 2
sigma2 <- 1
sigma2_error <- sigma2/stn
nu <- 1.5
reps <- 10

results <- array(dim = c(6, length(alg_names), 2, 2, reps))
dimnames(results)[[1]] <- c('error_var', 'GP_var', 'range', 'NLL', 'time', 'convergence')
dimnames(results)[[2]] <- alg_names
dimnames(results)[[3]] <- n_l
dimnames(results)[[4]] <- gp_approx_s


set.seed(2024)
for(rep in 1:reps){
  for(i in 1:2){
    n <- n_l[i]
    
    coords <- matrix(runif(n*d), ncol=d)
    y = gen_data(coords, n, rho, sigma2, sigma2_error, nu)
    
    ###### default initial values 
    init_pars <- gen_init(coords, y, nu)
    
    for(j in 1:2){
      gp_approx <- gp_approx_s[j]
      
      results[,1,i,j,rep] = fit_matern(y, coords, optimizer_cov='gradient_descent',
                                       init_cov_pars = init_pars, nu=nu, trace = FALSE,
                                       gp_approx = gp_approx)[1:6]
      print(paste0('1/',length(alg_names),' algorithms finished!'))
      
      results[,2,i,j,rep] = fit_matern(y, coords, optimizer_cov='gradient_descent',
                                       use_nesterov_acc=TRUE,
                                       init_cov_pars = init_pars, nu=nu, trace = FALSE,
                                       gp_approx = gp_approx)[1:6]
      print(paste0('2/',length(alg_names),' algorithms finished!'))
      
      results[,3,i,j,rep] = fit_matern(y, coords, optimizer_cov='fisher_scoring',
                                       init_cov_pars = init_pars, nu=nu, trace = FALSE,
                                       gp_approx = gp_approx)[1:6]
      print(paste0('3/',length(alg_names),' algorithms finished!'))
      
      results[,4,i,j,rep] = fit_matern(y, coords, optimizer_cov='nelder_mead',
                                       init_cov_pars = init_pars, nu=nu, trace = FALSE,
                                       gp_approx = gp_approx)[1:6]
      print(paste0('4/',length(alg_names),' algorithms finished!'))
      
      results[,5,i,j,rep] = fit_matern(y, coords, optimizer_cov='newton',
                                       init_cov_pars = init_pars, nu=nu, trace = FALSE,
                                       gp_approx = gp_approx)[1:6]
      print(paste0('5/',length(alg_names),' algorithms finished!'))
      
      results[,6,i,j,rep] = fit_matern(y, coords, optimizer_cov='lbfgs',
                                       init_cov_pars = init_pars, nu=nu, trace = FALSE,
                                       gp_approx = gp_approx)[1:6]
      print(paste0('6/',length(alg_names),' algorithms finished!'))
      
      
      # results[,1,i,j,rep] = fit_matern(y, coords, optimizer_cov='lbfgs_not_profile_out_nugget',
      #                                init_cov_pars = init_pars, nu=nu, trace = FALSE,
      #                                gp_approx = gp_approx)[1:6]
      # print(paste0('1/',length(alg_names),' algorithms finished!'))
      # 
      # results[,2,i,j,rep] = fit_matern(y, coords, optimizer_cov='lbfgs',
      #                                init_cov_pars = init_pars, nu=nu, trace = FALSE,
      #                                gp_approx = gp_approx)[1:6]
      # print(paste0('2/',length(alg_names),' algorithms finished!'))
      # 
      # results[,3,i,j,rep] = fit_matern_optim(y, coords, optimizer='L-BFGS',
      #                                      init_cov_pars = init_pars, nu=nu, trace = 0,
      #                                      gp_approx = gp_approx)
      # print(paste0('3/',length(alg_names),' algorithms finished!'))
      # 
      # results[,4,i,j,rep] = fit_matern_optim(y, coords, optimizer='L-BFGS-B',
      #                                      init_cov_pars = init_pars, nu=nu, trace = 0,
      #                                      gp_approx = gp_approx)
      # print(paste0('4/',length(alg_names),' algorithms finished!'))
      # 
      # results[,5,i,j,rep] = fit_matern_optim(y, coords, optimizer='BFGS',
      #                                      init_cov_pars = init_pars, nu=nu, trace = 0,
      #                                      gp_approx = gp_approx)
      # print(paste0('5/',length(alg_names),' algorithms finished!'))
      
      print(paste0("Round ",rep," with n=",n," and ", gp_approx, " approximation finished!"))
    }
  }
  print(paste0("Round ",rep," finished!"))
}

# save(d, alg_names, results, file='results/bfgs_approx_0606.RData')
# save(d, alg_names, results, file='results/matern_approx_0606.RData')

##### load saved results
# load(file='results/bfgs_approx_0528.RData', verbose = TRUE)



####### check convergence with approximation #######
source('check_convergence.R')
for(rep in 1:reps){
  converge <- check_converge(results["NLL",,,,rep], alg_names)
  for(i in 1:dim(converge)[2]){
    for(j in 1:dim(converge)[3]){
      if(any(!converge[,i,j], na.rm = TRUE)){
        print(paste0(alg_names[which(!converge[,i,j])]," algorithm fails to converge in repetition ", 
                     rep, " of n=", n_l[i], " with ", gp_approx_s[j], " approximation"))
      }
    }
  }
}



