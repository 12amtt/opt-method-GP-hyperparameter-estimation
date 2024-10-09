
source('generate_fit_matern.R')

library(dplyr)

alg_names <- c("Gradient Descent", "GD with Nesterov",
               "Nelder-Mead", "Newton", "LBFGS") 
#five algorithms - Fisher-scoring is not supported for covariance parameters for non-Gaussian likelihoods


##### parameter settings #####
d=2
n_s <- c(200, 500, 1000, 2000)  #sample size
rho_s <- c(0, 0.01, 0.05, 0.1, 0.5, sqrt(d))  #range 
sigma2_s <- c(0.1, 0.5, 1, 2, 10)
nu_s <- c(0.5, 1.5, 2.5, Inf)  #smoothness

comb_bi <- data.frame() #combinations of parameters
comb_bi <- rbind(comb_bi, data.frame(sample_size = n_s, range = rep(0.1, length(n_s)), 
                                     sigma2 = rep(1, length(n_s)), smoothness = rep(1.5, length(n_s))))
comb_bi <- rbind(comb_bi, data.frame(sample_size = rep(500, length(rho_s)), range = rho_s, 
                                     sigma2 = rep(1, length(rho_s)), smoothness = rep(1.5, length(rho_s))))
comb_bi <- rbind(comb_bi, data.frame(sample_size = rep(500, length(sigma2_s)), range = rep(0.1, length(sigma2_s)),
                                     sigma2 = sigma2_s, smoothness = rep(1.5, length(sigma2_s))))
comb_bi <- rbind(comb_bi, data.frame(sample_size = rep(500, length(nu_s)), range = rep(0.1, length(nu_s)), 
                                     sigma2 = rep(1, length(nu_s)), smoothness = nu_s))
comb_bi <- distinct(comb_bi)

likelihood <- "bernoulli_logit"



reps <- 10
results = array(dim = c(5, length(alg_names), dim(comb_bi)[1], reps))
dimnames(results)[[1]] <- c('GP_var', 'GP_range', 'NLL', 'time', 'convergence')
dimnames(results)[[2]] <- alg_names

set.seed(2024)

for(rep in 1:reps){
  for(i in 1:dim(comb_bi)[1]){
    n <- as.numeric(comb_bi[i,][1])
    rho <- as.numeric(comb_bi[i,][2])
    sigma2_gp <- as.numeric(comb_bi[i,][3])
    nu <- as.numeric(comb_bi[i,][4])
    
    coords <- matrix(runif(n*d), ncol=d)
    gp_eff <- gen_data(coords, n, rho, sigma2_gp, 0, nu)
    # gp_eff <- gp_eff - mean(gp_eff)
    y <- simulate_response_variable(lp=0, rand_eff=gp_eff, likelihood=likelihood)
    
    init_pars <- gen_init(coords, y, nu)[-2]
    
    results[,1,i,rep] = fit_matern_binary(y, coords, optimizer_cov='gradient_descent', 
                                          init_cov_pars = init_pars, 
                                          likelihood = likelihood, nu=nu, trace = FALSE)[1:5]
    print(paste0('1/',length(alg_names),' algorithms finished!'))
    
    results[,2,i,rep] = fit_matern_binary(y, coords, optimizer_cov='gradient_descent',
                                          use_nesterov_acc=TRUE, init_cov_pars = init_pars, 
                                          likelihood = likelihood, nu=nu, trace = FALSE)[1:5]
    print(paste0('2/',length(alg_names),' algorithms finished!'))
    
    results[,3,i,rep] = fit_matern_binary(y, coords, optimizer_cov='nelder_mead',
                                          init_cov_pars = init_pars, 
                                          likelihood = likelihood, nu=nu, trace = FALSE)[1:5]
    print(paste0('4/',length(alg_names),' algorithms finished!'))
    
    results[,4,i,rep] = fit_matern_binary(y, coords, optimizer_cov='newton',
                                          init_cov_pars = init_pars, 
                                          likelihood = likelihood, nu=nu, trace = FALSE)[1:5]
    print(paste0('5/',length(alg_names),' algorithms finished!'))
    
    results[,5,i,rep] = fit_matern_binary(y, coords, optimizer_cov='lbfgs',
                                          init_cov_pars = init_pars, 
                                          likelihood = likelihood, nu=nu, trace = FALSE)[1:5]
    print(paste0('6/',length(alg_names),' algorithms finished!'))
    
    print(paste0("Round ",rep," with comb ",i," finished!"))
  }
  print(paste0("Round ",rep," finished!"))
}


# save(d, comb_bi, alg_names, results, file='results/matern_binary.RData')

load(file='results/matern_binary.RData', verbose = TRUE)


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

## packages (for plotting)
library(ggplot2)
library(viridis)
library(ggpubr)


times = results["time",,,]

time_results = cbind(algorithm=rep(alg_names, 16), comb_bi[rep(1:16, each=length(alg_names)),])

row.names(time_results) = 1:(16*length(alg_names))
time_results$mean_time = c(apply(times, 1:2, mean, na.rm = TRUE))
time_results$min_time = c(apply(times, 1:2, min, na.rm = TRUE))
time_results$max_time = c(apply(times, 1:2, max, na.rm = TRUE))

### convert parameters into factor variables
time_results$sample_size <- as.factor(time_results$sample_size)
time_results$range <- as.factor(time_results$range)
time_results$sigma2 <- as.factor(time_results$sigma2)
time_results$smoothness <- as.factor(time_results$smoothness)

attach(time_results)

### y-axis time is plotted on a log scale
colors = viridis_pal(option = "H")(6)
colors = colors[c(1:2,4:6)]

## remove scientific notation
options(scipen=10000)

remove(sigma2)

plot1 <- 
  ggplot(data = time_results[range==0.1 & sigma2==1 & smoothness==1.5 & alg_names != "Fisher Scoring",], 
         aes(x = sample_size, y = mean_time, group = algorithm, color = algorithm)) + 
  geom_line(position=position_dodge(0.2)) +
  geom_errorbar(aes(ymin=min_time, ymax=max_time), 
                width=.3, position=position_dodge(0.2), alpha=.5) +
  geom_point(position=position_dodge(0.2)) +
  scale_color_manual(values = colors,limits=alg_names, drop = TRUE) +
  scale_y_log10(breaks = c(0.01, 0.02,0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1, 2, 5, 10, 20, 30, 50, 100,200,300,500,seq(1000,3000,1000))) +
  # coord_cartesian(ylim = c(0.01, 50)) +
  labs(x="sample size", y = "time(s)") +
  theme_pubclean()




plot2 <- 
  ggplot(data = time_results[sample_size==500 & sigma2==1 & smoothness==1.5 & alg_names != "Fisher Scoring",], 
         aes(x = range, y = mean_time, group = algorithm, color = algorithm)) + 
  geom_line(alpha=1,position=position_dodge(0.2)) +
  geom_errorbar(alpha=.5, aes(ymin=min_time, ymax=max_time), 
                width=.3, position=position_dodge(0.2)) +
  geom_point(position=position_dodge(0.2)) +
  scale_x_discrete(labels=c('0', '0.01', '0.05', '0.1', '0.5', 'sqrt(d)')) +
  scale_color_manual(values = colors,limits=alg_names, drop = TRUE) +
  labs(x="effective range", y = "time(s)") +
  scale_y_log10(breaks = c(0.05,0.1, 0.2, 0.4,0.6, 1, 1.5, 2, seq(4,10,2),15,20,30,seq(20,60,20),100,seq(200,600,200),seq(1000,5000,1000))) +
  # coord_cartesian(ylim = c(0, 2)) +
  theme_pubclean()


plot3 <- 
  ggplot(data = time_results[sample_size==500 & range==0.1 & smoothness==1.5 & alg_names != "Fisher Scoring",], 
         aes(x = sigma2, y = mean_time, group = algorithm, color = algorithm)) + 
  geom_line(alpha=1,position=position_dodge(0.2)) +
  geom_errorbar(alpha=.5, aes(ymin=min_time, ymax=max_time), 
                width=.3, position=position_dodge(0.2)) +
  geom_point(position=position_dodge(0.2)) +
  scale_color_manual(values = colors,limits=alg_names, drop = TRUE) +
  labs(x="GP variance", y = "time(s)") +
  scale_y_log10(breaks = c(0.05, 0.1, 0.15, 0.2,0.3,0.4,0.6,1,2,4,6,10, seq(20,60,20),seq(100,300,100),500,1000,2000)) +
  # coord_cartesian(ylim = c(0, 2)) +
  theme_pubclean()


plot4 <- 
  ggplot(data = time_results[sample_size==500 & range==0.1 & sigma2==1 & alg_names != "Fisher Scoring",], 
         aes(x = smoothness, y = mean_time, group = algorithm, color = algorithm)) + 
  geom_line(alpha=1,position=position_dodge(0.2)) +
  geom_errorbar(alpha=.3, aes(ymin=min_time, ymax=max_time), 
                width=.3, position=position_dodge(0.2)) +
  geom_point(position=position_dodge(0.2)) +  
  scale_color_manual(values = colors,limits=alg_names, drop = TRUE) +
  labs(x="smoothness", y = "time(s)") +
  # scale_y_log10(breaks = c(0.05,0.1,0.2,0.3,0.5,seq(0.6,1,0.2),0.15,0.4,1, seq(2,6,2), 10,seq(20,60,20),100,200,300,500,1000,1500)) +
  scale_y_continuous(breaks = seq(0,1.4,0.2)) +
  coord_cartesian(ylim = c(0, 1.4)) +
  theme_pubclean()

detach(time_results)

ggarrange(plot1,plot2,plot3,plot4, ncol = 2, nrow = 2, common.legend=TRUE)




##### large sample size (>2000) experiment #####
gp_approx_s = c("vecchia", "fitc_stable")
n_l = c(5000, 20000)

likelihood <- "bernoulli_logit"

d <- 2
rho <- 0.1
stn <- 2
sigma2_gp <- 1
nu <- 1.5
reps <- 10

results <- array(dim = c(5, length(alg_names), 2, 2, reps))
dimnames(results)[[1]] <- c('GP_var', 'GP_range', 'NLL', 'time', 'convergence')
dimnames(results)[[2]] <- alg_names
dimnames(results)[[3]] <- n_l
dimnames(results)[[4]] <- gp_approx_s

set.seed(2024)
for(rep in 1:reps){
  for(i in 1:2){
    n <- n_l[i]
    coords <- matrix(runif(n*d), ncol=d)
    gp_eff <- gen_data(coords, n, rho, sigma2_gp, 0, nu)
    # gp_eff <- gp_eff - mean(gp_eff)
    y <- simulate_response_variable(lp=0, rand_eff=gp_eff, likelihood=likelihood)
    ###### default initial values #######
    init_pars <- gen_init(coords, y, nu)[-2]
    
    for(j in 1:2){
      gp_approx <- gp_approx_s[j]
      results[,1,i,j,rep] = fit_matern_binary(y, coords, optimizer_cov='gradient_descent',
                                              init_cov_pars = init_pars, nu=nu, trace = TRUE,
                                              likelihood = likelihood, gp_approx = gp_approx)[1:5]
      print(paste0('1/',length(alg_names),' algorithms finished!'))
      
      results[,2,i,j,rep] = fit_matern_binary(y, coords, optimizer_cov='gradient_descent',
                                              use_nesterov_acc=TRUE, likelihood = likelihood, 
                                              init_cov_pars = init_pars, nu=nu, trace = TRUE,
                                              gp_approx = gp_approx)[1:5]
      print(paste0('2/',length(alg_names),' algorithms finished!'))
      
      results[,3,i,j,rep] = fit_matern_binary(y, coords, optimizer_cov='nelder_mead',
                                              init_cov_pars = init_pars, nu=nu, trace = TRUE,
                                              likelihood = likelihood, gp_approx = gp_approx)[1:5]
      print(paste0('3/',length(alg_names),' algorithms finished!'))
      
      results[,4,i,j,rep] = fit_matern_binary(y, coords, optimizer_cov='newton',
                                              init_cov_pars = init_pars, nu=nu, trace = TRUE,
                                              likelihood = likelihood, gp_approx = gp_approx)[1:5]
      print(paste0('4/',length(alg_names),' algorithms finished!'))
      
      results[,5,i,j,rep] = fit_matern_binary(y, coords, optimizer_cov='lbfgs',
                                              init_cov_pars = init_pars, nu=nu, trace = TRUE,
                                              likelihood = likelihood, gp_approx = gp_approx)[1:5]
      print(paste0('5/',length(alg_names),' algorithms finished!'))
      
      print(paste0("Round ",rep," with n=",n," and ", gp_approx, " approximation finished!"))
    }
  }
  print(paste0("Round ",rep," finished!"))
}

# save(d, alg_names, results, file='results/binary_approx.RData')
load('results/binary_approx.RData', verbose = TRUE)

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


