######## Experiments of model mis-specification based on GP-matern ########

## source the function code of generating samples and fitting models
source('generate_fit_matern.R')


### set parameters
d <- 2
sigma2 <- 1 # marginal variance
rho=0.1
sigma2_error = 0.5

#### generating smoothness
gen_nu <- 0.5
# gen_nu <- Inf

n_s <- c(200, 500, 1000, 2000)  #sample size
nu_s <- c(0.5, 1.5, 2.5, Inf)  #smoothness

alg_names <- c("Gradient Descent", "GD with Nesterov", "Fisher scoring",
               "Nelder-Mead", "Newton", "LBFGS")  #six algorithms
reps=10
results = array(dim = c(6,length(alg_names), 3, length(n_s), reps))


#### fitting models
set.seed(2024)
for(rep in 1:reps){
  for(i in 1:length(n_s)){
    n <- n_s[i]
    coords <- matrix(runif(n*d), ncol=d)
    y <- gen_data(coords, n, rho, sigma2, sigma2_error, gen_nu)
    print(paste0("dataset with ", n, " points in rep ",rep ," is generated."))
    j <- 1
    for(nu in nu_s){
      if(nu != gen_nu){
        init_pars <- gen_init(coords, y, nu)
        if(n==5000 && nu == Inf){
          results[,1,j,i,rep] = fit_matern(y, coords, optimizer_cov='gradient_descent',
                                           init_cov_pars = init_pars, nu=nu, trace = TRUE)[1:6]
          print(paste0('1/',length(alg_names),' algorithms finished!'))
          
          results[,2,j,i,rep] = fit_matern(y, coords, optimizer_cov='gradient_descent',
                                           use_nesterov_acc=TRUE,
                                           init_cov_pars = init_pars, nu=nu, trace = TRUE)[1:6]
          print(paste0('2/',length(alg_names),' algorithms finished!'))
          
          results[,3,j,i,rep] = fit_matern(y, coords, optimizer_cov='fisher_scoring',
                                           init_cov_pars = init_pars, nu=nu, trace = TRUE)[1:6]
          print(paste0('3/',length(alg_names),' algorithms finished!'))
          
          results[,4,j,i,rep] = fit_matern(y, coords, optimizer_cov='nelder_mead',
                                           init_cov_pars = init_pars, nu=nu, trace = TRUE)[1:6]
          print(paste0('4/',length(alg_names),' algorithms finished!'))
          
          results[,5,j,i,rep] = fit_matern(y, coords, optimizer_cov='newton',
                                           init_cov_pars = init_pars, nu=nu, trace = TRUE)[1:6]
          print(paste0('5/',length(alg_names),' algorithms finished!'))
          
          results[,6,j,i,rep] = fit_matern(y, coords, optimizer_cov='lbfgs',
                                           init_cov_pars = init_pars, nu=nu, trace = TRUE)[1:6]
          print(paste0('6/',length(alg_names),' algorithms finished!'))
        }
        j = j + 1
        print(paste0("Repetition ",rep,", n=",n,", nu=",nu," finished!"))
      }
    }
  }
}
dimnames(results)[[1]] <- c('error_var', 'GP_var', 'range', 'NLL', 'time', 'convergence')
dimnames(results)[[2]] <- alg_names
dimnames(results)[[3]] <- c(1.5,2.5,Inf)
# dimnames(results)[[3]] <- paste0("nu=",c(0.5,1.5,2.5))
dimnames(results)[[4]] <- paste0("n=",n_s)


# save(gen_nu, alg_names, results, file = "results/matern_misspec_nu=0.5.RData")
# save(gen_nu, alg_names, results, file = "results/matern_misspec_nu=Inf.RData")

# load("results/matern_misspec_nu=0.5.RData", verbose = TRUE)
# load("results/matern_misspec_nu=Inf.RData", verbose = TRUE)

####### check convergence #######
source('check_convergence.R')
for(rep in 1:reps){
  converge <- check_converge(results["NLL",,,,rep], alg_names)
  for(i in 1:dim(converge)[2]){
    for(j in 1:dim(converge)[1]){
      if(any(!converge[j,i,], na.rm = TRUE)){
        print(paste0(alg_names[j]," algorithm fails to converge when n=",n_s[which(!converge[j,i,])],
                     " and nu=",dimnames(results)[[3]][i], " in repetition ", rep))
        
      }
    }
  }
}


####### make plots #######
library(ggplot2)
library(ggpubr)
library(viridis)

times = results["time",,,,]

time_results = as.data.frame(cbind(algorithm=rep(alg_names, 12), 
                                   sample_size = rep(n_s[1:4], each=length(alg_names)*3),
                                   smoothness = rep(rep(dimnames(results)[[3]],each=length(alg_names)),length(n_s)-1)))

# smoothness = rep(rep(c(0.5,1.5,2.5),each=length(alg_names)),length(n_s))))
# row.names(time_results) = 1:90
time_results$mean_time = c(apply(times, 1:3, mean))
time_results$min_time = c(apply(times, 1:3, min))
time_results$max_time = c(apply(times, 1:3, max))

time_results$algorithm <- as.factor(time_results$algorithm)
time_results$sample_size <- as.factor(as.numeric(time_results$sample_size))
time_results$smoothness <- as.factor(time_results$smoothness)

colors = viridis_pal(option = "H")(6)

plot <- 
  ggplot(data = time_results, aes(x = sample_size, y = mean_time, group = interaction(algorithm, smoothness), 
                                  color = algorithm, linetype = smoothness)) +
  geom_line(position=position_dodge(0.2)) +
  geom_point(position=position_dodge(0.2)) +
  geom_errorbar(aes(ymin=min_time, ymax=max_time), 
                width=.3, position=position_dodge(0.2), alpha=.5) +
  # scale_color_viridis(discrete = TRUE, option = "H") +
  scale_color_manual(values = colors, limits = alg_names, drop = FALSE) +
  scale_y_log10(breaks = c(0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,seq(1000,5000,1000),10000,40000)) +
  # coord_cartesian(ylim = c(0.005, 1800)) +
  labs(x="sample size", y = "time(s)") +
  theme_pubclean()

pdf(file=paste0('plots_l/matern_misspec_nu=0.5.pdf'),width=10.0,height=7.0)
par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.3,.1)) # bltr
plot
dev.off()

pdf(file=paste0('plots_l/matern_misspec_nu=Inf.pdf'),width=10.0,height=7.0)
par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.3,.1)) # bltr
plot
dev.off()



####### large sample size (>2000) experiment #######
gp_approx_s = c("vecchia", "fitc_stable")

d <- 2
sigma2 <- 1 # marginal variance
rho=0.1
sigma2_error = 0.5
# gen_nu <- 0.5
gen_nu <- Inf

n <- 5000
nu_s <- c(0.5, 1.5, 2.5, Inf)  #smoothness

reps <- 10
results = array(dim = c(6,length(alg_names), 3, 2, reps))
dimnames(results)[[1]] <- c('error_var', 'GP_var', 'range', 'NLL', 'time', 'convergence')
dimnames(results)[[2]] <- alg_names
# dimnames(results)[[3]] <- c(1.5, 2.5, Inf) ## gen_nu=0.5
dimnames(results)[[3]] <- c(0.5, 1.5, 2.5) ## gen_nu=Inf
dimnames(results)[[4]] <- gp_approx_s

set.seed(2024)
for(rep in 1:reps){
  coords <- matrix(runif(n*d), ncol=d)
  y <- gen_data(coords, n, rho, sigma2, sigma2_error, gen_nu)
  print(paste0("dataset with ", n, " points in rep ",rep ," is generated."))
  j <- 1
  for(nu in nu_s){
    if(nu != gen_nu){
      init_pars <- gen_init(coords, y, nu)
      for(i in 1:2){
        gp_approx = gp_approx_s[i]
        results[,1,j,i,rep] = fit_matern(y, coords, optimizer_cov='gradient_descent',
                                         init_cov_pars = init_pars, nu=nu, trace = FALSE,
                                         gp_approx = gp_approx)[1:6]
        print(paste0('1/',length(alg_names),' algorithms finished!'))
        
        results[,2,j,i,rep] = fit_matern(y, coords, optimizer_cov='gradient_descent',
                                         use_nesterov_acc=TRUE,
                                         init_cov_pars = init_pars, nu=nu, trace = FALSE,
                                         gp_approx = gp_approx)[1:6]
        print(paste0('2/',length(alg_names),' algorithms finished!'))
        
        results[,3,j,i,rep] = fit_matern(y, coords, optimizer_cov='fisher_scoring',
                                         init_cov_pars = init_pars, nu=nu, trace = FALSE,
                                         gp_approx = gp_approx)[1:6]
        print(paste0('3/',length(alg_names),' algorithms finished!'))
        
        results[,4,j,i,rep] = fit_matern(y, coords, optimizer_cov='nelder_mead',
                                         init_cov_pars = init_pars, nu=nu, trace = FALSE,
                                         gp_approx = gp_approx)[1:6]
        print(paste0('4/',length(alg_names),' algorithms finished!'))
        
        results[,5,j,i,rep] = fit_matern(y, coords, optimizer_cov='newton',
                                         init_cov_pars = init_pars, nu=nu, trace = FALSE,
                                         gp_approx = gp_approx)[1:6]
        print(paste0('5/',length(alg_names),' algorithms finished!'))
        
        results[,6,j,i,rep] = fit_matern(y, coords, optimizer_cov='lbfgs',
                                         init_cov_pars = init_pars, nu=nu, trace = FALSE,
                                         gp_approx = gp_approx)[1:6]
        print(paste0('6/',length(alg_names),' algorithms finished!'))
      }
      j = j + 1
      print(paste0("Repetition ",rep,", n=5000, nu=",nu," with approximation finished!"))
    }
  }
}

# save(d, alg_names, results, file='results/misspec0.5_approx.RData')
# save(d, alg_names, results, file='results/misspecinf_approx.RData')

load('results/misspec0.5_approx.RData', verbose = TRUE)
load('results/misspecinf_approx.RData', verbose = TRUE)

####### check convergence with approximation #######
source('check_convergence.R')
for(rep in 1:reps){
  converge <- check_converge(results["NLL",,,,rep], alg_names)
  for(i in 1:dim(converge)[2]){
    for(j in 1:dim(converge)[3]){
      if(any(!converge[,i,j], na.rm = TRUE)){
        print(paste0(alg_names[which(!converge[,i,j])]," algorithm fails to converge in repetition ", 
                     rep, " of nu=", dimnames(results)[[3]][i], " with ", gp_approx_s[j], " approximation"))
      }
    }
  }
}

####### make plots with approximation #######
library(ggplot2)
library(ggpubr)
library(viridis)

times = results["time",,,,]

time_results = data.frame(algorithm=rep(alg_names, 6), smoothness = rep(dimnames(results)[[3]], times = 2, each = 6),
                          gp_approx = rep(gp_approx_s, each = 18))

time_results$mean_time = c(apply(times, 1:3, mean))
time_results$min_time = c(apply(times, 1:3, min))
time_results$max_time = c(apply(times, 1:3, max))

time_results$algorithm <- as.factor(time_results$algorithm)
time_results$gp_approx <- as.factor(time_results$gp_approx)

colors = viridis_pal(option = "H")(6)

plot1 <- 
  ggplot(data = time_results, aes(x = gp_approx, y = mean_time, group = interaction(algorithm, smoothness), 
                                  color = algorithm, linetype = smoothness)) +
  geom_line(position=position_dodge(0.2)) +
  geom_point(position=position_dodge(0.2)) +
  geom_errorbar(aes(ymin=min_time, ymax=max_time), 
                width=.3, position=position_dodge(0.2), alpha=.5) +
  scale_color_manual(values = colors, limits = alg_names, drop = FALSE) +
  scale_y_log10(breaks = c(0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000,2000,4000,10000,40000)) +
  # coord_cartesian(ylim = c(0.2, 2500)) +
  labs(x="gp_approximation", y = "time(s)") +
  theme_pubclean()+ 
  guides(color = guide_legend(nrow = 2))

pdf(file=paste0('plots_l/misspec_0.5_approx.pdf'),width=10,height=5)
par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.3,.1)) # bltr
ggarrange(plot, plot1, ncol=2, common.legend = TRUE)
dev.off()



pdf(file=paste0('plots_l/misspec_inf_approx.pdf'),width=10,height=5)
par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.3,.1)) # bltr
ggarrange(plot, plot1, ncol=2, common.legend = TRUE)
dev.off()





