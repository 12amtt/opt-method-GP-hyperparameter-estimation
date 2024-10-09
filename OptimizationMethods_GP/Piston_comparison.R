####### comparison on the Piston function   #######

source('generate_fit_matern.R')
source('Piston.R')


####### repeated Experiments #######

###### dimensions
d = 7
n = 1000
nu_s <- c(0.5, 1.5, 2.5, Inf) 
alg_names <- c("Gradient Descent", "GD with Nesterov", "Fisher scoring",
               "Nelder-Mead", "newton", "lbfgs")  #six algorithms

sigma2_error <- 0.01

######## comparison of different algorithms (rel_tol=1e-8)
reps = 10

results <- array(dim = c(6, length(alg_names), length(nu_s),reps))

set.seed(2024)
for(rep in 1:reps){
  coords <- matrix(runif(n*d), ncol=d)
  # y <- apply(coords, 1, piston)
  # var(y)
  y <- apply(coords, 1, piston) + rnorm(n,0,sqrt(sigma2_error))
  y <- y - mean(y)
  print("Dataset is generated.")
  
  for (i in 1:length(nu_s)){
    
    nu <- nu_s[i]
    init_pars <- gen_init(coords, y, nu)
    
    results[,1,i,rep] = fit_matern(y, coords, optimizer_cov='gradient_descent',
                                   init_cov_pars = init_pars, nu=nu, trace = FALSE)[1:6]
    print(paste0('1/',length(alg_names),' algorithms finished!'))
    
    results[,2,i,rep] = fit_matern(y, coords, optimizer_cov='gradient_descent',
                                   use_nesterov_acc=TRUE,
                                   init_cov_pars = init_pars, nu=nu, trace = FALSE)[1:6]
    print(paste0('2/',length(alg_names),' algorithms finished!'))
    
    results[,3,i,rep] = fit_matern(y, coords, optimizer_cov='fisher_scoring',
                                   init_cov_pars = init_pars, nu=nu, trace = FALSE)[1:6]
    print(paste0('3/',length(alg_names),' algorithms finished!'))
    
    results[,4,i,rep] = fit_matern(y, coords, optimizer_cov='nelder_mead',
                                   init_cov_pars = init_pars, nu=nu, trace = FALSE)[1:6]
    print(paste0('4/',length(alg_names),' algorithms finished!'))
    
    results[,5,i,rep] = fit_matern(y, coords, optimizer_cov='newton',
                                   init_cov_pars = init_pars, nu=nu, trace = FALSE)[1:6]
    print(paste0('5/',length(alg_names),' algorithms finished!'))
    
    results[,6,i,rep] = fit_matern(y, coords, optimizer_cov='lbfgs',
                                   init_cov_pars = init_pars, nu=nu, trace = FALSE)[1:6]
    print(paste0('6/',length(alg_names),' algorithms finished!'))
    
    print(paste0("Round ",rep," with nu = ",nu," finished..."))
  }
  print(paste0("Round ",rep," finished!"))
}

dimnames(results)[[1]] <- c('error_var', 'GP_var', 'range', 'NLL', 'time', 'convergence')
dimnames(results)[[2]] <- alg_names
dimnames(results)[[3]] <- nu_s

####### display optimization results (time, estimated parameters)

#load saved results
load('results/piston_comparison_0426.RData', verbose = TRUE)

####### check convergence #######
source('check_convergence.R')
converge <- check_converge(results["NLL",,,], alg_names)
for(i in 1:dim(converge)[2]){
  for(j in 1:dim(converge)[1]){
    if(any(!converge[j,i,], na.rm = TRUE)){
      print(paste0(alg_names[j]," algorithm fails to converge in repetition ",toString(which(!converge[j,i,]))," when \nu=",nu_s[i]))
    }
  }
}


####### plotting comparison of times #######
### load necessary packages:
library(ggplot2)
library(viridis)
library(ggpubr)


times_mean <- apply(results["time",,,], 1:2, mean)

times <- as.data.frame(as.table(times_mean))
times$min_time = c(apply(results["time",,,], 1:2, min, na.rm = TRUE))
times$max_time = c(apply(results["time",,,], 1:2, max, na.rm = TRUE))
colnames(times) <- c('algorithm', 'smoothness', 'mean_time',"min_time","max_time")

colors = viridis_pal(option = "H")(6)

plot <- 
  ggplot(data = times, aes(x = smoothness, y = mean_time, group = algorithm, color = algorithm)) + 
  geom_line(position=position_dodge(0.2)) +
  geom_errorbar(aes(ymin=min_time, ymax=max_time), width=.3, position=position_dodge(0.2), alpha=.5) +
  geom_point(position=position_dodge(0.2)) +
  scale_color_manual(values = colors, limits = alg_names, drop = FALSE) +
  scale_y_log10(breaks = c(0.5,seq(1,6,1),1.5,7.5,15, seq(10,50,10),100,200,500,700,1000)) +
  # coord_cartesian(ylim = c(1, 500)) +
  labs(x="smoothness", y = "time") +
  theme_pubclean()



########  plot of number of iterations  ########

num_mean <- apply(results["convergence",,,], 1:2, mean)

num_it <- as.data.frame(as.table(num_mean))
num_it$min_num <-  c(apply(results["convergence",,,], 1:2, min, na.rm = TRUE))
num_it$max_num <-  c(apply(results["convergence",,,], 1:2, max, na.rm = TRUE))
colnames(num_it) <- c('algorithm', 'smoothness', 'mean_num', 'min_num', 'max_num')

plot_num <- 
  ggplot(data = num_it, aes(x = smoothness, y = mean_num, group = algorithm, color = algorithm)) + 
  geom_line(position=position_dodge(0.2)) +
  geom_errorbar(aes(ymin=min_num, ymax=max_num), width=.3, position=position_dodge(0.2), alpha=.5) +
  geom_point(position=position_dodge(0.2)) +
  #  scale_color_viridis(discrete = TRUE, option = "H") +  
  scale_color_manual(values = colors, limits = alg_names, drop = FALSE) +
  scale_y_log10(breaks = c(1,seq(2,12,2),15,20,30,40,50,75,seq(100,400,100),150, 500,600,750,1000)) +
  # coord_cartesian(ylim = c(6, 1000)) +
  labs(x="smoothness", y = "number of iterations") +
  theme_pubclean()


pdf(file=paste0('plots_l/Piston_comparison.pdf'),width=10,height=5)
par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.3,.1)) # bltr
ggarrange(plot, plot_num, ncol=2, common.legend = TRUE)
dev.off()

pdf(file=paste0('plots_l/Piston_20000.pdf'),width=10,height=5)
par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.3,.1)) # bltr
ggarrange(plot, plot_num, ncol=2, common.legend = TRUE)
dev.off()





# ggarrange(plot1,plot2,common.legend=TRUE)



####### irrelevant part #######
#### turn off scientific notation for numbers
options(scipen = 999)
#### revert it back to the default
options(scipen = 0, "digits"= 4)

