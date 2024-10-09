source('generate_fit_matern.R')
source('Branin.R')
source('Borehole.R')
source('Piston.R')

gp_approx_s = c("vecchia", "fitc_stable")

n <- 20000

nu_s <- c(0.5, 1.5, 2.5, Inf) 

alg_names <- c("Gradient Descent", "GD with Nesterov", "Fisher scoring",
               "Nelder-Mead", "Newton", "LBFGS")  #six algorithms

reps = 10
results = array(dim = c(6, length(alg_names), length(nu_s), 2, reps))
dimnames(results)[[1]] <- c('error_var', 'GP_var', 'range', 'NLL', 'time', 'convergence')
dimnames(results)[[2]] <- alg_names
dimnames(results)[[3]] <- nu_s
dimnames(results)[[4]] <- gp_approx_s

set.seed(2024)
for(rep in 1:reps){
  ######## Branin function ########
  ### dimensions
  # d <- 2
  # sigma2_error <- 260 # variance of error (=nugget) = 1/10 * var(y)
  # 
  # coords <- matrix(c(runif(n, -5, 10), runif(n, 0, 15)), ncol=2)
  # # y <- branin(xx = coords, nugget = 0)
  # # var(y)
  # y <- branin(xx = coords, nugget = sigma2_error)
  # y <- y - mean(y)
  # print(paste0("Branin dataset of round ", rep, " is generated."))
  
  ######## Borehole function #########
  d = 8
  sigma2_error = 205
  
  coords <- matrix(runif(n*d), ncol=d)
  y <- apply(coords, 1, bore) + rnorm(n,0,sqrt(sigma2_error))
  y <- y - mean(y)
  print(paste0("Borehole dataset of round ", rep, " is generated."))
  
  ####### Piston function #######
  # d = 7
  # sigma2_error <- 0.01 ## 1/2*var(y)
  # 
  # coords <- matrix(runif(n*d), ncol=d)
  # # y <- apply(coords, 1, piston)
  # # var(y)
  # y <- apply(coords, 1, piston) + rnorm(n,0,sqrt(sigma2_error))
  # y <- y - mean(y)
  # print(paste0("Piston dataset of round ", rep, " is generated."))
  
  # #### intercept term
  # intercept = matrix(1, nrow = n, ncol = 1)
  
  for (i in 1:length(nu_s)){  
    nu <- nu_s[i]
    init_pars <- gen_init(coords, y, nu)
    
    for (j in 1:1){
      gp_approx <- gp_approx_s[j]
      
      results[,1,i,j,rep] = fit_matern(y, coords, optimizer_cov='gradient_descent',
                                       nu=nu, trace = TRUE, gp_approx = gp_approx,
                                       init_cov_pars = init_pars)[1:6]
      print(paste0('1/',length(alg_names),' algorithms finished!'))
      
      results[,2,i,j,rep] = fit_matern(y, coords, optimizer_cov='gradient_descent',
                                       use_nesterov_acc=TRUE, gp_approx = gp_approx,
                                       nu=nu, trace = TRUE,
                                       init_cov_pars = init_pars)[1:6]
      print(paste0('2/',length(alg_names),' algorithms finished!'))
      
      results[,3,i,j,rep] = fit_matern(y, coords, optimizer_cov='fisher_scoring',
                                       nu=nu, trace = TRUE, gp_approx = gp_approx,
                                       init_cov_pars = init_pars)[1:6]
      print(paste0('3/',length(alg_names),' algorithms finished!'))
      
      results[,4,i,j,rep] = fit_matern(y, coords, optimizer_cov='nelder_mead',
                                       nu=nu, trace = TRUE, gp_approx = gp_approx,
                                       init_cov_pars = init_pars)[1:6]
      print(paste0('4/',length(alg_names),' algorithms finished!'))
      
      results[,5,i,j,rep] = fit_matern(y, coords, optimizer_cov='newton',
                                       nu=nu, trace = TRUE, gp_approx = gp_approx,
                                       init_cov_pars = init_pars)[1:6]
      print(paste0('5/',length(alg_names),' algorithms finished!'))
      
      results[,6,i,j,rep] = fit_matern(y, coords, optimizer_cov='lbfgs',
                                       nu=nu, trace = TRUE, gp_approx = gp_approx,
                                       init_cov_pars = init_pars)[1:6]
      print(paste0('6/',length(alg_names),' algorithms with nu=', nu, ' and ', gp_approx, ' approximation finished!'))
    }
    print(paste0("Round ", rep, ", experiments with nu = ",nu," finished..."))
  }
}

#### save results locally

# save(d, n, nu_s, alg_names, results, file = 'results/Branin_20000_0604.RData')
# save(d, n, nu_s, alg_names, results, file = 'results/Borehole_20000_0604.RData')
# save(d, n, nu_s, alg_names, results, file = 'results/Piston_20000_0604.RData')

##### load saved results

# load('results/Branin_20000_0604.RData', verbose = TRUE)
# load('results/Borehole_20000_0604.RData', verbose = TRUE)
# load('results/Piston_20000_0604.RData', verbose = TRUE)


####### check convergence #######
source('check_convergence.R')
converge <- check_converge_10000(results["NLL",,,], alg_names)
for(i in 1:dim(converge)[2]){
  for(j in 1:dim(converge)[1]){
    if(any(!converge[j,i,], na.rm = TRUE)){
      print(paste0(alg_names[j]," algorithm fails to converge in repetition ",
                   toString(which(!converge[j,i,]))," when \nu=",nu_s[i]))
    }
  }
}


####### plots #######

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
  ggplot(data = times, 
         aes(x = smoothness, y = mean_time, group = algorithm, color = algorithm)) + 
  geom_line(position=position_dodge(0.2)) +
  geom_errorbar(aes(ymin=min_time, ymax=max_time), width=.3, position=position_dodge(0.2), alpha=.5) +
  geom_point(position=position_dodge(0.2)) +
  scale_color_manual(values = colors, limits = alg_names, drop = FALSE) +
  scale_y_log10(breaks = c(0.5,seq(1,5,1), 1.5,7, 10,20,30,50,70,100,200,300,500, 150,1000, 1500,2000,3000,4000)) +
  # coord_cartesian(ylim = c(1, 1000)) +
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
  scale_color_manual(values = colors, drop = FALSE) +
  scale_y_log10(breaks = c(seq(6,10,2),15,20,30,50,75,seq(100,400,100),150, 500,700,1000)) +
  # coord_cartesian(ylim = c(6, 1000)) +
  labs(x="smoothness", y = "number of iterations") +
  theme_pubclean()

ggarrange(plot, plot_num, ncol=2, common.legend = TRUE)






