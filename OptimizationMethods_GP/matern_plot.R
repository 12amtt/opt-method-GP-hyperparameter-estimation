####### plot results #######


#### load saved data
# load(file='results/matern_d=2.RData', verbose = TRUE)


## packages (for plotting)
library(ggplot2)
library(viridis)
library(ggpubr)

#### comparison of 6 algorithms


times = results["time",,,]
times = results["time",,,]

time_results = cbind(algorithm=rep(alg_names, 17), comb[rep(1:17, each=length(alg_names)),])

row.names(time_results) = 1:(17*length(alg_names))
time_results$mean_time = c(apply(times, 1:2, mean, na.rm = TRUE))
#time_results$se_time = c(apply(times, 1:2, sd, na.rm = TRUE)/sqrt(reps))
time_results$min_time = c(apply(times, 1:2, min, na.rm = TRUE))
time_results$max_time = c(apply(times, 1:2, max, na.rm = TRUE))

### convert parameters into factor variables
time_results$sample_size <- as.factor(time_results$sample_size)
time_results$range <- as.factor(time_results$range)
time_results$signal_to_noise <- as.factor(time_results$signal_to_noise)
time_results$smoothness <- as.factor(time_results$smoothness)

attach(time_results)

### y-axis time is plotted on a log scale
colors = viridis_pal(option = "H")(6)

# colors = colors[c(1:3,5:6)] #for linear regression term experiments

# colors = viridis_pal(option = "D")(5) # for matern with bfgs plots
# colors[1] = "#440156FF"; colors[2] = "#7A0403FF"; colors[3] = "#61908CFF"

options(scipen=10000)

plot1 <- 
  ggplot(data = time_results[range==0.1 & signal_to_noise==2 & smoothness==1.5,], 
         aes(x = sample_size, y = mean_time, group = algorithm, color = algorithm)) + 
  geom_line(position=position_dodge(0.2)) +
  geom_errorbar(aes(ymin=min_time, ymax=max_time), 
                width=.3, position=position_dodge(0.2), alpha=.5) +
  geom_point(position=position_dodge(0.2)) +
  scale_color_manual(values = colors,limits=alg_names, drop = FALSE) +
  scale_y_log10(breaks = c(0.01,0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100,200,500,1000,3000)) +
  # coord_cartesian(ylim = c(0.01, 50)) +
  labs(x="sample size", y = "time(s)") +
  theme_pubclean()




plot2 <- 
  ggplot(data = time_results[sample_size==500 & signal_to_noise==2 & smoothness==1.5,], 
         aes(x = range, y = mean_time, group = algorithm, color = algorithm)) + 
  geom_line(alpha=1,position=position_dodge(0.2)) +
  geom_errorbar(alpha=.5, aes(ymin=min_time, ymax=max_time), 
                width=.3, position=position_dodge(0.2)) +
  geom_point(position=position_dodge(0.2)) +
  scale_x_discrete(labels=c('0', '0.01', '0.05', '0.1', '0.5', 'sqrt(d)')) +
  scale_color_manual(values = colors,limits=alg_names, drop = FALSE) +
  labs(x="effective range", y = "time(s)") +
  scale_y_log10(breaks = c(0.05,0.1, 0.2, 0.4,0.6, 1,1.5,2,3,5,10,20,50,100,200,500,seq(1000,5000,1000))) +
  # coord_cartesian(ylim = c(0, 2)) +
  theme_pubclean()


plot3 <- 
  ggplot(data = time_results[sample_size==500 & range==0.1 & smoothness==1.5,], 
         aes(x = signal_to_noise, y = mean_time, group = algorithm, color = algorithm)) + 
  geom_line(alpha=1,position=position_dodge(0.2)) +
  geom_errorbar(alpha=.5, aes(ymin=min_time, ymax=max_time), 
                width=.3, position=position_dodge(0.2)) +
  geom_point(position=position_dodge(0.2)) +
  scale_color_manual(values = colors,limits=alg_names, drop = FALSE) +
  labs(x="signal-to-noise ratio", y = "time(s)") +
  scale_y_log10(breaks = c(0.05, 0.1, 0.2,0.3,0.5,1,2,5,10,20,50,100,300,500,1000,2000)) +
  # scale_y_continuous(breaks = seq(0,1,0.1))+
  # coord_cartesian(ylim = c(0, 1)) +
  theme_pubclean()
# theme_classic2() 


plot4 <- 
  ggplot(data = time_results[sample_size==500 & range==0.1 & signal_to_noise==2,], 
         aes(x = smoothness, y = mean_time, group = algorithm, color = algorithm)) + 
  geom_line(alpha=1,position=position_dodge(0.2)) +
  geom_errorbar(alpha=.3, aes(ymin=min_time, ymax=max_time), 
                width=.3, position=position_dodge(0.2)) +
  geom_point(position=position_dodge(0.2)) +  
  scale_color_manual(values = colors,limits=alg_names, drop = FALSE) +
  labs(x="smoothness", y = "time(s)") +
  scale_y_log10(breaks = c(0.05,0.1,0.2,0.3,0.5,1,seq(2,6,2), 10,20,50,100,200,500,1000,2000)) +
  # scale_y_continuous(breaks = seq(0,2,0.1))+
  # coord_cartesian(ylim = c(0.05, 0.8)) +
  theme_pubclean()

detach(time_results)

ggarrange(plot1,plot2,plot3,plot4, ncol = 2, nrow = 2, common.legend=TRUE)


########  plot of number of iterations  ########

num_it = cbind(algorithm=rep(alg_names, 17), comb[rep(1:17, each=length(alg_names)),])

row.names(num_it) <- 1:(17*length(alg_names))
num_it$mean_num <- c(apply(results["convergence",,,], 1:2, mean))
# num_it$se_num <- c(apply(results["convergence",,,], 1:2, sd)/sqrt(reps))
num_it$min_num <- c(apply(results["convergence",,,], 1:2, min))
num_it$max_num <- c(apply(results["convergence",,,], 1:2, max))

### convert parameters into factor variables
num_it$sample_size <- as.factor(num_it$sample_size)
num_it$range <- as.factor(num_it$range)
num_it$signal_to_noise <- as.factor(num_it$signal_to_noise)
num_it$smoothness <- as.factor(num_it$smoothness)

attach(num_it)
#### make plots
plot_num1 <- 
  ggplot(data = num_it[range==0.1 & signal_to_noise==2 & smoothness==1.5,], 
         aes(x = sample_size, y = mean_num, group = algorithm, color = algorithm)) + 
  geom_line(position=position_dodge(0.2)) +
  geom_errorbar(aes(ymin=min_num, ymax=max_num), 
                width=.3, position=position_dodge(0.2), alpha=.5) +
  geom_point(position=position_dodge(0.2)) +
  scale_color_manual(values = colors,limits=alg_names, drop = FALSE) +
  scale_y_log10(breaks = c(4,5,7,9,12,15,20,30,40,50,70,100,150,200,300,400,600,1000)) +
  # coord_cartesian(ylim = c(4, 1000)) +
  labs(x="sample size", y = "number of iterations") +
  theme_pubclean()


plot_num2 <- 
  ggplot(data = num_it[sample_size==500 & signal_to_noise==2 & smoothness==1.5,], 
         aes(x = range, y = mean_num, group = algorithm, color = algorithm)) + 
  geom_line(alpha=1,position=position_dodge(0.2)) +
  geom_errorbar(alpha=.5, aes(ymin=min_num, ymax=max_num), 
                width=.3, position=position_dodge(0.2)) +
  geom_point(position=position_dodge(0.2)) +
  scale_x_discrete(labels=c('0', '0.01', '0.05', '0.1', '0.5', 'sqrt(d)')) +
  scale_color_manual(values = colors,limits=alg_names, drop = FALSE) +
  labs(x="effective range", y = "number of iterations") +
  scale_y_log10(breaks = c(3,4,5,7,15,seq(10,50,10),70,100,150,200,300,500,750,1000)) +
  # coord_cartesian(ylim = c(4, 1000)) +
  theme_pubclean()

plot_num3 <- 
  ggplot(data = num_it[sample_size==500 & range==0.1 & smoothness==1.5,], 
         aes(x = signal_to_noise, y = mean_num, group = algorithm, color = algorithm)) + 
  geom_line(alpha=1,position=position_dodge(0.2)) +
  geom_errorbar(alpha=.5, aes(ymin=min_num, ymax=max_num), 
                width=.3, position=position_dodge(0.2)) +
  geom_point(position=position_dodge(0.2)) +
  scale_color_manual(values = colors,limits=alg_names, drop = FALSE) +
  labs(x="signal-to-noise ratio", y = "number of iterations") +
  scale_y_log10(breaks = c(4,6,8,15, seq(10,50,10),70,100,150,200,300,500,750,1000)) +
  # coord_cartesian(ylim = c(5, 1000)) +
  theme_pubclean()


plot_num4 <- 
  ggplot(data = num_it[sample_size==500 & range==0.1 & signal_to_noise==2,], 
         aes(x = smoothness, y = mean_num, group = algorithm, color = algorithm)) + 
  geom_line(alpha=1,position=position_dodge(0.2)) +
  geom_errorbar(alpha=.5, aes(ymin=min_num, ymax=max_num), 
                width=.3, position=position_dodge(0.2)) +
  geom_point(position=position_dodge(0.2)) +
  scale_color_manual(values = colors,limits=alg_names, drop = FALSE) +
  labs(x="smoothness", y = "number of iterations") +
  scale_y_log10(breaks = c(3,4,5,7,10,15,20,30,50,70,seq(100,300,100),500,750,1000)) +
  # coord_cartesian(ylim = c(5, 300)) +
  theme_pubclean()

detach(num_it)

ggarrange(plot_num1,plot_num2,plot_num3,plot_num4, ncol = 2, nrow = 2, common.legend=TRUE)



####### GP-matern with vecchia and fitc approximation #######
gp_approx_s = c("vecchia", "fitc")

times = results["time",,,,]
# times = results["time",,,,1:5] # binary model with linear terms

time_results = data.frame(algorithm=rep(alg_names, 4), sample_size = rep(rep(n_l, each = length(alg_names)), length(gp_approx_s)), 
                          gp_approx = rep(gp_approx_s, each = length(gp_approx_s)*length(alg_names)))
time_results$mean_time = c(apply(times, 1:3, mean, na.rm = TRUE))
time_results$min_time = c(apply(times, 1:3, min, na.rm = TRUE))
time_results$max_time = c(apply(times, 1:3, max, na.rm = TRUE))

### convert parameters into factor variables
time_results$algorithm <- as.factor(time_results$algorithm)
time_results$sample_size <- as.factor(time_results$sample_size)
time_results$gp_approx <- as.factor(time_results$gp_approx)
time_results$x_var <- paste0(time_results$sample_size, ",",time_results$gp_approx)
time_results$x_var <- factor(time_results$x_var, 
                             levels = c("5000,vecchia", "5000,fitc", "20000,vecchia", "20000,fitc"))
# attach(time_results)
colors = viridis_pal(option = "H")(6)

# colors = colors[c(1:3,5:6)] # gp with linear regression term
# colors = colors[c(1:2,4:6)] # binary model
# colors = colors[c(1:2,5:6)] # binary model with linear terms

# colors = viridis_pal(option = "D")(5) # bfgs variants plots
# colors[1] = "#440156FF"; colors[2] = "#7A0403FF"; colors[3] = "#61908CFF"

# options(scipen=10000)

plot <- 
  ggplot(data = time_results, 
         aes(x = x_var, y = mean_time, group = algorithm, color = algorithm)) + 
  geom_line(position=position_dodge(0.2)) +
  geom_errorbar(aes(ymin=min_time, ymax=max_time), 
                width=.3, position=position_dodge(0.2), alpha=.5) +
  geom_point(position=position_dodge(0.2)) +
  scale_color_manual(values = colors, limits=alg_names, drop = FALSE) +
  scale_y_log10(breaks = c(0.01,0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100,200,500,1000,2000,4000)) +
  # coord_cartesian(ylim = c(0.01, 50)) +
  labs(x="sample size + gp approximation", y = "time(s)") +
  theme_pubclean() + 
  guides(color = guide_legend(nrow = 3))

detach(time_results)

plot




