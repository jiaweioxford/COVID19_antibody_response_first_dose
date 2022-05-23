
####################################################################################################
# Example of code written by Jia Wei (Nuffield Department of Medicine, University of Oxford)       #
# and Koen Pouwels (Health Economics Research Centre, University of Oxford)                        #
# for simulation analyses of antibody waning following the first vaccination.                      #
# Accompanying paper: SARS-CoV-2 antibody trajectories after a single COVID-19 vaccination with    # 
# and without prior infection.                                                                     #
####################################################################################################

library(MASS)
library(parallel)
library(brms)

visit_schedule <- data.frame(id=1:1000)
visit_schedule$time <- sample(0:28,1000,TRUE)
# create the extra rows
visit_schedule$time1 <- visit_schedule$time + 14
visit_schedule$time2 <- visit_schedule$time + 28
visit_schedule$time3 <- visit_schedule$time + 42
# reshape from wide to long: 
visit_schedule <- gather(visit_schedule, condition, time, time:time3)



get_model <- function(x){
  n <- 1000
  ### average intercept and slope
  beta0 <- 10
  beta1 <- -0.01
  ### true error
  err.val <- 0.1
  ### true error SD, intercept SD, slope SD, and intercept-slope cor
  sigma <- 0.05
  tau0  <- 0.5
  tau1  <- 0.01
  tau01 <- 0.5
  ### maximum number of possible observations
  m <- 70
  ### get number of possible observations for each individual (if measured daily)
  p <- rep(m,n)
  ### set up data frame
  dat <- data.frame(id=rep(1:n, times=p), time=rep(1:m, times=n))
  ### simulate (correlated) random effects for intercepts and slopes
  mu  <- c(0,0)
  S   <- matrix(c(1, tau01, tau01, 1), nrow=2)
  tau <- c(tau0, tau1)
  S   <- diag(tau) %*% S %*% diag(tau)
  U   <- mvrnorm(n, mu=mu, Sigma=S)
  ### simulate errors and then the actual outcomes
  dat$eij <- round(unlist(rnorm(n*m, mean=0, sd=err.val)),3)

  dat$yij <- 2^((beta0 + rep(U[,1], times=p)) + (beta1 + rep(U[,2], times=p))*dat$time + dat$eij)
  
  # only keep observations in line with visit schedule: 
  dat_new <- merge(dat, visit_schedule, by=c("id","time"))
  
  dat_new$censy <- dat_new$y
  dat_new$censy[dat_new$censy >=450] <- 450
  dat_new=dat_new %>% mutate(cen=ifelse(censy==450,1,0))
  
  dat_new$id=as.factor(dat_new$id)
  
  a=prop.table(table(dat_new$cen))[[2]]
  
  priors=c(
    set_prior("normal(10,2)",class="Intercept"),
    set_prior("normal(0,0.1)",coef="time",class="b"),
    
    set_prior("cauchy(0,0.1)",coef="time",class="sd",group="id"),
    set_prior("cauchy(0,1)",coef="Intercept",class="sd",group="id"),
    set_prior("cauchy(0,0.5)",class="sigma"),
    set_prior("lkj_corr_cholesky(1)",class="L")
  )
  
  
  fit <- brm(formula=log2(censy)|cens(cen)~1+time+(1+time|id),
             data=dat_new,cores=4,family=gaussian(),
             prior = priors,
             chains = 4, iter=4000, warmup = 2000, seed=42,init_r=10,
             control = list(adapt_delta=0.95))
  
  sum=summary(fit)
  b=sum$fixed[,1]
  c=sum$random$id[,1]
  
  table=data.frame(prop=numeric(1),intercept=numeric(1),inter_low=numeric(1),inter_high=numeric(1),
                   slope=numeric(1),slope_low=numeric(1),slope_high=numeric(1),
                   sd_inter=numeric(1),sd_slope=numeric(1),correlation=numeric(1),cor_low=numeric(1),cor_high=numeric(1))
  table[,1]=a
  table[,2]=b[1]
  table[,3]=sum$fixed[1,3]
  table[,4]=sum$fixed[1,4]
  table[,5]=b[2]
  table[,6]=sum$fixed[2,3]
  table[,7]=sum$fixed[2,4]
  table[,8]=c[1]
  table[,9]=c[2]
  table[,10]=c[3]
  table[,11]=sum$random$id[3,3]
  table[,12]=sum$random$id[3,4]
  
  return(table)
}


x.list <- sapply(1:100, list)


system.time({
  clust <- makeCluster(24)
  clusterExport(clust, c("mvrnorm","visit_schedule","%>%","mutate","set_prior","brm"))
  a <- parLapply(clust, x.list, get_model)})

stopCluster(clust)
