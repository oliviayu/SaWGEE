library(geeM)
library(dplyr)
library(SimCorMultRes)
library(purrr)
library(glmmML)
library(CorBin)

source("sawgee_code.r")

parallel = F; mc.cores=1
gee_check = c(T, F)[2] # run regular GEE?

####################################################
T_n=50  # number of repeated sub-sampling
n_id=100 # number of clusters in a simulated dataset

# Specify the scenario to run, where
# scenarios[1] controls cluster size, c(500, 5000)
# scenarios[2] controls intercept, c(-2, -4, -6)
# scenarios[3] controls rho, c(2, 5, 0.01)
# scenarios[4] controls proportion of sampled data, by_percent = c(F, T, T, T), frac=c(NA, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9)
scenarios = c(1, 2, 3, 1)
## Parameters for generating cluster sizes from negative binomial 
## Fix cluster size at cluster_mean for small/moderate-sized case
cluster_mean=c(500, 5000)[scenarios[1]]
dispersion=c(NA, 50)[scenarios[1]]
# print(c("standard deviation", (sqrt(cluster_mean+cluster_mean**2/dispersion))))
## variance: cluster_mean+cluster_mean**2/dispersion

## Model para
fm = yvar~x1+x2
intercept_par = c(-2, -4, -6)[scenarios[2]]
true_par = c(intercept_par, 1, 1)
rho = c(0.2, 0.5, 0.01)[scenarios[3]] #controling the strength of within-cluster correlation
## Remove clusters with number of cases less than cut
remove_lowcase = T; cut = 3 
## Specifying number of obs to be sampled
## By percentage if event rate is moderate to high, e.g. when beta0 = -2 and E(y)>0.1
## Otherwise, 2*ncases by default
by_percent = c(F, rep(T, 6))[scenarios[4]] ; 
frac = c(NA, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9)[scenarios[4]]

####################################################
output_os2 = output_os = output_rs = output_cc = output_ccc = c()
output_glm = output_glm2 = output_gee = c()
output_os2_hybrid = output_os_hybrid = c()


B = 300
for(b in 1:B){
  set.seed(b)
  
  if(cluster_mean <= 500){
    n_time = rep(cluster_mean, n_id)
    
    # Generate binary data with exchangeable within-cluster correlation = rho
    DF=lapply(1:n_id, function(k){
      x1=rbinom(n_time[k], 1, prob=0.5)
      x2=rnorm(n_time[k], 0, 1)
      true_mu = true_par[1] + x1*true_par[2] + x2*true_par[3]
      true_prob = exp(true_mu)/(1+exp(true_mu))
      
      success = F
      ss=0
      while(success == F & ss < 10){
        if(ss>0) cat("attempt:", ss, '\n')
        simy = try(cBern(1, true_prob, rho=rho, type="exchange"), silent = T)
        if(all(class(simy) != 'try-error') & (!anyNA(simy))){
          res = data.frame(yvar=as.integer(simy), x1=x1, x2=x2, time=1:n_time[k],
                           idvar=k)
          if(remove_lowcase==T){
            success = (sum(res$y) >= cut)
          } else {
            success =T
          }
        } else {
          success = F
        }
        ss=ss+1
      }
      if(success) res else NULL
    }) %>% do.call(rbind, .)
    
  } else {
    n_time=rnbinom(n_id, size=dispersion, mu=cluster_mean)
    intercept = rnorm(n_id, true_par[1], rho)
    beta1 = rep(true_par[2], n_id)
    beta2= rep(true_par[3], n_id)
    idvar = rep(1:n_id, n_time)
    n_sample=length(idvar)
    
    x1=rbinom(n_sample, 1, prob=rep(runif(n_id, 0.3, 0.7), n_time))
    x2=rnorm(n_sample, rep(rnorm(n_id, 0, 0.5), n_time), 1)
    true_mu = rep(intercept, n_time) + x1*rep(beta1, n_time) + x2*rep(beta2, n_time)
    true_prob = exp(true_mu)/(1+exp(true_mu))
    yvar=rbinom(n_sample, 1, prob=true_prob)
    DF = data.frame(idvar, x1,x2, yvar)
    rm(x1); rm(x2); rm(true_prob); rm(idvar); rm(true_mu); rm(yvar)
    
    if(remove_lowcase == T){
      rmid= group_by(DF, idvar) %>% summarise(case=sum(yvar)) %>%
        filter(case < cut)
      DF = filter(DF, !idvar %in% rmid$idvar)
      if(length(unique(DF$idvar)) < n_id){
        cat(paste("Only", length(unique(DF$idvar)), "clusters remained. Consider increase n_id."), '\n')
      }
    }
  }
  count_cases = group_by(DF, idvar) %>% 
    summarise(ncase=sum(yvar), ncontrol=sum(1-yvar), nk=n())
  
  ## GLM model
  md0 = glm(fm, data=DF, family=binomial())
  output_glm = rbind(output_glm, c(b, coef(md0), summary(md0)$coef[,2]))
  mu = predict.glm(md0, DF, type="response")
  ini.beta = coef(md0)
  rho_ini  = 0
  rm(md0)
  
  ## GEE
  if(gee_check){
    gee_md = gee::gee(fm, id=idvar, data=DF, family='binomial', corstr='exchangeable')
    output_gee = rbind(output_gee, 
                       c(b, c(summary(gee_md)$coef[,1], 
                              summary(gee_md)$coef[,4],
                              gee_md$scale, gee_md$working.correlation[2,1])))
    rho_ini = gee_md$working.correlation[2,1]
    rm(gee_md)
  }
  
  ## Optimal sub-sampling
  if(by_percent==T){
    M = round(frac*nrow(DF)/2)*2
    Mccc = round(count_cases$ncase/sum(count_cases$ncase)*M/2)*2
  } else {
    M <- 2*sum(DF$yvar==1)
    Mccc = count_cases$ncase*2
  }
  
  
  # Hybrid + MSE_minimization
  OS.result2_hybrid = subWGEE(mu, DF, M, T_n, ini.beta, fm, strategy = 'os',
                              hybrid = T,
                              parallel = parallel, mc.cores = mc.cores,
                              method=2, rho_ini = rho_ini)
  output_os2_hybrid = rbind(output_os2_hybrid, cbind(b, OS.result2_hybrid$output))
  
  # Hybrid + Lopt
  OS.result_hybrid = subWGEE(mu, DF, M, T_n, ini.beta, fm, strategy = 'os',
                             hybrid = T,
                             parallel = parallel, mc.cores = mc.cores, track=2,
                             method = 1,  rho_ini=rho_ini)
  output_os_hybrid = rbind(output_os_hybrid, cbind(b, OS.result_hybrid$output))
  
  # MSE_mini
  OS.result2 = subWGEE(mu, DF, M, T_n, ini.beta, fm, strategy = 'os',
                       parallel = parallel, mc.cores = mc.cores,
                       method=2, rho_ini = rho_ini)
  output_os2 = rbind(output_os2, cbind(b, OS.result2$output))
  
  # Lopt      
  OS.result = subWGEE(mu, DF, M, T_n, ini.beta, fm, strategy = 'os',
                      parallel = parallel, mc.cores = mc.cores, track=0,
                      method = 1,  rho_ini=rho_ini)
  output_os = rbind(output_os, cbind(b, OS.result$output))
  
  # Random sampling
  RS.result = subWGEE(mu, DF, M, T_n, ini.beta, fm, strategy = 'random',
                      parallel = parallel, mc.cores = mc.cores)
  output_rs = rbind(output_rs, cbind(b, RS.result$output))
  
  # 1:1 case-control
  CC.result = subWGEE(mu, DF, M, T_n, ini.beta, fm, strategy = 'cc',
                      parallel = parallel, mc.cores = mc.cores, track=2)
  output_cc = rbind(output_cc, cbind(b, CC.result$output))
  
  # within-cluster 1:1 case-control
  CCC.result = subWGEE(mu, DF, Mccc, T_n, ini.beta, fm, strategy = 'ccc',
                       parallel = parallel, mc.cores = mc.cores)
  if(all(class(CCC.result)!='try-error') & !is.null(CCC.result)){
    output_ccc = rbind(output_ccc, cbind(b, CCC.result$output))
  }
  
  cat("###########  Simulation b=", b, 'is done. ########### \n')
}

