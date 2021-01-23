############################################################################
#####    codes for Information_incorporated_prevalence_clustering  ##########
##############################################################################


library("trcadmm")
library("foreach")
library("doParallel")
source("gen_cluster_y.R")
source("fun.R") 

error = 0.2
path <- paste0("simdata/500Error=",error,".RData")
load(path)

################################parameter setting######################
eps = 0.001
eta.seq<- seq(0,1,0.1)
#lambda
lambda.seq <- exp(seq(log(10),log(0.001),len= 30))
lambda.p.seq <- exp(seq(log(10),log(0.001),len= 30))
ntimes <- length(times)

ptm <- proc.time()
cl <- makeCluster(4)
registerDoParallel(cl) 
select <- sample(1:500, 4)

result <- foreach(r= select, .combine = rbind, .packages = c("fda","igraph"))%dopar%{
  ##load data
  ydat <- simdata[[r]]
  y = ydat[,-(ntimes+1)]
  ##center data
  y_center = t(apply(y, 1, function(v){ v- mean(v)}))
  beta_ols=matrix(0,nrow= N, ncol=ncol(B) )
  for( pop in 1:N)
  {
    beta_ols[pop,] = lm(y_center[pop,]~ 0 + B)$coefficients
  }
  res = NULL
  compare.res <- rep(NA,9)
  
  ###############################comparsion methods####################
  tryCatch({
    compare <- compare.method(beta_ols, y_center = y_center, times = times, groups = groups,
                                  lambda.seq = lambda.seq, eps = eps)
    compare.res <-compare$res
  })

  res <- c(res, compare.res)
  ######################3proposed method : with the prior##############################
  ##weight setting
  w <- setweight(N, groups, sigma = 0.5)
  #w2 <- setweight(N, groups, sigma = 0.2)
  #w  <- c(w1, w2)
  nweight <- length(w)
  fit <- list()
  for (k in seq_len(nweight)){
    weight <- w[[k]]
    fit.cv = cv.eta.lambda(y = y_center, B=B, lambda.p.seq = lambda.p.seq, lambda.seq = lambda.seq,
                           eps= eps, weight= weight, beta0 = beta_ols, eta.seq = eta.seq, maxit = 1000)
    Prior = evaluta.result(fit.cv$beta.p, groups = groups,  y = y_center, B = B)
    res <- c(res, P = Prior$eval.result, Prss = Prior$rss)
    pF = evaluta.result(fit.cv$beta, groups = groups,  y = y_center, B = B)
    res <- c(res, pF = pF$eval.result, pFrss = pF$rss,  eta= fit.cv$eta.min)
    fit[[k]] <- fit.cv
  }
  save(y_center, B, groups, res, compare, fit, file = paste0("Error=",error,"/",r,".RData"))
  res
}
stopCluster(cl)
save(result, file = paste0("Error=",error,".RData"))
proc.time() - ptm


