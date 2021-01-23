
#@parm y 
#@B basic function
##nPair, the number of pair i<j
ADMM_b <- function(y, B, v,  z, track, beta, lambda, rho, N, K, nPair){
  # K <- ncol(B)
  # nPair <- nrow(z)
  
  beta_up <- t (sapply(1:N, function(i){
    ###update beta
    tmp1 <- crossprod(B) + diag((N-1)*rho, K)
    
    sm <- which(track[,1]==i)
    lr <- which(track[,2]==i)
    ztmp1=vtmp1=rep(0, K)
    if(length(sm)>0)
    {
      vtmp1 <- -apply(v[sm,,drop=FALSE],2,sum)
      ztmp1 <- apply(z[sm,,drop=FALSE],2,sum)
    }
    ztmp2=vtmp2=rep(0,K)
    if(length(lr)>0)
    {
      vtmp2 <- apply(v[lr,,drop=FALSE],2,sum)
      ztmp2<- -apply(z[lr,,drop=FALSE],2,sum)
    }
    tmp2 <- t(B)%*%y[i,] + rho*(vtmp1+vtmp2)+rho*(colSums(beta[-i,]))+rho*(ztmp1+ztmp2)
    
    return(solve(tmp1)%*%tmp2)
  }))
  
  ###update z
  z_up <- z
  v_up <- v
  for (d in seq_len(nPair)){
    i <- track[d,1]
    j <- track[d,2]
    r <- beta_up[i,] - beta_up[j,] + v[d,1:K]
    norm <- sqrt(sum(r^2))
    z_up[d,] <- as.numeric(max(0, norm - lambda/ rho) * r /norm)
    v_up[d,] <- v[d,1:K] + rho*(beta_up[i,] - beta_up[j,] - z_up[d, 1:K])
  }
  res<- list(beta = beta_up , v = v_up, z=z_up)
  return(res)
}

#@parm y the matrix N*ntime
#@parm B the basic fucntion, matrix, ntimes*K
#$parm beta0 the initial parmater. matrix, N*K
#@parm lambda: the tuning parmater. double
estimate.beta <- function(y, B, beta0, lambda, eps, rho = 1, maxit = 1000){
  N <- nrow(y)
  # D <- ncol(y)
  K <- ncol(B)
  track <- expand.grid(1:N, 1:N)
  track <- as.matrix(track[track[,1] < track[,2],])
  nPair <- nrow(track) 
  z <- matrix(0, nrow = nPair, ncol = K)
  v <- z
  beta <- beta0
  iter <- 0
  
  maxdiff <- 1
  ## update by Iteration
  while (maxdiff > eps){
    # admm_up <- ADMM_b(y, B, v, z, track, beta, lambda, rho ,N, K, nPair )
    admm_up <- trcadmm::ADMMcpp(y, B, v, z, track, beta, lambda, rho ,N, K, nPair )
    
    beta_up <- admm_up$beta
    z_up <- admm_up$z
    v_up <- admm_up$v
    distdiff <- rep(NA, N)
    for (i in 1:N){
      distdiff[i] <- sum(((beta[i, ] - beta_up[i, ])/beta[i, ])^2)
    }
    maxdiff <- sqrt(max(distdiff))
    
    # maxdiff <- sqrt(max(apply( (beta - beta_up)^2, 1, sum)))
    iter = iter + 1
    if (iter > maxit){break}
    beta =  beta_up
    z = z_up
    v = v_up
  }
  cat(paste0("the iter is:",iter, "; The max diff:", maxdiff, "\n"))
  res <- list(beta = beta_up, z = z_up, iter = iter, maxdiff = maxdiff)
  return(res)
}

##choose the tuning parameter
cv.estimate <- function(y, B, beta0, lambda.seq, eps, rho = 1, maxit = 1000){
  # lambda.seq <- exp(seq(log(1.5),log(0.001),len= 10))
  D <- ncol(y)
  split <- seq(1, D, 2) 
  train.y <- y[,split]
  test.y <- y[,-split]
  L <- length(lambda.seq)
  obj <- numeric(L)
  for (i in seq_len(L)){
    fit <-  estimate.beta(y = train.y, B = B[split,], beta0 = beta0, eps = eps, lambda= lambda.seq[i], maxit = maxit)
    dis <- as.matrix(dist(fit$beta))
    cutoff = quantile(as.numeric(dist(fit$beta)),0.1)
    # groups_dist(dis =  dis, cutoff = cutoff)
    ####igraph
    g <- graph.adjacency(dis < cutoff)
    # subcomponent(g, 4, mode="out")
    myc = igraph::clusters(g, mode="weak")
    obj[i] <- obj_fun(y = test.y, B = B[-split, ], beta = fit$beta, groups = myc)$L
  }
  lambda.min <- lambda.seq[which.min(obj)][1]
  
  fit.opt <-  estimate.beta(y = y, B = B, beta0 = beta0, eps = eps, lambda= lambda.min, maxit = maxit)
  res <- list(beta = fit.opt$beta, lambda = lambda.seq, lambda.min = lambda.min, obj = obj)
  return (res)
}

### minimizes the function for tuning paramater
obj_fun <- function(y, B, beta, groups, weight = NULL){

  N <- nrow(y)
  if (is.null(weight)){
    weight = matrix(1,N,N)
  }
  loss <- sum( sapply(1:N, function(i) { (y[i,] - B%*%beta[i,])^2 }))
  gc <- groups$no
  gcd <- 0
  for (i in seq_len(gc)){
    if (groups$csize[i] > 1){
      gcd  <- gcd + sum(as.numeric(dist(beta[groups$membership == i, ])))
    }
    # else{
    #   gcd <- gcd + sqrt(sum(beta[i,]^2))
    # }
  }
  
  # for(i in seq_len(N-1)){
  #   for (j in seq(i+1, N, 1)){
  #     # if (groups$membership[i] == groups$membership[j]){
  #       gcd <- gcd + weight[i,j]*sum((beta[i,] - beta[j,])^2)
  #     # }
  #   }
  # }
  
  L <- loss + gcd
  res <- list(L = L, loss = loss,  gcd = gcd)
  return (res)
}

##estimate the beta based on prior information
prior.beta <- function(y, B, beta0, lambda, weight = NULL, eps, alpha = 0.2, maxit = 1000){
  # pset <- sample(1:N, 6)
  # yp <- y[pset,]
  # beta0p <- beta0[pset,]
  # weight <- as.matrix(dist(beta0p)^(-1))
  # alpha <- .1
  # eps <- 0.01
  if (is.null(weight)){
    weight = matrix(1, N, N)
  }
  N <- nrow(y)
  K <- ncol(B)
  maxdiff <- 1
  beta.up <- beta0
  iter <- 0
  while(maxdiff > eps){
    beta.old <- beta.up
    for (i in seq_len(N)){
      w.neqi = lambda* (sum(weight[i,]) - weight[i,i])
      tmp1 <- crossprod(B) + diag(2*w.neqi,K)
      w.neqij <- t(beta.up) %*% weight[i,] - weight[i,i]*beta.up[i,]
      tmp2 <- t(B)%*%y[i,] + w.neqij*lambda
      update <- solve(tmp1) %*% tmp2 - beta.old[i,]
      beta.up[i,] <- beta.old[i,] + alpha*update
    }
    distdiff <- rep(NA, N)
    for (i in 1:N){
      distdiff[i] <- sum(((beta.old[i, ] - beta.up[i, ])/beta.old[i, ])^2)
    }
    maxdiff <- sqrt(max(distdiff))
    # maxdiff <- max(apply( (beta.up - beta.old)^2, 1, sum)) 
    iter <- iter + 1
    if (iter > maxit) {break}
  }
  cat(paste0("the iter is:",iter, "; The max diff:", maxdiff, "\n"))
  res <- list(iter = iter, maxdiff = maxdiff, beta = beta.up)
  return(res)
}

prior.analytical.est <- function(y, B, lambda, weight=NULL){
  
  N <- nrow(y)
  
  if (is.null(weight)){
    weight = matrix(1, N, N)
  }
  K <- ncol(B)
  ## vectorize y
  yk <- as.vector(t(y))
  ## Bk
  IN <- diag(N)
  Bk <- kronecker(IN, B)
  ## D, A
  D <- matrix(0, N*(N-1)/2, N)
  l = 0
  for ( i in seq_len(N-1)){
    for (j in (i+1):N){
      l = l+ 1
      w = weight[i,j]
      D[l, i] = w
      D[l, j] = -w
    }
  }
  A <- kronecker(D, diag(K))
  betak <- solve(crossprod(Bk) + lambda*(crossprod(A))) %*% t(Bk) %*% yk
  beta <- matrix(betak, nrow = N, ncol = K, byrow = T)
  res <- list(beta = beta, lambda = lambda)
  return(res)
}

#cv to choose lambda.prior for prior model
cv.prior <- function(y, B, beta0, lambda.seq, eps, maxit = 1000, 
                     weight = NULL, method = "iter"){
  # lambda.seq <- exp(seq(log(1.5),log(0.001),len= 10))
  N = nrow(y)
  D <- ncol(y)
  split <- seq(1, D, 2) 
  train.y <- y[,split]
  test.y <- y[,-split]
  L <- length(lambda.seq)
  obj <- numeric(L)
  obj2 <- numeric(L)
  for (i in seq_len(L)){
    if (method == "iter"){
      fit.prior <- prior.beta(y = train.y, B = B[split,], beta0 = beta0, lambda = lambda.seq[i], 
                              weight= weight, eps = eps, alpha = 0.2)
    }else{
      fit.prior <- prior.analytical.est(y= train.y, B =B[split, ], lambda = lambda.seq[i], weight = weight)
    }
    dis = as.matrix(dist(fit.prior$beta))
    cutoff = quantile(as.numeric(dist(fit.prior$beta)),0.1)
    # groups_dist(dis =  dis, cutoff = cutoff)
    ####igraph
    g <- graph.adjacency(dis < cutoff)
    # subcomponent(g, 4, mode="out")
    myc = igraph::clusters(g, mode="weak")
    obj[i] <- obj_fun(y = test.y, B = B[-split, ], beta = fit.prior$beta, groups = myc)$L
    obj2[i] <-  obj_fun(y = test.y, B = B[-split, ], beta = fit.prior$beta, groups = myc)$loss
  }
  lambda.min <- lambda.seq[which.min(obj)][1]

  fit.opt <- prior.analytical.est(y= y, B =B, lambda = lambda.min, weight = weight)
  # fit.opt <- prior.beta(y = train.y, B = B[split,], beta0 = beta0, lambda = lambda.seq[i], 
  #                       weight= weight, eps = eps, alpha = 0.3)
  res <- list(beta = fit.opt$beta, lambda = lambda.seq, lambda.min = lambda.min, obj = obj, loss = obj2)
  return (res)
}


##cross validation to choose the tuning parameter: eta and lambda
cv.eta.lambda <- function(y, B, beta0, lambda.p.seq, lambda.seq, eps, weight, eta.seq, maxit = 1000){
  
  D <- ncol(y)
  N <- nrow(y)
  ##prior information
  fit.prior.cv <- cv.prior(y = y, B =B, beta0 = beta0, lambda.seq = lambda.p.seq, 
                           weight= weight, eps = eps, maxit = maxit)
  
  beta.prior=  fit.prior.cv$beta
  yhat = y
  ytilde = y
  for (i in seq_len(N)){
    yhat[i,] = B %*% beta.prior[i,]
    yhat[i,] = yhat[i,] - mean(yhat[i,])
  }
  Le <- length(eta.seq)
  L <- length(lambda.seq)
  obj <- matrix(NA, Le, L)
  ###cv choose the best eta and lambda
  split <- seq(1, D, 2) 
  for (j in seq_len(Le)){
    eta <- eta.seq[j]
    for (i in seq_len(N)){
      ytilde[i,] = (1-eta)* y[i,] + eta*yhat[i,]
      # ytilde[i,] = ytilde[i,] - mean(ytilde[i,])
    }
    for (k in seq_len(L)){
      lambda <- lambda.seq[k]
      ##train.fit
      fit <- estimate.beta(y = ytilde[,-split], B = B[-split,], beta0 = beta0, eps = eps, lambda= lambda*(1-eta), maxit = maxit)
      dis = as.matrix(dist(fit$beta))
      cutoff = quantile(as.numeric(dist(fit$beta)),0.1)
      ####igraph
      g <- graph.adjacency(dis < cutoff)
      myc = igraph::clusters(g, mode="weak")
      ##test
      obj[j,k] <- obj_fun(y = y[,split] , B = B[split,], beta = fit$beta, groups = myc)$L
    }
  }
  ##find the optimal tuning paramater
  index = which.min(obj)[1]
  eta.min = eta.seq[ifelse(index %% Le ==0, Le, index %% Le)]
  lambda.min = lambda.seq[ceiling(index /Le)]
  for (i in seq_len(N)){
    ytilde[i,] = (1-eta.min)* y[i,] + eta.min*yhat[i,]
    # ytilde[i,] = ytilde[i,] - mean(ytilde[i,])
  }
  fit.opt <- estimate.beta(y = ytilde, B = B, beta0 = beta0, eps = eps, lambda= lambda.min*(1-eta), maxit = maxit)
  
  res <- list(beta= fit.opt$beta, lambda.min = lambda.min, eta.min = eta.min,
              lambda.seq= lambda.seq, eta.seq = eta.seq, beta.p = beta.prior, 
              lambda.p.min = fit.prior.cv$lambda.min, weight= weight, obj = obj)
  
  return (res)
}

################################################################
################## Unity Function #############################
###############################################################3
evaluta.result <- function(beta, groups, y, B){
  
  q.seq <- seq(0.05, 0.2, len = 10)
  len.cut <- length(q.seq)
  diff <- numeric(len.cut)
  dis = as.matrix(dist(beta))
  
  for (i in seq_len(len.cut)){
    q = q.seq[i]
    cutoff = quantile(as.numeric(dist(beta)), q)
    ####igraph
    g <- graph.adjacency(dis < cutoff)
    myc = igraph::clusters(g, mode="weak")
    pred.group <- myc$membership
    pred.mat <- as.numeric(dist(pred.group) == 0)
    true.mat <- as.numeric(dist(groups) == 0)
    diff[i] = sum(abs(pred.mat - true.mat))
  }
  ind <- which.min(diff)
  q = q.seq[ind]
  cutoff = quantile(as.numeric(dist(beta)), q)
  ####igraph
  g <- graph.adjacency(dis < cutoff)
  myc = igraph::clusters(g, mode="weak")
  
  ###sum of residual square
  rss <- 0
  for (i in nrow(y)){
    rss <-rss+ sum((y[i,] - B %*% beta[i,])^2)
  }
  res <- list(cutoff = cutoff, eval.result = diff[ind], group = myc$membership, myc = myc, rss = rss)
  return (res)
}

evaluta.result.obj <- function(beta, y, B){
  
  q.seq <- seq(0.05, 0.2, len = 10)
  len.cut <- length(q.seq)
  obj <- numeric(len.cut)
  dis = as.matrix(dist(beta))
  
  for (i in seq_len(len.cut)){
    q = q.seq[i]
    cutoff = quantile(as.numeric(dist(beta)), q)
    ####igraph
    g <- graph.adjacency(dis < cutoff)
    myc = igraph::clusters(g, mode="weak")
    obj[i] <- obj_fun(y = y , B = B, beta = beta, groups = myc)$L
  }
  ind <- which.min(obj)
  q = q.seq[ind]
  cutoff = quantile(as.numeric(dist(beta)), q)
  ####igraph
  g <- graph.adjacency(dis < cutoff)
  myc = igraph::clusters(g, mode="weak")
  ###sum of residual square
  rss <- 0
  for (i in nrow(y)){
    rss <-rss+ sum((y[i,] - B %*% beta[i,])^2)
  }
  res <- list(cutoff = cutoff, myc = myc, rss = rss, group = myc$membership)
  return (res)
}

diffcluster <- function(pred, true){
  pred.mat <- as.numeric(dist(pred) == 0)
  true.mat <- as.numeric(dist(true) == 0)
  diff = sum(abs(pred.mat - true.mat))
  return(diff)
}


dist2mat <- function(vec,N){
  mat = matrix(0, N, N)
  l = 0
  for (i in seq_len(N-1)){
    for (j in (i+1):N){
      l = l + 1
      mat[i,j]= mat[j,i]= vec[l]
    }
  }
  return(mat)
}


#######################################################################3
###############compare Function #####################################
####################################################################
#@parm dis, the distance matrix, N x N
#@cutoff, if dis< cutoff, then i <-> j
##find a connected subgragh in the graph
groups_dist <- function(dis, cutoff){
###deep find search, graph 
  dfs <- function(dis, vertex){
    set <<- c(set, vertex)
    pset <- which(dis[vertex,] < cutoff)
    for (i in pset){
      # print(i)
      if(!(i %in% set)){
        cat(i, file = "cluster.txt", sep ="\n",append = T)
        # print(i)
        dfs(dis, i)
      }
    }
    return(set)
  }
  cluster_list <- list()
  n <-  nrow(dis)
  tag <- numeric(length = n)
  
  while (all(tag ==1) == FALSE){
    unlabel_set <- which(!(tag==1))
    if (length(unlabel_set) == 1){
      vertex = unlabel_set
    }else{
      vertex <- sample(unlabel_set, 1)
    }
    set <- NULL
    set  <- dfs(dis, vertex)
    cluster_list = c(cluster_list, list(set))
    tag[set] <- 1
  }
  return(cluster_list)
}

##
compare.method <- function(beta_ols, y_center, times, groups, lambda.seq, eps){
  k <- length(table(groups))
  comfun <- list()
  ################################comparsion method ########################3
  ##method1 ols 
  fun1grp <- evaluta.result.obj(beta_ols, y_center, B)
  m1diff = diffcluster(fun1grp$group, groups)
  comfun[[1]] <- fun1grp
  # evaluta.result(beta_ols, groups = groups, y, B)$eval.result
  
  ###method 2, kmeans
  fun2grp <- kmeans(beta_ols, centers = k)
  m2diff = diffcluster(fun2grp$cluster, groups)
  comfun[[2]] <- fun2grp
  
  ##method3. funFEM, Bouveryron et al. 2015
  yfd <- smoothingf(geno = y_center, pos = times)
  fun3grp = funFEM::funFEM(yfd, K= 2:10, init = "kmeans")
  m3diff = diffcluster(fun3grp$cls, groups)
  comfun[[3]] <- fun3grp
  
  
  ### method4 fdapac
  Ly=Lt=list()
  for(i in 1: nrow(y_center))
  {
    Ly[[i]]=y_center[i,]
    Lt[[i]]= times
  }
  fun4grp= fdapace::FClust(Ly,Lt,k=k)
  m4diff <- diffcluster(fun4grp$cluster, groups)
  comfun[[4]] <- fun4grp
  
  ###method 5 funHDDC.  (Bouveyron and Jacques, 2011, <doi:10.1007/s11634-011-0095-6>)
  ## or multivariate data (Schmutz et al., 2018)
  fun5grp=funHDDC::funHDDC(yfd,K=2:10)
  m5diff <- diffcluster(fun5grp$class, groups)
  comfun[[5]] <- fun5grp$class
  
  ##method 6  waveclust: 
  #M. Giacofci, S. Lambert-Lacroix, G. Marot, and F. Picard. Wavelet-based clustering for
  # mixed-effects functional models in high dimension. Biometrics, in press, 2012.
  data = t(y_center)
  fun6grp = funcy::funcit(data  = data, k=k, methods = "waveclust")
  m6diff <- diffcluster(funcy::Cluster(fun6grp), groups)
  comfun[[6]] <- fun6grp
  
  ##method 7, funcy, Gareth James and Catherine A. Sugar. Clustering for Sparsely Sampled Functional Data. 
  # Journal of the American Statistical Association. 98 (462). 297–408. 2003
  fun7grp = funcy::funcit(data  = data, k=k, methods = "fitfclust")
  m7diff <- diffcluster(funcy::Cluster(fun7grp), groups)
  comfun[[7]] <- fun7grp
  
  ##method 8 
  fit.est.cv <- cv.estimate(y_center,B, beta0 = beta_ols, lambda.seq = lambda.seq, eps = eps,maxit = 1000)
  m8diff <- evaluta.result(fit.est.cv$beta, groups = groups, y_center, B)$eval.result
  fun81grp <- evaluta.result.obj(fit.est.cv$beta, y_center, B)
  m81diff <- diffcluster(fun81grp$group, groups)
  comfun[[8]] <- fit.est.cv
  
  
  res <- c(m1diff, m2diff, m3diff, m4diff, m5diff, m6diff, m7diff, m8diff, m81diff, fun1grp$rss, fun81grp$rss)
  names(res) <- c("ols","kmeans", "funFEM", "FClust", "funHDDC", "wavelust", "fitclust", "Flasso","Flasso2","olsrss", "Flassorss")
  return(list(res = res, com = comfun))
}


smoothingf<-function(geno,pos)
{
  nN=nrow(geno)
  pos <- (pos-pos[1])/(pos[length(pos)]-pos[1])
  #pos2 = seq(min(pos), max(pos), length=2*length(pos))
  intrvllen<-0.0001
  pos3<-seq(0,1,intrvllen)
  shouldContinue=TRUE;
  
  if(shouldContinue==TRUE)
  {
    lambda<-rep(0,nN)
    for(j in 1:nN)
    {
      genospline<-try(smooth.spline(pos,geno[j,],nknots=length(pos)))
      if("try-error" %in% class(genospline)){
        lambda[j]<-NA}
      else {lambda[j]<-genospline$lambda}
    }
  }
  
  lambda_all<-mean(lambda,na.rm=TRUE)
  Knots = pos
  norder = 4
  nbasis=length(Knots) + norder - 2
  ybasis=create.bspline.basis(range(Knots), nbasis, norder)
  Lfdobj= 2
  yfdPar = fdPar(ybasis, Lfdobj, lambda=lambda_all)
  yfd = smooth.basis(pos,t(geno),yfdPar)$fd
  yfd<-center.fd(yfd) 
  # yij.f = t(eval.fd(pos,yfd))
  return(yfd)
  
}


