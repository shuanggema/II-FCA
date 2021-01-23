
gen.y <- function(group, g.sizes, error, time){
  number.g = length(group)
  y <- NULL
  for (i in seq_len(number.g)){
    yg <- gen.ycluster(group = group[i], g.size = g.sizes[i], error, time)
    
    y <- rbind(y, cbind(yg, group = group[i]))
  }
  return(y)
}


gen.ycluster.old <- function(group, g.size, error, time){
  ntime <- length(time)
  bias <- 0
  if(group == 1){
    coef <- c(-1,-1.5)
    Bsp <- create.monomial.basis(rangeval=c(0,2),nbasis= 2)
    B <- eval.basis(time, Bsp, returnMatrix=T)
    bias <- 2
  }
  if(group == 2){
    coef <- c(1,1.5)
    Bsp <- create.monomial.basis(rangeval=c(0,2), nbasis= 2)
    B <- eval.basis(time, Bsp, returnMatrix=T)
  }
  if(group == 3){
    coef <- c(1, 1)
    B1 <- (time <= 1)*(1 - abs(time - 1))
    B2 <- as.numeric(time > 1) 
    B <- cbind(B1, B2)
  }
  if(group == 4){
    coef <- 2
    B <- (time <= 1)*(1 - abs(time - 1))
  }
  if(group == 5){
    bias = 3
    coef <- c(-1.974947, -.8, -1, -1.1, -.5)
    Bsp <- create.fourier.basis(rangeval=c(0,2),nbasis= 5)
    B <- eval.basis(time, Bsp, returnMatrix=T)
  }
  if(group == 6){
    coef <- c(1.974947, .8, 1, 1.1, .5)
    Bsp <- create.fourier.basis(rangeval=c(0,2),nbasis= 5)
    B <- eval.basis(time, Bsp, returnMatrix=T)
  }
  if(group == 7){
    coef <- 2
    B = 1 - abs(time - 1)
  }
  
  if(group == 8){
    coef <- 2
    B = abs(time - 1)
  }
  
  if(group == 9)
  {
    coef=-1
    Bsp=create.constant.basis(rangeval=c(0,2))
    B =eval.basis(time, Bsp, returnMatrix=T)
  }
  if(group ==10)
  {
    coef=1
    Bsp=create.constant.basis(rangeval=c(0,2))
    B =eval.basis(time, Bsp, returnMatrix=T)
  }
  yg <- NULL
  B <- as.matrix(B)
  for (i in seq_len(g.size)){
    ntime <- nrow(B)
    ys <- B %*% coef + rnorm(ntime, mean = 0, sd = error) + bias
    yg <- cbind(yg, ys)
  }  
  return(t(yg))
}



gen.ycluster <- function(group, g.size, error, times){
  ntimes <- length(times)
  bias <- 0
  if(group == 1){
    Bsp <- create.monomial.basis(rangeval=c(0,1), nbasis= 2)
    B <- eval.basis(times, Bsp, returnMatrix=T)
    coef <- c(3,-1)
  }
  if(group == 2){
    Bsp <- create.monomial.basis(rangeval=c(0,1), nbasis= 2)
    B <- eval.basis(times, Bsp, returnMatrix=T)
    coef <- c(1, 1)
  }
  if(group == 3){
    Bsp <- create.monomial.basis(rangeval=c(0,1), nbasis= 4)
    B <- eval.basis(times, Bsp, returnMatrix=T)
    coef <- c(1, 0.5, 2, 1)
  }
  if(group == 6){
    # B = (times < 0.8)*(times) + (times <= 1.2) *(times >= 0.8) +  (times > 1.2) *(2 - times)
    # coef <- 1
    Bsp <- create.fourier.basis(rangeval=c(0,2),nbasis= 3)
    B <- eval.basis(times, Bsp, returnMatrix=T)
    coef <- c(0.5, 0.1, -1)
  }
  if(group == 7){
    # B = (times < 0.8)*(times) + (times <= 1.2) *(times >= 0.8) +  (times > 1.2) *(2 - times)
    # coef <- -1
    # bias <- 1
    Bsp <- create.fourier.basis(rangeval=c(0,2),nbasis= 3)
    B <- eval.basis(times, Bsp, returnMatrix=T)
    coef <- c(0.5, 0.1, 1)
  }
  if(group == 4){
    B1 <- (times <= 1.2)*(times)
    B2 <- 1.2*as.numeric(times > 1.2) 
    B <- cbind(B1, B2)
    coef <- c(1, 1)
  }
  if(group == 5){
    B = (times <= 1.2)*(1.2 - times)
    coef <- 1
  }
  
  if(group == 8){
    Bsp <- create.fourier.basis(rangeval=c(0,2),nbasis= 5)
    B <- eval.basis(times, Bsp, returnMatrix=T)
    coef <- c(1, 0.5,  1, 0.8, 1)
  }
  
  if(group == 9)
  {
    Bsp <- create.fourier.basis(rangeval=c(0,2),nbasis= 5)
    B <- eval.basis(times, Bsp, returnMatrix=T)
    coef <- c(1, 0.5, - 1, 0.8, -1)
  }
  
  yg <- NULL
  B <- as.matrix(B)
  for (i in seq_len(g.size)){
    ntime <- nrow(B)
    ys <- B %*% coef + rnorm(ntime, mean = 0, sd = error) + bias
    yg <- cbind(yg, ys)
  }  
  return(t(yg))
}

################################unity fucntion###############3
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
# a = as.numeric(dist(mat))
# dist2mat(a, nrow(mat))


setweight <- function(N, groups, sigma){
  w <- list()
  ###type 1
  weight.vec = as.numeric(dist(groups)==0) * pmax(rnorm(N*(N-1)/2,1, sigma),0)
  w[[1]] <- dist2mat(vec = weight.vec, N)
  # w[[1]][1:12,1:12]
  ##weight type 2
  weight.vec2 = weight.vec*rbinom(N*(N-1)/2,1,0.8)
  w[[2]] <- dist2mat(vec = weight.vec2, N)
  ##weight type 3
  weight.vec3 = weight.vec*rbinom(N*(N-1)/2,1,0.5)
  w[[3]] <- dist2mat(vec = weight.vec3, N)
  ## type 4
  weight.vec4 <- weight.vec
  noindex <- weight.vec == 0
  weight.vec4[noindex] <- pmax(rnorm(sum(noindex),1, sigma),0)*rbinom(sum(noindex),1,0.1)
  w[[4]] <- dist2mat(vec = weight.vec4, N)
  ## type 5
  weight.vec5 <- weight.vec
  noindex <- weight.vec == 0
  weight.vec5[noindex] <- pmax(rnorm(sum(noindex),1, sigma),0)*rbinom(sum(noindex),1,0.3)
  w[[5]] <- dist2mat(vec = weight.vec5, N)
  
  
  ## type 6
  weight.vec6 <- weight.vec2
  noindex <- weight.vec == 0
  weight.vec6[noindex] <- pmax(rnorm(sum(noindex),1, sigma),0)*rbinom(sum(noindex),1,0.1)
  w[[6]] <- dist2mat(vec = weight.vec6, N)
  ## type 7
  weight.vec7 <- weight.vec2
  noindex <- weight.vec == 0
  weight.vec7[noindex] <- pmax(rnorm(sum(noindex),1, sigma),0)*rbinom(sum(noindex),1,0.3)
  w[[7]] <- dist2mat(vec = weight.vec7, N)
  
  ## type 8
  weight.vec8 <- weight.vec3
  noindex <- weight.vec == 0
  weight.vec8[noindex] <- pmax(rnorm(sum(noindex),1, sigma),0)*rbinom(sum(noindex),1,0.1)
  w[[8]] <- dist2mat(vec = weight.vec8, N)
  ## type 9
  weight.vec9 <- weight.vec3
  noindex <- weight.vec == 0
  weight.vec9[noindex] <- pmax(rnorm(sum(noindex),1, sigma),0)*rbinom(sum(noindex),1,0.3)
  w[[9]] <- dist2mat(vec = weight.vec9, N)
  
  ## type 10
  weight.vec10 <- pmax(rnorm(N*(N-1)/2, 1, sigma),0) *rbinom(N*(N-1)/2,1,0.1)
  w[[10]] <- dist2mat(vec = weight.vec10, N)
  ## type 11
  weight.vec11 <- pmax(rnorm(N*(N-1)/2, 1, sigma),0) *rbinom(N*(N-1)/2,1,0.3)
  w[[11]] <- dist2mat(vec = weight.vec11, N)
  ## type 12
  weight.vec12 <- pmax(rnorm(N*(N-1)/2, 1, sigma),0) *rbinom(N*(N-1)/2,1,0.5)
  w[[12]] <- dist2mat(vec = weight.vec12, N)
  return(w)
}


