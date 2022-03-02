######################################################################################################
############## For the greedy algorithm
############## Theorem 1 & 2
######################################################################################################
library(glmnet)
library(RGCCA)
library(cluster)
library(kernlab)
library(CCA)
library(mclust)
library(class)
library("R.matlab")
library(R1magic)
library(psych)
library(fpc)
library(caret) ##### confusionMatrix
library(RSNNS)
library(MASS)  ###lda
library(e1071) ###nbc
require(nnet)
library(mda) 

##### convariance matrix generation
pd.matrix <- function(rrho, N){
  pd <- matrix(0, N, N)
  for (i in 2:N){
    for (j in 1:(i-1))
      pd[i,j] <- rrho^(abs(i-j))
  }
  pd.final <- pd + t(pd) + diag(N)
  result <- list( pd = pd.final)
  return(result)
}

PD1 <- pd.matrix(0.2, 1000)
PD2 <- pd.matrix(0.8, 1000)
Sigma1 <- PD1$pd
Sigma2 <- PD2$pd
Sigma <- matrix(0,2000,2000)
Sigma[1:1000, 1:1000] <- Sigma1
Sigma[1001:2000, 1001:2000] <- Sigma2


######################################################################################################
############## Experiment 1: Using simulation to verify
######################################################################################################

set.seed(20)
##### to let them have the same X, to have consistent performance of the objective value
Random.sample <- mvrnorm(n = 1000, rep(0,2000), Sigma)
##### The final function for the experiment
Convergence.ver1 <- function(N, Random.sample, sparsity, Number.repetition){
  ##### some preparation for initial parameter
  Numb.Cluster <- length(sparsity)
  lambda <- 0.08  #### this can be tuned 
  s <- max(sparsity)
  
  ##### final output variable
  Running.time <- rep(0, Number.repetition)
  Accuracy <- rep(0, Number.repetition)
  Objective.final <- rep(0, Number.repetition)
  
  ##### generation of X
 
  X <- as.matrix(Random.sample[1:N,])
  DATA <- as.matrix(Random.sample[1:N,])
  Numb.feature <- ncol(X)
  Numb.sample <- nrow(X)
  
  ##### generation of the position and value of nonzero \beta 
  index1 <- 501:(500 + sparsity[1])
  index2 <- 1301:(1300 + sparsity[2])
  set.seed(10)
  aa <- matrix(runif(2 * sparsity[1], min = 2.9, max = 3), sparsity[1], 2)
  set.seed(10)
  bb <- matrix(runif(2 * sparsity[2], min = -1.6, max = -1.5), sparsity[2], 2)
  Beta1 <- matrix(0,2000, 2) 
  Beta2 <- matrix(0,2000, 2) 
  Beta1[index1, ] <- aa
  Beta2[index2, ] <- bb
  
  ##### setup for Y
  Y1 <- matrix(0, N, 2)
  Y2 <- matrix(0, N, 2)
  
  ##### repetition of experiments
  set.seed(43)
  for (jjj in 1:Number.repetition){
  # set.seed(43)

  Y1[,1] <- X %*% Beta1[,1] + rnorm(N, mean = 0, sd =  sqrt(var(X %*% Beta1[,1])/1))   # sqrt(var(X %*% Beta1[,1])/9)
  # set.seed(43)
  Y1[,2] <- X %*% Beta1[,2] + rnorm(N, mean = 0, sd =  sqrt(var(X %*% Beta1[,2])/1))   # sqrt(var(X %*% Beta1[,2])/9)
  # set.seed(43)
  Y2[,1] <- X %*% Beta2[,1] + rnorm(N, mean = 0, sd =  sqrt(var(X %*% Beta2[,1])/1))  # sqrt(var(X %*% Beta2[,1])/9)
  # set.seed(43)
  Y2[,2] <- X %*% Beta2[,2] + rnorm(N, mean = 0, sd =  sqrt(var(X %*% Beta2[,2])/1) )   # sqrt(var(X %*% Beta2[,2])/9)

  Y <- NULL
  Y[[1]] <- Y1
  Y[[2]] <- Y2

  A.inverse <- NULL
  for (j in 1:Numb.Cluster)
    A.inverse[[j]] <- 1/(lambda * N) * diag(Numb.sample)
  G <- NULL
  for (j in 1:Numb.Cluster)
    G[[j]] <- rep(0, sparsity[j])
  
  feature <- 1:Numb.feature
  zzz.val <-  (norm(Y1,"F" )^2 +  norm(Y2,"F" )^2)/N
  
  start_time <- Sys.time()

  for(iter in 1:s){
    for(iter.j in 1:Numb.Cluster){
      if(sparsity[iter.j] < iter) {
        zzz.val[(iter-1)*Numb.Cluster + iter.j + 1 ] <- zzz.val[(iter-1)*Numb.Cluster + iter.j] 
        next}
      
      BB <- as.vector(diag(t(DATA[,feature])%*% A.inverse[[iter.j]] %*% DATA[,feature])) + 1
      obj.val <- - colSums ( (t(Y[[iter.j]]) %*% A.inverse[[iter.j]] %*% DATA[,feature])^2 )/BB
      index.i <- which.min(obj.val)
      index.feature <- feature[index.i]
      G[[iter.j]][iter] <- index.feature  
      feature <- feature[-index.i]
      AA <- as.matrix(A.inverse[[iter.j]] %*% DATA[, index.feature] %*% t(DATA[, index.feature])%*%  A.inverse[[iter.j]])
      
      A.med <-  A.inverse[[iter.j]] - AA/BB[index.i]
      A.inverse[[iter.j]] <- A.med
      kkk <- zzz.val[(iter-1)*Numb.Cluster + iter.j]
      zzz.val[(iter-1)*Numb.Cluster + iter.j + 1 ] <- kkk   +   lambda*min(obj.val)
    }
  }
  
  end_time <- Sys.time()
  Running.time[jjj] <- end_time - start_time
  Objective.final[jjj] <- min(zzz.val)
  Accuracy[jjj] <- (length(intersect(G[[1]], index1)) + length(intersect(G[[2]], index2)) )/(sparsity[1]+sparsity[2])
  
  cat("Number of Repetition", jjj, "\n")
  }
  
  result <- list( fvalue= Objective.final,time = Running.time, Accuracy = Accuracy)
  return(result)
  
}

ds <- Convergence.ver1(1000, Random.sample, c(20,20), 5)

qds$Accuracy
ds$time
ds$fvalue

mean(ds$Accuracy)
mean(ds$time)
mean(ds$fvalue)

######################################################################################################
############## Experiment 2: Using simulation to verify
######################################################################################################
pd.matrix <- function(rrho, N){
  pd <- matrix(0, N, N)
  for (i in 2:N){
    for (j in 1:(i-1))
      pd[i,j] <- rrho^(abs(i-j))
  }
  pd.final <- pd + t(pd) + diag(N)
  result <- list( pd = pd.final)
  return(result)
}

PD1 <- pd.matrix(0.2, 500)
PD2 <- pd.matrix(0.8, 500)
Sigma1 <- PD1$pd
Sigma2 <- PD2$pd
Sigma <- matrix(0,1000,1000)
Sigma[1:500, 1:500] <- Sigma1
Sigma[501:1000, 501:1000] <- Sigma2

set.seed(20)
DaD <- mvrnorm(n = 1000,rep(0,1000) , Sigma)
set.seed(20)
label.index1 <- sample(c(1,2),1000, replace=T)

set.seed(600)
label.index2 <- sample(c(1,2),1000, replace=T)
##### The final function for the experiment
Convergence.ver2 <- function(N, DaD, sparsity, Number.repetition, lamb){
  ##### some preparation for initial parameter
  
  Numb.Cluster <- length(sparsity)
  lambda <- lamb  #### this can be tuned 
  s <- max(sparsity)
  
  ##### final output variable
  Running.time <- rep(0, Number.repetition)
  Accuracy <- rep(0, Number.repetition)
  Objective.final <- rep(0, Number.repetition)
  
  ##### generation of X
  set.seed(20)
  Random.sample <- DaD[1:N,]
  X <- as.matrix(Random.sample)
  DATA <- as.matrix(Random.sample)
  Numb.feature <- ncol(X)
  Numb.sample <- nrow(X)
  
  ##### generation of the position and value of nonzero \beta 
  index1 <- 251:(250 + sparsity[1])
  index2 <- 701:(700 + sparsity[2])

  ##### setup for Y
  Y1 <- matrix(0, N, 2)
  Y2 <- matrix(0, N, 2)
  for (jj in 1:N){
    Y1[jj,label.index1[jj]] <- 1
    Y2[jj,label.index2[jj]] <- 1
  }

  for (iii in 1:3){
    aa <- solve(t(DaD[1:N,index1]) %*% DaD[1:N,index1]) %*% t(DaD[1:N,index1])  %*% Y1
    bb <- solve(t(DaD[1:N,index2]) %*% DaD[1:N,index2]) %*% t(DaD[1:N,index2])  %*% Y2
    Y.prime1 <- matrix(0, N, 2)
    Y.prime2 <- matrix(0, N, 2)
    for (jj in 1:N){
      Y.prime1[jj, apply(DaD[1:N,index1] %*% aa, 1, function(x) which.max(x))[jj] ]  <- 1
      Y.prime2[jj, apply(DaD[1:N,index2] %*% bb, 1, function(x) which.max(x))[jj] ]  <- 1
    }
    Y11 <- Y1
    Y22 <- Y2
    Y1 <- Y.prime1
    Y2 <- Y.prime2
  }

  ##### repetition of experiments
  set.seed(43)
  for (jjj in 1:Number.repetition){
    Y111 <- Y11 + matrix(rnorm(2 * N, mean = 0, sd = 0.01), N, 2)
    Y222 <- Y22 + matrix(rnorm(2 * N, mean = 0, sd = 0.01), N, 2)
    
    Y <- NULL
    Y[[1]] <- Y111
    Y[[2]] <- Y222
    
    A.inverse <- NULL
    for (j in 1:Numb.Cluster)
      A.inverse[[j]] <- 1/(lambda * N) * diag(Numb.sample)
    G <- NULL
    for (j in 1:Numb.Cluster)
      G[[j]] <- rep(0, sparsity[j])
    
    feature <- 1:Numb.feature
    zzz.val <-  (norm(Y111,"F" )^2 +  norm(Y222,"F" )^2)/N
    
    start_time <- Sys.time()
    
    for(iter in 1:s){
      for(iter.j in 1:Numb.Cluster){
        if(sparsity[iter.j] < iter) {
          zzz.val[(iter-1)*Numb.Cluster + iter.j + 1 ] <- zzz.val[(iter-1)*Numb.Cluster + iter.j] 
          next}
        
        BB <- as.vector(diag(t(DATA[,feature])%*% A.inverse[[iter.j]] %*% DATA[,feature])) + 1
        obj.val <- - colSums ( (t(Y[[iter.j]]) %*% A.inverse[[iter.j]] %*% DATA[,feature])^2 )/BB
        index.i <- which.min(obj.val)
        index.feature <- feature[index.i]
        G[[iter.j]][iter] <- index.feature  
        feature <- feature[-index.i]
        AA <- as.matrix(A.inverse[[iter.j]] %*% DATA[, index.feature] %*% t(DATA[, index.feature])%*%  A.inverse[[iter.j]])
        
        A.med <-  A.inverse[[iter.j]] - AA/BB[index.i]
        A.inverse[[iter.j]] <- A.med
        kkk <- zzz.val[(iter-1)*Numb.Cluster + iter.j]
        zzz.val[(iter-1)*Numb.Cluster + iter.j + 1 ] <- kkk   +   lambda*min(obj.val)
      }
    }
    
    end_time <- Sys.time()
    Running.time[jjj] <- end_time - start_time
    Objective.final[jjj] <- min(zzz.val)
    Accuracy[jjj] <- (length(intersect(G[[1]], index1)) + length(intersect(G[[2]], index2)) )/(sparsity[1]+sparsity[2])
    
    cat("Number of Repetition", jjj, "\n")
  }
  
  result <- list( fvalue = Objective.final,time = Running.time, Accuracy = Accuracy)
  return(result)
  
}

ds <- Convergence.ver2(500, DaD, c(5,5),5, 0.02)
ds$Accuracy
ds$time
ds$fvalue

mean(ds$Accuracy)
mean(ds$time)
mean(ds$fvalue)

##### generation of the position and value of nonzero \beta 
# index1 <- 251:(250 + sparsity[1])
# index2 <- 701:(700 + sparsity[2])
# set.seed(10)
# aa <- matrix(runif(2 * sparsity[1], min = 2.9, max = 3), sparsity[1], 2)
# set.seed(10)
# bb <- matrix(runif(2 * sparsity[2], min = -3, max = -0), sparsity[2], 2)
# Beta1 <- matrix(0,1000, 2) 
# Beta2 <- matrix(0,1000, 2) 
# Beta1[index1, 1] <- aa[,1]
# Beta2[index2, ] <- bb
# 
# ##### setup for Y
# Y <- matrix(0, 500, 2)




  

