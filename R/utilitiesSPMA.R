# R package functions

library(MASS)
library(grpreg)
library(gglasso)
library(survival)
library(glmnet)
library(ncvreg)

# Interval Censor package for Group lasso, group SCAD, and group MCP
library(grplasso)
library(simsurv)
library(pracma)
# --------------- Section One SPMA utilities ------------------

# ADD print result
# BETA_0 is from norm distribution now

#Rcpp::sourceCpp("OneDrive - The University of Texas at Dallas/metadata analysis/code/spmaC.cpp")
#Rcpp::sourceCpp("metaADMM/spmaC_ADMM.cpp")

# generate one set of simulation data
generate_data <- function(N = c(500,500,500),p = 4, correlate = FALSE, rightCensor = FALSE,intervalCensor = FALSE, M =3, q = 0.50,seed){
  # q controls the percentage of significant coefficient
  p1 = floor(p*q)
  p2 = p - p1
  set.seed(seed)
  #beta_01 <- c(runif(p1,-1,1), rep(0,p2))#c(0.8,0.1,0,0)
  #beta_02 <- c(runif(p1,-1,1), rep(0,p2))
  #beta_03 <- c(runif(p1,-1,1), rep(0,p2))
  beta_01 <- c((-1)^c(1:p1)*rnorm(p1,0.5,0.25), rep(0,p2))#c(0.8,0.1,0,0)
  beta_02 <- c((-1)^c(1:p1)*rnorm(p1,0.5,0.25), rep(0,p2))
  beta_03 <- c((-1)^c(1:p1)*rnorm(p1,0.5,0.25), rep(0,p2))

  #zero_idx <- sample(1:p,p,replace = FALSE)
  #beta_01 <- beta_01[zero_idx]
  #beta_02 <- beta_02[zero_idx]
  #beta_03 <- beta_03[zero_idx]
  augN = N*2
  if(correlate){
    # consider scenarios with autoregressive structure
    sigma  = matrix(ncol = p, nrow = p)
    for(i in 1:p){
      for(j in 1:p){
        sigma[i,j] = 0.5^(abs(i-j)) # debug for the correlation
        # rho = 0.2 works well, try rho = 0.3
      }
    }

    X1 <- mvrnorm(n=augN[1],mu =rep(0,p),Sigma = sigma )
    X2 <- mvrnorm(n=augN[2],mu =rep(0,p),Sigma = sigma )
    X3 <- mvrnorm(n=augN[3],mu =rep(0,p),Sigma = sigma )

  }
  else{ # block-correlation structure
    # try r_11 = 0.25, r_22 = 0.125
    sigma = matrix(0,ncol =p, nrow = p )
    sigma[1:p1,1:p1] <- 0.25
    sigma[(p1+1):p,(p1+1):p] = 0.125
    diag(sigma) <- 1

    # sigma = diag(p), identity matrix
    X1 <- mvrnorm(n=augN[1],mu =rep(0,p),Sigma = sigma )
    X2 <- mvrnorm(n=augN[2],mu =rep(0,p),Sigma = sigma )
    X3<-  mvrnorm(n=augN[3],mu =rep(0,p),Sigma = sigma )
  }

  # RMSE = (beta*-beta0)^T E(x^Tx) (beta*-beta0)
  exTx1 = t(X1)%*%X1/augN[1]
  exTx2 = t(X2)%*%X2/augN[2]
  exTx3 = t(X3)%*%X3/augN[3]

  if(rightCensor){
    # right censored data
    eventU <- runif(augN[1],0,1)
    event1 <- -log(eventU)/(exp(X1%*%beta_01))
    #event1 <- rexp(augN[1],exp(X1%*%beta_01)) # event time
    u = runif(1,2,3)
    censor1 <- rexp(augN[1], rate = exp(u*X1%*%beta_01))
    # if(length(censor1)!=length(event1)){
    #   print("Length of censor time and event time should be the same!")
    #   break;
    # }
    status1 <- ifelse(censor1>event1,1,0)
    event1 = ifelse(censor1>event1,event1,censor1)



    #event2 <- rexp(augN[2],exp(X2%*%beta_02)) # event time
    eventU <- runif(augN[2],0,1)
    event2 <- -log(eventU)/(exp(X2%*%beta_02))
    u = runif(1,2,3)
    print(paste("u is ", u))
    censor2 <- rexp(augN[2], exp(u*X2%*%beta_02))
    status2 <-  ifelse(censor2>event2,1,0)
    event2 = ifelse(censor2>event2,event2,censor2)



    #event3 <- rexp(augN[3],exp(X3%*%beta_03)) # event time
    eventU <- runif(augN[3],0,1)
    event3 <- -log(eventU)/(exp(X3%*%beta_03))
    u = runif(1,2,3)
    censor3 <- rexp(augN[3], exp(u*X3%*%beta_03))
    status3 <- ifelse(censor3>event3,1,0)
    event3 = ifelse(censor3>event3,event3,censor3)

    # return data
    sample1 <- cbind(X1, event1, status1)
    colnames(sample1)[c(p+1,p+2)] <- c("event", "status")
    sample2 <- cbind(X2, event2,status2)
    colnames(sample2)[c(p+1,p+2)] <- c("event", "status")
    sample3 <- cbind(X3, event3, status3)
    colnames(sample3)[c(p+1,p+2)] <- c("event", "status")
    beta_0 <- rbind(beta_01, beta_02, beta_03)

  } else if(intervalCensor){

    eventU <- runif(augN[1],0,1)
    event1 <- -log(eventU)/(exp(X1%*%beta_01))
    eventU <- runif(augN[2],0,1)
    event2 <- -log(eventU)/(exp(X2%*%beta_02))
    eventU <- runif(augN[3],0,1)
    event3 <- -log(eventU)/(exp(X3%*%beta_03))

    # Five inspection time
    inspection_time1 = matrix(ncol = 5, nrow = augN[1])
    interval_time1 = matrix(ncol = 3, nrow = augN[1])
    colnames(interval_time1) <- c("Left", "Right", "Status")
    for(i in 1:augN[1]){
      t1 <- runif(1,0,1)
      t2 <- t1 + 0.1+runif(1,0,1)
      t3 <- t2 + 0.1+runif(1,0,1)
      t4 <- t3 + 0.1+runif(1,0,1)
      t5 <- t4 + 0.1+runif(1,0,1)
      inspection_time1[i,] <- c(t1,t2,t3,t4,t5)
      tmp <- c(t1,t2,t3,t4,t5)

      # find interval for each subject
      if(event1[i] <= t1) {interval_time1[i,1:3] = c(0,t1,2)}
      else if(event1[i] >t1 & event1[i]<=t2){interval_time1[i,1:3] = c(t1,t2,3)}
      else if(event1[i] >t2 & event1[i]<=t3){interval_time1[i,1:3] = c(t2,t3,3)}
      else if(event1[i] >t3 & event1[i]<=t4){interval_time1[i,1:3] = c(t3,t4,3)}
      else if(event1[i] >t4 & event1[i]<=t5){interval_time1[i,1:3] = c(t4,t5,3)}
      else {interval_time1[i,1:3] = c(t5,Inf,0)}
    }

    inspection_time2 = matrix(ncol = 5, nrow = augN[2])
    interval_time2 = matrix(ncol = 3, nrow = augN[2])
    colnames(interval_time2) <- c("Left", "Right", "Status")
    for(i in 1:augN[2]){
      t1 <- runif(1,0,1)
      t2 <- t1 + 0.1+runif(1,0,1)
      t3 <- t2 + 0.1+runif(1,0,1)
      t4 <- t3 + 0.1+runif(1,0,1)
      t5 <- t4 + 0.1+runif(1,0,1)
      inspection_time2[i,] <- c(t1,t2,t3,t4,t5)
      tmp <- c(t1,t2,t3,t4,t5)

      # find interval for each subject
      if(event2[i] <= t1) {interval_time2[i,1:3] = c(0,t1,2)}
      else if(event2[i] >t1 & event2[i]<=t2){interval_time2[i,1:3] = c(t1,t2,3)}
      else if(event2[i] >t2 & event2[i]<=t3){interval_time2[i,1:3] = c(t2,t3,3)}
      else if(event2[i] >t3 & event2[i]<=t4){interval_time2[i,1:3] = c(t3,t4,3)}
      else if(event2[i] >t4 & event2[i]<=t5){interval_time2[i,1:3] = c(t4,t5,3)}
      else {interval_time2[i,1:3] = c(t5,Inf,0)}
    }

    inspection_time3 = matrix(ncol = 5, nrow = augN[3])
    interval_time3 = matrix(ncol = 3, nrow = augN[3])
    colnames(interval_time3) <- c("Left", "Right", "Status")
    for(i in 1:augN[3]){
      t1 <- runif(1,0,1)
      t2 <- t1 + 0.1+runif(1,0,1)
      t3 <- t2 + 0.1+runif(1,0,1)
      t4 <- t3 + 0.1+runif(1,0,1)
      t5 <- t4 + 0.1+runif(1,0,1)
      inspection_time3[i,] <- c(t1,t2,t3,t4,t5)
      tmp <- c(t1,t2,t3,t4,t5)

      # find interval for each subject
      if(event3[i] <= t1) {interval_time3[i,1:3] = c(0,t1,2)}
      else if(event3[i] >t1 & event3[i]<=t2){interval_time3[i,1:3] = c(t1,t2,3)}
      else if(event3[i] >t2 & event3[i]<=t3){interval_time3[i,1:3] = c(t2,t3,3)}
      else if(event3[i] >t3 & event3[i]<=t4){interval_time3[i,1:3] = c(t3,t4,3)}
      else if(event3[i] >t4 & event3[i]<=t5){interval_time3[i,1:3] = c(t4,t5,3)}
      else {interval_time3[i,1:3] = c(t5,Inf,0)}
    }

    # return data
    sample1 <- cbind(X1, interval_time1[,1:3] )
    colnames(sample1)[c(p+1,p+2)] <- c("L", "R")
    sample2 <- cbind(X2, interval_time2[,1:3] )
    colnames(sample2)[c(p+1,p+2)] <- c("L", "R")
    sample3 <- cbind(X3, interval_time3[,1:3] )
    colnames(sample3)[c(p+1,p+2)] <- c("L", "R")

    beta_0 <- rbind(beta_01, beta_02, beta_03)

  }else {
    #while(TRUE){
    eventU <- runif(augN[1],0,1)
    event1 <- -log(eventU)/(exp(X1%*%beta_01))
    # if(sum((log(2)/exp(X1%*%beta_01) - event1)^2)/N[1]/2 <500){break;}
    #  }

    #event1 <- rexp(augN[1],exp(X1%*%beta_01)) # event time
    status1 <- rep(1,length(event1))

    # while(TRUE){
    eventU <- runif(augN[2],0,1)
    event2 <- -log(eventU)/(exp(X2%*%beta_02))
    #  if(sum((log(2)/exp(X2%*%beta_02) - event2)^2)/N[2]/2 <500){break;}
    #  }

    #event2 <- rexp(augN[2],exp(X2%*%beta_02)) # event time
    status2 <-  rep(1,length(event2))


    # while(TRUE){
    eventU <- runif(augN[3],0,1)
    event3 <- -log(eventU)/(exp(X3%*%beta_03))
    #  if(sum((log(2)/exp(X3%*%beta_03) - event3)^2)/N[3]/2 <500){break;}
    #}

    #event3 <- rexp(augN[3],exp(X3%*%beta_03)) # event time
    status3<- rep(1,length(event3))

    # return data
    sample1 <- cbind(X1, event1, status1)
    colnames(sample1)[c(p+1,p+2)] <- c("event", "status")
    sample2 <- cbind(X2, event2,status2)
    colnames(sample2)[c(p+1,p+2)] <- c("event", "status")
    sample3 <- cbind(X3, event3, status3)
    colnames(sample3)[c(p+1,p+2)] <- c("event", "status")
    beta_0 <- rbind(beta_01, beta_02, beta_03)
  }




  return(list(sample1 = sample1, sample2 = sample2, sample3 = sample3,beta_0 = beta_0, p = p,N = N, M = M, q=q,
              correlate = correlate, rightCensor = rightCensor, intervalCensor = intervalCensor, exTx1 = exTx1, exTx2 = exTx2, exTx3=exTx3))
}

# pipeline of SPMA
# Aug 24, 2023 ADD SCAD and MCP
trainSPMA <- function(simulData,iter = 1000, m = 5){
  # fit model
  sample1 = simulData$sample1;
  sample2 = simulData$sample2;
  sample3 = simulData$sample3;
  N = simulData$N;
  p = simulData$p;
  M = simulData$M
  rightCensor = simulData$rightCensor
  intervalCensor = simulData$intervalCensor

  # N is for training data
  #N = augN/2;
  if(rightCensor){

    colnames(sample1) <- c(paste0("V",1:p),"time","status")
    sample1 <- data.frame(sample1)

    X1 = as.matrix(sample1[1:N[1],1:p])
    y1 = cbind(sample1[1:N[1],c(p+1)],sample1[1:N[1],c(p+2)] )
    colnames(y1) <- c("time","status")


    colnames(sample2) <- c(paste0("V",1:p),"time","status")
    sample2 <- data.frame(sample2)
    X2 = as.matrix(sample2[1:N[2],1:p])
    y2 = cbind(sample2[1:N[2],c(p+1)],sample2[1:N[2],c(p+2)] )
    colnames(y2) <- c("time","status")



    colnames(sample3) <- c(paste0("V",1:p),"time","status")
    sample3 <- data.frame(sample3)
    X3 = as.matrix(sample3[1:N[3],1:p])
    y3 = cbind(sample3[1:N[3],c(p+1)],sample3[1:N[3],c(p+2)] )
    colnames(y3) <- c("time","status")

    # unpenalized estimation Cox proportional model for study one
    covF1 = coxph(Surv(time,status)~.,sample1)$var
    theta1F = coxph(Surv(time,status)~.,sample1)$coefficients

    covF2 = coxph(Surv(time,status)~.,sample2)$var
    theta2F = coxph(Surv(time,status)~.,sample2)$coefficients

    covF3 = coxph(Surv(time,status)~.,sample3)$var
    theta3F = coxph(Surv(time,status)~.,sample3)$coefficients

    #View(covF1)
    #print(eigen(covF1))
    covF1 <- solve(covF1)/N[1]
    covF2 <- solve(covF2)/N[2]
    covF3 <- solve(covF3)/N[3]


    # LASSO estimation
    # add cross validation
    #cv_model1 <- cv.glmnet(X1,y1,family = "cox", alpha = 1)
    lasso1<- glmnet(X1,y1, alpha = 1, family = "cox",nlambda =100, keep = TRUE)
    #theta1 = as.vector(coef(best_model1))

    # LASSO Cox proportional model for study two
    #cv_model2 <- cv.glmnet(X2,y2,family = "cox", alpha = 1)
    lasso2<- glmnet(X2,y2, alpha = 1, family = "cox",nlambda =100,keep = TRUE)
    #theta2 = as.vector(coef(best_model2))

    # LASSO Cox proportional model for study three
    #cv_model3 <- cv.glmnet(X3,y3,family = "cox", alpha = 1)
    lasso3<- glmnet(X3,y3, alpha = 1, family = "cox",nlambda = 100, keep =TRUE)
    #theta3 = as.vector(coef(best_model3))

    # SCAD estimation
    scad1 <- ncvsurv(X1,y1,alpha = 1, penalty = "SCAD")
    scad2 <- ncvsurv(X2,y2,alpha = 1, penalty = "SCAD")
    scad3 <- ncvsurv(X3,y3,alpha = 1, penalty = "SCAD")


    # MCP estimation
    mcp1 <- ncvsurv(X1,y1,alpha = 1, penalty = "MCP")
    mcp2 <- ncvsurv(X2,y2,alpha = 1, penalty = "MCP")
    mcp3 <- ncvsurv(X3,y3,alpha = 1, penalty = "MCP")

    # grouplasso
    yy = cbind(c(y1[,1], y2[,1], y3[,1]), c(y1[,2], y2[,2],y3[,2]))
    xx = matrix(0,nrow = (N[1]+N[2]+N[3]), ncol = p*3)
    xx[1:N[1], 1:p] = X1
    xx[(N[1]+1):(N[1]+N[2]), (p+1):(2*p)] =X2
    xx[(N[1]+N[2]+1):(N[1]+N[2]+N[3]), (2*p+1):(3*p)] =X3
    group = c(c(1:p), c(1:p),c(1:p))


    grlasso <- grpsurv(xx,yy,group = group, penalty = "grLasso",alpha = 1)

    grscad <- grpsurv(xx,yy,group = group, penalty = "grSCAD",alpha = 1)

    grmcp <- grpsurv(xx,yy,group = group, penalty = "grMCP",alpha = 1)


    # SPMA
    lambda <- 2^(seq(-30,10,1))
    rho = 10;# range of rho
    red = list()
    for (la in 1:length(lambda)){
      augIMat = matrix(M*p, M*p);
      augIMat = as.matrix(Matrix::bdiag(covF1, covF2, covF3))
      # augNMat = bdiag(200*diag(p),500*diag(p),800*diag(p))
      #print(augNMat%*%augIMat)

      you = spmaC(alpha = 1.2, epi_abs = 10e-5,epi_rel = 10e-3,rho, lambda[la],iter = iter,covF1 = covF1,  covF2 = covF2, covF3 = covF3,augIMat = augIMat, p =p,M=M,setN = N,thetaF = rbind(theta1F,theta2F,theta3F))

      #write.csv(you, "tmp.csv")
      #print(paste0(dim(you$thetaM)," lllll "))
      you$thetaM = you$thetaM[,-c(1:p)]
      # this is not an error, bbk is it+1,
      # but for the one hit the limit, if we want to keep the last one.
      #print(paste0(dim(you$thetaM)," nnnnnn "))
      you$thetaM = array(you$thetaM, dim = c(M,p,iter))
      #print(paste0(dim(you$thetaM)," okok "))
      #print(you$thetaM)
      you$alphaM = you$alphaM[,-c(1:p)]
      you$alphaM = array(you$alphaM, dim = c(M,p,iter))

      you$wM = you$wM[,-c(1:p)]
      red[[la]] = you
    }

    #lasso = rbind(theta1, theta2,theta3)
    #lasso.lambda = c(best_model1$lambda, best_model2$lambda, best_model3$lambda)
    unpel = rbind(theta1F, theta2F, theta3F)
    beta_0 = simulData$beta_0
    return(list(sample1 = sample1, sample2 = sample2, sample3 = sample3, lasso1 = lasso1, lasso2 = lasso2, lasso3 = lasso3,
                scad1 = scad1, scad2 = scad2, scad3 = scad3, mcp1 = mcp1, mcp2 = mcp2, mcp3 = mcp3,
                unpel = unpel, beta_0 = beta_0,
                grlasso = grlasso, grscad = grscad, grmcp = grmcp, spma = red, spma.lambda = lambda,
                N = N, p = p, covF1 =covF1, covF2 = covF2, covF3 = covF3,correlate = simulData$correlate,intervalCensor = simulData$intervalCensor))


  }
  else if(intervalCensor){
    # source two files fun.R, Wu & Cook
    X1 = sample1[1:N[1],1:p]
    X2 = sample2[1:N[2],1:p]
    X3 = sample3[1:N[3],1:p]
    ic1 = sample1[1:N[1],c(p+1,p+2)]
    ic2 = sample2[1:N[2],c(p+1,p+2)]
    ic3 = sample3[1:N[3],c(p+1,p+2)]

    # Unpenalized
    # Wu & Sun
    dataList = list(observation = data.frame(ic1), covdat = data.frame(X1))
    unpel1 = oracle_regression_case2(initial_val=list(beta=rep(0,p),phi=rep(0,m+1)),n = N[1],J = p,p = rep(1,p),dataList,m = m,epsilon = 0.0001,maxiter = 200)
    covF1 = -1*ddl_n_case2(unpel1$phi_hat,dataList,unpel1$beta_hat,m=m,n=N[1])
    theta1F = unpel1$beta_hat

    dataList = list(observation = data.frame(ic2), covdat = data.frame(X2))
    unpel2 = oracle_regression_case2(initial_val=list(beta=rep(0,p),phi=rep(0,m+1)),n = N[2],J = p,p = rep(1,p),dataList,m = m,epsilon = 0.0001,maxiter = 200)
    covF2 = -1*ddl_n_case2(unpel2$phi_hat,dataList,unpel2$beta_hat,m=m,n=N[2])
    theta2F = unpel2$beta_hat


    dataList = list(observation = data.frame(ic3), covdat = data.frame(X3))
    unpel3 = oracle_regression_case2(initial_val=list(beta=rep(0,p),phi=rep(0,m+1)),n = N[3],J = p,p = rep(1,p),dataList,m = m,epsilon = 0.0001,maxiter = 200)
    covF3 = -1*ddl_n_case2(unpel3$phi_hat,dataList,unpel3$beta_hat,m=m,n=N[3])
    theta3F = unpel3$beta_hat

    # fit = ic_sp(cbind(l, u) ~ x1 + x2,
    #             data,
    #             # Need bootstrap samples to compute covariance matrix
    #             bs_samples = 1000)
    # vcov(fit)
    # data = cc$sample1[,-(cc$p+3)]
    # unpel1 = ic_sp(Surv(L,R, type = "interval2")~ ., data = data.frame(data) , model = "ph",
    #                # Need bootstrap samples to compute covariance matrix
    #                bs_samples = 1000)
    # covF1 = vcov(unpel1)
    # covF1 = solve(covF1)/N[1]
    # theta1F = unpel1$coefficients
    #
    # data = cc$sample2[,-(cc$p+3)]
    # unpel2 = ic_sp(Surv(L,R, type = "interval2")~ ., data = data.frame(data) , model = "ph",
    #                # Need bootstrap samples to compute covariance matrix
    #                bs_samples = 1000)
    # covF2 = vcov(unpel2)
    # covF2 = solve(covF2)/N[2]
    # theta2F = unpel2$coefficients
    #
    # data = cc$sample3[,-(cc$p+3)]
    # unpel3 = ic_sp(Surv(L,R, type = "interval2")~ ., data = data.frame(data) , model = "ph",
    #                # Need bootstrap samples to compute covariance matrix
    #                bs_samples = 1000)
    # covF3 = vcov(unpel3)
    # covF3 = solve(covF3)/N[3]
    # theta3F = unpel3$coefficients
    ### Note: Sep 28, estimation of covF1 covF2, covF3 has large impact on the estimation of theta.spma
    ### If covF1, covF2, covF3 are wrong, ADMM wont converge.
    ### information matrix I(theta) = - second derivative of loglikelihood
    ### var(theta_hat) = I^-1(theta)
    ### SPMA asks for (var(theta_hat))^-1 = I(theta) = - second derivative of loglikelihood

    print("unpel is done")
    # LASSO
    # Wu & Cook
    if(FALSE){
      # study one
      indata = data.frame(cbind( 1:nrow(X1),  ic1[,1], ic1[,2], X1,  sample1[1:N[1],p+3] ))
      colnames(indata) = c("id", "timeL", "timeR", paste0("x",1:ncol(X1)), "status")

      # decide break point
      fit0 <- survfit(Surv(timeL, timeR, status, type='interval')~1, data = indata)
      surv.prob <- summary(fit0)$surv;
      npieces <- 1
      ncov <- p
      probs <- seq(from = 1-max(surv.prob), to = 1-min(surv.prob), length.out = npieces + 1)
      probs <- probs[-c(1, npieces+1)]
      cpoints <- quantile(fit0, probs = probs, conf.int = FALSE)
      cutpoints <- data.frame(start=c(0,cpoints),
                              stop=c(cpoints,9999),
                              piece=1:(length(cpoints)+1))

      lambda.lasso <- 10^seq(-2,1,length.out = 5)
      ss = 1
      lasso1 = list()
      scad1 = list()
      bic.lasso1 = c()
      bic.scad1 = c()
      for( lam in lambda.lasso){
        lasso1.tmp <- EM.f(indata = indata, lam0 = rep(1, npieces), beta0 = rep(0, ncov),
                           lasso.lam = lam, ncov = ncov, npieces = npieces,
                           cutpoints = cutpoints, penalty.function = "lasso")

        scad1.tmp <- EM.f(indata = indata, lam0 = rep(1, npieces), beta0 = rep(0, ncov),
                          lasso.lam = lam, ncov = ncov, npieces = npieces,
                          cutpoints = cutpoints, penalty.function = "scad")

        lasso1[[ss]] = lasso1.tmp
        scad1[[ss]] = scad1.tmp
        # Calculate bic
        bic.lasso1[ss] = -2*logL.f(ic1,Z = X1, lam = lasso1.tmp$lam, beta= lasso1.tmp$beta,ncov = p, cutpoints = cutpoints )/log(N[1]) + log(N[1])*sum(lasso1.tmp$beta!=0)
        bic.scad1[ss] = -2*logL.f(ic1,Z = X1, lam = scad1.tmp$lam, beta= scad1.tmp$beta,ncov = p, cutpoints = cutpoints )/log(N[1]) + log(N[1])*sum(scad1.tmp$beta!=0)

        ss =ss +1
      }

      # Calculate bic
      # bic.lasso = -2*logL.f(ic1,Z = X1, lam = lasso1$lam, beta= lasso1$beta,ncov = p, cutpoints = cutpoints )/log(N[1]) + log(N[1])*sum(lasso1$beta!=0)

      # study 2
      indata = data.frame(cbind( 1:nrow(X2),  ic2[,1], ic2[,2], X2,  sample1[1:N[2],p+3] ))
      colnames(indata) = c("id", "timeL", "timeR", paste0("x",1:ncol(X2)), "status")

      # decide break point
      fit0 <- survfit(Surv(timeL, timeR, status, type='interval')~1, data = indata)
      surv.prob <- summary(fit0)$surv;
      npieces <- 1
      ncov <- p
      probs <- seq(from = 1-max(surv.prob), to = 1-min(surv.prob), length.out = npieces + 1)
      probs <- probs[-c(1, npieces+1)]
      cpoints <- quantile(fit0, probs = probs, conf.int = FALSE)
      cutpoints <- data.frame(start=c(0,cpoints),
                              stop=c(cpoints,9999),
                              piece=1:(length(cpoints)+1))

      ss = 1
      lasso2 = list()
      scad2 = list()
      bic.lasso2 = c()
      bic.scad2 = c()
      for( lam in lambda.lasso){
        lasso2.tmp <- EM.f(indata = indata, lam0 = rep(1, npieces), beta0 = rep(0, ncov),
                           lasso.lam = lam, ncov = ncov, npieces = npieces,
                           cutpoints = cutpoints, penalty.function = "lasso")
        scad2.tmp <- EM.f(indata = indata, lam0 = rep(1, npieces), beta0 = rep(0, ncov),
                          lasso.lam = lam, ncov = ncov, npieces = npieces,
                          cutpoints = cutpoints, penalty.function = "scad")
        lasso2[[ss]] = lasso1.tmp
        scad2[[ss]] = scad1.tmp
        # Calculate bic
        bic.lasso2[ss] = -2*logL.f(ic2,Z = X2, lam = lasso2.tmp$lam, beta= lasso2.tmp$beta,ncov = p, cutpoints = cutpoints )/log(N[2]) + log(N[2])*sum(lasso2.tmp$beta!=0)
        bic.scad2[ss] = -2*logL.f(ic2,Z = X2, lam = scad2.tmp$lam, beta= scad2.tmp$beta,ncov = p, cutpoints = cutpoints )/log(N[2]) + log(N[2])*sum(scad2.tmp$beta!=0)

        ss =ss +1
      }
      # study 3
      indata = data.frame(cbind( 1:nrow(X2),  ic2[,1], ic2[,2], X2,  sample1[1:N[2],p+3] ))
      colnames(indata) = c("id", "timeL", "timeR", paste0("x",1:ncol(X2)), "status")

      # decide break point
      fit0 <- survfit(Surv(timeL, timeR, status, type='interval')~1, data = indata)
      surv.prob <- summary(fit0)$surv;
      npieces <- 1
      ncov <- p
      probs <- seq(from = 1-max(surv.prob), to = 1-min(surv.prob), length.out = npieces + 1)
      probs <- probs[-c(1, npieces+1)]
      cpoints <- quantile(fit0, probs = probs, conf.int = FALSE)
      cutpoints <- data.frame(start=c(0,cpoints),
                              stop=c(cpoints,9999),
                              piece=1:(length(cpoints)+1))
      ss = 1
      lasso3 = list()
      scad3 = list()
      bic.lasso3 = c()
      bic.scad3 = c()
      for (lam in lambda.lasso){
        lasso3.tmp <- EM.f(indata = indata, lam0 = rep(1, npieces), beta0 = rep(0, ncov),
                           lasso.lam = lam, ncov = ncov, npieces = npieces,
                           cutpoints = cutpoints, penalty.function = "lasso")
        scad3.tmp <- EM.f(indata = indata, lam0 = rep(1, npieces), beta0 = rep(0, ncov),
                          lasso.lam = lam, ncov = ncov, npieces = npieces,
                          cutpoints = cutpoints, penalty.function = "scad")
        lasso3[[ss]] = lasso3.tmp
        scad3[[ss]] = scad3.tmp
        # Calculate bic
        bic.lasso3[ss] = -2*logL.f(ic3,Z = X3, lam = lasso3.tmp$lam, beta= lasso3.tmp$beta,ncov = p, cutpoints = cutpoints )/log(N[3]) + log(N[3])*sum(lasso3.tmp$beta!=0)
        bic.scad3[ss] = -2*logL.f(ic3,Z = X3, lam = scad3.tmp$lam, beta= scad3.tmp$beta,ncov = p, cutpoints = cutpoints )/log(N[3]) + log(N[3])*sum(scad3.tmp$beta!=0)
        ss =ss +1
      }
      print("lasso scad are done")
      # Data preparation for Group Lasso, Group SCAD, Group MCP
      xx = matrix(0,nrow = (N[1]+N[2]+N[3]), ncol = p*3)
      xx[1:N[1], 1+3*c(0:(p-1))] = X1
      xx[(N[1]+1):(N[1]+N[2]), 2+3*c(0:(p-1))] =X2
      xx[(N[1]+N[2]+1):(N[1]+N[2]+N[3]), 3+3*c(0:(p-1))] =X3

      # intervals
      ic = rbind(ic1,ic2,ic3)
      colnames(ic) <- c("L", "R")

      dataList = list(observation = data.frame(ic), covdat = data.frame(xx))
      # lambda = seq(0,50,length.out = 10)
      lambda = 2^seq(-1,6,length.out = 10)
      #m = 5 # tuning parameter

      # Group Lasso
      grlasso = list()
      grlasso.bic = c()
      it = 1
      for(la in lambda){
        tmp = lasso_regression_group_case2(initial_val=list(beta=rep(0,p*3),phi=rep(0,m+1)),n = sum(N),J = p,p = rep(M,p),dataList,m = m,lambda = la,epsilon = 0.0001,maxiter = 200)
        grlasso[[it]] = tmp
        grlasso.bic[it] = -2*l_n_case2(tmp$phi_hat,dataList,tmp$beta_hat,m,sum(N))+sum(abs(tmp$beta_hat)>0)*log(sum(N))
        it = it+1
      }
      print("group lasso is done")

      # Group SCAD
      grscad = list()
      grscad.bic = c()
      it = 1
      for(la in lambda){
        tmp = scad_regression_group_case2(initial_val=list(beta=rep(0,p*3),phi=rep(0,m+1)),n = sum(N),J = p,p = rep(M,p),dataList,m = m,gamma = 4,lambda = la,epsilon = 0.0001,maxiter = 200)
        grscad[[it]] = tmp
        grscad.bic[it] = -2*l_n_case2(tmp$phi_hat,dataList,tmp$beta_hat,m,sum(N))+sum(abs(tmp$beta_hat)>0)*log(sum(N))
        it = it+1
      }
      print("group scad is done")

      # Group MCP
      grmcp = list()
      grmcp.bic = c()
      it = 1 # iter
      for(la in lambda){
        tmp = mcp_regression_group_case2(initial_val=list(beta=rep(0,p*3),phi=rep(0,m+1)),n = sum(N),J = p,p = rep(M,p),dataList,m = m,gamma = 3,lambda = la,epsilon = 0.0001,maxiter = 200)
        grmcp[[it]] = tmp
        grmcp.bic[it] = -2*l_n_case2(tmp$phi_hat,dataList,tmp$beta_hat,m,sum(N))+sum(abs(tmp$beta_hat)>0)*log(sum(N))
        it = it+1
      }
      print("group mcp is done")
    }
    # SPMA
    lambda.spma <- 2^(seq(-30,10,1))
    red = list()
    #print(rbind(theta1F,theta2F,theta3F))
    for (la in 1:length(lambda.spma)){
      augIMat = matrix(M*p, M*p);
      augIMat = as.matrix(Matrix::bdiag(covF1, covF2, covF3))
      # augNMat = bdiag(200*diag(p),500*diag(p),800*diag(p))
      #print(augNMat%*%augIMat)
      rho = la;# range of rho
      you = spmaC(alpha = 1.5, epi_abs = 10e-3,epi_rel = 10e-2,rho, lambda.spma[la],iter = iter,covF1 = covF1,  covF2 = covF2, covF3 = covF3,augIMat = augIMat, p =p,M=M,setN = N,thetaF = rbind(theta1F,theta2F,theta3F))
      # epi_abs = 10e-2, epi_rel = 10e-1 achieves 100% selection accuracy.
      # tried with 10e-5 10e-3, selection accuracy is not good
      # tried 10e-4 10e-2, estimation is improved, selection is not
      # try 10e-3 10e-2, the 10 rounds performance on my laptop is good, but it is not on the server.
      # try rho = lambda, alpha = 1.5, 10e-3, 10e-2 seems to be good with 10 replicates
      #write.csv(you, "tmp.csv")
      #print(paste0(dim(you$thetaM)," lllll "))
      you$thetaM = you$thetaM[,-c(1:p)]
      # this is not an error, bbk is it+1,
      # but for the one hit the limit, if we want to keep the last one.
      #print(paste0(dim(you$thetaM)," nnnnnn "))
      you$thetaM = array(you$thetaM, dim = c(M,p,iter))
      #print(paste0(dim(you$thetaM)," okok "))
      #print(you$thetaM)
      you$alphaM = you$alphaM[,-c(1:p)]
      you$alphaM = array(you$alphaM, dim = c(M,p,iter))

      you$wM = you$wM[,-c(1:p)]
      red[[la]] = you
    }
    print("spma is done")

    unpel = rbind(theta1F, theta2F, theta3F)
    # return(list(sample1 = sample1, sample2 = sample2, sample3 = sample3,
    #             unpel = unpel, beta_0 =simulData$beta_0,
    #             lasso1 = lasso1, lasso2 = lasso2, lasso3 = lasso3,
    #             scad1 = scad1, scad2 = scad2, scad3 = scad3,
    #             bic.lasso1 = bic.lasso1, bic.lasso2 = bic.lasso2, bic.lasso3 = bic.lasso3,
    #             bic.scad1 = bic.scad1, bic.scad2 = bic.scad2, bic.scad3 = bic.scad3,
    #             grlasso = grlasso, grlasso.bic = grlasso.bic,
    #             grscad = grscad, grscad.bic = grscad.bic,
    #             grmcp = grmcp, grmcp.bic = grmcp.bic, spma = red, lambda.group = lambda, spma.lambda = lambda.spma,
    #             N=N, p = p, covF1 = covF1, covF2 = covF2, covF3 = covF3, correlate = simulData$correlate,intervalCensor = simulData$intervalCensor))
    # #
    # When to test SPMA
    return(list(sample1 = sample1, sample2 = sample2, sample3 = sample3,
                unpel = unpel, beta_0 =simulData$beta_0,intervalCensor = simulData$intervalCensor,spma = red,
                spma.lambda = lambda.spma,N=N, p = p, covF1 = covF1, covF2 = covF2, covF3 = covF3,
                correlate = simulData$correlate))
  }
  else{
    colnames(sample1) <- c(paste0("V",1:p),"time","status")
    sample1 <- data.frame(sample1)
    X1 = as.matrix(sample1[1:N[1],1:p])
    y1 = cbind(sample1[1:N[1],c(p+1)],sample1[1:N[1],c(p+2)] )
    colnames(y1) <- c("time","status")


    colnames(sample2) <- c(paste0("V",1:p),"time","status")
    sample2 <- data.frame(sample2)
    X2 = as.matrix(sample2[1:N[2],1:p])
    y2 = cbind(sample2[1:N[2],c(p+1)],sample2[1:N[2],c(p+2)] )
    colnames(y2) <- c("time","status")



    colnames(sample3) <- c(paste0("V",1:p),"time","status")
    sample3 <- data.frame(sample3)
    X3 = as.matrix(sample3[1:N[3],1:p])
    y3 = cbind(sample3[1:N[3],c(p+1)],sample3[1:N[3],c(p+2)] )
    colnames(y3) <- c("time","status")

    # unpenalized estimation Cox proportional model for study one
    covF1 = coxph(Surv(time,status)~.,sample1)$var
    theta1F = coxph(Surv(time,status)~.,sample1)$coefficients

    covF2 = coxph(Surv(time,status)~.,sample2)$var
    theta2F = coxph(Surv(time,status)~.,sample2)$coefficients

    covF3 = coxph(Surv(time,status)~.,sample3)$var
    theta3F = coxph(Surv(time,status)~.,sample3)$coefficients


    covF1 <- solve(covF1)/N[1]
    covF2 <- solve(covF2)/N[2]
    covF3 <- solve(covF3)/N[3]


    # LASSO estimation
    # add cross validation
    #cv_model1 <- cv.glmnet(X1,y1,family = "cox", alpha = 1)
    lasso1<- glmnet(X1,y1, alpha = 1, family = "cox",nlambda =100, keep = TRUE)
    #theta1 = as.vector(coef(best_model1))

    # LASSO Cox proportional model for study two
    #cv_model2 <- cv.glmnet(X2,y2,family = "cox", alpha = 1)
    lasso2<- glmnet(X2,y2, alpha = 1, family = "cox",nlambda =100,keep = TRUE)
    #theta2 = as.vector(coef(best_model2))

    # LASSO Cox proportional model for study three
    #cv_model3 <- cv.glmnet(X3,y3,family = "cox", alpha = 1)
    lasso3<- glmnet(X3,y3, alpha = 1, family = "cox",nlambda = 100, keep =TRUE)
    #theta3 = as.vector(coef(best_model3))

    # SCAD estimation
    scad1 <- ncvsurv(X1,y1,alpha = 1, penalty = "SCAD")
    scad2 <- ncvsurv(X2,y2,alpha = 1, penalty = "SCAD")
    scad3 <- ncvsurv(X3,y3,alpha = 1, penalty = "SCAD")


    # MCP estimation
    mcp1 <- ncvsurv(X1,y1,alpha = 1, penalty = "MCP")
    mcp2 <- ncvsurv(X2,y2,alpha = 1, penalty = "MCP")
    mcp3 <- ncvsurv(X3,y3,alpha = 1, penalty = "MCP")

    # grouplasso
    yy = cbind(c(y1[,1], y2[,1], y3[,1]), c(y1[,2], y2[,2],y3[,2]))
    xx = matrix(0,nrow = (N[1]+N[2]+N[3]), ncol = p*3)
    xx[1:N[1], 1:p] = X1
    xx[(N[1]+1):(N[1]+N[2]), (p+1):(2*p)] =X2
    xx[(N[1]+N[2]+1):(N[1]+N[2]+N[3]), (2*p+1):(3*p)] =X3
    group = c(c(1:p), c(1:p),c(1:p))


    grlasso <- grpsurv(xx,yy,group = group, penalty = "grLasso",alpha = 1)

    grscad <- grpsurv(xx,yy,group = group, penalty = "grSCAD",alpha = 1)

    grmcp <- grpsurv(xx,yy,group = group, penalty = "grMCP",alpha = 1)

    # SPMA
    lambda <- 2^(seq(-30,10,1))
    rho = 10;# range of rho
    red = list()
    for (la in 1:length(lambda)){
      augIMat = matrix(M*p, M*p);
      augIMat = as.matrix(Matrix::bdiag(covF1, covF2, covF3))
      # augNMat = bdiag(200*diag(p),500*diag(p),800*diag(p))
      #print(augNMat%*%augIMat)

      you = spmaC(alpha = 1.2, epi_abs = 10e-5,epi_rel = 10e-3,rho, lambda[la],iter = iter,covF1 = covF1,  covF2 = covF2, covF3 = covF3,augIMat = augIMat, p =p,M=M,setN = N,thetaF = rbind(theta1F,theta2F,theta3F))

      #write.csv(you, "tmp.csv")
      #print(paste0(dim(you$thetaM)," lllll "))
      you$thetaM = you$thetaM[,-c(1:p)]
      # this is not an error, bbk is it+1,
      # but for the one hit the limit, if we want to keep the last one.
      #print(paste0(dim(you$thetaM)," nnnnnn "))
      you$thetaM = array(you$thetaM, dim = c(M,p,iter))
      #print(paste0(dim(you$thetaM)," okok "))
      #print(you$thetaM)
      you$alphaM = you$alphaM[,-c(1:p)]
      you$alphaM = array(you$alphaM, dim = c(M,p,iter))

      you$wM = you$wM[,-c(1:p)]
      red[[la]] = you
    }

    #lasso = rbind(theta1, theta2,theta3)
    #lasso.lambda = c(best_model1$lambda, best_model2$lambda, best_model3$lambda)
    unpel = rbind(theta1F, theta2F, theta3F)
    beta_0 = simulData$beta_0
    return(list(sample1 = sample1, sample2 = sample2, sample3 = sample3, lasso1 = lasso1, lasso2 = lasso2, lasso3 = lasso3,
                scad1 = scad1, scad2 = scad2, scad3 = scad3, mcp1 = mcp1, mcp2 = mcp2, mcp3 = mcp3,
                unpel = unpel, beta_0 = beta_0,
                grlasso = grlasso, grscad = grscad, grmcp = grmcp, spma = red, spma.lambda = lambda,
                N = N, p = p, covF1 =covF1, covF2 = covF2, covF3 = covF3,correlate = simulData$correlate,intervalCensor = simulData$intervalCensor))

  }


}

# T_i is event time
bayes_bic_spma <- function(thetaF, thetahat, covF1, covF2, covF3, N, p ){
  if(nrow(thetahat)!=nrow(thetaF)){print("Parameter dimension should be the same"); break; }
  #print(dim( as.matrix(thetahat[1,] - thetaF[1,], nrow =1,ncol = p)))
  #print(dim(covF1))
  p1 = t(as.matrix(thetahat[1,] - thetaF[1,], nrow =1,ncol = p))%*% covF1 %*% as.matrix(thetahat[1,] - thetaF[1,])
  p1 = p1+t(as.matrix(thetahat[2,] - thetaF[2,], nrow =1,ncol = p))%*% covF2 %*% as.matrix(thetahat[2,] - thetaF[2,])
  p1 = p1+t(as.matrix(thetahat[3,] - thetaF[3,], nrow =1,ncol = p))%*% covF3 %*% as.matrix(thetahat[3,] - thetaF[3,])
  p2 = log(sum(N))*sum(thetahat!=0) # decide the form of n_m or n
  return (list(p1 = p1,p2=p2))
}

bic_non_spma <- function(set.n, thetahat, sample1, sample2, sample3, p, M = 3,intervalCensor = FALSE){

  # first component is still 500*1 by setting nrow =1, to transpose
  p1 = t(as.matrix(sample1[,p+2], nrow = 1,ncol = set.n[1]))%*%(as.matrix(sample1[,1:p])%*%as.matrix(thetahat[1:p], ncol =1))
  p1 = p1+ t(as.matrix(sample2[,p+2], nrow = 1, ncol = set.n[2]))%*%(as.matrix(sample2[,1:p])%*%as.matrix(thetahat[(p+1):(2*p)],ncol =1))
  p1 = p1+ t(as.matrix(sample3[,p+2], nrow = 1, ncol = set.n[3]))%*%(as.matrix(sample3[,1:p])%*%as.matrix(thetahat[(2*p+1):(3*p)],ncol =1))

  sample1 = sample1[order(sample1[,p+1], decreasing = FALSE), ]
  delta1 = sample1[,p+2]
  prod = exp(as.matrix(sample1[,1:p])%*%as.matrix(thetahat[1:p], ncol =1))
  for(i in 1:(set.n[1])){
    tt =  -1*log(sum(prod[i:set.n[1]]))
    p1 = p1 + delta1[i]*tt
  }

  sample2 = sample2[order(sample2[,p+1], decreasing = FALSE), ]
  delta2 = sample2[,p+2]
  prod = exp(as.matrix(sample2[,1:p])%*%as.matrix(thetahat[(p+1):(2*p)], ncol =1))
  for(i in 1:(set.n[2])){
    tt =  -1*log(sum(prod[i:set.n[2]]))
    p1 = p1 + delta2[i]*tt
  }

  sample3 = sample3[order(sample3[,p+1], decreasing = FALSE), ]
  delta3 = sample3[,p+2]
  prod = exp(as.matrix(sample3[,1:p])%*%as.matrix(thetahat[(2*p+1):(3*p)], ncol =1))
  for(i in 1:(set.n[3])){
    tt =  -1*log(sum(prod[i:set.n[3]]))
    p1 = p1 + delta3[i]*tt
  }


  p1 = (-2)*p1
  #p2 = log(sum(set.n)) / sum(set.n) * log(sum(thetahat!=0) )
  p2 = log(sum(set.n)) *(sum(thetahat!=0)) # ?
  return (list(p1 = p1,p2=p2))
  #return(p1+p2)


}

bic_lasso <- function(ni, thetahat, saple, p){
  # first coponent is still 500*1 by setting nrow =1, to transpose
  p1 = t(as.atrix(saple[,p+2], nrow = 1,ncol = ni))%*%
    (as.matrix(sample[,1:p])%*%as.matrix(thetahat[1:p], ncol =1))

  sample = sample[order(sample[,p+1], decreasing = FALSE), ]
  delta = sample[,p+2]
  prod = exp(as.matrix(sample[,1:p])%*%as.matrix(thetahat[1:p], ncol =1))
  for(i in 1:ni){
    tt =  -1*log(sum(prod[i:ni]))
    p1 = p1 + delta[i]*tt
  }

  p1 = (-2)*p1
  #p2 = log(sum(set.n)) / sum(set.n) * log(sum(thetahat!=0) )
  p2 = log(ni)* (sum(thetahat!=0)) # ?
  return (list(p1 = p1,p2=p2))
  #return(p1+p2)
}

bic_model_selection <- function(trainSPMA.object){

  intervalCensor = trainSPMA.object$intervalCensor
  if(intervalCensor){
    # unpel, lasso, scad is free of selecting
    idx.lasso1= which.min(trainSPMA.object$bic.lasso1)
    bic.lasso1 = trainSPMA.object$bic.lasso1
    idx.lasso2= which.min(trainSPMA.object$bic.lasso2)
    bic.lasso2 = trainSPMA.object$bic.lasso2
    idx.lasso3= which.min(trainSPMA.object$bic.lasso3)
    bic.lasso3 = trainSPMA.object$bic.lasso3

    idx.scad1 = which.min(trainSPMA.object$bic.scad1)
    bic.scad1 = trainSPMA.object$bic.scad1
    idx.scad2 = which.min(trainSPMA.object$bic.scad2)
    bic.scad2 = trainSPMA.object$bic.scad2
    idx.scad3 = which.min(trainSPMA.object$bic.scad3)
    bic.scad3 = trainSPMA.object$bic.scad3
    # Group Lasso
    idx.grlasso = which.min(trainSPMA.object$grlasso.bic)
    bic.grlasso = trainSPMA.object$grlasso.bic

    idx.grscad = which.min(trainSPMA.object$grscad.bic)
    bic.grscad = trainSPMA.object$grscad.bic

    idx.grmcp = which.min(trainSPMA.object$grmcp.bic)
    bic.grmcp = trainSPMA.object$grmcp.bic


    return (list(idx.lasso1 = idx.lasso1, idx.lasso2 = idx.lasso2, idx.lasso3 = idx.lasso3,
                 idx.scad1 = idx.scad1, idx.scad2 = idx.scad2, idx.scad3 = idx.scad3,
                 #idx.mcp1 = idx.mcp1, idx.mcp2 = idx.mcp2, idx.mcp3 = idx.mcp3,
                 idx.grlasso = idx.grlasso, idx.grscad = idx.grscad,idx.grmcp = idx.grmcp,
                 bic.lasso1 = bic.lasso1, bic.lasso2 = bic.lasso2, bic.lasso3 = bic.lasso3,
                 bic.scad1 = bic.scad1, bic.scad2 = bic.scad2, bic.scad3 = bic.scad3,
                 bic.grlasso = bic.grlasso,bic.grscad = bic.grscad, bic.grmcp = bic.grmcp))
  }
  else{
    # 100 is numebr of lambda in group
    # parameter extract
    N = trainSPMA.object$N
    sample1 = trainSPMA.object$sample1[1:N[1],];
    sample2 = trainSPMA.object$sample2[1:N[2],];
    sample3 = trainSPMA.object$sample3[1:N[3],];
    p = trainSPMA.object$p

    # LASSO sample 1
    thetahat.lasso1 =as.matrix(trainSPMA.object$lasso1$beta)
    num.col = ncol(thetahat.lasso1)
    bic.lasso1 = cbind(rep(0,num.col), rep(0,num.col))

    for(i in 1:num.col){
      bb = bic_lasso(N[1],thetahat = thetahat.lasso1[,i],sample1,p, M = 3)
      bic.lasso1[i,1:2] = c(bb$p1,bb$p2)
    }
    bic.lasso1 = bic.lasso1[-1,]
    idx.lasso1 = which.min(bic.lasso1[,1]+bic.lasso1[,2])+1

    # LASSO sample 2
    thetahat.lasso2 =as.matrix(trainSPMA.object$lasso2$beta)
    num.col = ncol(thetahat.lasso2)
    bic.lasso2 = cbind(rep(0,num.col), rep(0,num.col))

    for(i in 1:num.col){
      bb = bic_lasso(N[2],thetahat = thetahat.lasso2[,i],sample2,p, M = 3)
      bic.lasso2[i,1:2] = c(bb$p1,bb$p2)
    }
    bic.lasso2 = bic.lasso2[-1,]
    idx.lasso2 = which.min(bic.lasso2[,1]+bic.lasso2[,2])+1


    # LASSO sample 3
    thetahat.lasso3 =as.matrix(trainSPMA.object$lasso3$beta)
    num.col = ncol(thetahat.lasso3)
    bic.lasso3 = cbind(rep(0,num.col), rep(0,num.col))

    for(i in 1:num.col){
      bb = bic_lasso(N[3],thetahat = thetahat.lasso3[,i],sample3,p, M = 3)
      bic.lasso3[i,1:2] = c(bb$p1,bb$p2)
    }
    bic.lasso3 = bic.lasso3[-1,]
    idx.lasso3 = which.min(bic.lasso3[,1]+bic.lasso3[,2])+1

    # SCAD
    thetahat.scad1 =as.matrix(trainSPMA.object$scad1$beta)
    num.col = ncol(thetahat.scad1)
    bic.scad1 = cbind(rep(0,num.col), rep(0,num.col))

    for(i in 1:num.col){
      bb = bic_lasso(N[1],thetahat = thetahat.scad1[,i],sample1,p, M = 3)
      bic.scad1[i,1:2] = c(bb$p1,bb$p2)
    }
    bic.scad1 = bic.scad1[-1,]
    idx.scad1 = which.min(bic.scad1[,1]+bic.scad1[,2])+1

    # group 2
    thetahat.scad2 =as.matrix(trainSPMA.object$scad2$beta)
    num.col = ncol(thetahat.scad2)
    bic.scad2 = cbind(rep(0,num.col), rep(0,num.col))

    for(i in 1:num.col){
      bb = bic_lasso(N[2],thetahat = thetahat.scad2[,i],sample2,p, M = 3)
      bic.scad2[i,1:2] = c(bb$p1,bb$p2)
    }
    bic.scad2 = bic.scad2[-1,]
    idx.scad2 = which.min(bic.scad2[,1]+bic.scad2[,2])+1

    # group 3
    thetahat.scad3 =as.matrix(trainSPMA.object$scad3$beta)
    num.col = ncol(thetahat.scad3)
    bic.scad3 = cbind(rep(0,num.col), rep(0,num.col))

    for(i in 1:num.col){
      bb = bic_lasso(N[3],thetahat = thetahat.scad3[,i],sample3,p, M = 3)
      bic.scad3[i,1:2] = c(bb$p1,bb$p2)
    }
    bic.scad3 = bic.scad3[-1,]
    idx.scad3 = which.min(bic.scad3[,1]+bic.scad3[,2])+1

    # MCP
    thetahat.mcp1 =as.matrix(trainSPMA.object$mcp1$beta)
    num.col = ncol(thetahat.mcp1)
    bic.mcp1 = cbind(rep(0,num.col), rep(0,num.col))

    for(i in 1:num.col){
      bb = bic_lasso(N[1],thetahat = thetahat.mcp1[,i],sample1,p, M = 3)
      bic.mcp1[i,1:2] = c(bb$p1,bb$p2)
    }
    bic.mcp1 = bic.mcp1[-1,]
    idx.mcp1 = which.min(bic.mcp1[,1]+bic.mcp1[,2])+1

    # group 2
    thetahat.mcp2 =as.matrix(trainSPMA.object$mcp2$beta)
    num.col = ncol(thetahat.mcp2)
    bic.mcp2 = cbind(rep(0,num.col), rep(0,num.col))

    for(i in 1:num.col){
      bb = bic_lasso(N[2],thetahat = thetahat.mcp2[,i],sample2,p, M = 3)
      bic.mcp2[i,1:2] = c(bb$p1,bb$p2)
    }
    bic.mcp2 = bic.mcp2[-1,]
    idx.mcp2 = which.min(bic.mcp2[,1]+bic.mcp2[,2])+1

    # group 3
    thetahat.mcp3 =as.matrix(trainSPMA.object$mcp3$beta)
    num.col = ncol(thetahat.mcp3)
    bic.mcp3 = cbind(rep(0,num.col), rep(0,num.col))

    for(i in 1:num.col){
      bb = bic_lasso(N[3],thetahat = thetahat.mcp3[,i],sample3,p, M = 3)
      bic.mcp3[i,1:2] = c(bb$p1,bb$p2)
    }
    bic.mcp3 = bic.mcp3[-1,]
    idx.mcp3 = which.min(bic.mcp3[,1]+bic.mcp3[,2])+1


    # GROUP LASSO
    #
    thetahat.grlasso = trainSPMA.object$grlasso$beta
    ncolb = ncol(thetahat.grlasso)
    bic.grlasso = cbind(rep(0,ncolb),rep(0,ncolb))

    for(i in 1:ncolb){
      bb = bic_non_spma(set.n = N,thetahat = thetahat.grlasso[,i],sample1,sample2,sample3,p, M = 3)
      bic.grlasso[i,1:2] = c(bb$p1,bb$p2)
    }
    bic.grlasso = bic.grlasso[-1,]
    idx.grlasso = which.min(bic.grlasso[,1]+bic.grlasso[,2])+1

    rm(ncolb)
    # GROUP SCAD
    thetahat.grscad = trainSPMA.object$grscad$beta
    ncolb = ncol(thetahat.grscad)
    bic.grscad = cbind(rep(0,ncolb),rep(0,ncolb))

    for(i in 1:ncolb){
      bb = bic_non_spma(N,thetahat = thetahat.grscad[,i],sample1,sample2,sample3,p, M = 3)
      bic.grscad[i,1:2] = c(bb$p1,bb$p2)
    }
    bic.grscad = bic.grscad[-1,]
    idx.grscad = which.min(bic.grscad[,1]+bic.grscad[,2])+1
    rm(ncolb)
    # GROUP MCP
    thetahat.grmcp = trainSPMA.object$grmcp$beta
    ncolb =ncol(thetahat.grmcp)
    bic.grmcp = cbind(rep(0,ncolb),rep(0,ncolb))

    for(i in 1:ncolb){
      bb = bic_non_spma(N,thetahat = thetahat.grmcp[,i],sample1,sample2,sample3,p, M = 3)
      bic.grmcp[i,1:2] = c(bb$p1,bb$p2)
    }
    bic.grmcp = bic.grmcp[-1,]
    idx.grmcp = which.min(bic.grmcp[,1]+bic.grmcp[,2])+1

    return (list(idx.lasso1 = idx.lasso1, idx.lasso2 = idx.lasso2, idx.lasso3 = idx.lasso3,
                 idx.scad1 = idx.scad1, idx.scad2 = idx.scad2, idx.scad3 = idx.scad3,
                 idx.mcp1 = idx.mcp1, idx.mcp2 = idx.mcp2, idx.mcp3 = idx.mcp3,
                 idx.grmcp = idx.grmcp, idx.grlasso = idx.grlasso, idx.grscad = idx.grscad,
                 bic.lasso1 = bic.lasso1, bic.lasso2 = bic.lasso2, bic.lasso3 = bic.lasso3,
                 bic.grlasso = bic.grlasso,bic.grscad = bic.grscad, bic.grmcp = bic.grmcp))
  }
}


select_model_spma <- function(trainSPMA.object){
  n = length(trainSPMA.object$spma)
  bic.spma = cbind(rep(0,n), rep(0,n))
  for(i in 1:n){
    #  print(i);
    # bb = bayes_bic_spma(trainSPMA.object$unpel,trainSPMA.object$spma[[i]]$thetaM[,,trainSPMA.object$spma[[i]]$bbk],
    #                     trainSPMA.object$covF1,trainSPMA.object$covF2,trainSPMA.object$covF3,trainSPMA.object$N)
    bb = bayes_bic_spma(trainSPMA.object$unpel,trainSPMA.object$spma[[i]]$alphaM[,,trainSPMA.object$spma[[i]]$bbk],
                        trainSPMA.object$covF1,trainSPMA.object$covF2,trainSPMA.object$covF3,trainSPMA.object$N)
    bic.spma[i,1:2] = c(bb$p1, bb$p2)
  }
  idx.spma = which.min(bic.spma[,1]+bic.spma[,2])
  return(list(idx.spma = idx.spma, bic.spma = bic.spma))
}


dist2beta0 <- function(beta0, betahat,p, M ){
  #print(betahat)
  betahat = matrix(betahat, nrow = M, ncol = p, byrow = TRUE)
  if(nrow(beta0)!= nrow(betahat)){ print("Dimension must be same"); break; }
  tmp = sum((beta0 - betahat)^2)
  return (tmp)
}

# 07/10/2023 Update confuse matrix 0 and 1
rmse <- function(trainSPMA.object,best.idx, spma.bic,Data){

  # 'Data' provide exTx1 exTx2 exTx3
  intervalCensor = trainSPMA.object$intervalCensor
  beta0 = trainSPMA.object$beta_0
  binbeta0 = ifelse(beta0==0,0,1)
  N = trainSPMA.object$N
  # use test set
  sample1 = trainSPMA.object$sample1[(N[1]+1):(2*N[1]), ]
  sample2 = trainSPMA.object$sample2[(N[2]+1):(2*N[2]), ]
  sample3 = trainSPMA.object$sample3[(N[3]+1):(2*N[3]), ]

  #sample1 = trainSPMA.object$sample1[(1):(N[1]), ]
  #sample2 = trainSPMA.object$sample2[(1):(N[2]), ]
  #sample3 = trainSPMA.object$sample3[(1):(N[3]), ]
  p = trainSPMA.object$p
  exTx1 = Data$exTx1
  exTx2 = Data$exTx2
  exTx3 = Data$exTx3



  if(!intervalCensor){
    # unpel

    # rmse.unpel1 = sum(( log(2)/(exp(matrix(as.matrix(sample1[,1:p]),ncol = p)%*%trainSPMA.object$unpel[1,])) - sample1[,1+p])^2)/N[1]
    # rmse.unpel2 = sum(( log(2)/(exp(matrix(as.matrix(sample2[,1:p]),ncol = p)%*%trainSPMA.object$unpel[2,])) - sample2[,1+p])^2)/N[2]
    # rmse.unpel3 = sum(( log(2)/(exp(matrix(as.matrix(sample3[,1:p]),ncol = p)%*%trainSPMA.object$unpel[3,])) - sample3[,1+p])^2)/N[3]
    # rmse.unpel = c(rmse.unpel1, rmse.unpel2, rmse.unpel3)
    # rmse.unpel <- c(rmse.unpel, sum(rmse.unpel) )

    dist2beta.unpel = sum((beta0 - trainSPMA.object$unpel)^2)

    # test for the alternative rmse
    # rmse.unpel = c(0,0,0,sum((matrix(as.matrix(sample1[,1:p]),ncol = p)%*%trainSPMA.object$unpel[1,] - matrix(as.matrix(sample1[,1:p]),ncol = p)%*%trainSPMA.object$beta_0[1,])^2/N[1] +
    #                   (matrix(as.matrix(sample2[,1:p]),ncol = p)%*%trainSPMA.object$unpel[2,] - matrix(as.matrix(sample2[,1:p]),ncol = p)%*%trainSPMA.object$beta_0[2,])^2/N[2] +
    #                    (matrix(as.matrix(sample3[,1:p]),ncol = p)%*%trainSPMA.object$unpel[3,] - matrix(as.matrix(sample3[,1:p]),ncol = p)%*%trainSPMA.object$beta_0[3,])^2/N[3]
    # ))


    # use exTx1 exTx2 exTx3 instead
    # if(trainSPMA.object$correlate){
    #   sigma  = matrix(ncol = p, nrow = p)
    #   for(i in 1:p){
    #     for(j in 1:p){
    #       sigma[i,j] = 0.5^(abs(i-j)) # debug for the correlation
    #       # rho = 0.2 works well, try rho = 0.3
    #     }
    #   }
    # }
    # else{
    #   q = 0.2
    #   p1 = floor(p*q) # q = 0.2
    #   sigma = matrix(ncol =p, nrow = p )
    #   sigma[1:p1,1:p1] <- 0.25
    #   sigma[(p1+1):p,(p1+1):p] = 0.125
    # }

    # rmse.unpel = c(0,0,0,sum(matrix(rbind(trainSPMA.object$unpel[1,] -trainSPMA.object$beta_0[1,],
    #                                       trainSPMA.object$unpel[2,] -trainSPMA.object$beta_0[2,],
    #                                       trainSPMA.object$unpel[3,] -trainSPMA.object$beta_0[3,]),ncol = p)%*%sigma%*%
    #                            t(matrix(rbind(trainSPMA.object$unpel[1,] -trainSPMA.object$beta_0[1,],
    #                                           trainSPMA.object$unpel[2,] -trainSPMA.object$beta_0[2,],
    #                                           trainSPMA.object$unpel[3,] -trainSPMA.object$beta_0[3,]),ncol = p)) ))

    rmse.unpel = c(0,0,0,sum(matrix(trainSPMA.object$unpel[1,] -trainSPMA.object$beta_0[1,], ncol =p)%*%exTx1%*%t(matrix(trainSPMA.object$unpel[1,] -trainSPMA.object$beta_0[1,],ncol =p))) +
                     sum(matrix(trainSPMA.object$unpel[2,] -trainSPMA.object$beta_0[2,], ncol =p)%*%exTx2%*%t(matrix(trainSPMA.object$unpel[2,] -trainSPMA.object$beta_0[2,], ncol = p))) +
                     sum(matrix(trainSPMA.object$unpel[3,] -trainSPMA.object$beta_0[3,], ncol =p)%*%exTx3%*%t(matrix(trainSPMA.object$unpel[3,] -trainSPMA.object$beta_0[3,], ncol = p)))
    )

    #print(rmse.unpel)
    confuse.unpel = table(binbeta0, ifelse(trainSPMA.object$unpel==0,0,1) )

    # LASSO

    # we didnt adopt this rmse
    # rmse.lasso1 = sum(( log(2)/(exp(matrix(as.matrix(sample1[,1:p]),ncol = p)%*%matrix(trainSPMA.object$lasso1$beta[,best.idx$idx.lasso1], nrow = p))) -sample1[,1+p]  )^2)/N[1]
    # rmse.lasso2 = sum(( log(2)/(exp(matrix(as.matrix(sample2[,1:p]),ncol = p)%*%matrix(trainSPMA.object$lasso2$beta[,best.idx$idx.lasso2], nrow = p))) - sample2[,1+p] )^2)/N[2]
    # rmse.lasso3 = sum(( log(2)/(exp(matrix(as.matrix(sample3[,1:p]),ncol = p)%*%matrix(trainSPMA.object$lasso3$beta[,best.idx$idx.lasso3], nrow = p))) - sample3[,1+p] )^2)/N[3]
    # rmse.lasso = c(rmse.lasso1, rmse.lasso2,rmse.lasso3)
    # rmse.lasso <- c(rmse.lasso, sum(rmse.lasso) )

    dist2beta.lasso = dist2beta0(beta0, c(trainSPMA.object$lasso1$beta[,best.idx$idx.lasso1],
                                          trainSPMA.object$lasso2$beta[,best.idx$idx.lasso2],
                                          trainSPMA.object$lasso3$beta[,best.idx$idx.lasso3]), p, M = 3)


    # test new idea
    # rmse.lasso <- c(0,0,0,sum((matrix(as.matrix(sample1[,1:p]),ncol = p)%*%matrix(trainSPMA.object$lasso1$beta[,best.idx$idx.lasso1], nrow = p)
    #                    - matrix(as.matrix(sample1[,1:p]),ncol = p)%*%matrix(trainSPMA.object$beta_0[1,]))^2/N[1]
    #                   + (matrix(as.matrix(sample2[,1:p]),ncol = p)%*%matrix(trainSPMA.object$lasso2$beta[,best.idx$idx.lasso2], nrow = p)
    #                      - matrix(as.matrix(sample2[,1:p]),ncol = p)%*%matrix(trainSPMA.object$beta_0[2,]))^2/N[2]
    #                   +(matrix(as.matrix(sample3[,1:p]),ncol = p)%*%matrix(trainSPMA.object$lasso3$beta[,best.idx$idx.lasso3], nrow = p)
    #                     - matrix(as.matrix(sample3[,1:p]),ncol = p)%*%matrix(trainSPMA.object$beta_0[3,]))^2/N[3] ))
    #

    rmse.lasso = c(0,0,0,sum(matrix(trainSPMA.object$lasso1$beta[,best.idx$idx.lasso1] -trainSPMA.object$beta_0[1,],ncol = p) %*%exTx1 %*% t(matrix(trainSPMA.object$lasso1$beta[,best.idx$idx.lasso1] -trainSPMA.object$beta_0[1,],ncol = p)))+
                     sum(matrix(trainSPMA.object$lasso2$beta[,best.idx$idx.lasso2] -trainSPMA.object$beta_0[2,],ncol = p) %*%exTx2 %*% t(matrix(trainSPMA.object$lasso2$beta[,best.idx$idx.lasso2] -trainSPMA.object$beta_0[2,],ncol = p)))+
                     sum(matrix(trainSPMA.object$lasso3$beta[,best.idx$idx.lasso3] -trainSPMA.object$beta_0[3,],ncol = p) %*%exTx3 %*% t(matrix(trainSPMA.object$lasso3$beta[,best.idx$idx.lasso3] -trainSPMA.object$beta_0[3,],ncol = p))))





    lambda.lasso =  c(trainSPMA.object$lasso1$lambda[best.idx$idx.lasso1], trainSPMA.object$lasso2$lambda[best.idx$idx.lasso2],
                      trainSPMA.object$lasso3$lambda[best.idx$idx.lasso3])
    theta = rbind(trainSPMA.object$lasso1$beta[,best.idx$idx.lasso1],trainSPMA.object$lasso2$beta[,best.idx$idx.lasso2],trainSPMA.object$lasso3$beta[,best.idx$idx.lasso3] )
    confuse.lasso = table(binbeta0, ifelse(theta==0,0,1))

    # SCAD
    # rmse.scad1 = sum(( log(2)/(exp(matrix(as.matrix(sample1[,1:p]),ncol = p)%*%matrix(trainSPMA.object$scad1$beta[,best.idx$idx.scad1], nrow = p))) -sample1[,1+p]  )^2)/N[1]
    # rmse.scad2 = sum(( log(2)/(exp(matrix(as.matrix(sample2[,1:p]),ncol = p)%*%matrix(trainSPMA.object$scad2$beta[,best.idx$idx.scad2], nrow = p))) - sample2[,1+p] )^2)/N[2]
    # rmse.scad3 = sum(( log(2)/(exp(matrix(as.matrix(sample3[,1:p]),ncol = p)%*%matrix(trainSPMA.object$scad3$beta[,best.idx$idx.scad3], nrow = p))) - sample3[,1+p] )^2)/N[3]
    # rmse.scad = c(rmse.scad1, rmse.scad2,rmse.scad3)
    # rmse.scad <- c(rmse.scad, sum(rmse.scad) )

    dist2beta.scad = dist2beta0(beta0, c(trainSPMA.object$scad1$beta[,best.idx$idx.scad1],
                                         trainSPMA.object$scad2$beta[,best.idx$idx.scad2],
                                         trainSPMA.object$scad3$beta[,best.idx$idx.scad3]), p, M = 3)
    lambda.scad =  c(trainSPMA.object$scad1$lambda[best.idx$idx.scad1], trainSPMA.object$scad2$lambda[best.idx$idx.scad2],
                     trainSPMA.object$scad3$lambda[best.idx$idx.scad3])
    theta = rbind(trainSPMA.object$scad1$beta[,best.idx$idx.scad1],trainSPMA.object$scad2$beta[,best.idx$idx.scad2],trainSPMA.object$scad3$beta[,best.idx$idx.scad3] )
    confuse.scad = table(binbeta0, ifelse(theta==0,0,1))

    # test new idea
    # rmse.scad <- c(0,0,0,sum((matrix(as.matrix(sample1[,1:p]),ncol = p)%*%matrix(trainSPMA.object$scad1$beta[,best.idx$idx.scad1], nrow = p)
    #                    - matrix(as.matrix(sample1[,1:p]),ncol = p)%*%matrix(trainSPMA.object$beta_0[1,]))^2/N[1]
    #                   + (matrix(as.matrix(sample2[,1:p]),ncol = p)%*%matrix(trainSPMA.object$scad2$beta[,best.idx$idx.scad2], nrow = p)
    #                      - matrix(as.matrix(sample2[,1:p]),ncol = p)%*%matrix(trainSPMA.object$beta_0[2,]))^2/N[2]
    #                   +(matrix(as.matrix(sample3[,1:p]),ncol = p)%*%matrix(trainSPMA.object$scad3$beta[,best.idx$idx.scad3], nrow = p)
    #                     - matrix(as.matrix(sample3[,1:p]),ncol = p)%*%matrix(trainSPMA.object$beta_0[3,]))^2/N[3] ))
    #
    # rmse.scad = c(0,0,0,sum(matrix(rbind(trainSPMA.object$scad1$beta[,best.idx$idx.scad1] -trainSPMA.object$beta_0[1,],
    #                                      trainSPMA.object$scad2$beta[,best.idx$idx.scad2]-trainSPMA.object$beta_0[2,],
    #                                      trainSPMA.object$scad3$beta[,best.idx$idx.scad3] -trainSPMA.object$beta_0[3,]),ncol = p)%*%sigma%*%
    #                           t(matrix(rbind(trainSPMA.object$scad1$beta[,best.idx$idx.scad1] -trainSPMA.object$beta_0[1,],
    #                                          trainSPMA.object$scad2$beta[,best.idx$idx.scad2]-trainSPMA.object$beta_0[2,],
    #                                          trainSPMA.object$scad3$beta[,best.idx$idx.scad3] -trainSPMA.object$beta_0[3,]),ncol = p)) ))
    #


    rmse.scad = c(0,0,0,sum(matrix(trainSPMA.object$scad1$beta[,best.idx$idx.scad1] -trainSPMA.object$beta_0[1,],ncol = p) %*%exTx1 %*% t(matrix(trainSPMA.object$scad1$beta[,best.idx$idx.scad1] -trainSPMA.object$beta_0[1,],ncol = p)))+
                    sum(matrix(trainSPMA.object$scad2$beta[,best.idx$idx.scad2] -trainSPMA.object$beta_0[2,],ncol = p) %*%exTx2 %*% t(matrix(trainSPMA.object$scad2$beta[,best.idx$idx.scad2] -trainSPMA.object$beta_0[2,],ncol = p)))+
                    sum(matrix(trainSPMA.object$scad3$beta[,best.idx$idx.scad3] -trainSPMA.object$beta_0[3,],ncol = p) %*%exTx3 %*% t(matrix(trainSPMA.object$scad3$beta[,best.idx$idx.scad3] -trainSPMA.object$beta_0[3,],ncol = p))))




    # MCP
    # rmse.mcp1 = sum(( log(2)/(exp(matrix(as.matrix(sample1[,1:p]),ncol = p)%*%matrix(trainSPMA.object$mcp1$beta[,best.idx$idx.mcp1], nrow = p))) -sample1[,1+p]  )^2)/N[1]
    # rmse.mcp2 = sum(( log(2)/(exp(matrix(as.matrix(sample2[,1:p]),ncol = p)%*%matrix(trainSPMA.object$mcp2$beta[,best.idx$idx.mcp2], nrow = p))) - sample2[,1+p] )^2)/N[2]
    # rmse.mcp3 = sum(( log(2)/(exp(matrix(as.matrix(sample3[,1:p]),ncol = p)%*%matrix(trainSPMA.object$mcp3$beta[,best.idx$idx.mcp3], nrow = p))) - sample3[,1+p] )^2)/N[3]
    # rmse.mcp = c(rmse.mcp1, rmse.mcp2,rmse.mcp3)
    # rmse.mcp <- c(rmse.mcp, sum(rmse.mcp) )

    dist2beta.mcp = dist2beta0(beta0, c(trainSPMA.object$mcp1$beta[,best.idx$idx.mcp1],
                                        trainSPMA.object$mcp2$beta[,best.idx$idx.mcp2],
                                        trainSPMA.object$mcp3$beta[,best.idx$idx.mcp3]), p, M = 3)
    lambda.mcp =  c(trainSPMA.object$mcp1$lambda[best.idx$idx.mcp1], trainSPMA.object$mcp2$lambda[best.idx$idx.mcp2],
                    trainSPMA.object$mcp3$lambda[best.idx$idx.mcp3])
    theta = rbind(trainSPMA.object$mcp1$beta[,best.idx$idx.mcp1],trainSPMA.object$mcp2$beta[,best.idx$idx.mcp2],trainSPMA.object$mcp3$beta[,best.idx$idx.mcp3] )
    confuse.mcp = table(binbeta0, ifelse(theta==0,0,1))

    # test new idea
    # rmse.mcp <- c(0,0,0,sum((matrix(as.matrix(sample1[,1:p]),ncol = p)%*%matrix(trainSPMA.object$mcp1$beta[,best.idx$idx.mcp1], nrow = p)
    #                    - matrix(as.matrix(sample1[,1:p]),ncol = p)%*%matrix(trainSPMA.object$beta_0[1,]))^2/N[1]
    #                   + (matrix(as.matrix(sample2[,1:p]),ncol = p)%*%matrix(trainSPMA.object$mcp2$beta[,best.idx$idx.mcp2], nrow = p)
    #                      - matrix(as.matrix(sample2[,1:p]),ncol = p)%*%matrix(trainSPMA.object$beta_0[2,]))^2/N[2]
    #                   +(matrix(as.matrix(sample3[,1:p]),ncol = p)%*%matrix(trainSPMA.object$mcp3$beta[,best.idx$idx.mcp3], nrow = p)
    #                     - matrix(as.matrix(sample3[,1:p]),ncol = p)%*%matrix(trainSPMA.object$beta_0[3,]))^2 /N[3]))
    #
    # rmse.mcp = c(0,0,0,sum(matrix(rbind(trainSPMA.object$mcp1$beta[,best.idx$idx.mcp1] -trainSPMA.object$beta_0[1,],
    #                                     trainSPMA.object$mcp2$beta[,best.idx$idx.mcp2]-trainSPMA.object$beta_0[2,],
    #                                     trainSPMA.object$mcp3$beta[,best.idx$idx.mcp3] -trainSPMA.object$beta_0[3,]),ncol = p)%*%sigma%*%
    #                          t(matrix(rbind(trainSPMA.object$mcp1$beta[,best.idx$idx.mcp1] -trainSPMA.object$beta_0[1,],
    #                                         trainSPMA.object$mcp2$beta[,best.idx$idx.mcp2]-trainSPMA.object$beta_0[2,],
    #                                         trainSPMA.object$mcp3$beta[,best.idx$idx.mcp3] -trainSPMA.object$beta_0[3,]),ncol = p)) ))
    #

    rmse.mcp = c(0,0,0,sum(matrix(trainSPMA.object$mcp1$beta[,best.idx$idx.mcp1] -trainSPMA.object$beta_0[1,],ncol = p) %*%exTx1 %*% t(matrix(trainSPMA.object$mcp1$beta[,best.idx$idx.mcp1] -trainSPMA.object$beta_0[1,],ncol = p)))+
                   sum(matrix(trainSPMA.object$mcp2$beta[,best.idx$idx.mcp2] -trainSPMA.object$beta_0[2,],ncol = p) %*%exTx2 %*% t(matrix(trainSPMA.object$mcp2$beta[,best.idx$idx.mcp2] -trainSPMA.object$beta_0[2,],ncol = p)))+
                   sum(matrix(trainSPMA.object$mcp3$beta[,best.idx$idx.mcp3] -trainSPMA.object$beta_0[3,],ncol = p) %*%exTx3 %*% t(matrix(trainSPMA.object$mcp3$beta[,best.idx$idx.mcp3] -trainSPMA.object$beta_0[3,],ncol = p))))




    # GROUP LASSO
    theta = trainSPMA.object$grlasso$beta[,best.idx$idx.grlasso]

    # rmse.grlasso1 = sum((log(2)/(exp(matrix(as.matrix(sample1[,1:p]),ncol = p)%*%matrix(theta[1:p], nrow = p))) -sample1[,1+p]   )^2)/N[1]
    # rmse.grlasso2 = sum((log(2)/(exp( matrix(as.matrix(sample2[,1:p]),ncol = p)%*%matrix(theta[(p+1):(2*p)], nrow = p))) - sample2[,1+p]  )^2)/N[2]
    # rmse.grlasso3 = sum(( log(2)/(exp(matrix(as.matrix(sample3[,1:p]),ncol = p)%*%matrix(theta[(2*p+1):(3*p)], nrow = p))) - sample3[,1+p] )^2)/N[3]
    # rmse.grlasso = c(rmse.grlasso1, rmse.grlasso2,rmse.grlasso3)
    # rmse.grlasso <- c(rmse.grlasso, sum(rmse.grlasso))

    dist2beta.grlasso = dist2beta0(beta0,theta,p, M=3)
    lambda.grlasso <-  trainSPMA.object$grlasso$lambda[best.idx$idx.grlasso]
    ok = matrix(theta, ncol = p, nrow = 3, byrow = TRUE)
    confuse.grlasso = table(binbeta0, ifelse(ok==0,0,1))

    # test new idea
    # rmse.grlasso = c(0,0,0,sum(
    #   (matrix(as.matrix(sample1[,1:p]),ncol = p)%*%matrix(theta[1:p], nrow = p) - matrix(as.matrix(sample1[,1:p]),ncol = p)%*%matrix(trainSPMA.object$beta_0[1,], nrow = p) )^2/N[1]+
    #   (matrix(as.matrix(sample2[,1:p]),ncol = p)%*%matrix(theta[(p+1):(2*p)], nrow = p) - matrix(as.matrix(sample2[,1:p]),ncol = p)%*%matrix(trainSPMA.object$beta_0[2,], nrow = p) )^2/N[2]+
    #   (matrix(as.matrix(sample3[,1:p]),ncol = p)%*%matrix(theta[(2*p+1):(3*p)], nrow = p) - matrix(as.matrix(sample3[,1:p]),ncol = p)%*%matrix(trainSPMA.object$beta_0[3,], nrow = p) )^2/N[3]
    #
    # ))



    rmse.grlasso = c(0,0,0,sum(matrix(theta[1:p] -trainSPMA.object$beta_0[1,], ncol = p) %*%exTx1 %*% t(matrix(theta[1:p] -trainSPMA.object$beta_0[1,], ncol = p)))+
                       sum(matrix(theta[(p+1):(2*p)]-trainSPMA.object$beta_0[2,], ncol = p) %*%exTx2 %*% t(matrix( theta[(p+1):(2*p)]-trainSPMA.object$beta_0[2,], ncol = p)))+
                       sum(matrix(theta[(2*p+1):(3*p)] -trainSPMA.object$beta_0[3,],ncol = p) %*% exTx3 %*% t(matrix(theta[(2*p+1):(3*p)] -trainSPMA.object$beta_0[3,],ncol = p))))





    # GROUP SCAD
    theta = trainSPMA.object$grscad$beta[,best.idx$idx.grscad]

    # rmse.grscad1 = sum((log(2)/(exp(matrix(as.matrix(sample1[,1:p]),ncol = p)%*%matrix(theta[1:p], nrow = p))) -sample1[,1+p])^2)/N[1]
    # rmse.grscad2 = sum((log(2)/(exp( matrix(as.matrix(sample2[,1:p]),ncol = p)%*%matrix(theta[(p+1):(2*p)], nrow = p))) - sample2[,1+p])^2)/N[2]
    # rmse.grscad3 = sum(( log(2)/(exp(matrix(as.matrix(sample3[,1:p]),ncol = p)%*%matrix(theta[(2*p+1):(3*p)], nrow = p))) - sample3[,1+p])^2)/N[3]
    # rmse.grscad = c(rmse.grscad1, rmse.grscad2,rmse.grscad3)
    # rmse.grscad <- c(rmse.grscad, sum(rmse.grscad) )

    dist2beta.grscad = dist2beta0(beta0,theta,p,M = 3)
    lambda.grscad <-  trainSPMA.object$grscad$lambda[best.idx$idx.grscad]
    ok = matrix(theta, ncol = p, nrow = 3, byrow = TRUE)
    confuse.grscad = table(binbeta0, ifelse(ok==0,0,1))

    # test new idea
    # rmse.grscad = c(0,0,0,sum(
    #   (matrix(as.matrix(sample1[,1:p]),ncol = p)%*%matrix(theta[1:p], nrow = p) - matrix(as.matrix(sample1[,1:p]),ncol = p)%*%matrix(trainSPMA.object$beta_0[1,], nrow = p) )^2/N[1]+
    #     (matrix(as.matrix(sample2[,1:p]),ncol = p)%*%matrix(theta[(p+1):(2*p)], nrow = p) - matrix(as.matrix(sample2[,1:p]),ncol = p)%*%matrix(trainSPMA.object$beta_0[2,], nrow = p) )^2/N[2]+
    #     (matrix(as.matrix(sample3[,1:p]),ncol = p)%*%matrix(theta[(2*p+1):(3*p)], nrow = p) - matrix(as.matrix(sample3[,1:p]),ncol = p)%*%matrix(trainSPMA.object$beta_0[3,], nrow = p) )^2/N[3]
    #
    # ))

    # rmse.grscad = c(0,0,0,sum(matrix(rbind(theta[1:p] -trainSPMA.object$beta_0[1,],
    #                                        theta[(p+1):(2*p)]-trainSPMA.object$beta_0[2,],
    #                                        theta[(2*p+1):(3*p)] -trainSPMA.object$beta_0[3,]),ncol = p)%*%sigma%*%
    #                             t(matrix(rbind(theta[1:p] -trainSPMA.object$beta_0[1,],
    #                                            theta[(p+1):(2*p)]-trainSPMA.object$beta_0[2,],
    #                                            theta[(2*p+1):(3*p)] -trainSPMA.object$beta_0[3,]),ncol = p)) ))

    rmse.grscad = c(0,0,0,sum(matrix(theta[1:p] -trainSPMA.object$beta_0[1,], ncol = p) %*%exTx1 %*% t(matrix(theta[1:p] -trainSPMA.object$beta_0[1,], ncol = p)))+
                      sum(matrix(theta[(p+1):(2*p)]-trainSPMA.object$beta_0[2,], ncol = p) %*%exTx2 %*% t(matrix( theta[(p+1):(2*p)]-trainSPMA.object$beta_0[2,], ncol = p)))+
                      sum(matrix(theta[(2*p+1):(3*p)] -trainSPMA.object$beta_0[3,],ncol = p) %*% exTx3 %*% t(matrix(theta[(2*p+1):(3*p)] -trainSPMA.object$beta_0[3,],ncol = p))))



    # GROUP MCP
    theta = trainSPMA.object$grmcp$beta[,best.idx$idx.grmcp]

    # rmse.grmcp1 = sum((log(2)/(exp(matrix(as.matrix(sample1[,1:p]),ncol = p)%*%matrix(theta[1:p], nrow = p))) -sample1[,1+p])^2)/N[1]
    # rmse.grmcp2 = sum((log(2)/(exp( matrix(as.matrix(sample2[,1:p]),ncol = p)%*%matrix(theta[(p+1):(2*p)], nrow = p))) - sample2[,1+p])^2)/N[2]
    # rmse.grmcp3 = sum(( log(2)/(exp(matrix(as.matrix(sample3[,1:p]),ncol = p)%*%matrix(theta[(2*p+1):(3*p)], nrow = p))) - sample3[,1+p])^2)/N[3]
    # rmse.grmcp = c(rmse.grmcp1, rmse.grmcp2,rmse.grmcp3)
    # rmse.grmcp <- c(rmse.grmcp, sum(rmse.grmcp) )

    dist2beta.grmcp = dist2beta0(beta0,theta,p,M =3)
    lambda.grmcp <-  trainSPMA.object$grmcp$lambda[best.idx$idx.grmcp]
    ok = matrix(theta, ncol = p, nrow = 3, byrow = TRUE)
    confuse.grmcp = table(binbeta0, ifelse(ok==0,0,1))

    # test new idea
    # rmse.grmcp = c(0,0,0,sum(
    #   (matrix(as.matrix(sample1[,1:p]),ncol = p)%*%matrix(theta[1:p], nrow = p) - matrix(as.matrix(sample1[,1:p]),ncol = p)%*%matrix(trainSPMA.object$beta_0[1,], nrow = p) )^2/N[1]+
    #     (matrix(as.matrix(sample2[,1:p]),ncol = p)%*%matrix(theta[(p+1):(2*p)], nrow = p) - matrix(as.matrix(sample2[,1:p]),ncol = p)%*%matrix(trainSPMA.object$beta_0[2,], nrow = p) )^2/N[2]+
    #     (matrix(as.matrix(sample3[,1:p]),ncol = p)%*%matrix(theta[(2*p+1):(3*p)], nrow = p) - matrix(as.matrix(sample3[,1:p]),ncol = p)%*%matrix(trainSPMA.object$beta_0[3,], nrow = p) )^2/N[3]
    #
    # ))

    rmse.grmcp = c(0,0,0,sum(matrix(theta[1:p] -trainSPMA.object$beta_0[1,], ncol = p) %*%exTx1 %*% t(matrix(theta[1:p] -trainSPMA.object$beta_0[1,], ncol = p)))+
                     sum(matrix(theta[(p+1):(2*p)]-trainSPMA.object$beta_0[2,], ncol = p) %*%exTx2 %*% t(matrix( theta[(p+1):(2*p)]-trainSPMA.object$beta_0[2,], ncol = p)))+
                     sum(matrix(theta[(2*p+1):(3*p)] -trainSPMA.object$beta_0[3,],ncol = p) %*% exTx3 %*% t(matrix(theta[(2*p+1):(3*p)] -trainSPMA.object$beta_0[3,],ncol = p))))

    # SPMA
    # theta = trainSPMA.object$spma[[spma.bic$idx.spma]]$thetaM[,,trainSPMA.object$spma[[spma.bic$idx.spma]]$bbk]
    theta = trainSPMA.object$spma[[spma.bic$idx.spma]]$alphaM[,,trainSPMA.object$spma[[spma.bic$idx.spma]]$bbk]
    # rmse.spma1 = sum((log(2)/(exp(matrix(as.matrix(sample1[,1:p]),ncol = p)%*%matrix(theta[1,1:p], nrow = p))) - sample1[,1+p])^2)/N[1]
    # rmse.spma2 = sum((log(2)/(exp( matrix(as.matrix(sample2[,1:p]),ncol = p)%*%matrix(theta[2,1:p], nrow = p))) - sample2[,1+p])^2)/N[2]
    # rmse.spma3 = sum((log(2)/(exp(matrix(as.matrix(sample3[,1:p]),ncol = p)%*%matrix(theta[3,1:p], nrow = p))) - sample3[,1+p])^2)/N[3]
    # rmse.spma = c(rmse.spma1, rmse.spma2,rmse.spma3)
    # rmse.spma = c(rmse.spma, sum(rmse.spma))

    dist2beta.spma = sum((beta0-theta)^2)
    lambda.spma = trainSPMA.object$spma.lambda[spma.bic$idx.spma]
    confuse.spma = table(binbeta0,ifelse(theta==0,0,1))

    # test new idea
    # rmse.spma = c(0,0,0,sum(
    #   (matrix(as.matrix(sample1[,1:p]),ncol = p)%*%matrix(theta[1,1:p], nrow = p) - matrix(as.matrix(sample1[,1:p]),ncol = p)%*%matrix(trainSPMA.object$beta_0[1,], nrow = p) )^2/N[1]+
    #     (matrix(as.matrix(sample2[,1:p]),ncol = p)%*%matrix(theta[2,1:p], nrow = p) - matrix(as.matrix(sample2[,1:p]),ncol = p)%*%matrix(trainSPMA.object$beta_0[2,], nrow = p) )^2/N[2]+
    #     (matrix(as.matrix(sample3[,1:p]),ncol = p)%*%matrix(theta[3,1:p], nrow = p) - matrix(as.matrix(sample3[,1:p]),ncol = p)%*%matrix(trainSPMA.object$beta_0[3,], nrow = p) )^2/N[3]
    #
    # ))

    rmse.spma = c(0,0,0,sum(matrix(theta[1,1:p] -trainSPMA.object$beta_0[1,],ncol = p)%*% exTx1 %*% t(matrix(theta[1,1:p] -trainSPMA.object$beta_0[1,],ncol = p))) +
                    sum(matrix(theta[2,1:p]-trainSPMA.object$beta_0[2,], ncol = p)%*% exTx2 %*% t(matrix(theta[2,1:p]-trainSPMA.object$beta_0[2,], ncol = p))) +
                    sum(matrix(theta[3,1:p] -trainSPMA.object$beta_0[3,],ncol = p) %*% exTx3 %*% t(matrix(theta[3,1:p] -trainSPMA.object$beta_0[3,],ncol = p))))



    return(list(rmse.unpel = rmse.unpel,rmse.lasso = rmse.lasso,rmse.scad = rmse.scad, rmse.mcp = rmse.mcp, rmse.grlasso = rmse.grlasso,
                rmse.grscad = rmse.grscad, rmse.grmcp = rmse.grmcp, rmse.spma = rmse.spma,
                dist2beta.unpel = dist2beta.unpel, dist2beta.lasso = dist2beta.lasso,dist2beta.scad = dist2beta.scad,
                dist2beta.mcp = dist2beta.mcp,
                dist2beta.grlasso = dist2beta.grlasso,dist2beta.grscad = dist2beta.grscad,dist2beta.grmcp = dist2beta.grmcp,
                dist2beta.spma = dist2beta.spma,
                lambda.lasso = lambda.lasso, lambda.scad = lambda.scad,
                lambda.mcp = lambda.mcp, lambda.grlasso = lambda.grlasso,  lambda.grscad = lambda.grscad, lambda.grmcp = lambda.grmcp,
                lambda.spma = lambda.spma,
                confuse.unpel =   confuse.unpel, confuse.lasso =   confuse.lasso, confuse.scad = confuse.scad,confuse.mcp = confuse.mcp,
                confuse.grlasso =   confuse.grlasso,
                confuse.grscad =   confuse.grscad,  confuse.grmcp =   confuse.grmcp,
                confuse.spma = confuse.spma))


  }
  else{

    # unpel
    dist2beta.unpel = sum((beta0 - trainSPMA.object$unpel)^2)
    rmse.unpel = c(0,0,0,sum(matrix(trainSPMA.object$unpel[1,] -trainSPMA.object$beta_0[1,], ncol =p)%*%exTx1%*%t(matrix(trainSPMA.object$unpel[1,] -trainSPMA.object$beta_0[1,],ncol =p))) +
                     sum(matrix(trainSPMA.object$unpel[2,] -trainSPMA.object$beta_0[2,], ncol =p)%*%exTx2%*%t(matrix(trainSPMA.object$unpel[2,] -trainSPMA.object$beta_0[2,], ncol = p))) +
                     sum(matrix(trainSPMA.object$unpel[3,] -trainSPMA.object$beta_0[3,], ncol =p)%*%exTx3%*%t(matrix(trainSPMA.object$unpel[3,] -trainSPMA.object$beta_0[3,], ncol = p))))
    confuse.unpel = table(binbeta0, ifelse(trainSPMA.object$unpel==0,0,1) )
    # lasso
    dist2beta.lasso = dist2beta0(beta0, c(trainSPMA.object$lasso1[[best.idx$idx.lasso1]]$beta,
                                          trainSPMA.object$lasso2[[best.idx$idx.lasso2]]$beta,
                                          trainSPMA.object$lasso3[[best.idx$idx.lasso3]]$beta), p, M = 3)
    theta = rbind(trainSPMA.object$lasso1[[best.idx$idx.lasso1]]$beta,trainSPMA.object$lasso2[[best.idx$idx.lasso2]]$beta,trainSPMA.object$lasso3[[best.idx$idx.lasso3]]$beta)
    rmse.lasso = c(0,0,0,sum(matrix(trainSPMA.object$lasso1[[best.idx$idx.lasso1]]$beta-trainSPMA.object$beta_0[1,],ncol = p) %*%exTx1 %*% t(matrix(trainSPMA.object$lasso1[[best.idx$idx.lasso1]]$beta -trainSPMA.object$beta_0[1,],ncol = p)))+
                     sum(matrix(trainSPMA.object$lasso2[[best.idx$idx.lasso2]]$beta -trainSPMA.object$beta_0[2,],ncol = p) %*%exTx2 %*% t(matrix(trainSPMA.object$lasso2[[best.idx$idx.lasso2]]$beta -trainSPMA.object$beta_0[2,],ncol = p)))+
                     sum(matrix(trainSPMA.object$lasso3[[best.idx$idx.lasso3]]$beta -trainSPMA.object$beta_0[3,],ncol = p) %*%exTx3 %*% t(matrix(trainSPMA.object$lasso3[[best.idx$idx.lasso3]]$beta -trainSPMA.object$beta_0[3,],ncol = p))))

    confuse.lasso = table(binbeta0, ifelse(theta==0,0,1))

    # SCAD
    dist2beta.scad = dist2beta0(beta0, c(trainSPMA.object$scad1[[best.idx$idx.scad1]]$beta,
                                         trainSPMA.object$scad2[[best.idx$idx.scad2]]$beta,
                                         trainSPMA.object$scad3[[best.idx$idx.scad3]]$beta), p, M = 3)
    theta = rbind(trainSPMA.object$scad1[[best.idx$idx.scad1]]$beta,trainSPMA.object$scad2[[best.idx$idx.scad2]]$beta,trainSPMA.object$scad3[[best.idx$idx.scad3]]$beta)
    rmse.scad = c(0,0,0,sum(matrix(trainSPMA.object$scad1[[best.idx$idx.scad1]]$beta -trainSPMA.object$beta_0[1,],ncol = p) %*%exTx1 %*% t(matrix(trainSPMA.object$scad1[[best.idx$idx.scad1]]$beta -trainSPMA.object$beta_0[1,],ncol = p)))+
                    sum(matrix(trainSPMA.object$scad2[[best.idx$idx.scad2]]$beta -trainSPMA.object$beta_0[2,],ncol = p) %*%exTx2 %*% t(matrix(trainSPMA.object$scad2[[best.idx$idx.scad2]]$beta -trainSPMA.object$beta_0[2,],ncol = p)))+
                    sum(matrix(trainSPMA.object$scad3[[best.idx$idx.scad3]]$beta -trainSPMA.object$beta_0[3,],ncol = p) %*%exTx3 %*% t(matrix(trainSPMA.object$scad3[[best.idx$idx.scad3]]$beta -trainSPMA.object$beta_0[3,],ncol = p))))
    confuse.scad = table(binbeta0, ifelse(theta==0,0,1))

    # Group Lasso
    # beta from group lasso, group scad, group mcp is different
    # Group lasso
    thetaH = trainSPMA.object$grlasso[[best.idx$idx.grlasso]]$beta_hat
    theta = matrix(thetaH, nrow = 3, byrow = FALSE)
    dist2beta.grlasso = sum((beta0-theta)^2)
    lambda.grlasso <-  trainSPMA.object$lambda.group[best.idx$idx.grlasso]
    confuse.grlasso = table(binbeta0,ifelse(theta==0,0,1))
    rmse.grlasso = c(0,0,0,sum(matrix(theta[1,1:p] -trainSPMA.object$beta_0[1,],ncol = p)%*% exTx1 %*% t(matrix(theta[1,1:p] -trainSPMA.object$beta_0[1,],ncol = p))) +
                       sum(matrix(theta[2,1:p]-trainSPMA.object$beta_0[2,], ncol = p)%*% exTx2 %*% t(matrix(theta[2,1:p]-trainSPMA.object$beta_0[2,], ncol = p))) +
                       sum(matrix(theta[3,1:p] -trainSPMA.object$beta_0[3,],ncol = p) %*% exTx3 %*% t(matrix(theta[3,1:p] -trainSPMA.object$beta_0[3,],ncol = p))))


    # Group SCAD
    thetaH = trainSPMA.object$grscad[[best.idx$idx.grscad]]$beta_hat
    theta = matrix(thetaH, nrow = 3, byrow = FALSE)
    dist2beta.grscad = sum((beta0-theta)^2)
    lambda.grscad <-  trainSPMA.object$lambda.group[best.idx$idx.grscad]
    confuse.grscad = table(binbeta0,ifelse(theta==0,0,1))
    rmse.grscad = c(0,0,0,sum(matrix(theta[1,1:p] -trainSPMA.object$beta_0[1,],ncol = p)%*% exTx1 %*% t(matrix(theta[1,1:p] -trainSPMA.object$beta_0[1,],ncol = p))) +
                      sum(matrix(theta[2,1:p]-trainSPMA.object$beta_0[2,], ncol = p)%*% exTx2 %*% t(matrix(theta[2,1:p]-trainSPMA.object$beta_0[2,], ncol = p))) +
                      sum(matrix(theta[3,1:p] -trainSPMA.object$beta_0[3,],ncol = p) %*% exTx3 %*% t(matrix(theta[3,1:p] -trainSPMA.object$beta_0[3,],ncol = p))))


    # Group MCP
    thetaH = trainSPMA.object$grmcp[[best.idx$idx.grmcp]]$beta_hat
    theta = matrix(thetaH, nrow = 3, byrow = FALSE)
    dist2beta.grmcp = sum((beta0-theta)^2)
    lambda.grmcp <-  trainSPMA.object$lambda.group[best.idx$idx.grmcp]
    confuse.grmcp = table(binbeta0,ifelse(theta==0,0,1))
    rmse.grmcp = c(0,0,0,sum(matrix(theta[1,1:p] -trainSPMA.object$beta_0[1,],ncol = p)%*% exTx1 %*% t(matrix(theta[1,1:p] -trainSPMA.object$beta_0[1,],ncol = p))) +
                     sum(matrix(theta[2,1:p]-trainSPMA.object$beta_0[2,], ncol = p)%*% exTx2 %*% t(matrix(theta[2,1:p]-trainSPMA.object$beta_0[2,], ncol = p))) +
                     sum(matrix(theta[3,1:p] -trainSPMA.object$beta_0[3,],ncol = p) %*% exTx3 %*% t(matrix(theta[3,1:p] -trainSPMA.object$beta_0[3,],ncol = p))))


    # SPMA
    # theta = trainSPMA.object$spma[[spma.bic$idx.spma]]$thetaM[,,trainSPMA.object$spma[[spma.bic$idx.spma]]$bbk]
    theta = trainSPMA.object$spma[[spma.bic$idx.spma]]$alphaM[,,trainSPMA.object$spma[[spma.bic$idx.spma]]$bbk]
    # rmse.spma1 = sum((log(2)/(exp(matrix(as.matrix(sample1[,1:p]),ncol = p)%*%matrix(theta[1,1:p], nrow = p))) - sample1[,1+p])^2)/N[1]
    # rmse.spma2 = sum((log(2)/(exp( matrix(as.matrix(sample2[,1:p]),ncol = p)%*%matrix(theta[2,1:p], nrow = p))) - sample2[,1+p])^2)/N[2]
    # rmse.spma3 = sum((log(2)/(exp(matrix(as.matrix(sample3[,1:p]),ncol = p)%*%matrix(theta[3,1:p], nrow = p))) - sample3[,1+p])^2)/N[3]
    # rmse.spma = c(rmse.spma1, rmse.spma2,rmse.spma3)
    # rmse.spma = c(rmse.spma, sum(rmse.spma))

    dist2beta.spma = sum((beta0-theta)^2)
    lambda.spma = trainSPMA.object$spma.lambda[spma.bic$idx.spma]
    confuse.spma = table(binbeta0,ifelse(theta==0,0,1))

    # test new idea
    # rmse.spma = c(0,0,0,sum(
    #   (matrix(as.matrix(sample1[,1:p]),ncol = p)%*%matrix(theta[1,1:p], nrow = p) - matrix(as.matrix(sample1[,1:p]),ncol = p)%*%matrix(trainSPMA.object$beta_0[1,], nrow = p) )^2/N[1]+
    #     (matrix(as.matrix(sample2[,1:p]),ncol = p)%*%matrix(theta[2,1:p], nrow = p) - matrix(as.matrix(sample2[,1:p]),ncol = p)%*%matrix(trainSPMA.object$beta_0[2,], nrow = p) )^2/N[2]+
    #     (matrix(as.matrix(sample3[,1:p]),ncol = p)%*%matrix(theta[3,1:p], nrow = p) - matrix(as.matrix(sample3[,1:p]),ncol = p)%*%matrix(trainSPMA.object$beta_0[3,], nrow = p) )^2/N[3]
    #
    # ))

    rmse.spma = c(0,0,0,sum(matrix(theta[1,1:p] -trainSPMA.object$beta_0[1,],ncol = p)%*% exTx1 %*% t(matrix(theta[1,1:p] -trainSPMA.object$beta_0[1,],ncol = p))) +
                    sum(matrix(theta[2,1:p]-trainSPMA.object$beta_0[2,], ncol = p)%*% exTx2 %*% t(matrix(theta[2,1:p]-trainSPMA.object$beta_0[2,], ncol = p))) +
                    sum(matrix(theta[3,1:p] -trainSPMA.object$beta_0[3,],ncol = p) %*% exTx3 %*% t(matrix(theta[3,1:p] -trainSPMA.object$beta_0[3,],ncol = p))))


    return(list(rmse.unpel = rmse.unpel,rmse.lasso = rmse.lasso,rmse.scad = rmse.scad,
                rmse.grlasso = rmse.grlasso,
                rmse.grscad = rmse.grscad, rmse.grmcp = rmse.grmcp, rmse.spma = rmse.spma,
                dist2beta.unpel = dist2beta.unpel, dist2beta.lasso = dist2beta.lasso,dist2beta.scad = dist2beta.scad,
                # dist2beta.mcp = dist2beta.mcp,
                dist2beta.grlasso = dist2beta.grlasso,dist2beta.grscad = dist2beta.grscad,dist2beta.grmcp = dist2beta.grmcp,
                dist2beta.spma = dist2beta.spma,
                #lambda.lasso = lambda.lasso, lambda.scad = lambda.scad,
                #lambda.mcp = lambda.mcp,
                lambda.grlasso = lambda.grlasso,  lambda.grscad = lambda.grscad, lambda.grmcp = lambda.grmcp,
                lambda.spma = lambda.spma,
                confuse.unpel =   confuse.unpel, confuse.lasso =   confuse.lasso, confuse.scad = confuse.scad,
                # confuse.mcp = confuse.mcp,
                confuse.grlasso =   confuse.grlasso,
                confuse.grscad =   confuse.grscad,  confuse.grmcp =   confuse.grmcp,
                confuse.spma = confuse.spma))


  }


}


printResult<-function(rre,intervalCensor = FALSE){

  if(!intervalCensor){
    rep = length(rre)
    #rmseFinal <- rep(0,8)
    # report median
    rmseFinal <- matrix(nrow = rep,ncol = 8)
    dist2beta <- rep(0,8)
    tp<- rep(0,8)
    fp<- rep(0,8)
    sdd <- matrix(nrow = rep, ncol = 8)
    std.rmse <- rep(0,8)

    sdd.beta <- matrix(nrow = rep, ncol = 8)
    std.beta <- rep(0,8)
    for (i in 1:length(rre)){
      #print(i)
      # RMSE
      sdd[i,1] = rre[[i]]$rmse.unpel[4]
      sdd[i,2] = rre[[i]]$rmse.lasso[4]
      sdd[i,3] = rre[[i]]$rmse.grlasso[4]
      sdd[i,4] = rre[[i]]$rmse.scad[4]
      sdd[i,5] = rre[[i]]$rmse.grscad[4]
      sdd[i,6] = rre[[i]]$rmse.mcp[4]
      sdd[i,7] = rre[[i]]$rmse.grmcp[4]
      sdd[i,8] = rre[[i]]$rmse.spma[4]


      sdd.beta[i,1] = rre[[i]]$dist2beta.unpel
      sdd.beta[i,2] = rre[[i]]$dist2beta.lasso
      sdd.beta[i,3] = rre[[i]]$dist2beta.grlasso
      sdd.beta[i,4] = rre[[i]]$dist2beta.scad
      sdd.beta[i,5] = rre[[i]]$dist2beta.grscad
      sdd.beta[i,6] = rre[[i]]$dist2beta.mcp
      sdd.beta[i,7] = rre[[i]]$dist2beta.grmcp
      sdd.beta[i,8] = rre[[i]]$dist2beta.spma


      # rmseFinal[1]  = rmseFinal[1] + rre[[i]]$rmse.unpel[4]
      # rmseFinal[2]  = rmseFinal[2] + rre[[i]]$rmse.lasso[4]
      # rmseFinal[3]  = rmseFinal[3] + rre[[i]]$rmse.grlasso[4]
      # rmseFinal[4]  = rmseFinal[4] + rre[[i]]$rmse.scad[4]
      # rmseFinal[5]  = rmseFinal[5] + rre[[i]]$rmse.grscad[4]
      # rmseFinal[6]  = rmseFinal[6] + rre[[i]]$rmse.mcp[4]
      # rmseFinal[7]  = rmseFinal[7] + rre[[i]]$rmse.grmcp[4]
      # rmseFinal[8]  = rmseFinal[8] + rre[[i]]$rmse.spma[4]

      rmseFinal[i,1]  = rre[[i]]$rmse.unpel[4]
      rmseFinal[i,2]  = rre[[i]]$rmse.lasso[4]
      rmseFinal[i,3]  = rre[[i]]$rmse.grlasso[4]
      rmseFinal[i,4]  = rre[[i]]$rmse.scad[4]
      rmseFinal[i,5]  = rre[[i]]$rmse.grscad[4]
      rmseFinal[i,6]  = rre[[i]]$rmse.mcp[4]
      rmseFinal[i,7]  = rre[[i]]$rmse.grmcp[4]
      rmseFinal[i,8]  = rre[[i]]$rmse.spma[4]

      # Distance to beta_0
      dist2beta[1]  = dist2beta[1] + rre[[i]]$dist2beta.unpel
      dist2beta[2]  = dist2beta[2] + rre[[i]]$dist2beta.lasso
      dist2beta[3]  = dist2beta[3] + rre[[i]]$dist2beta.grlasso
      dist2beta[4]  = dist2beta[4] + rre[[i]]$dist2beta.scad
      dist2beta[5]  = dist2beta[5] + rre[[i]]$dist2beta.grscad
      dist2beta[6]  = dist2beta[6] + rre[[i]]$dist2beta.mcp
      dist2beta[7]   = dist2beta[7] + rre[[i]]$dist2beta.grmcp
      dist2beta[8]  = dist2beta[8] + rre[[i]]$dist2beta.spma


      # TRUE POSITIVE
      # 07/11/2023
      # update tp, fp formula
      tp[1] = tp[1]+ true_positive(rre[[i]]$confuse.unpel)
      tp[2] = tp[2]+ true_positive(rre[[i]]$confuse.lasso)
      tp[3] = tp[3]+ true_positive(rre[[i]]$confuse.grlasso)
      tp[4] = tp[4]+ true_positive(rre[[i]]$confuse.scad)
      tp[5] = tp[5]+ true_positive(rre[[i]]$confuse.grscad)
      tp[6] = tp[6]+ true_positive(rre[[i]]$confuse.mcp)
      tp[7] = tp[7]+ true_positive(rre[[i]]$confuse.grmcp)
      tp[8] = tp[8]+ true_positive(rre[[i]]$confuse.spma)


      # FALSE POSITIVE

      fp[1] = fp[1]+ false_positive(rre[[i]]$confuse.unpel)
      fp[2] = fp[2]+ false_positive(rre[[i]]$confuse.lasso)
      fp[3] = fp[3]+ false_positive(rre[[i]]$confuse.grlasso)
      fp[4] = fp[4]+ false_positive(rre[[i]]$confuse.scad)
      fp[5] = fp[5]+ false_positive(rre[[i]]$confuse.grscad)
      fp[6] = fp[6]+ false_positive(rre[[i]]$confuse.mcp)
      fp[7] = fp[7]+ false_positive(rre[[i]]$confuse.grmcp)
      fp[8] = fp[8]+ false_positive(rre[[i]]$confuse.spma)


    }


    std.rmse[1]  = sd(sdd[,1])
    std.rmse[2]  = sd(sdd[,2])
    std.rmse[3]  = sd(sdd[,3])
    std.rmse[4]  = sd(sdd[,4])
    std.rmse[5]  = sd(sdd[,5])
    std.rmse[6]  = sd(sdd[,6])
    std.rmse[7]  = sd(sdd[,7])
    std.rmse[8]  = sd(sdd[,8])

    std.beta[1]  = sd(sdd.beta[,1])
    std.beta[2]  = sd(sdd.beta[,2])
    std.beta[3]  = sd(sdd.beta[,3])
    std.beta[4]  = sd(sdd.beta[,4])
    std.beta[5]  = sd(sdd.beta[,5])
    std.beta[6]  = sd(sdd.beta[,6])
    std.beta[7]  = sd(sdd.beta[,7])
    std.beta[8]  = sd(sdd.beta[,8])
    # return(list(rmseFinal = log(rmseFinal/rep), sdd = log(sdd), std.rmse = log(std.rmse), tp =tp/rep,
    #             fp =fp/rep,dist2beta =dist2beta/rep,std.beta = std.beta))

    # return(list(rmseFinal = round(log(rmseFinal/rep),3), std.rmse = round(log(std.rmse),3), tp = round(tp/rep*100,3),
    #             fp =round(fp/rep*100,3),dist2beta = round(dist2beta/rep,3),std.beta = round(std.beta,3)))

    # oops = cbind(round(dist2beta/rep,3),round(std.beta,3),
    #              round(tp/rep*100,4),round(fp/rep*100,4),round((rmseFinal/rep),3),round((std.rmse),3)
    #              )
    print(rmseFinal)
    oops = cbind(round(dist2beta/rep,3),round(std.beta,3),
                 round(tp/rep*100,4),round(fp/rep*100,4),round(apply(rmseFinal,2,median),3),round((std.rmse),3)
    )
    colnames(oops) <- c( "dist2beta","std.beta","tp", "fp","rmseFinal", "std.rmse")
    rownames(oops) <- c("unpel", "lasso", "group lasso","scad","group scad","mcp","group mcp","spma")
    return(oops)
  }
  else{
    rep = length(rre)
    #rmseFinal <- rep(0,8)
    # report median
    rmseFinal <- matrix(nrow = rep,ncol = 7)
    dist2beta <- rep(0,7)
    tp<- rep(0,7)
    fp<- rep(0,7)
    sdd <- matrix(nrow = rep, ncol = 7)
    std.rmse <- rep(0,7)

    sdd.beta <- matrix(nrow = rep, ncol = 7)
    std.beta <- rep(0,7)
    for (i in 1:length(rre)){
      #print(i)
      # RMSE
      sdd[i,1] = rre[[i]]$rmse.unpel[4]
      sdd[i,2] = rre[[i]]$rmse.lasso[4]
      sdd[i,3] = rre[[i]]$rmse.scad[4]
      sdd[i,4] = rre[[i]]$rmse.grlasso[4]

      sdd[i,5] = rre[[i]]$rmse.grscad[4]
      # sdd[i,6] = rre[[i]]$rmse.mcp[4]
      sdd[i,6] = rre[[i]]$rmse.grmcp[4]
      sdd[i,7] = rre[[i]]$rmse.spma[4]


      sdd.beta[i,1] = rre[[i]]$dist2beta.unpel
      sdd.beta[i,2] = rre[[i]]$dist2beta.lasso
      sdd.beta[i,3] = rre[[i]]$dist2beta.scad
      sdd.beta[i,4] = rre[[i]]$dist2beta.grlasso

      sdd.beta[i,5] = rre[[i]]$dist2beta.grscad
      #sdd.beta[i,6] = rre[[i]]$dist2beta.mcp
      sdd.beta[i,6] = rre[[i]]$dist2beta.grmcp
      sdd.beta[i,7] = rre[[i]]$dist2beta.spma


      # rmseFinal[1]  = rmseFinal[1] + rre[[i]]$rmse.unpel[4]
      # rmseFinal[2]  = rmseFinal[2] + rre[[i]]$rmse.lasso[4]
      # rmseFinal[3]  = rmseFinal[3] + rre[[i]]$rmse.grlasso[4]
      # rmseFinal[4]  = rmseFinal[4] + rre[[i]]$rmse.scad[4]
      # rmseFinal[5]  = rmseFinal[5] + rre[[i]]$rmse.grscad[4]
      # rmseFinal[6]  = rmseFinal[6] + rre[[i]]$rmse.mcp[4]
      # rmseFinal[7]  = rmseFinal[7] + rre[[i]]$rmse.grmcp[4]
      # rmseFinal[8]  = rmseFinal[8] + rre[[i]]$rmse.spma[4]

      rmseFinal[i,1]  = rre[[i]]$rmse.unpel[4]
      rmseFinal[i,2]  = rre[[i]]$rmse.lasso[4]
      rmseFinal[i,3]  = rre[[i]]$rmse.scad[4]
      rmseFinal[i,4]  = rre[[i]]$rmse.grlasso[4]

      rmseFinal[i,5]  = rre[[i]]$rmse.grscad[4]
      #rmseFinal[i,6]  = rre[[i]]$rmse.mcp[4]
      rmseFinal[i,6]  = rre[[i]]$rmse.grmcp[4]
      rmseFinal[i,7]  = rre[[i]]$rmse.spma[4]

      # Distance to beta_0
      dist2beta[1]  = dist2beta[1] + rre[[i]]$dist2beta.unpel
      dist2beta[2]  = dist2beta[2] + rre[[i]]$dist2beta.lasso
      dist2beta[3]  = dist2beta[3] + rre[[i]]$dist2beta.scad
      dist2beta[4]  = dist2beta[4] + rre[[i]]$dist2beta.grlasso

      dist2beta[5]  = dist2beta[5] + rre[[i]]$dist2beta.grscad
      #dist2beta[6]  = dist2beta[6] + rre[[i]]$dist2beta.mcp
      dist2beta[6]   = dist2beta[6] + rre[[i]]$dist2beta.grmcp
      dist2beta[7]  = dist2beta[7] + rre[[i]]$dist2beta.spma


      # TRUE POSITIVE
      # 07/11/2023
      # update tp, fp formula
      tp[1] = tp[1]+ true_positive(rre[[i]]$confuse.unpel)
      tp[2] = tp[2]+ true_positive(rre[[i]]$confuse.lasso)
      tp[3] = tp[3]+ true_positive(rre[[i]]$confuse.scad)
      tp[4] = tp[4]+ true_positive(rre[[i]]$confuse.grlasso)

      tp[5] = tp[5]+ true_positive(rre[[i]]$confuse.grscad)
      #tp[6] = tp[6]+ true_positive(rre[[i]]$confuse.mcp)
      tp[6] = tp[6]+ true_positive(rre[[i]]$confuse.grmcp)
      tp[7] = tp[7]+ true_positive(rre[[i]]$confuse.spma)


      # FALSE POSITIVE

      fp[1] = fp[1]+ false_positive(rre[[i]]$confuse.unpel)
      fp[2] = fp[2]+ false_positive(rre[[i]]$confuse.lasso)

      fp[3] = fp[3]+ false_positive(rre[[i]]$confuse.scad)
      fp[4] = fp[4]+ false_positive(rre[[i]]$confuse.grlasso)
      fp[5] = fp[5]+ false_positive(rre[[i]]$confuse.grscad)
      #fp[6] = fp[6]+ false_positive(rre[[i]]$confuse.mcp)
      fp[6] = fp[6]+ false_positive(rre[[i]]$confuse.grmcp)
      fp[7] = fp[7]+ false_positive(rre[[i]]$confuse.spma)


    }


    std.rmse[1]  = sd(sdd[,1])
    std.rmse[2]  = sd(sdd[,2])
    std.rmse[3]  = sd(sdd[,3])
    std.rmse[4]  = sd(sdd[,4])
    std.rmse[5]  = sd(sdd[,5])
    #std.rmse[6]  = sd(sdd[,6])
    std.rmse[6]  = sd(sdd[,6])
    std.rmse[7]  = sd(sdd[,7])

    std.beta[1]  = sd(sdd.beta[,1])
    std.beta[2]  = sd(sdd.beta[,2])
    std.beta[3]  = sd(sdd.beta[,3])
    std.beta[4]  = sd(sdd.beta[,4])
    std.beta[5]  = sd(sdd.beta[,5])
    #std.beta[6]  = sd(sdd.beta[,6])
    std.beta[6]  = sd(sdd.beta[,6])
    std.beta[7]  = sd(sdd.beta[,7])
    # return(list(rmseFinal = log(rmseFinal/rep), sdd = log(sdd), std.rmse = log(std.rmse), tp =tp/rep,
    #             fp =fp/rep,dist2beta =dist2beta/rep,std.beta = std.beta))

    # return(list(rmseFinal = round(log(rmseFinal/rep),3), std.rmse = round(log(std.rmse),3), tp = round(tp/rep*100,3),
    #             fp =round(fp/rep*100,3),dist2beta = round(dist2beta/rep,3),std.beta = round(std.beta,3)))

    # oops = cbind(round(dist2beta/rep,3),round(std.beta,3),
    #              round(tp/rep*100,4),round(fp/rep*100,4),round((rmseFinal/rep),3),round((std.rmse),3)
    #              )
    #print(rmseFinal)
    oops = cbind(round(dist2beta/rep,3),round(std.beta,3),
                 round(tp/rep*100,4),round(fp/rep*100,4),round(apply(rmseFinal,2,median),3),round((std.rmse),3)
    )
    colnames(oops) <- c( "dist2beta","std.beta","tp", "fp","rmseFinal", "std.rmse")
    #rownames(oops) <- c("unpel", "lasso", "group lasso","scad","group scad","mcp","group mcp","spma")
    rownames(oops) <- c("unpel", "lasso", "scad", "group lasso","group scad","group mcp","spma")
    return(oops)
  }
}


true_positive <- function(confuse){
  if(ncol(confuse) == 1){
    if(colnames(confuse)=="1"){
      return (confuse[2,1]/rowSums(confuse)[2])
    }
    else{
      return (0)
    }
  }

  else{
    return (confuse[2,2]/rowSums(confuse)[2])
  }

}


false_positive <- function(confuse){
  if(ncol(confuse)==1){
    if(colnames(confuse) == "1"){
      return(confuse[1,1]/rowSums(confuse)[1])
    }
    else{
      return (0)
    }
  }
  else{
    return(confuse[1,2]/rowSums(confuse)[1])
  }
}

# Calculate I_F^m for real data from unpenalized estimators and likelihood
# ifm_real <- function(sample, ni, thetahat, p , M =3  ){
#   # http://sfb649.wiwi.hu-berlin.de/fedc_homepage/xplore/tutorials/xaghtmlnode28.html
#   X <- sample[1:ni,1:p] # expected N x p
#   # thetahat is expected to be  p x1
#   y <- sample[1:ni, p+1] # time
#   status <- as.vector(sample[1:ni,p+2])
#   fisher <- 0
#   for(i in 1:ni){
#     # w_i(beta)
#     oneV = as.vector(ifelse(y>=y[i],1,0)) # N x 1
#     xb = X%*%thetahat
#     #print(paste("Dimension of xb", dim(xb)))
#     #print(paste("size of onev", length(oneV)))
#     wib = t(as.matrix(oneV)) %*% exp(xb) # 1x1, it is a scalar
#     # bar{Z}(i)
#     zbari =  diag(oneV)%*%X  # N x p
#     tmp = as.vector(exp(xb))
#     wib = as.double(wib) # convert wib into a scalar
#     mid = as.matrix(wib*as.matrix(diag(tmp)) - xb%*%t(xb), ncol = ni)
#     fisher = fisher + status[i]/wib/wib*t(zbari)%*%mid%*%zbari
#   }
#   # fisher is expecetd to be p x p
#   if(dim(fisher)!=c(p,p)){
#     print("Wrong dimension!")}
#   return (fisher)
# }


# For Wu and Cook loglikelihood


# --------------- Section Two Yuxiang Wu utilities -----------
#source("~/metaADMM/funICWu.R")

# We used function oracle_regression_case2 for unpenalized
# lasso_regression_group_case2() for group lasso
# scad_regression_group_case2() for group scad
# mcp_regression_group_case2() for group mcp

# Lambda_0n = sum_k=0^m Phi_k^*B_k(t,m,u,v)
Lambda_0n=function(t,phi,m,u,v){
  t=as.vector(t)
  res=rep(0,length(t))
  for (k in 0:m) {
    res[t!=Inf]=res[t!=Inf]+sum(exp(phi[1:(k+1)]))*bernsteinb(k=k,n = m,x = (t[t!=Inf]-v)/(u-v))
  }
  res[t==0]=0
  res[t==Inf]=Inf
  return(res)
}


# general interval censored data
l_n_case2=function(phi,data,beta,m,n){
  res=0
  O=data$observation
  u=ceiling(max(max(O$R[O$R!=Inf]),max(O$L)))
  v=floor(min(O$L))
  Z=data$covdat
  # for (i in 1:n) {
  #   res=res+log(exp(-Lambda_0n(O$L[i],phi,m,u,v)*exp(sum(beta*Z[i,])))-
  #               exp(-Lambda_0n(O$R[i],phi,m,u,v)*exp(sum(beta*Z[i,]))))
  # }
  res=sum(log(exp(-Lambda_0n(O$L,phi,m,u,v)*exp(as.matrix(Z)%*%as.matrix(beta)))-
                exp(-Lambda_0n(O$R,phi,m,u,v)*exp(as.matrix(Z)%*%as.matrix(beta)))))
  return(res)
}

# first derivative
dl_n_case2=function(phi,data,beta,m,n){
  res=rep(0,length(beta))
  O=data$observation
  u=ceiling(max(max(O$R[O$R!=Inf]),max(O$L)))
  v=floor(min(O$L))
  Z=data$covdat
  # for (i in 1:n) {
  #     temp_L=exp(-Lambda_0n(O$L[i],phi,m,u,v)*exp(sum(beta*Z[i,])))
  #     temp_R=exp(-Lambda_0n(O$R[i],phi,m,u,v)*exp(sum(beta*Z[i,])))
  #     res=res+((temp_L*log(temp_L)-temp_R*log(temp_R))/(temp_L-temp_R))*Z[i,]
  # }
  temp_L=exp(-Lambda_0n(O$L,phi,m,u,v)*exp(as.matrix(Z)%*%as.matrix(beta)))
  #print(temp_L)
  temp_R=exp(-Lambda_0n(O$R,phi,m,u,v)*exp(as.matrix(Z)%*%as.matrix(beta)))
  temp_L_log=-Lambda_0n(O$L,phi,m,u,v)*exp(as.matrix(Z)%*%as.matrix(beta))
  #print(paste(" log is ", temp_L_log))
  temp_R_log=-Lambda_0n(O$R,phi,m,u,v)*exp(as.matrix(Z)%*%as.matrix(beta))
  temp=((temp_L*temp_L_log-temp_R*temp_R_log)/(temp_L-temp_R))
  temp[O$R==Inf]=temp_L_log[O$R==Inf]
  res=t(as.matrix(Z))%*%as.matrix(temp)
  return(res)

}

# second order derivative
ddl_n_case2=function(phi,data,beta,m,n){
  res=rep(0,length(beta))
  O=data$observation
  u=ceiling(max(max(O$R[O$R!=Inf]),max(O$L)))
  #print(O)
  v=floor(min(O$L))
  Z=data$covdat
  # for (i in 1:n) {
  #   temp_L=exp(-Lambda_0n(O$L[i],phi,m,u,v)*exp(sum(beta*Z[i,])))
  #   temp_R=exp(-Lambdfa_0n(O$R[i],phi,m,u,v)*exp(sum(beta*Z[i,])))
  #   res=res+((temp_L*log(temp_L)*(1+log(temp_L))-temp_R*log(temp_R)*(1+log(temp_R)))/(temp_L-temp_R)-
  #             ((temp_L*log(temp_L)-temp_R*log(temp_R))/(temp_L-temp_R))^2)*(t(as.matrix(Z[i,]))%*%as.matrix(Z[i,]))
  # }
  temp_L=exp(-Lambda_0n(O$L,phi,m,u,v)*exp(as.matrix(Z)%*%as.matrix(beta)))
  temp_R=exp(-Lambda_0n(O$R,phi,m,u,v)*exp(as.matrix(Z)%*%as.matrix(beta)))
  temp_L_log=-Lambda_0n(O$L,phi,m,u,v)*exp(as.matrix(Z)%*%as.matrix(beta))
  temp_R_log=-Lambda_0n(O$R,phi,m,u,v)*exp(as.matrix(Z)%*%as.matrix(beta))
  temp=(
    (temp_L*temp_L_log*(1+temp_L_log)-temp_R*temp_R_log*(1+temp_R_log))/(temp_L-temp_R)-
      ((temp_L*temp_L_log-temp_R*temp_R_log)/(temp_L-temp_R))^2)
  temp[O$R==Inf]=temp_L_log[O$R==Inf]
  #print(temp_L_log[O$R==Inf])
  #print(temp[O$R==Inf])
  res=(t(as.matrix(Z))%*%diag(as.vector(temp))%*%as.matrix(Z)) # diagonal approximation
  return(res)
}

# Group Lasso
lasso_regression_group_case2=function(initial_val,n,J,p,data,m,lambda,epsilon,maxiter){
  penalty_oracle=function(beta_phi,data,p,m,n){
    beta=beta_phi[1:sum(p)]
    phi=beta_phi[(sum(p)+1):(sum(p)+m+1)]
    res=penalty_initial(beta,data,phi,m,50,n,2)#-2*l_n_case2(phi,data,beta,m,n)
    return(res)
  }

  Q=function(beta,phi_hat,data,m,n){
    res=-2*l_n_case2(phi_hat,data,beta,m,n)
    return(res)
  }
  dQ=function(beta,phi_hat,data,m,n){
    res=-2*dl_n_case2(phi_hat,data,beta,m,n)
    return(res)
  }
  condition=function(beta,beta_hat_temp,t){
    res=Q(beta_hat,phi_hat,data,m,n)<=(Q(beta_hat_temp,phi_hat,data,m,n)
                                       +t(beta_hat-beta_hat_temp)%*%dQ(beta_hat_temp,phi_hat,data,m,n)
                                       +1/(2*t)*t(beta_hat-beta_hat_temp)%*%(beta_hat-beta_hat_temp))
    return(res)
  }
  p_cum=c(0,cumsum(p))
  beta_hat=initial_val$beta#rep(0,sum(p))
  beta_hat_norm=rep(0,J)
  phi_hat=initial_val$phi#rep(0,m+1)
  beta_hat_temp=rep(0,sum(p))#beta_hat
  beta_hat_norm_temp=rep(0,J)
  phi_hat_temp=rep(0,m+1)#phi_hat
  beta_phi_hat=c(beta_hat,phi_hat)
  temp=optim(par=beta_phi_hat, penalty_oracle, gr = NULL,data=data,p=p,m=m,n=n,method = "BFGS",lower = -Inf, upper = Inf,hessian = FALSE)
  beta_hat=temp$par[1:sum(p)]
  phi_hat=temp$par[(sum(p)+1):(sum(p)+m+1)]

  for (j in 1:J) {
    beta_hat_norm[j]=norm(beta_hat[(p_cum[j]+1):p_cum[j+1]], type="2")
  }

  if(lambda==0){
    for (iterations in 1:maxiter) {
      beta_hat=optim(par=beta_hat,l_n_case2, gr = NULL,phi=phi_hat,data=data,m=m,n=n,method = "BFGS",lower = -Inf, upper = Inf,hessian = FALSE,control=list(fnscale=-1,maxit=10))$par
      phi_hat=optim(par=phi_hat,l_n_case2, gr = NULL,beta=beta_hat,data=data,m=m,n=n,method = "BFGS",lower = -Inf, upper = Inf,hessian = FALSE,control=list(fnscale=-1,maxit=10))$par
      if (sqrt((sum((beta_hat-beta_hat_temp)^2)+sum((phi_hat-phi_hat_temp)^2)))<epsilon) {
        break
      }else if(iterations==maxiter){
        warning("max iteration is attained")
      }
      beta_hat_temp=beta_hat
      phi_hat_temp=phi_hat
    }
    for (j in 1:J) {
      beta_hat_norm[j]=norm(beta_hat[(p_cum[j]+1):p_cum[j+1]], type="2")
    }
    res=list(beta_hat=beta_hat,beta_hat_norm=beta_hat_norm,phi_hat=phi_hat)
    return(res)
  }

  for (iterations in 1:maxiter) {

    beta_hat_temp=beta_hat
    beta_hat_norm_temp=beta_hat_norm
    phi_hat_temp=phi_hat


    for (j in 1:J) {
      t=0.1
      while (TRUE) {
        temp=beta_hat_temp[(p_cum[j]+1):p_cum[j+1]]-t*dQ(beta_hat_temp,phi_hat,data,m,n)[(p_cum[j]+1):p_cum[j+1]]
        temp_norm=sqrt(t(temp)%*%temp)
        if(temp_norm>=t*lambda*sqrt(p[j])){
          beta_hat[(p_cum[j]+1):p_cum[j+1]]=as.numeric((1-t*lambda*sqrt(p[j])/temp_norm))*temp
        }else{
          beta_hat[(p_cum[j]+1):p_cum[j+1]]=0
        }
        if(condition(beta_hat,beta_hat_temp,t)==TRUE){
          break
        }else{
          t=0.5*t
        }

      }
    }


    phi_hat=optim(phi_hat,l_n_case2, gr = NULL,data=data,beta=beta_hat,m=m,n=n,method = "BFGS",lower = -Inf, upper = Inf,hessian = FALSE,control=list(fnscale=-1,maxit=20))
    phi_hat=phi_hat$par




    for (j in 1:J) {

      beta_hat_norm[j]=norm(beta_hat[(p_cum[j]+1):p_cum[j+1]], type="2")
    }


    if(l_n_case2(phi_hat,data,beta_hat,m,n)==-Inf){
      warning('likelihood function is infinite')
      res=list(beta_hat=beta_hat,beta_hat_norm=beta_hat_norm,phi_hat=phi_hat)
      return(res)

    }

    #convergence condition
    if (sqrt((sum((beta_hat-beta_hat_temp)^2)+sum((phi_hat-phi_hat_temp)^2)))<epsilon) {
      break
    }else if(iterations==maxiter){
      warning("max iteration is attained")
    }
  }
  res=list(beta_hat=beta_hat,beta_hat_norm=beta_hat_norm,phi_hat=phi_hat)
  return(res)
}

# Group MCP
mcp_regression_group_case2=function(initial_val,n,J,p,data,m,gamma,lambda,epsilon,maxiter){
  # m = 3
  penalty_oracle=function(beta_phi,data,p,m,n){
    beta=beta_phi[1:sum(p)]
    phi=beta_phi[(sum(p)+1):(sum(p)+m+1)]
    res=penalty_initial(beta,data,phi,m,50,n,2)#-2*l_n_case2(phi,data,beta,m,n)
    return(res)
  }

  Q=function(beta,phi_hat,data,m,n){
    res=-2*l_n_case2(phi_hat,data,beta,m,n)
    return(res)
  }
  dQ=function(beta,phi_hat,data,m,n){
    res=-2*dl_n_case2(phi_hat,data,beta,m,n)
    return(res)
  }
  S=function(z,lambda){
    res=as.numeric((1-lambda/sqrt(t(z)%*%z))*(1>lambda/sqrt(t(z)%*%z)))*z
    return(res)
  }
  condition=function(beta,beta_hat_temp,t){
    res=Q(beta_hat,phi_hat,data,m,n)<=(Q(beta_hat_temp,phi_hat,data,m,n)
                                       +t(beta_hat-beta_hat_temp)%*%dQ(beta_hat_temp,phi_hat,data,m,n)
                                       +1/(2*t)*t(beta_hat-beta_hat_temp)%*%(beta_hat-beta_hat_temp))
    return(res)
  }


  p_cum=c(0,cumsum(p))
  beta_hat=initial_val$beta#rep(0,sum(p))
  beta_hat_norm=rep(0,J)
  phi_hat=initial_val$phi#rep(0,m+1)
  beta_hat_temp=rep(0,sum(p))#beta_hat
  beta_hat_norm_temp=rep(0,J)
  phi_hat_temp=rep(0,m+1)#phi_hat
  beta_phi_hat=c(beta_hat,phi_hat)
  temp=optim(par=beta_phi_hat, penalty_oracle, gr = NULL,data=data,p=p,m=m,n=n,method = "BFGS",lower = -Inf, upper = Inf,hessian = FALSE)
  beta_hat=temp$par[1:sum(p)]
  phi_hat=temp$par[(sum(p)+1):(sum(p)+m+1)]
  for (j in 1:J) {
    beta_hat_norm[j]=norm(beta_hat[(p_cum[j]+1):p_cum[j+1]], type="2")
  }



  for (iterations in 1:maxiter) {

    beta_hat_temp=beta_hat
    beta_hat_norm_temp=beta_hat_norm
    phi_hat_temp=phi_hat
    for (j in 1:J) {
      t=0.1
      while (TRUE) {
        temp=beta_hat_temp[(p_cum[j]+1):p_cum[j+1]]-t*dQ(beta_hat_temp,phi_hat,data,m,n)[(p_cum[j]+1):p_cum[j+1]]
        temp_norm=sqrt(t(temp)%*%temp)
        if(temp_norm>=t*gamma*lambda*sqrt(p[j])){
          beta_hat[(p_cum[j]+1):p_cum[j+1]]=temp
        }else{
          beta_hat[(p_cum[j]+1):p_cum[j+1]]=(gamma/(gamma-1))*S(temp,t*lambda*sqrt(p[j]))
        }
        if(condition(beta_hat,beta_hat_temp,t)==TRUE){
          break
        }else{
          t=0.5*t
        }
      }
    }
    phi_hat=optim(phi_hat,l_n_case2, gr = NULL,data=data,beta=beta_hat,m=m,n=n,method = "BFGS",lower = -Inf, upper = Inf,hessian = FALSE,control=list(fnscale=-1,maxit=20))
    phi_hat=phi_hat$par
    ## update beta_hat_norm
    for (j in 1:J) {
      beta_hat_norm[j]=norm(beta_hat[(p_cum[j]+1):p_cum[j+1]], type="2")
    }
    #step 4
    if(l_n_case2(phi_hat,data,beta_hat,m,n)==-Inf){
      warning('likelihood function is infinite')
      res=list(beta_hat=beta_hat,beta_hat_norm=beta_hat_norm,phi_hat=phi_hat)
      return(res)

    }

    if (sqrt((sum((beta_hat-beta_hat_temp)^2)+sum((phi_hat-phi_hat_temp)^2)))<epsilon) {
      break
    }else if(iterations==maxiter){
      warning("max iteration is attained")
    }
  }


  res=list(beta_hat=beta_hat,beta_hat_norm=beta_hat_norm,phi_hat=phi_hat)
  return(res)
}


# Group SCAD
scad_regression_group_case2=function(initial_val,n,J,p,data,m,gamma,lambda,epsilon,maxiter){
  penalty_oracle=function(beta_phi,data,p,m,n){
    beta=beta_phi[1:sum(p)]
    phi=beta_phi[(sum(p)+1):(sum(p)+m+1)]
    res=penalty_initial(beta,data,phi,m,50,n,2)#-2*l_n_case2(phi,data,beta,m,n)
    return(res)
  }

  Q=function(beta,phi_hat,data,m,n){
    res=-2*l_n_case2(phi_hat,data,beta,m,n)
    return(res)
  }
  dQ=function(beta,phi_hat,data,m,n){
    res=-2*dl_n_case2(phi_hat,data,beta,m,n)
    return(res)
  }
  S=function(z,lambda){
    res=as.numeric((1-lambda/sqrt(t(z)%*%z))*(1>lambda/sqrt(t(z)%*%z)))*z
    return(res)
  }
  condition=function(beta,beta_hat_temp,t){
    res=Q(beta_hat,phi_hat,data,m,n)<=(Q(beta_hat_temp,phi_hat,data,m,n)
                                       +t(beta_hat-beta_hat_temp)%*%dQ(beta_hat_temp,phi_hat,data,m,n)
                                       +1/(2*t)*t(beta_hat-beta_hat_temp)%*%(beta_hat-beta_hat_temp))
    return(res)
  }

  p_cum=c(0,cumsum(p))
  beta_hat=initial_val$beta#rep(0,sum(p))
  beta_hat_norm=rep(0,J)
  phi_hat=initial_val$phi#rep(0,m+1)
  beta_hat_temp=rep(0,sum(p))#beta_hat
  beta_hat_norm_temp=rep(0,J)
  phi_hat_temp=rep(0,m+1)#phi_hat
  beta_phi_hat=c(beta_hat,phi_hat)
  temp=optim(par=beta_phi_hat, penalty_oracle, gr = NULL,data=data,p=p,m=m,n=n,method = "BFGS",lower = -Inf, upper = Inf,hessian = FALSE)
  beta_hat=temp$par[1:sum(p)]
  phi_hat=temp$par[(sum(p)+1):(sum(p)+m+1)]
  for (j in 1:J) {
    beta_hat_norm[j]=norm(beta_hat[(p_cum[j]+1):p_cum[j+1]], type="2")
  }



  for (iterations in 1:maxiter) {

    beta_hat_temp=beta_hat
    beta_hat_norm_temp=beta_hat_norm
    phi_hat_temp=phi_hat
    for (j in 1:J) {
      t=0.1
      while (TRUE) {
        temp=beta_hat_temp[(p_cum[j]+1):p_cum[j+1]]-t*dQ(beta_hat_temp,phi_hat,data,m,n)[(p_cum[j]+1):p_cum[j+1]]
        temp_norm=sqrt(t(temp)%*%temp)
        if(temp_norm>=t*gamma*lambda*sqrt(p[j])){
          beta_hat[(p_cum[j]+1):p_cum[j+1]]=temp
        }else if((temp_norm<t*gamma*lambda*sqrt(p[j]))&(temp_norm>=2*t*lambda*sqrt(p[j]))){
          beta_hat[(p_cum[j]+1):p_cum[j+1]]=(gamma-1)/(gamma-2)*S(temp,gamma*t*lambda*sqrt(p[j])/(gamma-1))
        }else if(temp_norm<2*t*lambda*sqrt(p[j])){
          beta_hat[(p_cum[j]+1):p_cum[j+1]]=S(temp,t*lambda*sqrt(p[j]))
        }
        if(condition(beta_hat,beta_hat_temp,t)==TRUE){
          break
        }else{
          t=0.5*t
        }
      }
    }
    phi_hat=optim(phi_hat,l_n_case2, gr = NULL,data=data,beta=beta_hat,m=m,n=n,method = "BFGS",lower = -Inf, upper = Inf,hessian = FALSE,control=list(fnscale=-1,maxit=20))
    phi_hat=phi_hat$par
    ## update beta_hat_norm
    for (j in 1:J) {
      beta_hat_norm[j]=norm(beta_hat[(p_cum[j]+1):p_cum[j+1]], type="2")
    }
    #step 4
    if(l_n_case2(phi_hat,data,beta_hat,m,n)==-Inf){
      warning('likelihood function is infinite')
      res=list(beta_hat=beta_hat,beta_hat_norm=beta_hat_norm,phi_hat=phi_hat)
      return(res)

    }

    if (sqrt((sum((beta_hat-beta_hat_temp)^2)+sum((phi_hat-phi_hat_temp)^2)))<epsilon) {
      break
    }else if(iterations==maxiter){
      warning("max iteration is attained")
    }
  }


  res=list(beta_hat=beta_hat,beta_hat_norm=beta_hat_norm,phi_hat=phi_hat)
  return(res)
}

# Oracle/ Unpenalized
oracle_regression_case2=function(initial_val,n,J,p,data,m,epsilon,maxiter){
  # penalty_oracle=function(beta_phi,data,p,m,n){
  #   beta=beta_phi[1:sum(p)]
  #   phi=beta_phi[(sum(p)+1):(sum(p)+m+1)]
  #   res=-2*l_n_case2(phi,data,beta,m,n)
  #   return(res)
  # }
  # p_cum=c(0,cumsum(p))
  # beta_hat=initial_val$beta
  # beta_hat_norm=rep(0,J)
  # phi_hat=initial_val$phi
  # beta_phi_hat=c(beta_hat,phi_hat)
  # temp=optim(par=beta_phi_hat, penalty_oracle, gr = NULL,data=data,p=p,m=m,n=n,method = "BFGS",lower = -Inf, upper = Inf,hessian = FALSE)
  # beta_hat=temp$par[1:sum(p)]
  # phi_hat=temp$par[(sum(p)+1):(sum(p)+m+1)]
  p_cum=c(0,cumsum(p))
  beta_hat=initial_val$beta
  beta_hat_temp=rep(0,sum(p))
  beta_hat_norm=rep(0,J)
  phi_hat=initial_val$phi
  phi_hat_temp=rep(0,m+1)
  for (iterations in 1:maxiter) {
    beta_hat=optim(par=beta_hat,l_n_case2, gr = dl_n_case2,phi=phi_hat,data=data,m=m,n=n,method = "CG",lower = -Inf, upper = Inf,hessian = FALSE,control=list(fnscale=-1,maxit=10))$par
    phi_hat=optim(par=phi_hat,l_n_case2, gr = NULL,beta=beta_hat,data=data,m=m,n=n,method = "CG",lower = -Inf, upper = Inf,hessian = FALSE,control=list(fnscale=-1,maxit=10))$par
    if (sqrt((sum((beta_hat-beta_hat_temp)^2)+sum((phi_hat-phi_hat_temp)^2)))<epsilon) {
      break
    }else if(iterations==maxiter){
      warning("max iteration is attained")
    }
    beta_hat_temp=beta_hat
    phi_hat_temp=phi_hat
  }
  for (j in 1:J) {
    beta_hat_norm[j]=norm(beta_hat[(p_cum[j]+1):p_cum[j+1]], type="2")
  }
  res=list(beta_hat=beta_hat,beta_hat_norm=beta_hat_norm,phi_hat=phi_hat)
  return(res)
}





# ---------------- Section Three Wu & Cook utilities ---------------------

EM_f <- function(indata, lam0, beta0, lasso.lam, ncov, npieces, cutpoints, penalty.function, penalty.factor = NULL, nopenalty.index = NULL, thresh = 1e-6, maxit = 200) {

  cur.lam  <- lam0
  cur.beta <- beta0

  if(is.null(penalty.factor)){
    penalty.factor <- penalty.factor.f(beta0=cur.beta, lasso.lam=lasso.lam,
                                       ncov=ncov, npieces=npieces,
                                       penalty.function=penalty.function,
                                       nopenalty.index = nopenalty.index)
  }

  iter <- 0
  tol  <- 9999
  while ( tol > thresh ) {
    if(penalty.function=="alasso" | penalty.function == "scad"){
      penalty.factor <- penalty.factor.f(beta0=cur.beta, lasso.lam=lasso.lam,
                                         ncov=ncov, npieces=npieces,
                                         penalty.function=penalty.function,
                                         nopenalty.index = nopenalty.index)
    }
    iter <- iter + 1

    pre.lam  <- cur.lam
    pre.beta <- cur.beta

    pseudodata <- createdata.f(indata=indata, lam=pre.lam, beta=pre.beta, ncov=ncov, cutpoints=cutpoints)
    pseudodata <- pseudodata[!is.na(pseudodata$EIk),]
    pseudodata <- pseudodata[!is.na(pseudodata$logEwk),]

    if ( npieces == 1 ) {
      y <- pseudodata$EIk
      x <- as.matrix( pseudodata[,c(paste("x",1:ncov,sep=""))] )
      adjterm <- pseudodata$logEwk

      fit <- glmnet(x=x, y=y, family="poisson", offset=adjterm,
                    penalty.factor=penalty.factor,
                    lambda=lasso.lam, alpha=1)
      cf <- coef(fit)
      cur.lam  <- exp( cf[1] )
      cur.beta <- cf[-1]

      y <- NULL; x <- NULL; adjterm <- NULL; fit <- NULL; cf <- NULL
    }
    else {
      pieces <- matrix(0, nrow=nrow(pseudodata), ncol=(npieces-1))
      for (k in 2:npieces) {
        pieces[,k-1] <- ifelse(pseudodata$piece == k, 1, 0)
      }
      y <- pseudodata$EIk
      x <- as.matrix( cbind(pieces, as.matrix(pseudodata[,c(paste("x",1:ncov,sep=""))])) )
      adjterm <- pseudodata$logEwk

      fit <- glmnet(x=x, y=y, family="poisson", offset=adjterm,
                    penalty.factor=penalty.factor,
                    lambda=lasso.lam, alpha=1)

      cf <- as.vector(coef(fit))
      cur.beta <- cf[-c(1:npieces)]

      cf <- cf[c(1:npieces)]
      cur.lam <- exp( c(0, cf[-1]) + cf[1] )

      pieces <- NULL; y <- NULL; x <- NULL; adjterm <- NULL; fit <- NULL; cf <- NULL
    }
    cur.lam  <- ifelse(is.na(cur.lam), 0, cur.lam)
    cur.beta <- ifelse(is.na(cur.beta), 0, cur.beta)

    dif.lam  <- abs( (cur.lam - pre.lam) / pre.lam )
    dif.beta <- abs( (cur.beta - pre.beta) / pre.beta )
    dif.beta <- ifelse(pre.beta == 0, abs(cur.beta - pre.beta), dif.beta)

    tol <- max( c(dif.lam, dif.beta) )
    if ( iter > maxit ) { break }
  }

  out <- NULL
  out$tol  <- tol
  out$iter <- iter
  out$beta <- as.vector(cur.beta)
  out$lam  <- as.vector(cur.lam)
  return(out)

}


penalty.factor.f <- function(beta0 = NULL, lasso.lam = NULL, ncov, npieces, penalty.function, nopenalty.index = NULL){
  if(penalty.function == "lasso"){
    pf <- rep(1, ncov)
    pf[nopenalty.index] <- 0
    penalty.factor <- c(rep(0, npieces-1), pf)
  }
  if(penalty.function == "alasso"){
    beta0.na <- beta0
    beta0.na[beta0 == 0] <- 10^(-5)
    pf <- 1/abs(beta0.na)
    pf[nopenalty.index] <- 0
    penalty.factor <- c(rep(0, npieces-1), pf)
  }
  scad.pen.deriv <- function(beta0, lasso.lam){
    abs.beta0 <- abs(beta0)
    indicator <- ifelse(abs.beta0 > lasso.lam, 0, 1)
    deriv <- lasso.lam*(indicator + (1-indicator)*ifelse(3.7*lasso.lam - abs.beta0 >= 0, 3.7*lasso.lam - abs.beta0, 0)/2.7/lasso.lam)
    return(deriv)
  }
  if(penalty.function == "scad"){
    deriv <- scad.pen.deriv(beta0=beta0, lasso.lam=lasso.lam)
    deriv[nopenalty.index] <- 0
    penalty.factor <- c(rep(0, npieces-1),deriv)
  }
  return(penalty.factor)
}


# Create Pseudo Data Set
createdata.f <- function(indata, lam, beta, ncov, cutpoints) {
  outdata <- lapply(1:nrow(indata), function(ith, indata, lam, beta, ncov, cutpoints) {
    out <- Eterm.f(Li=indata$timeL[ith],
                   Ri=indata$timeR[ith],
                   zi=as.vector(unlist(indata[ith, paste("x",1:ncov,sep="")])),
                   lam=lam, beta=beta,
                   ncov=ncov,
                   cutpoints=cutpoints)
    nlen <- length(out$piece)

    out$id <- rep(ith, nlen)
    out <- data.frame(out[,c("id","piece","EIk","Ewk","logEwk")])

    outcov <- apply(indata[ith,paste("x",1:ncov,sep="")], 2, rep, times = nlen)
    outcov <- matrix(outcov, nrow=nlen, ncol=ncov)
    row.names(outcov) <- 1:nlen
    outcov <- data.frame(outcov)

    outdata <- cbind(out, outcov)
    return(outdata)
  }, indata=indata, lam=lam, beta=beta, ncov=ncov, cutpoints=cutpoints)
  outdata <- do.call("rbind", outdata)
  outdata <- data.frame(outdata)
  dimnames(outdata)[[2]] <- c("id","piece","EIk","Ewk","logEwk",paste("x",1:ncov,sep=""))
  outdata <- outdata[order(outdata$id),]
  return(outdata)
}

# Evaluate Fbar(t[i] | Z[i]; theta)
Fbar.f <- function(ti, nzi, zi, lam, beta, ncov, cutpoints) {
  #print(paste0(" ooo ", length(ti)))
  cstart <- cutpoints$start
  cstop  <- cutpoints$stop
  npiece <- length(cstart)
  nsubj  <- length(ti) # what is ti

  beta <- matrix(beta, nrow=ncov, ncol=1)
  zi   <- matrix(zi, nrow=nzi, ncol=ncov)
  zi.times.beta <- as.vector( zi %*% beta ) # zi is ith row in the Z matrix

  Ht <- rep(0,nsubj)
  for (k in 1:npiece) {
    wk <- rep(0,nsubj)
    wk <- ifelse( (ti <= cstop[k]) & (ti >= cstart[k]), ti - cstart[k], wk)
    wk <- ifelse( (ti > cstart[k]) & (ti > cstop[k]), rep(cstop[k] - cstart[k],nsubj), wk) # t_ij

    Ht <- Ht + ( (wk*lam[k])*exp( zi.times.beta ) )
  }

  Fbar <- exp( (-1)*Ht )
  return(Fbar)
}

logL.f <- function(T, Z, lam, beta, ncov, cutpoints){
  # T is the interval time
  len = nrow(T)
  logL = 0
  for(i in 1:len){
    tmpL = Fbar.f(T[i,1], 1, Z[i,], lam, beta, ncov, cutpoints)
    #print(tmpL)
    tmpR = Fbar.f(T[i,2], 1, Z[i,], lam, beta, ncov, cutpoints)
    #print(tmpR)
    logL = logL + tmpL - log(tmpR+0.0000000001)
  }

  return (logL)
}



# E-Step:
# Evaluate E(I[k](u[i]) | D[i]) and E(S[k](u[i]) | D[i]))

EIu.f <- function(ti, zi, lam, beta, ncov, cutpoints, Fterm) {
  Fbar <- Fbar.f(ti=ti, nzi=1, zi=zi, lam=lam, beta=beta, ncov=ncov, cutpoints=cutpoints)
  #print(paste0("OKOKOK ", Fbar, " tititit ", ti))
  term <- Fbar - Fterm
  term <- ifelse(is.na(term),  0, term)
  term <- ifelse(is.nan(term), 0, term)
  return( term )
}

Eterm.f <- function(Li, Ri, zi, lam, beta, ncov, cutpoints) {
  npiece <- max(cutpoints$piece)

  outi <- cutpoints
  outi$EIk <- rep(0, npiece)
  outi$Ewk <- rep(0, npiece)

  if(Ri > 9999) { # my upper bound is inf, In wu's example data, Ri = 9999 , if it is rightcensored case
    for (k in 1:npiece) {
      EIk <- 0
      Ewk <- min(c(Li,outi$stop[k])) - min( c(Li, outi$start[k]) )
      outi$EIk[k] <- EIk
      outi$Ewk[k] <- Ewk
    }
  }
  else{
    FbarLi <- Fbar.f(ti=Li, nzi=1, zi=zi, lam=lam, beta=beta, ncov=ncov, cutpoints=cutpoints)
    FbarRi <- Fbar.f(ti=Ri, nzi=1, zi=zi, lam=lam, beta=beta, ncov=ncov, cutpoints=cutpoints)

    for (k in 1:npiece) {
      EIk <- 0
      Ewk <- 0
      if (Ri < outi$start[k]) {
        EIk <- 0
        Ewk <- 0
      }
      else if (outi$stop[k] < Li) {
        EIk <- 0
        Ewk <- outi$stop[k] - outi$start[k]
      }
      else {
        FbarLik <- Fbar.f(ti=max(c(Li, outi$start[k])), nzi=1, zi=zi, lam=lam, beta=beta, ncov=ncov, cutpoints=cutpoints)
        FbarRik <- Fbar.f(ti=min(c(outi$stop[k], Ri)),  nzi=1, zi=zi, lam=lam, beta=beta, ncov=ncov, cutpoints=cutpoints)

        EIk <- (FbarLik - FbarRik)/(FbarLi - FbarRi)

        term1 <- max( c(Li - outi$start[k], 0) )

        lower.val <- max( c(Li, outi$start[k]) )
        upper.val <- min( c(outi$stop[k], Ri) )
        upper.val <- ifelse(upper.val >= 9999, Inf, upper.val)
        int.term <- integrate(EIu.f, lower=lower.val, upper=upper.val,
                              zi=zi, lam=lam, beta=beta, ncov=ncov, cutpoints=cutpoints, Fterm=FbarRi)
        term2 <- int.term$value/(FbarLi - FbarRi)
        Ewk <- term1 + term2
      }

      outi$EIk[k] <- EIk
      outi$Ewk[k] <- Ewk
    }
  }

  outi$logEwk <- log(outi$Ewk)
  outi$logEwk <- ifelse(outi$Ewk == 0, NA, outi$logEwk)

  return( outi[,c("piece","EIk","Ewk","logEwk")] )
}



