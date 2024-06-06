library(MASS)
library(grpreg)
library(gglasso)
library(survival)
library(glmnet)
library(ncvreg)


SPMAreal2 <- function(realData,iter = 1000, rightCensor = FALSE, intervalCensor = FALSE,m = 5){
  # fit model 
  sample1 = realData$sample1;
  sample2 = realData$sample2;
  N = realData$N; 
  p = realData$p;
  M = realData$M
  
  
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
    
    if(TRUE){  
      # unpenalized estimation Cox proportional model for study one 
      m1 = coxph(Surv(time,status)~.,sample1)
      covF1 = m1$var
      theta1F = m1$coefficients
      
      m2 = coxph(Surv(time,status)~.,sample2)
      covF2 = m2$var
      theta2F = m2$coefficients
      
    
      
      #View(covF1)
      #print(eigen(covF1))
      covF1 <- solve(covF1)/N[1]
      covF2 <- solve(covF2)/N[2]
      
      
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
      #theta3 = as.vector(coef(best_model3))
      
      # SCAD estimation
      scad1 <- ncvsurv(X1,y1,alpha = 1, penalty = "SCAD")
      scad2 <- ncvsurv(X2,y2,alpha = 1, penalty = "SCAD")
      
      
      # MCP estimation
      mcp1 <- ncvsurv(X1,y1,alpha = 1, penalty = "MCP")
      mcp2 <- ncvsurv(X2,y2,alpha = 1, penalty = "MCP")
    } 
    # grouplasso
    yy = cbind(c(y1[,1], y2[,1]), c(y1[,2], y2[,2]))
    xx = matrix(0,nrow = (N[1]+N[2]), ncol = p*2)
    xx[1:N[1], 1:p] = X1
    xx[(N[1]+1):(N[1]+N[2]), (p+1):(2*p)] =X2
    group = c(c(1:p), c(1:p))
    
    
    #grlasso <- grpsurv(xx,yy,group = group,nlambda = 4001, lambda = 2^(seq(-30,10,0.01)),gamma = 3, penalty = "grLasso",alpha = 1, keep = TRUE)
    # 
    #grscad <- grpsurv(xx,yy,group = group, nlambda = 201,lambda=2^(seq(-30,10,0.2)), gamma =4, penalty = "grSCAD",alpha = 1, keep = TRUE)
    # 
    #grmcp <- grpsurv(xx,yy,group = group, nlambda = 201, lambda=2^(seq(-30,10,0.2)),gamma =3, penalty = "grMCP",alpha = 1, keep = TRUE)
    

    grlasso <- grpsurv(xx,yy,group = group, penalty = "grLasso",alpha = 1)

    grscad <- grpsurv(xx,yy,group = group, penalty = "grSCAD",alpha = 1)

    grmcp <- grpsurv(xx,yy,group = group, penalty = "grMCP",alpha = 1)

        if(TRUE){   
      # SPMA 
      lambda <- 2^(seq(-30,10,1))
      #lambda = 2^(seq(-30,10,0.2))
      rho = 10;# range of rho
      red = list() # keep in mind that for simulatio study, rc = TRUE, rho = 10, alpha = 1.2, abs 10e-5 rel 10e-3
      for (la in 1:length(lambda)){
        #rho = la
        augIMat = matrix(M*p, M*p);
        augIMat = as.matrix(Matrix::bdiag(covF1, covF2))
        # augNMat = bdiag(200*diag(p),500*diag(p),800*diag(p))
        #print(augNMat%*%augIMat)
        # rightCensor simulation study
        # you = spmaC(alpha = 1.2, epi_abs = 10e-5,epi_rel = 10e-3,rho, lambda[la],iter = iter,covF1 = covF1,  covF2 = covF2, covF3 = covF3,augIMat = augIMat, p =p,M=M,setN = N,thetaF = rbind(theta1F,theta2F,theta3F))
        # test for REAL 
        you = spmaC(alpha = 1.5, epi_abs = 10e-5,epi_rel = 10e-3,rho, lambda[la],iter = iter,covF1 = covF1,  covF2 = covF2, augIMat = augIMat, p =p,M=M,setN = N,thetaF = rbind(theta1F,theta2F))
        
        #you = spmaC(alpha = 1.5, epi_abs = 10e-5,epi_rel = 10e-3,rho, lambda[la],iter = iter,covF1 = covF1,  covF2 = covF2, covF3 = covF3,augIMat = augIMat, p =p,M=M,setN = N,thetaF = rbind(theta1F,theta2F,theta3F))
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
      unpel = rbind(theta1F, theta2F)}
    beta_0 = realData$beta_0
    return(list(sample1 = sample1, sample2 = sample2, lasso1 = lasso1, lasso2 = lasso2, 
                scad1 = scad1, scad2 = scad2, mcp1 = mcp1, mcp2 = mcp2, 
                unpel = unpel, beta_0 = beta_0,
                grlasso = grlasso, grscad = grscad, grmcp = grmcp, spma = red, spma.lambda = lambda,
                N = N, p = p, covF1 =covF1, covF2 = covF2, correlate = realData$correlate,intervalCensor = realData$intervalCensor, m1= m1, m2 = m2))
    
    # return(list(sample1 = sample1, sample2 = sample2, sample3 = sample3,
    #             grlasso = grlasso, grscad = grscad, grmcp = grmcp, 
    #             N = N, p = p, correlate = realData$correlate,intervalCensor = realData$intervalCensor))
    
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
                unpel = unpel, beta_0 =realData$beta_0,intervalCensor = realData$intervalCensor,spma = red,
                spma.lambda = lambda.spma,N=N, p = p, covF1 = covF1, covF2 = covF2, covF3 = covF3,
                correlate = realData$correlate))
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

bic_model_selection <- function(trainSPMA.object, real = FALSE){
  
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
    p = trainSPMA.object$p
    
    if(TRUE){
      # LASSO sample 1 
      thetahat.lasso1 =as.matrix(trainSPMA.object$lasso1$beta)
      num.col = ncol(thetahat.lasso1)
      bic.lasso1 = cbind(rep(0,num.col), rep(0,num.col))
      
      for(i in 1:num.col){
        bb = bic_lasso(N[1],thetahat = thetahat.lasso1[,i],sample1,p, M = 2)
        bic.lasso1[i,1:2] = c(bb$p1,bb$p2)
      }
      bic.lasso1 = bic.lasso1[-1,]
      idx.lasso1 = which.min(bic.lasso1[,1]+bic.lasso1[,2])+1
      #idx.lasso1 = which.min(bic.lasso1[,1]+bic.lasso1[,2])
      
      # LASSO sample 2 
      thetahat.lasso2 =as.matrix(trainSPMA.object$lasso2$beta)
      num.col = ncol(thetahat.lasso2)
      bic.lasso2 = cbind(rep(0,num.col), rep(0,num.col))
      
      for(i in 1:num.col){
        bb = bic_lasso(N[2],thetahat = thetahat.lasso2[,i],sample2,p, M = 2)
        bic.lasso2[i,1:2] = c(bb$p1,bb$p2)
      }
      bic.lasso2 = bic.lasso2[-1,]
      idx.lasso2 = which.min(bic.lasso2[,1]+bic.lasso2[,2])+1
      #idx.lasso2 = which.min(bic.lasso2[,1]+bic.lasso2[,2])
      
      
    
      
      # SCAD
      thetahat.scad1 =as.matrix(trainSPMA.object$scad1$beta)
      num.col = ncol(thetahat.scad1)
      bic.scad1 = cbind(rep(0,num.col), rep(0,num.col))
      
      for(i in 1:num.col){
        bb = bic_lasso(N[1],thetahat = thetahat.scad1[,i],sample1,p, M = 2)
        bic.scad1[i,1:2] = c(bb$p1,bb$p2)
      }
      #bic.scad1 = bic.scad1[-1,]
      #idx.scad1 = which.min(bic.scad1[,1]+bic.scad1[,2])+1
      idx.scad1 = which.min(bic.scad1[,1]+bic.scad1[,2])
      
      # group 2
      thetahat.scad2 =as.matrix(trainSPMA.object$scad2$beta)
      num.col = ncol(thetahat.scad2)
      bic.scad2 = cbind(rep(0,num.col), rep(0,num.col))
      
      for(i in 1:num.col){
        bb = bic_lasso(N[2],thetahat = thetahat.scad2[,i],sample2,p, M = 2)
        bic.scad2[i,1:2] = c(bb$p1,bb$p2)
      }
      bic.scad2 = bic.scad2[-1,]
      idx.scad2 = which.min(bic.scad2[,1]+bic.scad2[,2])+1
      #idx.scad2 = which.min(bic.scad2[,1]+bic.scad2[,2])
      
     
      # MCP
      thetahat.mcp1 =as.matrix(trainSPMA.object$mcp1$beta)
      num.col = ncol(thetahat.mcp1)
      bic.mcp1 = cbind(rep(0,num.col), rep(0,num.col))
      
      for(i in 1:num.col){
        bb = bic_lasso(N[1],thetahat = thetahat.mcp1[,i],sample1,p, M = 2)
        bic.mcp1[i,1:2] = c(bb$p1,bb$p2)
      }
      bic.mcp1 = bic.mcp1[-1,]
      idx.mcp1 = which.min(bic.mcp1[,1]+bic.mcp1[,2])+1
      #idx.mcp1 = which.min(bic.mcp1[,1]+bic.mcp1[,2])
      
      # group 2
      thetahat.mcp2 =as.matrix(trainSPMA.object$mcp2$beta)
      num.col = ncol(thetahat.mcp2)
      bic.mcp2 = cbind(rep(0,num.col), rep(0,num.col))
      
      for(i in 1:num.col){
        bb = bic_lasso(N[2],thetahat = thetahat.mcp2[,i],sample2,p, M = 2)
        bic.mcp2[i,1:2] = c(bb$p1,bb$p2)
      }
      bic.mcp2 = bic.mcp2[-1,]
      idx.mcp2 = which.min(bic.mcp2[,1]+bic.mcp2[,2])+1
      #idx.mcp2 = which.min(bic.mcp2[,1]+bic.mcp2[,2])
      
     
    }
    
    # GROUP LASSO 
    # 
    thetahat.grlasso = trainSPMA.object$grlasso$beta
    ncolb = ncol(thetahat.grlasso)
    bic.grlasso = cbind(rep(0,ncolb),rep(0,ncolb))
    
    for(i in 1:ncolb){
      bb = bic_non_spma( N,thetahat = thetahat.grlasso[,i],sample1,sample2,p, M = 2, real = real,grlasso = TRUE)
      bic.grlasso[i,1:2] = c(bb$p1,bb$p2)
    }
    bic.grlasso = bic.grlasso[-1,]
    idx.grlasso = which.min(bic.grlasso[,1]+bic.grlasso[,2])+1
    #idx.grlasso = which.min(bic.grlasso[,1]+bic.grlasso[,2])
    
    rm(ncolb)
    # GROUP SCAD 
    thetahat.grscad = trainSPMA.object$grscad$beta
    ncolb = ncol(thetahat.grscad)
    bic.grscad = cbind(rep(0,ncolb),rep(0,ncolb))
    
    for(i in 1:ncolb){
      bb = bic_non_spma(N,thetahat = thetahat.grscad[,i],sample1,sample2,p, M = 2, real = real)
      bic.grscad[i,1:2] = c(bb$p1,bb$p2)
    }
    #bic.grscad = bic.grscad[-1,]
    #idx.grscad = which.min(bic.grscad[,1]+bic.grscad[,2])+1
    idx.grscad = which.min(bic.grscad[,1]+bic.grscad[,2])
    rm(ncolb)
    # GROUP MCP
    thetahat.grmcp = trainSPMA.object$grmcp$beta
    ncolb =ncol(thetahat.grmcp)
    bic.grmcp = cbind(rep(0,ncolb),rep(0,ncolb))
    
    for(i in 1:ncolb){
      bb = bic_non_spma(N,thetahat = thetahat.grmcp[,i],sample1,sample2,p, M = 2, real = real)
      bic.grmcp[i,1:2] = c(bb$p1,bb$p2)
    }
    #bic.grmcp = bic.grmcp[-1,]
    #idx.grmcp = which.min(bic.grmcp[,1]+bic.grmcp[,2])+1
    idx.grmcp = which.min(bic.grmcp[,1]+bic.grmcp[,2])
    
    return (list(idx.lasso1 = idx.lasso1, idx.lasso2 = idx.lasso2, 
                 idx.scad1 = idx.scad1, idx.scad2 = idx.scad2, 
                 idx.mcp1 = idx.mcp1, idx.mcp2 = idx.mcp2, 
                 idx.grmcp = idx.grmcp, idx.grlasso = idx.grlasso, idx.grscad = idx.grscad,
                 bic.lasso1 = bic.lasso1, bic.lasso2 = bic.lasso2, 
                 bic.grlasso = bic.grlasso,bic.grscad = bic.grscad, bic.grmcp = bic.grmcp))
    
    # return (list(
    #              idx.grmcp = idx.grmcp, idx.grlasso = idx.grlasso, idx.grscad = idx.grscad,
    #              bic.grlasso = bic.grlasso,bic.grscad = bic.grscad, bic.grmcp = bic.grmcp))
  }
}



select.model.spma <- function(trainSPMA.object,real =FALSE){
  n = length(trainSPMA.object$spma)
  bic.spma = cbind(rep(0,n), rep(0,n))
  for(i in 1:n){
    #  print(i);
    # bayes_bic_spma(thetaF, thetahat, covF1, covF2, covF3, set.n, p,real)
    # bb = bayes_bic_spma(trainSPMA.object$unpel,trainSPMA.object$spma[[i]]$thetaM[,,trainSPMA.object$spma[[i]]$bbk],
    #                     trainSPMA.object$covF1,trainSPMA.object$covF2,trainSPMA.object$covF3,trainSPMA.object$N)
    bb = bayes_bic_spma(trainSPMA.object$unpel,trainSPMA.object$spma[[i]]$alphaM[,,trainSPMA.object$spma[[i]]$bbk],
                        trainSPMA.object$covF1,trainSPMA.object$covF2,trainSPMA.object$N,
                        p= trainSPMA.object$p, real=real)
    bic.spma[i,1:2] = c(bb$p1, bb$p2)
  }
  idx.spma = which.min(bic.spma[,1]+bic.spma[,2])
  return(list(idx.spma = idx.spma, bic.spma = bic.spma))
}

bayes_bic_spma <- function(thetaF, thetahat, covF1, covF2, set.n, p,real ){
  if(nrow(thetahat)!=nrow(thetaF)){print("Parameter dimension should be the same"); break; }
  #print(dim( as.matrix(thetahat[1,] - thetaF[1,], nrow =1,ncol = p)))
  #print(dim(covF1))
  # p1 = t(as.matrix(thetahat[1,] - thetaF[1,], nrow =1,ncol = p))%*% covF1 %*% as.matrix(thetahat[1,] - thetaF[1,])
  # p1 = p1+t(as.matrix(thetahat[2,] - thetaF[2,], nrow =1,ncol = p))%*% covF2 %*% as.matrix(thetahat[2,] - thetaF[2,])
  # p1 = p1+t(as.matrix(thetahat[3,] - thetaF[3,], nrow =1,ncol = p))%*% covF3 %*% as.matrix(thetahat[3,] - thetaF[3,])
  # test for real
  
  if(real){ #%*% covF1 %*% 
    p1 = t(as.matrix(thetahat[1,] - thetaF[1,], nrow =1,ncol = p)) %*% covF1 %*% as.matrix(thetahat[1,] - thetaF[1,])#*(set.n[1])
    p1 = p1+t(as.matrix(thetahat[2,] - thetaF[2,], nrow =1,ncol = p))%*% covF2 %*% as.matrix(thetahat[2,] - thetaF[2,])#*(set.n[2])
    
    
    p2 = log(sum(set.n))/(sum(set.n))*(sum(thetahat!=0))
    #p2 = log(sum(set.n))*(sum(thetahat!=0))
    #p2 = (log(set.n[1])/set.n[1] + log(set.n[2])/set.n[2] + log(set.n[3])/set.n[3])*(sum(thetahat!=0))
  } # 
  else{
    p1 = t(as.matrix(thetahat[1,] - thetaF[1,], nrow =1,ncol = p))%*% covF1 %*% as.matrix(thetahat[1,] - thetaF[1,])*(set.n[1])
    p1 = p1+t(as.matrix(thetahat[2,] - thetaF[2,], nrow =1,ncol = p))%*% covF2 %*% as.matrix(thetahat[2,] - thetaF[2,])*(set.n[2])
   
    
    p2 = log(sum(set.n))*sum(thetahat!=0)}
  #p2 = log(sum(set.n))/sum(set.n)*sum(thetahat!=0)} # decide the form of n_m or n
  # if all thethat are 0, p2 is very small comparing to p1, so p2 cant use log(sum(set.n)) 
  #p1 = p1*log(sum(set.n))
  return (list(p1 = p1,p2=p2))
}



bic_lasso <- function(ni, thetahat, sample, p, M = 2){
  # first component is still 500*1 by setting nrow =1, to transpose
  p1 = t(as.matrix(sample[,p+2], nrow = 1,ncol = ni))%*%
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

bic_non_spma <- function(set.n, thetahat, sample1, sample2, p, M = 2,intervalCensor = FALSE,real = FALSE, grlasso = FALSE){
  
  # first component is still 500*1 by setting nrow =1, to transpose
  p1 = t(as.matrix(sample1[,p+2], nrow = 1,ncol = set.n[1]))%*%(as.matrix(sample1[,1:p])%*%as.matrix(thetahat[1:p], ncol =1))
  p1 = p1+ t(as.matrix(sample2[,p+2], nrow = 1, ncol = set.n[2]))%*%(as.matrix(sample2[,1:p])%*%as.matrix(thetahat[(p+1):(2*p)],ncol =1))
 
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
  

  
  if(real){
    if(grlasso){ p1 = (-2)*p1
    #p2 = log(sum(set.n))*log(sum(thetahat!=0))
    p2 = log(sum(set.n))*sum(thetahat!=0)
    #p2 = log(sum(set.n))/sum(set.n)*(sum(thetahat!=0)) 
    }
    else{
    p1 = (-2)*p1
    #p2 = log(sum(set.n))*log(sum(thetahat!=0))
    p2 = log(sum(set.n))*(sum(thetahat!=0))
    }
  }
  else{
    p1 = (-2)*p1
    #p2 = log(sum(set.n))/sum(set.n)*log(sum(thetahat!=0))
    p2 = log(sum(set.n)) *(sum(thetahat!=0)) # ? 
  }
  return (list(p1 = p1,p2=p2))
  #return(p1+p2)
  
  
}

