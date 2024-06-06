# time-dependent AUC 

library(Hmisc) 
library(timeROC)

# Data preparation 
#gene_id = intersect(g2,g3)

load("~/metaADMM/real/38genes.RData")
X1n = study1[,gene_id]
X2n = study2[,gene_id]
#X3n = study3[,gene_id]

#X1.scaled = scale(X1n)
#X2.scaled = scale(X2n)
#X3.scaled = scale(X3n)

y1<- study1[,c('event','status')]
colnames(y1) <- c("time","status")
y2<- study2[,c('event','status')]
colnames(y2) <- c("time","status")
#y3<- study3[,c('event','status')]
#colnames(y3) <- c("time","status")


# sample1.s = cbind(scale(X1n),y1)
# sample2.s = cbind(scale(X2n),y2)
# sample3.s = cbind(scale(X3n),y3)

sample1.s = cbind(X1n,y1)
sample2.s = cbind(X2n,y2)
#sample3.s = cbind(X3n,y3)

sample1o = cbind(X1n,y1)
sample2o = cbind(X2n,y2)
#sample3o = cbind(X3n,y3)

iter = 200
auc <- matrix(nrow = 8*iter , ncol =5  )
auc_time <-  matrix(nrow = iter , ncol =4 )
modelSize <- matrix(nrow = iter, ncol = 8)
tt = 1


for (i in 1:iter){
flag = TRUE  
while(flag){   
  print(i)
  N <- c(nrow(study1), nrow(study2), nrow(study3))
  N.train = floor(N*0.8)
  train.idx1 = sample(1:N[1], N.train[1],replace = FALSE )
  train_sample1 =  sample1.s[train.idx1,]
  test_sample1 = sample1o[-train.idx1,]
  train.idx2 = sample(1:N[2], N.train[2],replace = FALSE )
  train_sample2 =  sample2.s[train.idx2,]
  test_sample2 = sample2o[-train.idx2,]
  # train.idx3 = sample(1:N[3], N.train[3],replace = FALSE )
  # train_sample3 =  sample3.s[train.idx3,]
  # test_sample3 = sample3o[-train.idx3,]  
  
  exTx1 = t(train_sample1[,1:(ncol(train_sample1)-2)])%*%train_sample1[,1:(ncol(train_sample1)-2)]/N[1]
  exTx2 = t(train_sample2[,1:(ncol(train_sample2)-2)])%*%train_sample2[,1:(ncol(train_sample2)-2)]/N[2]
  #exTx3 = t(train_sample3[,1:(ncol(train_sample3)-2)])%*%train_sample3[,1:(ncol(train_sample3)-2)]/N[3]

  # model fitting   
  realData <- list(
    sample1 = train_sample1, sample2 = train_sample2,
    p=length(gene_id), M=2, N = N.train[c(1,2)], rightCensor = TRUE, intervalCensor = FALSE, exTx1 = exTx1, exTx2=exTx2
  )
  re2 = SPMAreal2(realData, rightCensor = TRUE)
# Model selection by Modified BIC 
  two_study <- TRUE
  if(two_study){ 
    best.idx = bic_model_selection(re2, real = TRUE)
    spma.bic = select.model.spma(re2,real = TRUE)
    
    unpel.beta = re2$unpel
    lasso1.beta = re2$lasso1$beta[,best.idx$idx.lasso1]
    lasso2.beta = re2$lasso2$beta[,best.idx$idx.lasso2]
    
    
    lasso.beta = rbind(lasso1.beta,lasso2.beta)
    
    scad1.beta = re2$scad1$beta[,best.idx$idx.scad1]
    scad2.beta = re2$scad2$beta[,best.idx$idx.scad2]
    
    
    scad.beta = rbind(scad1.beta,scad2.beta)
    
    
    mcp1.beta = re2$mcp1$beta[,best.idx$idx.mcp1]
    mcp2.beta = re2$mcp2$beta[,best.idx$idx.mcp2]
    
    
    mcp.beta = rbind(mcp1.beta,mcp2.beta)
    
    
    grlasso.beta = matrix(re2$grlasso$beta[,best.idx$idx.grlasso], nrow = 2, byrow = TRUE)
    grscad.beta = matrix(re2$grscad$beta[,best.idx$idx.grscad], nrow = 2, byrow = TRUE)
    grmcp.beta = matrix(re2$grmcp$beta[,best.idx$idx.grmcp], nrow = 2, byrow = TRUE)
    
    spma.beta = re2$spma[[spma.bic$idx.spma]]$alphaM[,,re2$spma[[spma.bic$idx.spma]]$bbk]
    
    colnames(unpel.beta) <- colnames(X1n)
    colnames(lasso.beta) <- colnames(X1n)
    colnames(scad.beta) <- colnames(X1n)
    colnames(mcp.beta) <- colnames(X1n)
    colnames(grlasso.beta) <- colnames(X1n)
    colnames(grscad.beta) <- colnames(X1n)
    colnames(grmcp.beta) <- colnames(X1n)
    colnames(spma.beta) <- colnames(X1n)
    
    modelSize[i,] = c(model_size(unpel.beta),model_size(lasso.beta), model_size(scad.beta), model_size(mcp.beta),
    model_size(grlasso.beta), model_size(grscad.beta),  model_size(grmcp.beta), model_size(spma.beta))
  }

  
  if (sum(unpel.beta!=0) & sum(lasso.beta!=0) & sum(scad.beta!=0) & sum(mcp.beta!=0) & sum(grlasso.beta!=0) & sum(grscad.beta!=0) & sum(grmcp.beta!=0) & sum(spma.beta!=0)){flag = FALSE}
  
} 
  
    
# Estimated hazadrd score   
estimate_score <- TRUE
if(estimate_score){
  unpel.pred.haz <- c(exp(test_sample1[,1:length(gene_id)]%*%as.matrix(unpel.beta[1,], nrow = length(gene_id))), 
                      exp(test_sample2[,1:length(gene_id)]%*%as.matrix(unpel.beta[2,], nrow = length(gene_id))))
  lasso.pred.haz <- c(exp(test_sample1[,1:length(gene_id)]%*%as.matrix(lasso.beta[1,], nrow = length(gene_id))), 
                      exp(test_sample2[,1:length(gene_id)]%*%as.matrix(lasso.beta[2,], nrow = length(gene_id))))
  
  scad.pred.haz <- c(exp(test_sample1[,1:length(gene_id)]%*%as.matrix(scad.beta[1,], nrow = length(gene_id))), 
                     exp(test_sample2[,1:length(gene_id)]%*%as.matrix(scad.beta[2,], nrow = length(gene_id))))
  mcp.pred.haz <- c(exp(test_sample1[,1:length(gene_id)]%*%as.matrix(mcp.beta[1,], nrow = length(gene_id))), 
                    exp(test_sample2[,1:length(gene_id)]%*%as.matrix(mcp.beta[2,], nrow = length(gene_id))))
  
  grlasso.pred.haz <- c(exp(test_sample1[,1:length(gene_id)]%*%as.matrix(grlasso.beta[1,], nrow = length(gene_id))), 
                        exp(test_sample2[,1:length(gene_id)]%*%as.matrix(grlasso.beta[2,], nrow = length(gene_id))))
  
  grscad.pred.haz <- c(exp(test_sample1[,1:length(gene_id)]%*%as.matrix(grscad.beta[1,], nrow = length(gene_id))), 
                       exp(test_sample2[,1:length(gene_id)]%*%as.matrix(grscad.beta[2,], nrow = length(gene_id))))
  
  grmcp.pred.haz <- c(exp(test_sample1[,1:length(gene_id)]%*%as.matrix(grmcp.beta[1,], nrow = length(gene_id))), 
                      exp(test_sample2[,1:length(gene_id)]%*%as.matrix(grmcp.beta[2,], nrow = length(gene_id))))
  
  spma.pred.haz <- c(exp(test_sample1[,1:length(gene_id)]%*%as.matrix(spma.beta[1,], nrow = length(gene_id))), 
                     exp(test_sample2[,1:length(gene_id)]%*%as.matrix(spma.beta[2,], nrow = length(gene_id))))
}

### Time dependent AUC 
T = c(as.vector(test_sample1[,'time']), as.vector(test_sample2[,'time']))
delta = c(as.vector(test_sample1[,'status']), as.vector(test_sample2[,'status']))
roc_cal <- TRUE
if(roc_cal){
  ROC.bili.cox.unpel1<-timeROC(T = T,
                               delta=delta,marker=unpel.pred.haz,
                               #other_markers=as.matrix(test_sample1[,2:38]),
                               cause=1,weighting="cox",
                               #times=quantile(T,probs=seq(0.2,0.8,0.1)))
                               times=c(6,12,18,24))
  ROC.bili.cox.lasso1<-timeROC(T = T,
                               delta=delta,marker=lasso.pred.haz,
                               #other_markers=as.matrix(test_sample1[,2:38]),
                               cause=1,weighting="cox",
                               #times=quantile(T,probs=seq(0.2,0.8,0.1)))
                               times=c(6,12,18,24))
  
  ROC.bili.cox.scad1<-timeROC(T = T,
                              delta=delta,marker=scad.pred.haz,
                              #other_markers=as.matrix(test_sample1[,2:38]),
                              cause=1,weighting="cox",
                              #times=quantile(T,probs=seq(0.2,0.8,0.1)))
                              times=c(6,12,18,24))
  
  ROC.bili.cox.mcp1<-timeROC(T = T,
                             delta=delta,marker=mcp.pred.haz,
                             #other_markers=as.matrix(test_sample1[,2:38]),
                             cause=1,weighting="cox",
                             #times=quantile(T,probs=seq(0.2,0.8,0.1)))
                             #times=c(0,10,20,50,100))
                             times=c(6,12,18,24))
  
  ROC.bili.cox.grlasso1<-timeROC(T = T,
                                 delta=delta,marker=grlasso.pred.haz,
                                 #other_markers=as.matrix(test_sample1[,2:38]),
                                 cause=1,weighting="cox",
                                 #times=quantile(T,probs=seq(0.2,0.8,0.1)))
                                 times=c(6,12,18,24))
  
  ROC.bili.cox.grscad1<-timeROC(T = T,
                                delta=delta,marker=grscad.pred.haz,
                                #other_markers=as.matrix(test_sample1[,2:38]),
                                cause=1,weighting="cox",
                                #times=c(0,10,20,50,100))
                                #times=quantile(T,probs=seq(0.2,0.8,0.1)))
                                times=c(6,12,18,24))
  
  ROC.bili.cox.grmcp1<-timeROC(T = T,
                               delta=delta,marker=grmcp.pred.haz,
                               #other_markers=as.matrix(test_sample1[,2:38]),
                               cause=1,weighting="cox",
                               #times=c(0,10,20,50,100))
                               #times=quantile(T,probs=seq(0.2,0.8,0.1)))
                               times=c(6,12,18,24))
  
  ROC.bili.cox.spma1<-timeROC(T = T,
                              delta=delta,marker=spma.pred.haz,
                              #other_markers=as.matrix(test_sample1[,2:38]),
                              cause=1,weighting="cox",
                              #times=c(0,10,20,50,100))
                              #times=quantile(T,probs=seq(0.2,0.8,0.1)))
                              times=c(6,12,18,24))
}
auc[tt:(tt+7), 2:5] <- rbind(ROC.bili.cox.unpel1$AUC, ROC.bili.cox.lasso1$AUC, ROC.bili.cox.scad1$AUC, ROC.bili.cox.mcp1$AUC,
             ROC.bili.cox.grlasso1$AUC, ROC.bili.cox.grscad1$AUC, ROC.bili.cox.grmcp1$AUC, ROC.bili.cox.spma1$AUC)
auc[tt:(tt+7),1] <- c("Unpenalized CoxPH", "Lasso", "SCAD", "MCP", "Group Lasso", "Group SCAD", "Group MCP", "SPMA")
#auc_time[tt, 1:8] <- ROC.bili.cox.unpel1$times
tt =tt + 8 
}

######### dataset for visualization ########## 
na_num <- apply(auc,1,function(x){ sum(is.na(x))})
ggdata=auc
ggdata <- auc[-which(na_num==5), ]
colnames(ggdata) <- c("Method", "6", "12", "18","24")
method <- rep(ggdata[,'Method'], 4)
metric <- as.vector(ggdata[,2:5])

time <- rep(colnames(ggdata)[2:5], each = nrow(ggdata))
ggd <- data.frame(cbind(method, metric,  time))
ggd$pp = pp = rep("AUC", nrow(ggd))
ggd$time <- factor(ggd$time, levels =c(6,12,18,24),
                   labels = c("6 Months","12 Months",
                              "18 Months","24 Months"))
ggd$method <- factor(ggd$method, levels = c("Unpenalized CoxPH", "Lasso", "SCAD", "MCP", "Group Lasso", "Group SCAD", "Group MCP", "MACS"))

save(ggd, file = "~/metaADMM/real/tdAUC.RData")

# library(ggplot2)
# g1 <- ggplot(ggd, aes( x = pp, y = as.numeric(metric)))+
#   #geom_violin(aes(fill = method))+
#   #geom_boxplot(width = 0.1)+
#   geom_boxplot(aes(fill = method),width = 1.3 )+ #"#D2AEAC","#797979","#7FBFBB",
#   facet_wrap(~time, ncol = 4)+ #"#76c893"
#   scale_fill_manual(values = c( "#FFCE54", "#f5ada9", "#AC92EB","#c19883","#EDE4D4","#A0D568", "#4FC1E8", "#ED5564"), name = "Method")+
#   theme(panel.grid.major = element_line(color = gray(.5), linetype = 'dashed', size = 0.1),
#         #panel.grid.major = element_blank(),
#         axis.title.x = element_blank(),
#         strip.text = element_text(size = 14),
#         panel.spacing.x = unit(1.5,"lines"),
#         legend.position = "bottom",
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         panel.background = element_rect(fill = "white"))+
#   
#   scale_y_continuous(breaks =c(0.5,0.6,0.7,0.8,0.9,1))+
#   guides(fill = guide_legend(nrow = 1))+
# xlab("")+ylab("Time-depedent AUC")
#   
# ggsave(file = "/Users/gilly/Library/CloudStorage/OneDrive-TheUniversityofTexasatDallas/metadata analysis/real/AML/time_depeddent_auc_43.pdf",
#        dpi = 300, width = 13, height = 6)
# 
# 
# colnames(modelSize) <- c("Unpenalized CoxPH", "Lasso", "SCAD", "MCP", "Group Lasso", "Group SCAD", "Group MCP", "SPMA")
# 
# boxplot(modelSize)
# 
# 
# par(mfrow = c(1,3))
# 
# plot(best.idx$bic.grlasso[,1], col ="red", ylim = c(12,2530), main = "grlasso")
# lines(best.idx$bic.grlasso[,2], col = "blue")
# lines(rowSums(best.idx$bic.grlasso), col = "green")
# 
# 
# plot(best.idx$bic.grscad[,1], col ="red", ylim = c(12,2530), main = "grscad")
# lines(best.idx$bic.grscad[,2], col = "blue")
# lines(rowSums(best.idx$bic.grscad), col = "green")
# 
# plot(best.idx$bic.grmcp[,1], col ="red", ylim = c(12,2530), main = "grmcp")
# lines(best.idx$bic.grmcp[,2], col = "blue")
# lines(rowSums(best.idx$bic.grmcp), col = "green")