# AML select gene 

load("~/OneDrive - The University of Texas at Dallas/metadata analysis/real/AML/processed_data_rpkm.RData")
all.equal(colnames(ohsu.2018.gene.rpkm)[-c(1,2)], ohsu.2018.sample$SAMPLE_ID)
#TRUE
all.equal(ohsu.2018.sample$PATIENT_ID,ohsu.2018.patient[,'PATIENT_ID']) # TRUE

all.equal(colnames(ohsu.2022.gene.rpkm)[-c(1)], ohsu.2022.sample$SAMPLE_ID)
#TRUE
all.equal(ohsu.2022.sample$PATIENT_ID,ohsu.2022.patient[,'PATIENT_ID']) # TRUE

all.equal(colnames(tcga.2018.gene.rpkm)[-c(1,2)], stringr::str_replace_all(tcga.2018.sample$SAMPLE_ID,"-","."))
#TRUE
all.equal(tcga.2018.sample$PATIENT_ID,tcga.2018.patient[,'PATIENT_ID']) # TRUE


# 351 x 13207
# OHSU 2018
study1 = cbind(t(ohsu.2018.gene.rpkm[,3:ncol(ohsu.2018.gene.rpkm)]), ohsu.2018.patient[,'OS_MONTHS'], ohsu.2018.patient[,'OS'])
colnames(study1)[c(13206,13207)] <- c("event", "status")
mode(study1) <- "numeric"

# OHSU 2022
# 165 x 13207
study2 = cbind(t(ohsu.2022.gene.rpkm[,2:ncol(ohsu.2022.gene.rpkm)]), ohsu.2022.patient[,'OS_MONTHS'], ohsu.2022.patient[,'OS'])
colnames(study2)[c(13206,13207)] <- c("event", "status")
mode(study2) <- "numeric"

# TCGA 2018
# 161 x 13205
study3 = cbind(t(tcga.2018.gene.rpkm[,3:ncol(tcga.2018.gene.rpkm)]),tcga.2018.patient[,'OS_MONTHS'], tcga.2018.patient[,'OS'])
colnames(study3)[c(13204,13205)] <- c("event", "status")
mode(study3) <- "numeric"

study3 <- study3[-which(is.na(study3[,'event'])), ]

# tcga.2018.os.month.na <- which(is.na(tcga.2018.patient$OS_MONTHS))
# tcga.2018.patient2 <- tcga.2018.patient[-tcga.2018.os.month.na,]
# tcga.2018.sample2 <- tcga.2018.sample[-tcga.2018.os.month.na,]
# tcga.2018.gene.rpkm2 <- tcga.2018.gene.rpkm[,-(tcga.2018.os.month.na+2)]
# 
# # TRUE
# all.equal(colnames(tcga.2018.gene.rpkm2)[-c(1,2)], stringr::str_replace_all(tcga.2018.sample2$SAMPLE_ID,"-","."))
# # TRUE
# all.equal(tcga.2018.patient2$PATIENT_ID, tcga.2018.sample2$PATIENT_ID)
          
          
ohsu_2018_hugo <- ohsu.2018.gene.rpkm$Hugo_Symbol
ohsu_2022_hugo <- ohsu.2022.gene.rpkm$Hugo_Symbol
tcga_2018_hugo <- tcga.2018.gene.rpkm$Hugo_Symbol

# step 1: variance; seems unreasonable
var_ohsu_2018 <- apply(study1[,1:(ncol(study1)-2)],2,var)
var_ohsu_2022 <- apply(study2[,1:(ncol(study2)-2)],2,var)
var_tcga_2018 <- apply(study3[,1:(ncol(study3)-2)],2,var)

# variance > 75% 
gene_ohsu_2018 <- ohsu_2018_hugo[which(var_ohsu_2018>summary(var_ohsu_2018)[5])] # 3301
gene_ohsu_2022 <- ohsu_2022_hugo[which(var_ohsu_2022>summary(var_ohsu_2022)[5])] # 3301
gene_tcga_2018 <- tcga_2018_hugo[which(var_tcga_2018>summary(var_tcga_2018)[5])] # 3301

ohsu_2018_idx <- which(var_ohsu_2018>summary(var_ohsu_2018)[5])
ohsu_2022_idx <- which(var_ohsu_2022>summary(var_ohsu_2022)[5])

univariate_selection = TRUE
#pval_z = matrix(nrow = ncol(study1)-2, ncol = 3)
pval_wald_ohsu = matrix(nrow = ncol(study1)-2, ncol = 2)
pval_score_ohsu = matrix(nrow = ncol(study1)-2, ncol = 2)
pval_wald_tcga = matrix(nrow = ncol(study3)-2, ncol = 1)
pval_score_tcga = matrix(nrow = ncol(study3)-2, ncol = 1)

set.seed(1000)
if(univariate_selection){
  library(survival)
  for (g in 1:(ncol(study1)-2)){
    data  = data.frame(cbind(study1[,g],study1[,c('event','status')]))
    m1 = coxph(Surv(event,status)~.,data)
    pval_score_ohsu[g,1] = summary(m1)$sctest[3]
    #pval_z[g,1] = summary(m1)$coefficients[5]
    pval_wald_ohsu[g,1] = summary(m1)$waldtest[3]
    #pval_like[g,1] = summary(m1)$logtest[3]
    
    data  = data.frame(cbind(study2[,g],study2[,c('event','status')]))
    m2 = coxph(Surv(event,status)~.,data)
    pval_score_ohsu[g,2] = summary(m2)$sctest[3]
    #pval_z[g,2] = summary(m2)$coefficients[5]
    pval_wald_ohsu[g,2] = summary(m2)$waldtest[3]
    #pval_like[g,2] = summary(m2)$logtest[3]
    
  }   
  for (g in 1:(ncol(study3)-2)){   
    data  = data.frame(cbind(study3[,g],study3[,c('event','status')]))
    m3 = coxph(Surv(event,status)~.,data)
    pval_score_tcga[g,1] = summary(m3)$sctest[3]
    #pval_z[g,3] = summary(m3)$coefficients[5]
    pval_wald_tcga[g,1] = summary(m3)$waldtest[3]
    #pval_like[g,3] = summary(m3)$logtest[3]
  }
  
  
}

pval_score_ohsu_sub_2018 <- pval_score_ohsu[ohsu_2018_idx,1]
pval_score_ohsu_sub_2022 <- pval_score_ohsu[ohsu_2022_idx,2]

pval_score_ohsu_sub_2018_adj <- p.adjust(pval_score_ohsu_sub_2018,method = "fdr")
pval_score_ohsu_sub_2022_adj <- p.adjust(pval_score_ohsu_sub_2022,method = "fdr")

pval_wald_ohsu2 <- apply(pval_wald_ohsu,2,function(x) p.adjust(x,method = "fdr"))
pval_score_ohsu2 <- apply(pval_score_ohsu,2,function(x) p.adjust(x,method = "fdr"))
pval_wald_tcga2 <- apply(pval_wald_tcga,2,function(x) p.adjust(x,method = "fdr"))
pval_score_tcga2 <- apply(pval_score_tcga,2,function(x) p.adjust(x,method = "fdr"))

apply(pval_wald_ohsu2,2,function(x) sum(x<0.1))
apply(pval_score_ohsu2,2,function(x) sum(x<0.1))

sum(pval_wald_tcga<0.1)
sum(pval_score_tcga<0.1)

apply(pval_wald_ohsu,2,function(x) sum(x<0.01))
apply(pval_score_ohsu,2,function(x) sum(x<0.01))

# three studies 
gene_ohsu_2018_2 <-ohsu_2018_hugo[which(pval_score_ohsu[,1]<0.01)]
gene_ohsu_2022_2 <-ohsu_2022_hugo[which(pval_score_ohsu[,2]<0.01)]
gene_tcga_2018_2 <-tcga_2018_hugo[which(pval_score_tcga[,1]<0.01)]

# gene_ohsu_2018_2 <-ohsu_2018_hugo[order(pval_score_ohsu[,1], decreasing = TRUE)[1:20]]
# gene_ohsu_2022_2 <-ohsu_2022_hugo[order(pval_score_ohsu[,2], decreasing = TRUE)[1:20]]
# gene_tcga_2018_2 <-tcga_2018_hugo[order(pval_score_tcga[,1], decreasing = TRUE)[1:20]]
# gene_id <- union(gene_ohsu_2018_2, gene_ohsu_2022_2)
# 
# gene_ohsu_2018_2 <-ohsu_2018_hugo[which(pval_score_ohsu2[,1]<0.01)]
# gene_ohsu_2022_2 <-ohsu_2022_hugo[which(pval_score_ohsu2[,2]<0.2)]
# gene_tcga_2018_2 <-tcga_2018_hugo[which(pval_score_tcga2[,1]<0.01)]
# 
# sum(gene_ohsu_2018_2%in%gene_tcga_2018_2)
# sum(gene_ohsu_2018_2%in%gene_ohsu_2018)
# sum(gene_tcga_2018_2%in%gene_tcga_2018)
# sum(gene_ohsu_2022_2%in%gene_ohsu_2022)
# sum(gene_ohsu_2022_2%in%gene_ohsu_2018_2)
# sum(gene_ohsu_2022_2%in%gene_tcga_2018_2)
#pval_like = matrix(nrow = p, ncol = 3)

gene_id <- union(gene_ohsu_2018_2, union(gene_ohsu_2022_2, gene_tcga_2018_2))

# g1 <- intersect(gene_tcga_2018_2,intersect(gene_tcga_2018 ,intersect(gene_ohsu_2018_2,gene_ohsu_2018)))
# g2 <- intersect(intersect(gene_ohsu_2022_2, gene_ohsu_2022),intersect(gene_ohsu_2018, gene_ohsu_2018_2))

g1 <- intersect(gene_tcga_2018_2,gene_tcga_2018)
g3 <- intersect(gene_ohsu_2018_2,gene_ohsu_2018)
g2 <- intersect(gene_ohsu_2022_2,gene_ohsu_2022)

length(intersect(g1,g2)) 
length(intersect(g1,g3))
length(intersect(g2,g3))

# top 10 significant, to narrow down the candidate genes, cant handle the case when 
# there are many candicate genes 

colnames(study1)[1:(ncol(study1)-2)] <- ohsu.2018.gene.rpkm$Hugo_Symbol
colnames(study2)[1:(ncol(study2)-2)] <- ohsu.2022.gene.rpkm$Hugo_Symbol
colnames(study3)[1:(ncol(study3)-2)] <- tcga.2018.gene.rpkm$Hugo_Symbol

#gene_id = union(intersect(g3,g2),union(intersect(g1,g2),intersect(g1,g3)))


# meta-analysis of OHSU 2018 and OHSU 2022
gene_id = intersect(g2,g3)

X1n = study1[,gene_id]
X2n = study2[,gene_id]
X3n = study3[,gene_id]


#X1.scaled = scale(X1n)
#X2.scaled = scale(X2n)
#X3.scaled = scale(X3n)

y1<- study1[,c('event','status')]
colnames(y1) <- c("time","status")
y2<- study2[,c('event','status')]
colnames(y2) <- c("time","status")
y3<- study3[,c('event','status')]
colnames(y3) <- c("time","status")


# sample1.s = cbind(scale(X1n),y1)
# sample2.s = cbind(scale(X2n),y2)
# sample3.s = cbind(scale(X3n),y3)

sample1.s = cbind(X1n,y1)
sample2.s = cbind(X2n,y2)
sample3.s = cbind(X3n,y3)

sample1o = cbind(X1n,y1)
sample2o = cbind(X2n,y2)
sample3o = cbind(X3n,y3)

#p = ncol(X1n)
N <- c(nrow(study1), nrow(study2), nrow(study3))
N.train = floor(N*1)
train.idx1 = sample(1:N[1], N.train[1],replace = FALSE )
train_sample1 =  sample1.s[train.idx1,]
test_sample1 = sample1o[-train.idx1,]
train.idx2 = sample(1:N[2], N.train[2],replace = FALSE )
train_sample2 =  sample2.s[train.idx2,]
test_sample2 = sample2o[-train.idx2,]
train.idx3 = sample(1:N[3], N.train[3],replace = FALSE )
train_sample3 =  sample3.s[train.idx3,]
test_sample3 = sample3o[-train.idx3,]

exTx1 = t(train_sample1[,1:(ncol(train_sample1)-2)])%*%train_sample1[,1:(ncol(train_sample1)-2)]/N[1]
exTx2 = t(train_sample2[,1:(ncol(train_sample2)-2)])%*%train_sample2[,1:(ncol(train_sample2)-2)]/N[2]
exTx3 = t(train_sample3[,1:(ncol(train_sample3)-2)])%*%train_sample3[,1:(ncol(train_sample3)-2)]/N[3]

# X1 = sample1[,1:p]
# X2 = sample2[,1:p]
# X3 = sample3[,1:p]
# train model

realData <- list(
  sample1 = train_sample1, sample2 = train_sample2,
  p=length(gene_id), M=2, N = N.train[c(1,2)], rightCensor = TRUE, intervalCensor = FALSE, exTx1 = exTx1, exTx2=exTx2
)

realData <- list(
  sample1 = train_sample1, sample2 = train_sample2, sample3 = train_sample3,
  p=length(gene_id), M=3, N = N.train, rightCensor = TRUE, intervalCensor = FALSE, exTx1 = exTx1, exTx2 = exTx2, exTx3=exTx3
)

#tt = generate_data(p=10, correlate=TRUE,rightCensor = TRUE,seed = 1)



realData <- list(
  sample1 = train_sample1, sample2 = train_sample3,
  p=length(gene.id), M=2, N = N.train[c(1,3)], rightCensor = TRUE, intervalCensor = FALSE, exTx1 = exTx1, exTx2 = exTx3
)

#tt = generate_data(p=10, correlate=TRUE,rightCensor = TRUE,seed = 1)
re2 = SPMAreal2(realData, rightCensor = TRUE)
re= trainSPMA(realData,rightCensor = T)



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
}


three_study<- FALSE
if(three_study){ 
  best.idx = bic_model_selection(re, real = TRUE)
  spma.bic = select.model.spma(re,real = TRUE)
  
  unpel.beta = re$unpel
  lasso1.beta = re$lasso1$beta[,best.idx$idx.lasso1]
  lasso2.beta = re$lasso2$beta[,best.idx$idx.lasso2]
  lasso3.beta = re$lasso3$beta[,best.idx$idx.lasso3]
  
  lasso.beta = rbind(lasso1.beta,lasso2.beta,lasso3.beta)
  
  scad1.beta = re$scad1$beta[,best.idx$idx.scad1]
  scad2.beta = re$scad2$beta[,best.idx$idx.scad2]
  scad3.beta = re$scad3$beta[,best.idx$idx.scad3]
  
  scad.beta = rbind(scad1.beta,scad2.beta,scad3.beta)
  
  
  mcp1.beta = re$mcp1$beta[,best.idx$idx.mcp1]
  mcp2.beta = re$mcp2$beta[,best.idx$idx.mcp2]
  mcp3.beta = re$mcp3$beta[,best.idx$idx.mcp3]
  
  mcp.beta = rbind(mcp1.beta,mcp2.beta,mcp3.beta)
  
  
  grlasso.beta = matrix(re$grlasso$beta[,best.idx$idx.grlasso], nrow = 3, byrow = TRUE)
  grscad.beta = matrix(re$grscad$beta[,best.idx$idx.grscad], nrow = 3, byrow = TRUE)
  grmcp.beta = matrix(re$grmcp$beta[,best.idx$idx.grmcp], nrow = 3, byrow = TRUE)
  
  spma.beta = re$spma[[spma.bic$idx.spma]]$alphaM[,,re$spma[[spma.bic$idx.spma]]$bbk]
  
  colnames(lasso.beta) <- colnames(X1n)
  colnames(scad.beta) <- colnames(X1n)
  colnames(mcp.beta) <- colnames(X1n)
  colnames(grlasso.beta) <- colnames(X1n)
  colnames(grscad.beta) <- colnames(X1n)
  colnames(grmcp.beta) <- colnames(X1n)
  colnames(spma.beta) <- colnames(X1n)
  colnames(unpel.beta) <- colnames(X1n)
}

############## TRAIN TEST ###############
# Somer D 


# ROC curve 
# t = 10

# t =30
# t = 50
# t = 100


ff1 = cbind(t(unpel.beta)[,1], t(lasso.beta)[,1], t(scad.beta)[,1], t(mcp.beta)[,1], t(grlasso.beta)[,1], t(grscad.beta)[,1], t(grmcp.beta)[,1],t(spma.beta)[,1])
ff2 = cbind(t(unpel.beta)[,2], t(lasso.beta)[,2], t(scad.beta)[,2], t(mcp.beta)[,2], t(grlasso.beta)[,2], t(grscad.beta)[,2], t(grmcp.beta)[,2], t(spma.beta)[,2])
ff3 = cbind(t(unpel.beta)[,3], t(lasso.beta)[,3], t(scad.beta)[,3], t(mcp.beta)[,3], t(grlasso.beta)[,3], t(grscad.beta)[,3], t(grmcp.beta)[,3], t(spma.beta)[,3])

colnames(ff1) <- c("unpel", "lasos", "scad", "mcp", "grlasso", "grscad", "grmcp", "spma")
colnames(ff2) <- c("unpel", "lasos", "scad", "mcp", "grlasso", "grscad", "grmcp", "spma")
colnames(ff3) <- c("unpel", "lasos", "scad", "mcp", "grlasso", "grscad", "grmcp", "spma")

#ff = matrix(nrow = nrow(ff1)*3, ncol = ncol(ff1)+1)
#ff[seq(1,nrow(ff)-2, 3),2:9 ] = ff1
#ff[seq(2,nrow(ff)-1,3),2:9 ] = ff2
#ff[seq(3,nrow(ff),3),2:9 ] = ff3

ff = matrix(nrow = nrow(ff1)*2, ncol = ncol(ff1)+1)
ff[seq(1,nrow(ff)-1, 2),2:9 ] = ff1
ff[seq(2,nrow(ff),2),2:9 ] = ff2


ff[,1] <- c(rep(colnames(spma.beta),each = 2))
ff[,2:9] <- round(as.numeric(ff[,2:9]),3)

tmp = matrix(nrow = nrow(ff1)*2, ncol = 1)
tmp[seq(1,nrow(tmp)-1, 2),1 ] = round(summary(re2$m1)$coefficient[,5],3)
tmp[seq(2,nrow(tmp),2),1] = round(summary(re2$m2)$coefficient[,5],3)

colnames(ff)<- c("gene","unpel", "lasso", "scad", "mcp", "grlasso", "grscad", "grmcp", "spma")


ff[,2] = paste0(ff[,2], " (",tmp, ")")


print(xtable::xtable(ff, type = "latex"), file = "filename2.tex")


model_size <- function(x){
  return(sum(x!=0))
}

model_size(unpel.beta)
model_size(lasso.beta)
model_size(scad.beta)
model_size(mcp.beta)
model_size(grlasso.beta)
model_size(grscad.beta)
model_size(grmcp.beta)
model_size(spma.beta)
