

# Ganymede Code for Meta Analysis
# Especially for interval censored data, p = 8/15  

library(optparse)
library(Rcpp)
option_list<- list(
  make_option("--rep", default =  500, help = "Number of replicate"),
  make_option("--p", default = 15, help = "Number of covariates"),
  make_option("--n1", default = 200, help = "sample size 1"),
  make_option("--n2", default = 200, help = "sample size 2"),
  make_option("--n3", default = 200, help = "sample size 3"),
  make_option("--rc", default = FALSE, help = "right censor"),
  make_option("--correlate", default = TRUE, help = "correlate covariates"),
  make_option("--ic", default = TRUE, help = "correlate covariates"),
  make_option("--n100",default=1,help = "Round of 100"),
make_option("--idxend",default=1,help = "idx of end"),
make_option("--idx",default=1,help = "idx of replicate to start")
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)


print(opt)
source("~/metaADMM/utilities_test_rcSep12_add_yuxiang_cook.R")
Rcpp::sourceCpp("~/metaADMM/spmaC_ADMM.cpp")
#Rcpp::sourceCpp("/Users/gilly/Library/CloudStorage/OneDrive-TheUniversityofTexasatDallas/metadata analysis/code/spmaC.cpp")
seed.save = c()
result = list()
seed = sample(1:1000000, opt$rep)
seed.save = c(seed.save, seed)
#write.table(seed.save, paste0("metaADMM/One1/seed_",opt$n1,"_",opt$n2,"_",opt$n3,"_ic_",opt$ic,"_cor_",opt$correlate,"_p_",opt$p,"_seed",".txt"))
seed = read.table(paste0("metaADMM/One1/seed_",opt$n1,"_",opt$n2,"_",opt$n3,"_ic_",opt$ic,"_cor_",opt$correlate,"_p_",opt$p,"_seed",".txt"))
# OnebyOne5 is for testing q = 0.2, 10e-3, 10e-2, rho = lambda, alpha = 1.5
# The performance is not good with 18 replicates
# try ceil(p*q) to generate more nonzero coefficients
#seed = as.vector(seed)
for(i in opt$idx:opt$idxend){
 print(i)
 print(seed[i,1]) # 222T needs seed[i,2]
  tt = generate_data(N = c(opt$n1,opt$n2,opt$n3), p= opt$p, seed = seed[i,1],rightCensor = FALSE, correlate = opt$correlate, intervalCensor =opt$ic)
  #tt = generate_data(N = c(opt$n1,opt$n2,opt$n3), p= opt$p, seed = seed[i],rightCensor = FALSE, correlate = opt$correlate, intervalCensor =opt$ic)
  re2 = trainSPMA(tt,iter = 1000,rightCensor = FALSE, intervalCensor = opt$ic,m = 5)
  best.idx = bic_model_selection(re2)
  spma.bic = select.model.spma(re2)
  print(paste0(best.idx$idx.grscad, " ", best.idx$idx.grmcp ))
  rmseO = rmse(trainSPMA.object=re2, best.idx = best.idx, spma.bic = spma.bic, Data = tt)
  result[[1]] = rmseO
  saveRDS(result, file = paste0("~/metaADMM/One1/sample_",opt$n1,"_",opt$n2,"_",opt$n3,"_ic_",opt$ic,"_cor_",opt$correlate,"_p_",opt$p,"_", i,".rds"))
}
