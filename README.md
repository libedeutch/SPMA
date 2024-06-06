# Sparse and heterogeneous meta-analysis with semiparametric regression models
Perform meta-analysis with semiparametric regression models based on summary-level statistics. 
## Usage

```R

tt = generate_data(N = c(opt$n1,opt$n2,opt$n3), p= opt$p, seed = seed[i,1],rightCensor = FALSE, correlate = opt$correlate, intervalCensor =opt$ic)
re2 = trainSPMA(tt,iter = 1000,rightCensor = FALSE, intervalCensor = opt$ic,m = 5)
best.idx = bic_model_selection(re2)
spma.bic = select.model.spma(re2)
print(paste0(best.idx$idx.grscad, " ", best.idx$idx.grmcp ))
rmseO = rmse(trainSPMA.object=re2, best.idx = best.idx, spma.bic = spma.bic, Data = tt)
result[[1]] = rmseO
