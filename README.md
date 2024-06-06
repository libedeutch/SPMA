# Sparse and heterogeneous meta-analysis with semiparametric regression models
Perform meta-analysis with semiparametric regression models based on summary-level statistics. 
## Usage

```R

data = generate_data(N = c(200,200,200), p= 10, seed = 2024,rightCensor = TRUE, correlate = TRUE, intervalCensor =FALSE

# Model fitting
re = trainSPMA(data,iter = 1000,rightCensor = TRUE, intervalCensor = FALSE,m = 5)

# Model selection
best.idx = bic_model_selection(re)
spma.bic = select.model.spma(re)

# Extract best estimates 
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

# Report MRME, TRR, FPR, ME
result = list()
rmseO = rmse(trainSPMA.object=re2, best.idx = best.idx, spma.bic = spma.bic, Data = tt)
result[[1]] = rmseO
printResult(result, intervalCensor = FALSE)
```
More function descriptions and examples are provided in [`\man`](https://github.com/libedeutch/SPMA/tree/main/man)
