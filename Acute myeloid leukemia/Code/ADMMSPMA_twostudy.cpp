#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
NumericVector arrayC(NumericVector input, IntegerVector dim);
double normf(NumericVector x);
//arma::vec spma(NumericVector lambda, int iter, NumericMatrix covF1, NumericMatrix covF2, NumericMatrix covF3, int p, int M, 
                 //               IntegerVector setN, NumericMatrix thetaF);
// [[Rcpp::export]]
Rcpp::List spmaC(double alpha, double epi_abs, double epi_rel, double rho, double lambda, int iter, NumericMatrix covF1, NumericMatrix covF2, arma::mat augIMat, int p, int M, IntegerVector setN, arma::mat thetaF){
  arma::mat augNMat(M*p, M*p);
  // NumericMatrix augIMat(M*p, M*p);
  NumericVector epi_pri_store(iter);
  NumericVector epi_dual_store(iter);
  NumericVector pri_res_store(iter);
  NumericVector dual_res_store(iter);
  
  NumericMatrix epi_store(iter,4);
  
  arma::mat UnityMatrix(M*p,M*p);
  arma::mat denom(M*p,M*p);
  arma::mat theta_s(2,p);
  theta_s = thetaF; 
  arma::mat thetaM(2,p*(iter+1));
  arma::mat alphaM(2,p*(iter+1));
  arma::mat wM(2,p*(iter+1));
  double epi_dual = 0.00; 
  int bbk=iter;
  arma::mat thetaFvec;
  thetaFvec = thetaF.as_row();
  // Rcout<< thetaF << "**********";
  //  Rcout<<thetaFvec << "&&&&&&&&&" ;
  // thetaM is array(M,p,iter+1) -> M*(p*(iter+1))
  // initialize thetaM 
  for(int rr =0;rr<2;rr++){
    for (int cc=0; cc<p*1;cc++){
      thetaM(rr,cc) = theta_s(rr,cc);
    }
  }
  
  //initialize alphaM and wm as 0 
  for(int rr=0;rr<2;rr++){
    for(int cc=0;cc<p;cc++){
      alphaM(rr,cc) = 0;
      wM(rr,cc) =0;
    }
  }
  
  for(int i=0;i<M*p;i++){
    UnityMatrix(i,i) = rho;
  }
  
  // Initialization
  //Rcout<< UnityMatrix;
  for(int i=0; i<M;i++ ){
    for(int j =0; j<p;j++){
      augNMat(j+(p*i),j+(p*i)) = setN(i);
    }
  }
  
  denom = augNMat*augIMat; //? 
    
    //Rcout<<augNMat.n_cols<<"^" << augNMat.n_rows <<"*"<< augIMat.n_rows<< "&"<<augIMat.n_cols<< "kk"<< denom.n_cols<<"oo"<< denom.n_rows;
  //Rcout<<denom;
  denom = denom*2; 
  arma::mat fullTmp = denom*(thetaFvec.t());
  denom = denom + UnityMatrix;
  denom = inv(denom); // 
    
    // estimation
  for(int it=1; it<iter; it++){
    // every iteration 
    // first update theta
    arma::mat flatten_alphaM = alphaM.submat(0,p*(it-1),1,(p*it)-1).as_row();
    arma::mat flatten_wM = wM.submat(0,p*(it-1),1,(p*it)-1).as_row();
    // Rcout<< flatten_alphaM << "********"<<flatten_alphaM.n_cols; 
    // Rcout<< alphaM <<"^^^^^^^^^^"; 
    arma::mat tmp  = denom*(fullTmp + rho*flatten_alphaM.t() - rho*flatten_wM.t());
    // Rcout<< tmp<<"****"; 
    tmp = reshape(tmp, p, M); // reshape(nrow, ncols)
    // dimension of tmp is M*p x 1 
    // Rcout<< tmp.n_rows<< "****";
    thetaM.submat(0,p*it,1,(p*(it+1))-1) = tmp.t() ;
    // Rcout<<thetaM<<"*******";
    // update alpha 
    
    for(int jj=0; jj<p;jj++){
      double norm_thetas_j = 0.00; 
      for(int kk=0; kk<M; kk++){
        norm_thetas_j += theta_s(kk,jj)*theta_s(kk,jj);
      }
      norm_thetas_j = sqrt(norm_thetas_j); 
      arma::mat thetaMhat = alpha*thetaM+ (1-alpha)*alphaM;
      // X.submat( first_row, first_col, last_row, last_col )
      //arma::mat AlphajTmp = thetaM.submat(0,it*p+jj,2,it*p+jj) + wM.submat(0,p*(it-1)+jj,2,p*(it-1)+jj) ; 
      arma::mat AlphajTmp = thetaMhat.submat(0,it*p+jj,1,it*p+jj) + wM.submat(0,p*(it-1)+jj,1,p*(it-1)+jj) ; 
      
      double denomConstNorm = 0.00;
      // double test = 0.00;
      for(int kk=0;kk<M;kk++){
        denomConstNorm+= AlphajTmp(kk)*AlphajTmp(kk);
        // test+= (thetaM(kk,it*p+jj)+wM(kk,p*(it-1)+jj))*(thetaM(kk,it*p+jj)+wM(kk,p*(it-1)+jj));
      }
      // Rcout<< denomConstNorm<<"******"<<test<<"%%%%%%%%";
      
      denomConstNorm = sqrt(denomConstNorm);
      double secTerm = 0.00; 
      secTerm = 1-lambda/rho/norm_thetas_j/denomConstNorm;
      // Rcout<<secTerm <<" secTerm " << lambda<< " lambda "; 
      // Rcout<<AlphajTmp<<" alphajtmp ";
      //alphaM.submat(0,it*p+jj,2,it*p+jj)  = secTerm*AlphajTmp;
      // Rcout << 1-secTerm << " 1-secterm " << it << " iter " << lambda << " Lambda ";
      // Rcout << 1/rho/norm_thetas_j/denomConstNorm << " test ";
      if((1-secTerm)>=1){
        alphaM.submat(0,it*p+jj,1,it*p+jj)  = 0*secTerm*AlphajTmp;
      }
      else{
        alphaM.submat(0,it*p+jj,1,it*p+jj)  = secTerm*AlphajTmp;
      }
      
      
      //Rcout<< alphaM<< "******";
      // Rcout<<wM <<"$$$$$$"; 
    }
    
    // update of wM 
    //arma::mat wMtmp = wM.submat(0,p*(it-1),2,(p*it)-1) + thetaM.submat(0,p*it,2,(p*(it+1))-1) - alphaM.submat(0,p*it,2,(p*(it+1))-1) ;
    arma::mat thetaMhat2 = alpha*thetaM.submat(0,p*it,1,(p*(it+1))-1)+ (1-alpha)*alphaM.submat(0,p*(it-1),1,(p*(it))-1);
    
    arma::mat wMtmp = wM.submat(0,p*(it-1),1,(p*it)-1) + thetaMhat2 - alphaM.submat(0,p*it,1,(p*(it+1))-1) ;
    
    wM.submat(0,p*it,1,(p*(it+1))-1) = wMtmp;
    
    // check convergence 
    double epi_primal = 0.00;
    double norm_thetaM_it = 0.00;
    double norm_alphaM_it = 0.00;
    double norm_pri_res = 0.00;
    double norm_wM_it = 0.00;
    double norm_dual_res = 0.00;
    
    for(int i = 0;i<M;i++){
      for(int j=0; j<p;j++){
        norm_thetaM_it += thetaM(i,p*it+j)*thetaM(i,p*it+j);
        norm_alphaM_it += alphaM(i,p*it+j)*alphaM(i,p*it+j);
        norm_pri_res   += (thetaM(i,p*it+j)-alphaM(i,p*it+j))*(thetaM(i,p*it+j)-alphaM(i,p*it+j));
        norm_wM_it += wM(i,p*it+j)*wM(i,p*it+j);
        norm_dual_res += rho*(alphaM(i,p*it+j)-alphaM(i,p*(it-1)+j))*(alphaM(i,p*it+j)-alphaM(i,p*(it-1)+j));
      }
    }
    norm_alphaM_it = sqrt(norm_alphaM_it); 
    norm_thetaM_it = sqrt(norm_thetaM_it);
    norm_wM_it = sqrt(norm_wM_it);
    
    norm_dual_res = sqrt(norm_dual_res);
    norm_pri_res = sqrt(norm_pri_res);
    
    NumericVector vTmp = {norm_alphaM_it,norm_thetaM_it};
    epi_primal = sqrt(p*M)*epi_abs + epi_rel*max(vTmp);
    // int sumN = sum(setN); may be too large
    epi_dual = sqrt(M*p)*epi_abs + epi_rel*rho*norm_wM_it;
    // Rcout << epi_primal << "epi_primal" << epi_dual << "epi_dual"; 
    // Rcout << norm_pri_res << " norm_pri_res " << norm_dual_res << " norm_dual_res ";
    
    epi_pri_store(it) = epi_primal;
    epi_dual_store(it) = epi_dual;
    pri_res_store(it) = norm_pri_res;
    dual_res_store(it) = norm_dual_res;
    
    epi_store(_,0) = epi_pri_store;
    epi_store(_,1) = epi_dual_store;
    epi_store(_,2) = pri_res_store;
    epi_store(_,3) = dual_res_store;
    
    if (norm_pri_res <= epi_primal & norm_dual_res <= epi_dual){
      bbk = it;
      return Rcpp::List::create(Rcpp::Named("epi_store" )= epi_store,Rcpp::Named("augNMat" )= augNMat, Rcpp::Named("thetaM") = thetaM, Rcpp::Named("alphaM") = alphaM, Rcpp::Named("wM") = wM, Rcpp::Named("bbk") = bbk,Rcpp::Named("lambda") = lambda,Rcpp::Named("epi_primal") = epi_pri_store,Rcpp::Named("epi_dual") = epi_dual_store,Rcpp::Named("primal_res") = pri_res_store,Rcpp::Named("dual_res") = dual_res_store);
      
    }
    // Rcout<< max(NumericVector::create(12.3,1.2,13.3,34,10,12.45)) << "maxxxxx";
    
  } 
  //thetaM.submat(0,0,2,3).ones();  
  //Rcout<< thetaM;
  
  Rcout<< " It is not converge "; 
  return Rcpp::List::create(Rcpp::Named("epi_store" )= epi_store,Rcpp::Named("augNMat" )= augNMat, Rcpp::Named("thetaM") = thetaM, Rcpp::Named("alphaM") = alphaM, Rcpp::Named("wM") = wM, Rcpp::Named("bbk") = bbk,Rcpp::Named("lambda") = lambda,Rcpp::Named("epi_primal") = epi_pri_store,Rcpp::Named("epi_dual") = epi_dual_store,Rcpp::Named("primal_res") = pri_res_store,Rcpp::Named("dual_res") = dual_res_store);
  
}

//return Rcpp::List::create(Rcpp::Named("mu") = mu_store, Rcpp::Named("phi") = phi_store, Rcpp::Named("z_summed_by_row") = z_summed_by_row, Rcpp::Named("z_colmeans") = z_colmeans,Rcpp::Named("omega") = omega_store, Rcpp::Named("Theta") = Theta,Rcpp::Named("theta") = theta_store, Rcpp::Named("lambda") = lambda_store, Rcpp::Named("accept_theta") = accept_theta,Rcpp::Named("accept_omega") = accept_omega, Rcpp::Named("accept_lambda") = accept_lambda, Rcpp::Named("accept_mu") = accept_mu_y, Rcpp::Named("accept_phi") = accept_phi_y );
// IntegerVector dim = {3,4,11};
//NumericVector s1 = arrayC(s,dim);
//Rcout<< "oioi " << s1[2+dim[0]];
//int nn = s1.length();
//return Rcpp::List::create(Rcpp::Named("theta_s" )= theta_s);

