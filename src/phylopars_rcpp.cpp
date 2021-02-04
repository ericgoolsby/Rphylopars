// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
List C_anc_recon(arma::mat Y,arma::vec anc,arma::vec des,arma::vec edge_vec,int nedge,int nvar,int nspecies)
{
  arma::vec p = arma::zeros(nedge+1);
  arma::mat Yhat = arma::zeros(nedge+1,nvar);

  int i=0;
  int anc_edge=0;
  int des_edge=0;
  double pA=0;
  double len=0;
  for(i=0;i<nedge;i++)
  {
    anc_edge = anc(i)-1;
    des_edge = des(i)-1;
    len = edge_vec(i);
    
    if(des_edge<nspecies)
    {
      p(des_edge) = 1/len;
      Yhat.row(des_edge) = Y.row(des_edge);
    } else
    {
      pA = p(des_edge);
      Yhat.row(des_edge) = Yhat.row(des_edge)/pA;
      p(des_edge) = pA/(1+len*pA);
    }
    p(anc_edge) += p(des_edge);
    Yhat.row(anc_edge) += Yhat.row(des_edge)*p(des_edge);
  }
  Yhat.row(anc_edge) = Yhat.row(anc_edge)/p(anc_edge);
  
  for(i=nedge-1;i>=0;i--)
  {
    anc_edge = anc(i)-1;
    des_edge = des(i)-1;
    len = edge_vec(i);
    
    if(des_edge>=nspecies)
    {
      Yhat.row(des_edge) = Yhat.row(des_edge)*p(des_edge)*len + Yhat.row(anc_edge) - Yhat.row(anc_edge)*p(des_edge)*len;
      p(des_edge) = p(des_edge)/(1-len*p(des_edge)) + (p(anc_edge)-p(des_edge))/(1+len*(p(anc_edge)-p(des_edge)));
    } else
    {
      p(des_edge) = 0;
    }
  }
  return List::create(_["Yhat"] = Yhat,_["p"] = p);
}


// [[Rcpp::export]]
List C_anc_recon_rates(arma::mat Y,arma::vec anc,arma::vec des,arma::vec edge_vec,int nedge,int nvar,int nspecies,int REML)
{
  arma::vec p = arma::zeros(nedge+1);
  arma::mat Yhat = arma::zeros(nedge+1,nvar);
  arma::mat XY = arma::zeros(nedge+1,nvar);
  arma::mat YY = arma::zeros(nvar*(nedge+1),nvar);
  
  int i=0;
  int anc_edge=0;
  int des_edge=0;
  double pA=0;
  double len=0;
  for(i=0;i<nedge;i++)
  {
    anc_edge = anc(i)-1;
    des_edge = des(i)-1;
    len = edge_vec(i);
    
    if(des_edge<nspecies)
    {
      p(des_edge) = 1/len;
      Yhat.row(des_edge) = Y.row(des_edge);
      XY.row(des_edge) = Y.row(des_edge)/len;
      YY.rows(des_edge*nvar,nvar+(des_edge*nvar)-1) = trans(Y.row(des_edge))*Y.row(des_edge)/len;
    } else
    {
      pA = p(des_edge);
      Yhat.row(des_edge) = Yhat.row(des_edge)/pA;
      YY.rows(des_edge*nvar,nvar+(des_edge*nvar)-1) = YY.rows(des_edge*nvar,nvar+(des_edge*nvar)-1) -
        trans(XY.row(des_edge))*XY.row(des_edge)*len/(1+len*pA);
      XY.row(des_edge) = XY.row(des_edge)/(1+len*pA);
      p(des_edge) = pA/(1+len*pA);
    }
    p(anc_edge) += p(des_edge);
    Yhat.row(anc_edge) += Yhat.row(des_edge)*p(des_edge);
    XY.row(anc_edge) += XY.row(des_edge);
    YY.rows(anc_edge*nvar,nvar+(anc_edge*nvar)-1) += YY.rows(des_edge*nvar,nvar+(des_edge*nvar)-1);
  }
  Yhat.row(anc_edge) = Yhat.row(anc_edge)/p(anc_edge);
  arma::mat Sigma = (YY.rows(anc_edge*nvar,nvar+(anc_edge*nvar)-1) - 
    2*(trans(Yhat.row(anc_edge))*XY.row(anc_edge)) + trans(Yhat.row(anc_edge))*Yhat.row(anc_edge)*p(anc_edge))/(nspecies-REML);
  for(i=nedge-1;i>=0;i--)
  {
    anc_edge = anc(i)-1;
    des_edge = des(i)-1;
    len = edge_vec(i);
    
    if(des_edge>=nspecies)
    {
      Yhat.row(des_edge) = Yhat.row(des_edge)*p(des_edge)*len + Yhat.row(anc_edge) - Yhat.row(anc_edge)*p(des_edge)*len;
      p(des_edge) = p(des_edge)/(1-len*p(des_edge)) + (p(anc_edge)-p(des_edge))/(1+len*(p(anc_edge)-p(des_edge)));
    } else
    {
      p(des_edge) = 0;
    }
  }
  
  return List::create(_["Yhat"] = Yhat,_["p"] = p,_["Sigma"] = Sigma);
}

// [[Rcpp::export]]
arma::mat try_inv(arma::mat M,int nvar)
{
  arma::mat Minv;
  try
  {
    std::ostream nullstream(0);
    arma::set_cerr_stream(nullstream);
    //Minv = pinv(M);
    //M = try_clip(M,nvar,1);
    Minv = inv(M); //Minv = inv(M,"std");
    return Minv;
  }
  catch(...)
  {
    return arma::eye(nvar,nvar)+.0123456789;
  }
  return arma::eye(nvar,nvar)+.0123456789;
}

arma::mat try_solve(arma::mat M,arma::mat V)
{
  arma::mat Msolve;
  try
  {
    std::ostream nullstream(0);
    arma::set_cerr_stream(nullstream);
    //M = try_clip(M,nvar,1);
    Msolve = solve(M,V); //Msolve = solve(M,V,"std");
    //Msolve = pinv(M)*V;
    return Msolve;
  }
  catch(...)
  {
    return arma::ones(M.n_rows,V.n_cols)+.0123456789;
  }
  return arma::ones(M.n_rows,V.n_cols)+.0123456789;
}

// [[Rcpp::export]]
List inv_subset(arma::mat mat_to_inv,List subset_list)
{
  int subset_list_size = subset_list.size();
  int i=0;
  List ret(subset_list_size);
  for(i=0;i<subset_list_size;i++)
  {
    arma::uvec temp_vec = subset_list(i);
    unsigned int temp_vec_size = temp_vec.size();
    ret(i) = try_inv(mat_to_inv.submat(temp_vec,temp_vec),temp_vec_size);
  }
  return ret;
}

double logdet(arma::mat A)
{
  double val = 0;
  double sign = 1;
  try
  {
    std::ostream nullstream(0);
    arma::set_cerr_stream(nullstream);
    log_det(val,sign,A);
    return val;
  }
  catch(...)
  {
    return 1.0123456789;
  }
  return 1.0123456789;
}

// [[Rcpp::export]]
arma::vec mat_to_pars2(arma::mat M,int nvar,int diag)
{
  int len = (nvar*nvar-nvar)/2+nvar;
  if(diag==1)
  {
    arma::vec ret = log(sqrt(M.diag()));
    return ret;
  } else
  {
    M = chol(M,"lower");
    M.diag() = log(M.diag());
  }
  if((diag!=1) & (nvar>1))
  {
    M(1,0) = M(1,0)/exp(M(0,0));
  }
  
  int vec_count = 0;
  int i = 0;
  int j = 0;
  
  arma::vec ret(len);
  
  for(i=0;i<nvar;i++)
  {
    for(j=0;j<=i;j++)
    {
      ret(vec_count) = M(i,j);
      vec_count += 1;
    }
  }
  // (eigen(logm(p$pars[[1]]))$vec %*% diag(exp(eigen(logm(p$pars[[1]]))$val)) %*% t(eigen(logm(p$pars[[1]]))$vec))
  // (eigen(p$pars[[1]])$vec %*% diag(log(eigen(p$pars[[1]])$val)) %*% t(eigen(p$pars[[1]])$vec))
  
  return ret;
}

// [[Rcpp::export]]
arma::mat pars_to_mat(arma::vec pars,int nvar,int diag,int log_chol=1,int mod_chol=1)//,int exp_mat=0)
{
  //int len = (nvar*nvar-nvar)/2+nvar;
  arma::mat M(nvar,nvar,arma::fill::zeros);
  int vec_count = 0;
  int i = 0;
  int j = 0;
  
  for(i=0;i<nvar;i++)
  {
    for(j=0;j<=i;j++)
    {
      if((i!=j) && ((diag==0) || (diag==2)))
      {
        M(i,j) = pars(vec_count);
        vec_count += 1;
      } else if(i==j)
      {
        M(i,j) = pars(vec_count);
        vec_count += 1;
      }
    }
  }
  
  /*if(exp_mat==1)
  {
  return expmat(symmatl(M));
  }*/
  
  if((nvar>1) && (diag!=1))
  {
    if(log_chol && mod_chol) M(1,0)=M(1,0)*(exp(M(0,0)));
  }
  if(log_chol) M.diag() = exp(M.diag());
  
  if(diag==1)
  {
    M.diag() = square(M.diag());
  } else
  {
    M = M * M.t();
  }
  
  // (eigen(logm(p$pars[[1]]))$vec %*% diag(exp(eigen(logm(p$pars[[1]]))$val)) %*% t(eigen(logm(p$pars[[1]]))$vec))
  // lp <- (eigen(p$pars[[1]])$vec %*% diag(log(eigen(p$pars[[1]])$val)) %*% t(eigen(p$pars[[1]])$vec))
  // (eigen(lp)$vec %*% diag(exp(eigen(lp)$val)) %*% t(eigen(lp)$vec))
  
  return M;
}

// [[Rcpp::export]]
List calc_OU_len(arma::vec heights,arma::mat edge_mat,arma::vec des_order,int nedge,arma::mat P,arma::vec lambda,arma::mat sigma,arma::vec anc,arma::vec des,int nvar,int nspecies)
{
  arma::mat tP = trans(P);
  arma::mat invP = try_inv(P,nvar);
  arma::mat tinvP = trans(invP);
  arma::mat lambda_mat(nvar,nvar);
  int i=0;
  int j=0;
  for(i=0;i<nvar;i++)
  {
    for(j=0;j<nvar;j++)
    {
      lambda_mat(i,j) = lambda(i)+lambda(j);
    }
  }
  
  double t1=0;
  double t2=0;
  double t12=0;
  
  arma::mat eAC1(nvar,nvar);
  arma::mat teAC2(nvar,nvar);
  arma::mat M(nvar,nvar);
  arma::mat height_mat(nvar*(nedge+1),nvar,arma::fill::zeros);
  arma::mat height_mat2(nvar*(nedge+1),nvar,arma::fill::zeros);
  List len_mat(nedge+1);
  
  for(i=0;i<nedge;i++)
  {
    t1 = heights(edge_mat(i,0)-1);
    t2 = heights(edge_mat(i,1)-1);
    t12 = heights(edge_mat(i,2)-1);
    if(nvar>1)
    {
      eAC1 = P*diagmat(exp(-lambda*(t1-t12)))*invP;
      teAC2 = P*diagmat(exp(-lambda*(t2-t12)))*invP;
    } else
    {
      eAC1 = P*(exp(-lambda*(t1-t12)))*invP;
      teAC2 = P*(exp(-lambda*(t2-t12)))*invP;
    }
    M = P*(((1/lambda_mat)%(1-exp(-lambda_mat*t12)))%(invP*sigma*tinvP))*tP;
  //M = P*((1/lambda_mat)%(1-exp(-lambda_mat*t12))%invP*sigma*tinvP)*tP;
    height_mat.rows(i*nvar,(i+1)*nvar-1) = eAC1*M*teAC2;
  }

  for(i=0;i<(nedge+1);i++)
  {
    height_mat2.rows(i*nvar,(i+1)*nvar-1) = height_mat.rows((des_order(i)-1)*nvar,des_order(i)*nvar-1);
  }
  //return List::create(_["height_mat"] = height_mat2);
  for(i=0;i<nedge;i++)
  {
    len_mat(i) = height_mat2.rows((des(i)-1)*nvar,des(i)*nvar-1) - height_mat2.rows((anc(i)-1)*nvar,anc(i)*nvar-1);
  }
  len_mat(nedge) = height_mat.rows(nedge*nvar,(nedge+1)*nvar-1);
  
  return len_mat;
}

// [[Rcpp::export]]
List tp(arma::mat L,arma::mat R,arma::mat Rmat,int mL,int mR,int pheno_error,arma::vec edge_vec,arma::vec edge_ind,arma::vec ind_edge,arma::vec parent_edges,arma::vec pars,unsigned int nvar,int phylocov_diag,int nind,int nob,int nspecies,int nedge,arma::vec anc,arma::vec des,int REML,List species_subset,List un_species_subset,List subset_list,List ind_list,arma::vec tip_combn,LogicalVector is_edge_ind,arma::mat fixed_mu,List OU_len,arma::mat phylocov_fixed,arma::mat phenocov_fixed,List phenocov_list,int is_phylocov_fixed=0,int is_phenocov_fixed=0,int OU_par=0,int ret_level=1,int use_LL=0,int is_phenocov_list=0)
{
  /*
  subset_list is a list of unique matrices to invert once
  species_subset is a list of length nspecies
  
  edge_vec <- tree$edge.length
  species_ind <- match(trait_data$species,tree$tip.label)
  edge_ind <- match(species_ind,tree$edge[,2])-1
  
  p = 1'(V^-1)1
  L1 = L'(V^-1)1
  LL = L'(V^-1)L
  R1 = R'(V^-1)1
  RR = R'(V^-1)R
  LR = L'(V^-1)R
  logd = log(det(V))
  L is a matrix of size n by mL
  R is a matrix of size n by mR
  
  mL <- ncol(L)
  LL <- rep(list(matrix(0,mL,mL)),nedge+1)
  L1 <- rep(list(matrix(0,mL,nvar)),nedge+1)
  RR <- rep(list(matrix(0,mR,mR)),nedge+1)
  R1 <- rep(list(matrix(0,nvar,mR)),nedge+1)
  LR <- rep(list(matrix(0,mL,mR)),nedge+1)
  
  */
  
  arma::mat p(nvar*(nind+nedge+1),nvar,arma::fill::zeros);
  //arma::mat pA(nvar*(nind+nedge+1),nvar,arma::fill::zeros);
  arma::mat LL(mL*(nind+nedge+1),mL,arma::fill::zeros);
  arma::mat L1(mL*(nind+nedge+1),nvar,arma::fill::zeros);
  arma::mat RR(mR*(nind+nedge+1),mR,arma::fill::zeros);
  arma::mat R1(nvar*(nind+nedge+1),mR,arma::fill::zeros);
  arma::mat LR(mL*(nind+nedge+1),mR,arma::fill::zeros);
  arma::vec lnW(nind+nedge+1,arma::fill::zeros);
  
  arma::mat len_mat(nvar,nvar,arma::fill::zeros);
  double len = 0;
  
  arma::mat phylocov(nvar,nvar,arma::fill::zeros);
  arma::mat phenocov(nvar,nvar,arma::fill::zeros);
  arma::mat pA(nvar,nvar,arma::fill::zeros);
  arma::mat itpa(nvar,nvar,arma::fill::zeros);
  arma::mat itpainv(nvar,nvar,arma::fill::zeros);
  int phylocov_par_length = (nvar*nvar-nvar)/2+nvar;
  if(phylocov_diag==1) phylocov_par_length = nvar;
  int phenocov_par_length = (nvar*nvar-nvar)/2+nvar;
  if(pheno_error==1) phenocov_par_length = nvar;
  
  if(is_phylocov_fixed==0)
  {
    phylocov = pars_to_mat(pars.subvec(0,phylocov_par_length-1),nvar,phylocov_diag);
  } else
  {
    phylocov = phylocov_fixed;
    phylocov_par_length = 0;
  }
  int pheno_diag = 0;
  if(pheno_error==0)
  {
    phenocov = phylocov;
  } else
  {
    if(is_phenocov_fixed==0)
    {
      if(pheno_error==1) pheno_diag = 1;
      phenocov = pars_to_mat(pars.subvec(phylocov_par_length,phylocov_par_length+phenocov_par_length-1),nvar,pheno_diag);
    } else
    {
      phenocov = phenocov_fixed;
    }
  }
  List sub_mats = inv_subset(phenocov,subset_list);
  int i = 0;
  int anc_edge = 0;
  int des_edge = 0;
  
  arma::uvec zero_vec(1);
  zero_vec(0) = 0;
  for(i=0;i<nind;i++)
  {
    anc_edge = nind + edge_ind(i);
    arma::uvec Ka = species_subset(i);
    if(Ka.size()>0)
    {
      arma::uvec Ia = ind_list(i);
      arma::mat Bainv = sub_mats(tip_combn(i));
      if(Bainv(0,0)==1.0123456789)
      {
        return List::create(_["logl"] = -1e128);
      }
      if(pheno_error==0)
      {
        if(OU_par==1)
        {
          phenocov = Rcpp::as<arma::mat>(OU_len(edge_ind(i)));
          if(is_phenocov_list>0) phenocov.submat(Ka,Ka) += Rcpp::as<arma::mat>(phenocov_list(i)).submat(Ka,Ka);
          Bainv = try_inv(phenocov.submat(Ka,Ka),Ka.size());
          if(Bainv(0,0)==1.0123456789)
          {
            return List::create(_["logl"] = -1e128);
          }
        } else
        {
          len = edge_vec(edge_ind(i));
          phenocov = phylocov*len;
          if(is_phenocov_list>0)
          {
            phenocov.submat(Ka,Ka) += Rcpp::as<arma::mat>(phenocov_list(i)).submat(Ka,Ka);
            Bainv = try_inv(phenocov.submat(Ka,Ka),Ka.size());
          } else Bainv = Bainv/len;
        }
      }
      arma::mat R_i = R.rows(Ia);
      arma::mat L_i = L.submat(Ia,Ka);
      p.submat((i*nvar)+Ka,Ka) = Bainv;
      
      if(use_LL==1)
      {
        L1.submat((i*nvar)+Ka,Ka) = trans(L_i)*Bainv;
        LL.submat((i*nvar)+Ka,Ka) = trans(L_i)*Bainv*L_i;
        LR.submat((i*nvar)+Ka,zero_vec) = trans(L_i)*Bainv*R_i;
        R1.submat((i*nvar)+Ka,zero_vec) = Bainv*R_i;
        
      } else
      {
        L1.submat((i*nvar)+Ka,Ka) = Bainv;
        LL.submat((i*nvar)+Ka,Ka) = Bainv;
        LR.submat((i*nvar)+Ka,zero_vec) = Bainv*R_i;
        R1.submat((i*nvar)+Ka,zero_vec) = LR.submat((i*nvar)+Ka,zero_vec);
      }
      RR.submat(zero_vec+i,zero_vec) = trans(R_i)*Bainv*R_i;
      lnW(i) = logdet(phenocov.submat(Ka,Ka));
      
      L1.submat((anc_edge*nvar)+Ka,Ka) += L1.submat((i*nvar)+Ka,Ka);
      LL.submat((anc_edge*nvar)+Ka,Ka) += LL.submat((i*nvar)+Ka,Ka);
      LR.submat((anc_edge*nvar)+Ka,zero_vec) += LR.submat((i*nvar)+Ka,zero_vec);
      R1.submat((anc_edge*nvar)+Ka,zero_vec) += R1.submat((i*nvar)+Ka,zero_vec);
      RR.submat(zero_vec+anc_edge,zero_vec) += RR.submat(zero_vec+i,zero_vec);
      p.submat((anc_edge*nvar)+Ka,Ka) += p.submat((i*nvar)+Ka,Ka);
      lnW(anc_edge) += lnW(i);
    }
  }
  
  for(i=0;i<=nedge;i++)
  {
    len = edge_vec(i);
    if((pheno_error==0) && is_edge_ind(i))
    {
      len = 0;
    }
    des_edge = i + nind;
    anc_edge = nind + parent_edges(i);
    pA = p.rows(des_edge*nvar,des_edge*nvar+nvar-1);
    if(OU_par==1)
    {
      if((pheno_error==0) && is_edge_ind(i))
      {
        itpa = arma::eye(nvar,nvar);
      } else
      {
        itpa = arma::eye(nvar,nvar) + Rcpp::as<arma::mat>(OU_len(i))*pA;
      }
    } else
    {
      itpa = arma::eye(nvar,nvar) + (len*phylocov)*pA;
    }
    itpainv = try_inv(itpa,nvar);
    if(itpainv(0,0)==1.0123456789)
    {
      return List::create(_["logl"] = -1e128);
    }
    lnW(des_edge) += logdet(itpa);
    if(des_edge!=anc_edge) lnW(anc_edge) += lnW(des_edge);
    
    p.rows(des_edge*nvar,des_edge*nvar+nvar-1) = pA*itpainv;
    if(des_edge!=anc_edge) p.rows(anc_edge*nvar,anc_edge*nvar+nvar-1) += p.rows(des_edge*nvar,des_edge*nvar+nvar-1);
    
    
    if(OU_par==1)
    {
      if(((pheno_error==0) &&  (is_edge_ind(i)==0)) || pheno_error) LR.rows(des_edge*nvar,des_edge*nvar+nvar-1) -= L1.rows(des_edge*nvar,des_edge*nvar+nvar-1)*itpainv*Rcpp::as<arma::mat>(OU_len(i))*R1.rows(des_edge*nvar,des_edge*nvar+nvar-1);
    } else
    {
      LR.rows(des_edge*nvar,des_edge*nvar+nvar-1) -= L1.rows(des_edge*nvar,des_edge*nvar+nvar-1)*itpainv*(len*phylocov)*R1.rows(des_edge*nvar,des_edge*nvar+nvar-1);
    }
    
    if(des_edge!=anc_edge)
    {
      LR.rows(anc_edge*nvar,anc_edge*nvar+nvar-1) += LR.rows(des_edge*nvar,des_edge*nvar+nvar-1);
    } else if(fixed_mu(0,0)!=0)
    {
      LR.rows(des_edge*nvar,des_edge*nvar+nvar-1) = LL.rows(des_edge*nvar,des_edge*nvar+nvar-1)*fixed_mu;
    }
    
    if(OU_par==1)
    {
      if(((pheno_error==0) && (is_edge_ind(i)==0)) || pheno_error) LL.rows(des_edge*nvar,des_edge*nvar+nvar-1) -= trans(L1.rows(des_edge*nvar,des_edge*nvar+nvar-1))*itpainv*Rcpp::as<arma::mat>(OU_len(i))*L1.rows(des_edge*nvar,des_edge*nvar+nvar-1);
    } else
    {
      LL.rows(des_edge*nvar,des_edge*nvar+nvar-1) -= trans(L1.rows(des_edge*nvar,des_edge*nvar+nvar-1))*itpainv*(len*phylocov)*L1.rows(des_edge*nvar,des_edge*nvar+nvar-1);
    }
    if(des_edge!=anc_edge) LL.rows(anc_edge*nvar,anc_edge*nvar+nvar-1) += LL.rows(des_edge*nvar,des_edge*nvar+nvar-1);
    
    if(OU_par==1)
    {
      if(((pheno_error==0) &&  (is_edge_ind(i)==0)) || pheno_error) RR.row(des_edge) -= trans(R1.rows(des_edge*nvar,des_edge*nvar+nvar-1))*itpainv*Rcpp::as<arma::mat>(OU_len(i))*R1.rows(des_edge*nvar,des_edge*nvar+nvar-1);
    } else
    {
      RR.row(des_edge) -= trans(R1.rows(des_edge*nvar,des_edge*nvar+nvar-1))*itpainv*(len*phylocov)*R1.rows(des_edge*nvar,des_edge*nvar+nvar-1);
    }
    
    if(des_edge!=anc_edge) RR.row(anc_edge) += RR.row(des_edge);
    
    L1.rows(des_edge*nvar,des_edge*nvar+nvar-1) = L1.rows(des_edge*nvar,des_edge*nvar+nvar-1)*itpainv;
    if(des_edge!=anc_edge) L1.rows(anc_edge*nvar,anc_edge*nvar+nvar-1) += L1.rows(des_edge*nvar,des_edge*nvar+nvar-1);
    
    R1.rows(des_edge*nvar,des_edge*nvar+nvar-1) = trans(trans(R1.rows(des_edge*nvar,des_edge*nvar+nvar-1))*itpainv);
    if(des_edge!=anc_edge) R1.rows(anc_edge*nvar,anc_edge*nvar+nvar-1) += R1.rows(des_edge*nvar,des_edge*nvar+nvar-1);
  }
  arma::mat anc_vec = try_solve(LL.rows(des_edge*nvar,des_edge*nvar+nvar-1),LR.rows(des_edge*nvar,des_edge*nvar+nvar-1));
  if(anc_vec(0,0)==1.0123456789)
  {
    return List::create(_["logl"] = -1e128);
  }
  
  /*
  if (method=="ML"){
  minus2ll <- nob*log(2*pi) + lnW[des_edge] +
  YY[des_edge] - 2 * t(anc) %*% XY[[des_edge]] + t(anc) %*% XX[[des_edge]] %*% anc
  }
  if (method=="REML"){
  minus2ll <- (nob-nvar)*log(2*pi) + log(det(XX[[des_edge]])) + lnW[des_edge] +
  YY[des_edge] - 2 * t(anc) %*% XY[[des_edge]] + t(anc) %*% XX[[des_edge]] %*% anc
  }
  -((nspecies*nvar-nvar)*log(2*pi) + determinant(r$LL)$modulus[[1]] + r$logd + r$RR - 2*t(r$theta) %*% r$LR + t(r$theta) %*% r$LL %*% r$theta)/2
  */
  arma::vec logl(1);
  if(REML==1)
  {
    logl = -.5*((nob-nvar)*log(2*arma::datum::pi) + logdet(LL.rows(des_edge*nvar,des_edge*nvar+nvar-1)) + lnW(des_edge) + RR.row(des_edge) - 2*trans(anc_vec) * LR.rows(des_edge*nvar,des_edge*nvar+nvar-1) + trans(anc_vec) *  LL.rows(des_edge*nvar,des_edge*nvar+nvar-1) * anc_vec);
  } else
  {
    logl = -.5*(nob*log(2*arma::datum::pi) + lnW(des_edge) + RR.row(des_edge) - 2*trans(anc_vec) * LR.rows(des_edge*nvar,des_edge*nvar+nvar-1) + trans(anc_vec) *  LL.rows(des_edge*nvar,des_edge*nvar+nvar-1) * anc_vec);
  }
  
  if(ret_level<2)
  {
    return List::create(_["logl"] = logl);
  }
  
  arma::mat ret_p = p.rows(des_edge*nvar,des_edge*nvar+nvar-1);
  arma::mat ret_L1 = L1.rows(des_edge*nvar,des_edge*nvar+nvar-1);
  arma::mat ret_LL = LL.rows(des_edge*nvar,des_edge*nvar+nvar-1);
  arma::mat ret_R1 = R1.rows(des_edge*nvar,des_edge*nvar+nvar-1);
  arma::mat ret_RR = RR.row(des_edge);
  arma::mat ret_LR = LR.rows(des_edge*nvar,des_edge*nvar+nvar-1);
  arma::mat ret_theta = anc_vec;
  
  if(ret_level<3)
  {
    return List::create(_["logl"] = logl,_["theta"]=ret_theta);
  }
  
  arma::mat des_p(nvar,nvar,arma::fill::zeros);
  arma::mat other_p(nvar,nvar,arma::fill::zeros);
  arma::mat other_pA(nvar,nvar,arma::fill::zeros);
  arma::mat other_LL(nvar,nvar,arma::fill::zeros);
  arma::mat other_LR(nvar,nvar,arma::fill::zeros);
  
  arma::mat anc_recon(des.size()+1,nvar,arma::fill::zeros);
  arma::mat recon_var(nvar*(des.size()+1),nvar,arma::fill::zeros);
  
  anc_recon.row(nspecies) = trans(anc_vec);
  recon_var.rows(nspecies*nvar,(nspecies+1)*nvar -1) = try_inv(LL.rows(des_edge*nvar,des_edge*nvar+nvar-1),nvar);
  arma::mat tip_uncertainty(nvar,nvar,arma::fill::zeros);
  
  for(i=nedge-1;i>=0;i--)
  {
    len = edge_vec(i);
    des_edge = i + nind;
    anc_edge = parent_edges(i) + nind;
    if((pheno_error==0) && is_edge_ind(i) && (is_phenocov_list==0))
    {
      /*
      XX[[des_edge]][-Ka,-Ka] <- other_p[-Ka,-Ka]
      XX[[des_edge]][Ka,Ka] <- 0
      
      imputed_data[match(i,edge_ind),(1:nvar)[-Ka]+1] <- 
      solve(XX[[anc_edge]],XY[[anc_edge]])[-Ka,] + swp(phylocov,Ka)[-Ka,Ka] %*% t(as.matrix(trait_data[match(i,edge_ind),1+Ka,drop=F])-solve(XX[[anc_edge]],XY[[anc_edge]])[Ka,])
      
      tip_uncertainty[-Ka,-Ka] <- tip_uncertainty[-Ka,-Ka] + swp(phylocov,(1:nvar)[Ka])[-Ka,-Ka,drop=FALSE]
      */
      arma::uvec Ka = species_subset(ind_edge(i));
      arma::mat temp_anc(nvar,1);
      if(Ka.size()<nvar)
      {
        arma::uvec un_Ka = un_species_subset(ind_edge(i));
        other_pA = p.rows(anc_edge*nvar,anc_edge*nvar+nvar-1) - p.rows(des_edge*nvar,des_edge*nvar+nvar-1);
        
        if(OU_par==1)
        {
          other_p = other_pA * try_inv(arma::eye(nvar,nvar) + Rcpp::as<arma::mat>(OU_len(i))*other_pA,nvar);
        } else
        {
          other_p = other_pA * try_inv(arma::eye(nvar,nvar) + (phylocov*len)*other_pA,nvar);
        }
        if(other_p(0,0)==1.0123456789)
        {
          return List::create(_["logl"] = -1e128);
        }
        p.submat((des_edge*nvar)+un_Ka,un_Ka) += other_p.submat(un_Ka,un_Ka);
        if(Ka.size()>0) p.submat((des_edge*nvar)+Ka,Ka) *= 0;
        
        temp_anc = try_solve(LL.rows(anc_edge*nvar,anc_edge*nvar+nvar-1),LR.rows(anc_edge*nvar,anc_edge*nvar+nvar-1));
        if(temp_anc(0,0)==1.0123456789)
        {
          return List::create(_["logl"] = -1e128);
        }
        arma::uvec Ia = ind_list(ind_edge(i));
        arma::mat R_i = R.rows(Ia);
        
        other_LL = LL.rows(anc_edge*nvar,anc_edge*nvar+nvar-1) - LL.rows(des_edge*nvar,des_edge*nvar+nvar-1);
        if(OU_par==1)
        {
          other_LL = trans(trans(other_LL)*try_inv(arma::eye(nvar,nvar) + Rcpp::as<arma::mat>(OU_len(i))*other_LL,nvar));
        } else
        {
          other_LL = trans(trans(other_LL)*try_inv(arma::eye(nvar,nvar) + (phylocov*len)*other_LL,nvar));
        }
        if(other_LL(0,0)==1.0123456789)
        {
          return List::create(_["logl"] = -1e128);
        }
        //LL.rows(des_edge*nvar,des_edge*nvar+nvar-1) = trans(trans(LL.rows(des_edge*nvar,des_edge*nvar+nvar-1))*itpainv) + other_LL;
        LL.submat(des_edge*nvar+un_Ka,un_Ka) += other_LL.submat(un_Ka,un_Ka);
        if(Ka.size()>0) LL.submat(des_edge*nvar+Ka,Ka) *= 0;
        
        //arma::mat L_i = L.submat(Ia,Ka);
        //row(des(i)).
        
        
        if(Ka.size()==0)
        {
          anc_recon.submat(zero_vec+des(i),un_Ka) = trans(temp_anc(un_Ka));
          Rmat.submat(zero_vec+ind_edge(i),un_Ka) = anc_recon.submat(zero_vec+des(i),un_Ka); //trans(temp_anc(un_Ka) + trans(try_solve(phenocov.submat(Ka,Ka),phenocov.submat(Ka,un_Ka)))*(R_i-temp_anc(Ka)));
          tip_uncertainty.submat(un_Ka,un_Ka) += phylocov.submat(un_Ka,un_Ka);
          recon_var.submat(des(i)*nvar+un_Ka,un_Ka) = try_inv(LL.submat(des_edge*nvar+un_Ka,un_Ka),un_Ka.size()); // try_inv(other_p,nvar);
        } else
        {
          anc_recon.submat(zero_vec+des(i),un_Ka) = trans(temp_anc(un_Ka) + trans(try_solve(phylocov.submat(Ka,Ka),phylocov.submat(Ka,un_Ka)))*(R_i-temp_anc(Ka)));
          Rmat.submat(zero_vec+ind_edge(i),un_Ka) = anc_recon.submat(zero_vec+des(i),un_Ka); //trans(temp_anc(un_Ka) + trans(try_solve(phenocov.submat(Ka,Ka),phenocov.submat(Ka,un_Ka)))*(R_i-temp_anc(Ka)));
          anc_recon.submat(zero_vec+des(i),Ka) = Rmat.submat(zero_vec+ind_edge(i),Ka);
          tip_uncertainty.submat(un_Ka,un_Ka) += phylocov.submat(un_Ka,un_Ka) - phylocov.submat(un_Ka,Ka)*try_inv(phylocov.submat(Ka,Ka),Ka.size())*phylocov.submat(Ka,un_Ka);
          recon_var.submat(des(i)*nvar+un_Ka,un_Ka) = try_inv(LL.submat(des_edge*nvar+un_Ka,un_Ka),un_Ka.size()); // try_inv(other_p,nvar);
        }
        
        
        //anc_recon.row(des(i)) = trans(try_solve(LL.rows(des_edge*nvar,des_edge*nvar+nvar-1),LR.rows(des_edge*nvar,des_edge*nvar+nvar-1)));
        //recon_var.rows(des(i)*nvar,des(i)*nvar+nvar-1) = try_inv(LL.rows(des_edge*nvar,des_edge*nvar+nvar-1),nvar);
        
        
        //imputed = temp_anc.rows(un_Ka) + trans(solve(phylocov.submat(Ka,Ka),phylocov.submat(Ka,un_Ka))) * trans(R_i-temp_anc.elem(Ka));
      } else
      {
        p.submat((des_edge*nvar)+Ka,Ka) *= 0;
        LL.submat(des_edge*nvar+Ka,Ka) *= 0;
        anc_recon.submat(zero_vec+des(i),Ka) = Rmat.submat(zero_vec+ind_edge(i),Ka);
        Rmat.submat(zero_vec+ind_edge(i),Ka) = anc_recon.submat(zero_vec+des(i),Ka);
        recon_var.submat(des(i)*nvar+Ka,Ka) *= 0;
      }
    } else
    {
      pA = p.rows(des_edge*nvar,des_edge*nvar+nvar-1);
      if(OU_par==1)
      {
        itpa = arma::eye(nvar,nvar) - Rcpp::as<arma::mat>(OU_len(i))*p.rows(des_edge*nvar,des_edge*nvar+nvar-1);
      } else
      {
        itpa = arma::eye(nvar,nvar) - (phylocov*len)*p.rows(des_edge*nvar,des_edge*nvar+nvar-1);
      }
      itpainv = try_inv(itpa,nvar);
      if(itpainv(0,0)==1.0123456789)
      {
        return List::create(_["logl"] = -1e128);
      }
      des_p = pA*itpainv;
      other_pA = p.rows(anc_edge*nvar,anc_edge*nvar+nvar-1) - p.rows(des_edge*nvar,des_edge*nvar+nvar-1);
      if(OU_par==1)
      {
        other_p = other_pA * try_inv(arma::eye(nvar,nvar) + Rcpp::as<arma::mat>(OU_len(i))*other_pA,nvar);
      } else
      {
        other_p = other_pA * try_inv(arma::eye(nvar,nvar) + (phylocov*len)*other_pA,nvar);
      }
      if(other_p(0,0)==1.0123456789)
      {
        return List::create(_["logl"] = -1e128);
      }
      p.rows(des_edge*nvar,des_edge*nvar+nvar-1) = des_p + other_p;
      
      other_LL = LL.rows(anc_edge*nvar,anc_edge*nvar+nvar-1) - LL.rows(des_edge*nvar,des_edge*nvar+nvar-1);
      if(OU_par==1)
      {
        other_LL = trans(trans(other_LL)*try_inv(arma::eye(nvar,nvar) + Rcpp::as<arma::mat>(OU_len(i))*other_pA,nvar));
      } else
      {
        other_LL = trans(trans(other_LL)*try_inv(arma::eye(nvar,nvar) + (phylocov*len)*other_pA,nvar));
      }
      LL.rows(des_edge*nvar,des_edge*nvar+nvar-1) = trans(trans(LL.rows(des_edge*nvar,des_edge*nvar+nvar-1))*itpainv) + other_LL;
      
      other_LR = LR.rows(anc_edge*nvar,anc_edge*nvar+nvar-1) - LR.rows(des_edge*nvar,des_edge*nvar+nvar-1);
      
      if(OU_par==1)
      {
        other_LR = trans(trans(other_LR)*try_inv(arma::eye(nvar,nvar) + Rcpp::as<arma::mat>(OU_len(i))*other_pA,nvar));
      } else
      {
        other_LR = trans(trans(other_LR)*try_inv(arma::eye(nvar,nvar) + (phylocov*len)*other_pA,nvar));
      }
      LR.rows(des_edge*nvar,des_edge*nvar+nvar-1) = trans(trans(LR.rows(des_edge*nvar,des_edge*nvar+nvar-1))*itpainv) + other_LR;
      
      anc_recon.row(des(i)) = trans(try_solve(LL.rows(des_edge*nvar,des_edge*nvar+nvar-1),LR.rows(des_edge*nvar,des_edge*nvar+nvar-1)));
      recon_var.rows(des(i)*nvar,des(i)*nvar+nvar-1) = try_inv(LL.rows(des_edge*nvar,des_edge*nvar+nvar-1),nvar);
      // node_var = try_inv(LL.rows(des_edge*nvar,des_edge*nvar+nvar-1),nvar);
      
    }
  }
  
  
  if(pheno_error!=0)
  {
    arma::mat temp_anc(nvar,1);
    for(i=0;i<nind;i++)
    {
      des_edge = i;
      anc_edge = nind + edge_ind(i);
      arma::uvec Ka = species_subset(i);
      arma::uvec Ia = ind_list(i);
      if(Ka.size()<nvar)
      {
        arma::uvec un_Ka = un_species_subset(i);
        temp_anc = trans(try_solve(LL.rows(anc_edge*nvar,anc_edge*nvar+nvar-1),LR.rows(anc_edge*nvar,anc_edge*nvar+nvar-1)));
        arma::mat Bainv = sub_mats(tip_combn(i));
        arma::mat R_i = R.rows(Ia);
        //arma::mat L_i = L.submat(Ia,Ka);
        if(Ka.size()>0)
        {
          Rmat.submat(zero_vec+i,un_Ka) = trans(temp_anc(un_Ka) + trans(try_solve(phenocov.submat(Ka,Ka),phenocov.submat(Ka,un_Ka)))*(R_i-temp_anc(Ka)));
          tip_uncertainty.submat(un_Ka,un_Ka) += phenocov.submat(un_Ka,un_Ka) - phenocov.submat(un_Ka,Ka)*try_inv(phenocov.submat(Ka,Ka),Ka.size())*phenocov.submat(Ka,un_Ka);
        } else
        {
          Rmat.submat(zero_vec+i,un_Ka) = trans(temp_anc(un_Ka));
          tip_uncertainty.submat(un_Ka,un_Ka) += phenocov.submat(un_Ka,un_Ka);
        }
      }
    }
  }
  return List::create(_["logl"] = logl,
                      _["theta"] = ret_theta,
                      _["p"] = ret_p,
                      _["L1"] = ret_L1,
                      _["LL"] = ret_LL,
                      _["R1"] = ret_R1,
                      _["RR"] = ret_RR,
                      _["LR"] = ret_LR,
                      _["anc_recon"] = anc_recon,
                      _["recon_var"] = recon_var,
                      _["recon_ind"] = Rmat,
                      _["tip_uncertainty"] = tip_uncertainty);
}

// [[Rcpp::export]]
arma::mat EM_Fels2008(arma::mat pics,arma::vec vars,arma::mat phylocov,arma::mat phenocov,int nvar,arma::mat phylocov_fixed,arma::mat phenocov_fixed,int is_phylocov_fixed=0,int is_phenocov_fixed=0,int diag_pheno=0,int EM_Fels_limit=5000,double tol=1e-6,int REML=1,int diag_phylo=0)
{
  arma::mat B_phylocov = arma::zeros(nvar,nvar);
  arma::mat B_phenocov = arma::zeros(nvar,nvar);
  arma::mat zer = arma::zeros(nvar,nvar);
  arma::mat I = arma::eye(nvar,nvar);
  arma::mat phylocov_accum = arma::zeros(nvar,nvar);
  arma::mat phenocov_accum = arma::zeros(nvar,nvar);
  arma::mat inv_phenocov = arma::zeros(nvar,nvar);
  arma::mat inv = arma::zeros(nvar,nvar);
  arma::mat tc = arma::zeros(nvar,nvar);
  arma::mat vap = arma::zeros(nvar,nvar);
  arma::mat max_phylocov = arma::zeros(nvar,nvar);
  arma::mat max_phenocov = arma::zeros(nvar,nvar);
  double det_phenocov = 0;
  double sum_det = 0;
  int i = 0;
  int j = 0;
  int len = vars.size();
  arma::vec res = arma::zeros(1);
  arma::vec logl = arma::zeros(1);
  double maxlogl = -1e128;

  if(is_phylocov_fixed==1) phylocov = phylocov_fixed;  
  if(is_phenocov_fixed==1) phenocov = phenocov_fixed;
  
  for(j=0;j<EM_Fels_limit;j++)
  {
    phylocov_accum = zer;
    phenocov_accum = zer;
    inv_phenocov = try_inv(phenocov,nvar);
    det_phenocov = logdet(phenocov);
    sum_det = 0;
    res(0) = 0;
    for(i=0;i<len;i++)
    {
      tc = trans(pics.row(i))*pics.row(i);
      if(vars(i)==0)
      {
        inv = inv_phenocov;
        phylocov_accum += phylocov;
        phenocov_accum += tc;
        sum_det += det_phenocov;
      } else
      {
        vap = vars(i)*phylocov + phenocov;
        inv = try_inv(vap,nvar);
        sum_det += logdet(vap);
        B_phylocov = sqrt(vars(i))*phylocov*inv;
        B_phenocov = phenocov*inv;
        phylocov_accum += (phylocov + (B_phylocov)*(tc-vap)*trans(B_phylocov));
        phenocov_accum += (phenocov + (B_phenocov)*(tc-vap)*trans(B_phenocov));
      }
      res += pics.row(i)*inv*trans(pics.row(i));
    }
    if(is_phylocov_fixed==0) phylocov = phylocov_accum / (len - (REML-1));
    if(is_phenocov_fixed==0) phenocov = phenocov_accum / (len - (REML-1));
    if(diag_pheno!=0) phenocov = diagmat(phenocov);
    if(diag_phylo!=0) phylocov = diagmat(phylocov);
    logl = -(res+sum_det);
    if((logl(0)-maxlogl)>tol)
    {
      maxlogl = logl(0);
      max_phylocov = phylocov;
      max_phenocov = phenocov;
    } else break;
  }
  arma::mat ret = arma::zeros(nvar*2,nvar);
  //ret.submat(0,0,nvar-1,nvar-1) = max_phylocov;
  //ret.submat(nvar,0,nvar*2-1,nvar-1) = max_phenocov;
  ret.submat(0,0,nvar-1,nvar-1) = phylocov;
  ret.submat(nvar,0,nvar*2-1,nvar-1) = phenocov;
  
  return ret;
}
