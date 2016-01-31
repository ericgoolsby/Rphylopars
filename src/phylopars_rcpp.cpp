// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat try_inv(arma::mat M,int nvar)
{
  arma::mat Minv;
  try
  {
    std::ostream nullstream(0);
    arma::set_stream_err2(nullstream);
    //Minv = pinv(M);
    //M = try_clip(M,nvar,1);
    Minv = inv(M,"std");    
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
    arma::set_stream_err2(nullstream);
    //M = try_clip(M,nvar,1);
    Msolve = solve(M,V,"std");
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
    arma::set_stream_err2(nullstream);
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
  if(diag!=1 & nvar>1)
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
arma::mat pars_to_mat(arma::vec pars,int nvar,int diag,int log_chol=0,int mod_chol=0)//,int exp_mat=0)
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
List tp(arma::mat L,arma::mat R,arma::mat Rmat,int mL,int mR,int pheno_error,arma::vec edge_vec,arma::vec edge_ind,arma::vec ind_edge,arma::vec parent_edges,arma::vec pars,int nvar,int phylocov_diag,int nind,int nob,int nspecies,int nedge,arma::vec anc,arma::vec des,int REML,List species_subset,List un_species_subset,List subset_list,List ind_list,arma::vec tip_combn,LogicalVector is_edge_ind,arma::mat fixed_mu,List OU_len,arma::mat phylocov_fixed,arma::mat phenocov_fixed,int is_phylocov_fixed=0,int is_phenocov_fixed=0,int OU_par=0,int ret_level=1,int use_LL=0)
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
          Bainv = try_inv(phenocov.submat(Ka,Ka),Ka.size());
          if(Bainv(0,0)==1.0123456789)
          {
            return List::create(_["logl"] = -1e128);
          }
        } else
        {
          len = edge_vec(edge_ind(i));
          phenocov = phylocov*len;
          Bainv = Bainv/len;
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
    if((pheno_error==0) && is_edge_ind(i))
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

// [[Rcpp::export]]
arma::mat cDot(arma::mat A,arma::mat B)
{
  return A*B;
}

// [[Rcpp::export]]
arma::mat fast_inv(int nspecies,int nedge,arma::vec len_vec,arma::vec anc,arma::vec des)
{
  arma::vec p = arma::zeros(nedge+1);
  arma::vec pA = arma::zeros(nedge+1);
  arma::vec detV = arma::zeros(nedge+1);
  int i = 0, des_node = 0, anc_node = 0;
  double len = 0, pA_des_value = 0;
  
  for(i=0;i<nedge;i++)
  {
    des_node = des(i);
    anc_node = anc(i);
    len = len_vec(i);
    
    if(des_node<nspecies)
    {
      p(des_node) = 1/len;
      detV(des_node) = log(len);
      pA(anc_node) = pA(anc_node) + p(des_node);
      detV(anc_node) = detV(anc_node) + detV(des_node);
      p(anc_node) = p(anc_node) + p(des_node);
      
    } else
    {
      pA_des_value = pA(des_node);
      detV(des_node) = detV(des_node) + log(1+len*pA_des_value);
      p(des_node) = pA_des_value / (1+len*pA_des_value);
      pA(anc_node) = pA(anc_node) + p(des_node);
      detV(anc_node) = detV(anc_node) + detV(des_node);
    }
  }
  
  len = len_vec(nedge);
  pA(nedge) = pA(anc_node);
  p(nedge) = pA(nedge) / (1+len*pA(nedge));
  detV(nedge) = detV(nspecies) + log(1+len*pA(nedge));
  
  //return detV(nedge);
  //return p(nedge);
  return pA;
}

// [[Rcpp::export]]
arma::mat testmat(int nrow)
{
  arma::mat matr = arma::zeros(nrow,nrow);
  return matr;
}

// [[Rcpp::export]]
arma::mat getInv(int nspecies,int nedge,arma::vec len_vec,arma::vec anc,arma::vec des,int inv=1)
{
  arma::vec p = arma::zeros(nedge+1);
  arma::vec pA = arma::zeros(nedge+1);
  arma::mat children = arma::zeros(nedge+1,2);
  arma::mat Vstar = arma::zeros(nspecies,nspecies);
  arma::mat Vinv = arma::zeros(nspecies,nspecies);
  children.col(0) = children.col(0) + nedge + 2;
  int i = 0;
  double des_node = 0, anc_node = 0;
  double len = 0, pA_des_value = 0;
  
  for(i=0;i<nedge;i++)
  {
    des_node = des(i);
    anc_node = anc(i);
    len = len_vec(i);
    
    if(des_node<nspecies)
    {
      p(des_node) = 1/len;
      pA(anc_node) = pA(anc_node) + p(des_node);
      p(anc_node) = p(anc_node) + p(des_node);
      children(des_node,0) = des_node;
      children(des_node,1) = 1;
      children(anc_node,0) = std::min(children(anc_node,0),des_node);
      children(anc_node,1) = children(anc_node,1) + 1;
      Vstar(des_node,des_node) = len;
      Vinv(des_node,des_node) = 1/len;
    } else
    {
      pA_des_value = pA(des_node);
      p(des_node) = pA_des_value / (1+len*pA_des_value);
      pA(anc_node) += pA(anc_node) + p(des_node);
      
      Vinv(children(des_node,0),children(des_node,0),arma::size(children(des_node,1),children(des_node,1)))
        -= (((len/(1+len*pA_des_value)) * Vinv(children(des_node,0),children(des_node,0),arma::size(children(des_node,1),children(des_node,1))))
              * arma::ones(children(des_node,1),children(des_node,1)) * Vinv(children(des_node,0),children(des_node,0),arma::size(children(des_node,1),children(des_node,1))));
              
              
              children(anc_node,0) = std::min(children(des_node,0),children(anc_node,0));
              children(anc_node,1) += children(des_node,1);
              
              
              Vstar(children(des_node,0),children(des_node,0),arma::size(children(des_node,1),children(des_node,1)))
                +=len;
              
    }
  }
  
  len = len_vec(nedge);
  pA(nedge) = pA(anc_node);
  pA_des_value = pA(nedge);
  p(nedge) = pA(nedge) / (1+len*pA(nedge));
  Vinv(children(nedge,0),children(nedge,0),arma::size(children(nedge,1),children(nedge,1)))
    -= (((len/(1+len*pA_des_value)) * Vinv(children(nedge,0),children(nedge,0),arma::size(children(nedge,1),children(nedge,1))))
          * arma::ones(children(nedge,1),children(nedge,1)) * Vinv(children(nedge,0),children(nedge,0),arma::size(children(nedge,1),children(nedge,1))));
          
          
          if(inv==1)  return Vinv; else if(inv==0) return Vstar;
          return pA;
}

// [[Rcpp::export]]
arma::mat fast_inv2(int nspecies,int nedge,arma::vec len_vec,arma::vec bool_vec,List des_list,arma::vec des_n,arma::vec is_tip,int ret_inv=1)
{
  arma::mat Vstar = arma::zeros(nspecies,nspecies);
  arma::mat Vinv = arma::zeros(nspecies,nspecies);
  arma::vec p = arma::zeros(nedge+1);
  arma::vec pA = arma::zeros(nedge+1);
  int i=0;
  int j=0;
  double len=0;
  double scal=0;
  for(i=0;i<=nedge;i++)
  {
    if(bool_vec(i)==1)
    {
      len = len_vec(i);
      arma::uvec des = des_list(i);
      if(is_tip[i]==1)
      {
        p(i) = 1/len;
      } else
      {
        for(j=0;j<des_n(i);j++)
        {
          pA(i) = pA(i) + p(des(j));
        }
        p(i) = pA(i)/(1+len*pA(i));
      }
      Vstar.submat(des,des) += len;
      arma::mat Ainv = diagmat(len_vec.elem(des));
      scal = arma::as_scalar(len/(1+len*pA(i)));
      if(i<nedge) Vinv.submat(des,des) = Ainv - (scal*Ainv)*(arma::ones<arma::colvec>(des_n(i)))*(arma::ones<arma::rowvec>(des_n(i)))*Ainv;
      //arma::mat Ainv = inv(arma::sympd(Vstar.submat(des,des)));
      Rcout << p(i) << std::endl;
      Rcout << pA(i) << std::endl;
      
      //arma::colvec temp_col = sum(Ainv,1);
      //arma::rowvec temp_row = sum(Ainv);
      //scal = arma::as_scalar(len/(1+(len*accu(Ainv))));
      //arma::mat temp_mat = scal*(temp_col*temp_row);
      //Vinv.submat(des,des) = Ainv - temp_mat;//(arma::as_scalar(len/(1+len*accu(Ainv))));// * (sum(Ainv)*trans(sum(Ainv))));// - (len/(1+len*((t_ones_vec*Ainv)*ones_vec)));// * (Ainv*ones_vec*t_ones_vec*Ainv);
    }
  }
  if(ret_inv==1) return Vinv;
  return Vstar;
}

// [[Rcpp::export]]
arma::mat Rcpp_chol(arma::mat M)
{
  return arma::chol(M);
}

// [[Rcpp::export]]
arma::mat try_clip(arma::mat M,int nvar,int verbose)
{
  arma::mat A(M);
  arma::vec eigval;
  arma::mat eigvec;
  try
  {
    std::ostream nullstream(0);
    arma::set_stream_err2(nullstream);
    
    arma::eig_sym(eigval, eigvec, A);
    unsigned int clip_i = 0;
    if(arma::min(eigval)<1e-12)
    {
      if(verbose==1)
      {
        Rcout << "clipping" << std::endl;
      }
      for(clip_i=0;clip_i<eigval.size();clip_i++)
      {
        if(eigval(clip_i)<1e-12)
        {
          eigval(clip_i)=1e-12;
        }
      }
      arma::mat eigvaldiag = arma::zeros(nvar,nvar);
      eigvaldiag.diag() = eigval;
      A=eigvec*eigvaldiag*eigvec.t();
    }
    return A;
  }
  catch(...)
  {
    return arma::eye(nvar,nvar)+.0123456789;
  }
  return M;
}

// [[Rcpp::export]]
arma::uvec which2(arma::vec x,int y)
{
  //IntegerVector v = Rcpp::seq(0, x.size()-1);
  //IntegerVector ret = v[x==y];
  arma::uvec ret = find(x==y);
  int uc_length = ret.size()+1;
  arma::uvec output(uc_length);
  output(0) = uc_length;
  output(arma::span(1,uc_length-1)) = ret;
  return output;
}

// [[Rcpp::export]]
arma::uvec condition_f(arma::vec iu,int i,int nspecies,arma::mat edge,int nob)
{
  if(i<nspecies)
  {
    arma::uvec uchildren = which2(iu,i);
    return uchildren;
  } else
  {
    //arma::uvec uchildren = which2(edge(_,0),i) + nob;
    arma::uvec uchildren = which2(edge.col(0),i) + nob;
    uchildren(0) = uchildren(0) - nob;
    return uchildren;
  }
  return arma::uvec(1);
}

// [[Rcpp::export]]
arma::mat convert_pars(arma::vec theta,arma::vec options,double T=1)
{
  // options
  // 0 = correlated
  // 1 = error
  // 2 = nvar
  // 3 = clip
  // 4 = verbose
  // 5 = nspecies
  // 6 = nedge
  // 7 = nob
  // 8 = nedge + nob + 1
  int correlated=options(0);
  int error=options(1);
  int nvar=options(2);
  int clip=options(3);
  int verbose=options(4);
  int i=0,j=0,counter=0;
  arma::mat A=arma::zeros(nvar,nvar);
  arma::mat B=arma::zeros(nvar,nvar);
  if(correlated==1)
  {
    for(i=0;i<nvar;i++)
    {
      for(j=0;j<=i;j++)
      {
        A(i,j)=theta(counter);
        counter++;
      }
    }
  } else
  {
    for(i=0;i<nvar;i++)
    {
      A(i,i)=theta(counter);
      counter++;
    }
  }
  if((nvar>1) && (correlated==1))
  {
    A(1,0)=A(1,0)*(exp(A(0,0)));
  }
  A.diag()=exp(A.diag());
  if(correlated==1)
  {
    A = A * A.t();
  }
  if(error==1)
  {
    for(i=0;i<nvar;i++)
    {
      B(i,i)=theta(counter);
      counter++;
    }
    B.diag()=exp(B.diag());
  } else if(error==0)
  {
    for(i=0;i<nvar;i++)
    {
      B(i,i)=0;
    }
  } else if(error==2)
  {
    for(i=0;i<nvar;i++)
    {
      for(j=0;j<=i;j++)
      {
        B(i,j)=theta(counter);
        counter++;
      }
    }
    if((nvar>1))
    {
      B(1,0)=B(1,0)*(exp(B(0,0)));
    }
    B.diag()=exp(B.diag());
    B = B * B.t();
    if(clip==1)
    {
      B = try_clip(B,nvar,verbose);
    }
  }
  if(clip==1)
  {
    A = try_clip(A,nvar,verbose);
  }
  if(arma::max(arma::max(A))>1e3)// || arma::min(A.diag()/T)<1e-6)
  {
    return arma::ones(nvar*2,nvar)+.0123456789;
  }
  /*if(error==1 || error==2)
  {
  if(arma::max(arma::max(B))>1e6 || arma::min(B.diag())<1e-32)
  {
  return arma::ones(nvar*2,nvar)+.0123456789;
  }
  }*/
  
  arma::mat ret=arma::zeros(nvar*2,nvar);
  ret(arma::span(0,nvar-1),arma::span(0,nvar-1)) = A;
  ret(arma::span(nvar,nvar*2-1),arma::span(0,nvar-1)) = B;
  return(ret);
}

// [[Rcpp::export]]
arma::mat convert_pars2(arma::vec theta,arma::vec options,arma::mat fixed_phylocov,arma::mat fixed_phenocov,double tree_height=1,int check1=1)
{
  // options
  // 0 = correlated
  // 1 = error
  // 2 = nvar
  // 3 = clip
  // 4 = verbose
  // 5 = nspecies
  // 6 = nedge
  // 7 = nob
  // 8 = nedge + nob + 1
  int correlated=options(0);
  int error=options(1);
  int nvar=options(2);
  int clip=options(3);
  int verbose=options(4);
  int estim_phylocov = options(12);
  int estim_phenocov = options(13);
  int i=0,j=0,counter=0;
  arma::mat A=arma::zeros(nvar,nvar);
  arma::mat B=arma::zeros(nvar,nvar);
  
  if(estim_phylocov==0)
  {
    A = fixed_phylocov;
  } else
  {
    if(correlated==1)
    {
      for(i=0;i<nvar;i++)
      {
        for(j=0;j<=i;j++)
        {
          A(i,j)=theta(counter);
          counter++;
        }
      }
    } else
    {
      for(i=0;i<nvar;i++)
      {
        A(i,i)=theta(counter);
        counter++;
      }
    }
    if((nvar>1) && (correlated==1))
    {
      A(1,0)=A(1,0)*(exp(A(0,0)));
    }
    A.diag()=exp(A.diag());
    if(correlated==1)
    {
      A = A * A.t();
    }
  }
  
  if(error==0)
  {
    for(i=0;i<nvar;i++)
    {
      B(i,i)=0;
    }
  } else
  {
    if(estim_phenocov==0)
    {
      B = fixed_phenocov;
    } else
    {
      if(error==1)
      {
        for(i=0;i<nvar;i++)
        {
          B(i,i)=theta(counter);
          counter++;
        }
        B.diag()=exp(B.diag());
      } else if(error==2)
      {
        for(i=0;i<nvar;i++)
        {
          for(j=0;j<=i;j++)
          {
            B(i,j)=theta(counter);
            counter++;
          }
        }
        if((nvar>1))
        {
          B(1,0)=B(1,0)*(exp(B(0,0)));
        }
        B.diag()=exp(B.diag());
        B = B * B.t();
      }
    }
  }
  if(clip==1 & error>0)
  {
    B = try_clip(B,nvar,verbose);
  }
  if(clip==1)
  {
    A = try_clip(A,nvar,verbose);
  }
  //if(estim_phylocov==1) A.diag()+=1e-6;
  //if((error>0) & (estim_phenocov==1)) B.diag()+=1e-6;
  if(check1==1)
  {
    if(((arma::max(arma::max(A))*tree_height)<.01) || ((arma::max(arma::max(A))*tree_height)>1e2) || (arma::max(arma::max(B))>1e3))// || ((arma::min(A.diag())*tree_height)<1e-6))
    {
      return arma::ones(nvar*2,nvar)+.0123456789;
    }
  }
  if(error>0)
  {
    if(arma::max(arma::max(B))>1e3 || arma::min(B.diag())<1e-6)
    {
      return arma::ones(nvar*2,nvar)+.0123456789;
    }
  }
  
  arma::mat ret=arma::zeros(nvar*2,nvar);
  ret(arma::span(0,nvar-1),arma::span(0,nvar-1)) = A;
  ret(arma::span(nvar,nvar*2-1),arma::span(0,nvar-1)) = B;
  return(ret);
}

// [[Rcpp::export]]
double threepoint(arma::vec theta,arma::vec options,arma::vec y,arma::vec ku,arma::vec iu,arma::mat edge,arma::vec edgelength,Rcpp::List uchildren_list,arma::vec models,arma::uvec externalEdge,arma::uvec not_externalEdge,arma::vec dist_anc,arma::vec dist_des,double Tmax,double Tmin,int nmodels,arma::vec lower_bounds,arma::vec upper_bounds,arma::vec OU_D,arma::vec times)
{
  // options
  // 0 = correlated
  // 1 = error
  // 2 = nvar
  // 3 = clip
  // 4 = verbose
  // 5 = nspecies
  // 6 = nedge
  // 7 = nob
  // 8 = nedge + nob + 1
  
  int nspecies = options(5);
  int nedge = options(6);
  
  double rootedge = 0;
  bool OU_mod = ((models(0) == 1) | (models(1) == 1));
  bool OUfixedRoot = models(0) == 1;
  bool lambda_mod = models(2) == 1;
  bool kappa_mod = models(3) == 1;
  bool delta_mod = models(4) == 1;
  bool EB_mod = models(5) == 1;
  bool pars_from_theta = true;
  double alpha = 0;
  double lambda = 0;
  double kappa = 0;
  double delta = 0;
  double rate = 0;
  int par_count = 0;
  int npars = options(11);
  if((unsigned int)(npars)==theta.size())
  {
    pars_from_theta = false;
  }
  if(OU_mod)
  {
    arma::ucolvec anc = arma::conv_to<arma::uvec>::from(edge.col(0));
    arma::ucolvec des = arma::conv_to<arma::uvec>::from(edge.col(1));
    if(pars_from_theta)
    {
      alpha = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      alpha = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    par_count++;
    if(alpha>1e-7/Tmax)
    {
      if(!OUfixedRoot)
      {
        arma::vec distFromRoot = exp(-2*alpha*times);
        arma::vec d1 = distFromRoot(anc-nspecies);
        arma::vec d2 = arma::zeros(nedge);
        d2(externalEdge) = exp(-2*alpha*OU_D(des(externalEdge)));
        d2(not_externalEdge) = distFromRoot(des(not_externalEdge)-nspecies);
        edgelength = d2 - d1;
        rootedge = arma::min(distFromRoot);
      } else
      {
        arma::vec distFromRoot = exp(-2*alpha*times)%(1 - exp(-2*alpha*(Tmax-times)));
        arma::vec d1 = distFromRoot(anc-nspecies);
        arma::vec d2 = arma::zeros(nedge);
        d2(externalEdge) = exp(-2*alpha*OU_D(des(externalEdge))) % (1-exp(-2*alpha*(Tmax-OU_D(des(externalEdge)))));
        d2(not_externalEdge) = distFromRoot(des(not_externalEdge)-nspecies);
        edgelength = d2 - d1;
        rootedge = arma::min(distFromRoot);
      }
    }
  }
  if(lambda_mod)
  {
    if(pars_from_theta)
    {
      lambda = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      lambda = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    if((lambda<1) & (lambda>1e-7))
    {
      edgelength = edgelength * lambda;
      edgelength(externalEdge) = edgelength(externalEdge) + (1-lambda)*dist_des(externalEdge);
    }
    par_count++;
  }
  if(kappa_mod)
  {
    if(pars_from_theta)
    {
      kappa = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      kappa = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    if((kappa<1) & (kappa>1e-6))
    {
      edgelength = exp(kappa*log(edgelength));
    }
    par_count++;
  }
  if(delta_mod)
  {
    if(pars_from_theta)
    {
      delta = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      delta = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    if(delta>1e-5)
    {
      edgelength = (exp(delta*log(dist_des)) - exp(delta*log(dist_anc)))*exp((1-delta)*log(Tmax));
    }
    par_count++;
  }
  if(EB_mod)
  {
    if(pars_from_theta)
    {
      rate = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      rate = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    if(rate!=0)
    {
      edgelength = (exp(rate*dist_des)-exp(rate*dist_anc))/rate;
    }
    par_count++;
  }
  
  int nvar = options(2);
  int nob = options(7);
  int nn = options(8);
  int verbose = options(4);
  int REML = options(10);
  int ret_YY = options(13);
  arma::mat AB = convert_pars(theta,options,Tmin);
  if(AB(0,0)==1.0123456789)
  {
    //Rcout << "pars conversion failed";
    return 1e128;
  }
  
  arma::mat phylocov = AB(0,0,arma::size(nvar,nvar));
  arma::mat phenocov = AB(nvar,0,arma::size(nvar,nvar));
  arma::mat itpainv = arma::zeros(nvar,nvar);
  arma::mat inverse = arma::zeros(nvar,nvar);
  
  arma::mat XU = arma::zeros(nvar*nn,nvar);
  arma::mat XY = arma::zeros(nvar*nn,1);
  arma::mat YU = arma::zeros(nn,nvar);
  arma::mat XX = arma::zeros(nvar*nn,nvar);
  arma::vec YY = arma::zeros(nn);
  arma::vec lnW = arma::zeros(nn);
  arma::mat YU_XY = arma::zeros(1,nvar);
  int u = 0;
  int k = 0;
  int e = 0;
  int uc_length = 0;
  int uci = 0;
  double len = 0;
  double determinant = 0;
  for(u=0;u<nob;u++)
  {
    k = ku(u);
    XU.row((u*nvar)+k) = phylocov.row(k) / phenocov(k,k);
    XY((u*nvar)+k,0) = y(u) / phenocov(k,k);
    XX((u*nvar)+k,k) = 1 / phenocov(k,k);
    YU.row(u) = y(u) * phylocov.row(k) / phenocov(k,k);
    YY(u) = y(u) * y(u) / phenocov(k,k);
    lnW(u) = log(phenocov(k,k));
  }
  for(e=0;e<(nedge+1);e++) // e=nedge+1 is to add the root edge
  {
    u = nob + e; // index of edge after observations were added as polytomies
    /*if(e==nedge)
    {
    i = nspecies;
    } else
    {
    i = edge(e,1);
    }*/
    arma::uvec uchildren = uchildren_list(e);
    uc_length = uchildren(0);
    
    for(uci=1;uci<uc_length;uci++)
    {
      XU((u*nvar),0,arma::size(nvar,nvar)) = XU((u*nvar),0,arma::size(nvar,nvar)) + XU((uchildren[uci]*nvar),0,arma::size(nvar,nvar));
      XY((u*nvar),0,arma::size(nvar,1)) = XY((u*nvar),0,arma::size(nvar,1)) + XY(uchildren[uci]*nvar,0,arma::size(nvar,1));
      XX((u*nvar),0,arma::size(nvar,nvar)) = XX((u*nvar),0,arma::size(nvar,nvar)) + XX((uchildren[uci]*nvar),0,arma::size(nvar,nvar));
      YU((u),0,arma::size(1,nvar)) = YU((u),0,arma::size(1,nvar)) + YU(uchildren[uci],0,arma::size(1,nvar));    
      YY(u) = YY(u) + YY(uchildren(uci));
      lnW(u) = lnW(u) + lnW(uchildren(uci));
    }
    if(e==nedge)
    {
      len = rootedge;
    } else
    {
      len = edgelength(e);
    }
    itpainv = arma::eye(nvar,nvar) + (XU((u*nvar),0,arma::size(nvar,nvar))*len);
    inverse = try_inv(itpainv,nvar);
    if(inverse(0,0)==1.0123456789)
    {
      return 1e128;
    }
    determinant = det(itpainv);
    lnW(u) = lnW(u) + log(determinant);
    XY((u*nvar),0,arma::size(nvar,1)) = inverse * XY((u*nvar),0,arma::size(nvar,1));
    XX((u*nvar),0,arma::size(nvar,nvar)) = inverse * XX((u*nvar),0,arma::size(nvar,nvar));
    arma::mat YU_XY = (YU((u),0,arma::size(1,nvar))*XY((u*nvar),0,arma::size(nvar,1)));
    YY(u) = YY(u) - (len * YU_XY(0,0));
    YU((u),0,arma::size(1,nvar)) = YU((u),0,arma::size(1,nvar)) * inverse;
    XU((u*nvar),0,arma::size(nvar,nvar)) = XU((u*nvar),0,arma::size(nvar,nvar)) * inverse;
  }
  if(ret_YY==1)
  {
    return YY(u);
  }
  
  arma::mat anc = try_solve(XX((u*nvar),0,arma::size(nvar,nvar)),XY((u*nvar),0,arma::size(nvar,1)));
  
  //int i=0;
  /*for(i=0;i<nvar;i++)
  {
  if(anc(i,0)>1e3 || anc(i,0)<-1e3)
  {
  //Rcout << "anc out of range" << std::endl;
  return 1e128;
  }
  }*/
  
  arma::mat anc_XY = anc.t() * XY((u*nvar),0,arma::size(nvar,1));
  arma::mat anc_XX = anc.t() * XX((u*nvar),0,arma::size(nvar,nvar)) * anc;
  
  double minus2ll = 0;
  if(REML==1)
  {
    minus2ll = (log(det(XX((u*nvar),0,arma::size(nvar,nvar)))) + lnW(u) + YY(u) - (2 * anc_XY(0,0)) + anc_XX(0,0)) / 2;
  } else
  {
    minus2ll = (lnW(u) + YY(u) - (2 * anc_XY(0,0)) + anc_XX(0,0)) / 2;
  }
  
  /*if(minus2ll<-1e3)
  {
  return 1e128;
  }*/
  
  if(verbose==1)
  {
    Rcout << minus2ll << std::endl;
  }
  return minus2ll;
  }

// [[Rcpp::export]]
double threepoint_nopheno(arma::vec theta,arma::vec options,arma::vec y,arma::vec ku,arma::vec iu,arma::mat edge,arma::vec edgelength,Rcpp::List uchildren_list,arma::vec edgevec,List subset_list,List species_subset,arma::vec tip_combn,arma::vec ymin,arma::vec ymax,arma::vec models,arma::uvec externalEdge,arma::uvec not_externalEdge,arma::vec dist_anc,arma::vec dist_des,double Tmax,double Tmin,int nmodels,arma::vec lower_bounds,arma::vec upper_bounds,arma::vec OU_D,arma::vec times)
{
  // edgevec is a vector of tip lengths
  // subset_list is a list of unique matrices to invert once
  // species_subset is a list of length nspecies
  // tip_combn is a vector signifyinc which matrix subset corresponds to each species
  // ymin, ymax
  
  // options
  // 0 = correlated
  // 1 = error
  // 2 = nvar
  // 3 = clip
  // 4 = verbose
  // 5 = nspecies
  // 6 = nedge
  // 7 = nob
  // 8 = nedge + nob + 1
  
  int nspecies = options(5);
  int nedge = options(6);
  
  double rootedge = 0;
  bool OU_mod = ((models(0) == 1) | (models(1) == 1));
  bool OUfixedRoot = models(0) == 1;
  bool lambda_mod = models(2) == 1;
  bool kappa_mod = models(3) == 1;
  bool delta_mod = models(4) == 1;
  bool EB_mod = models(5) == 1;
  bool pars_from_theta = true;
  double alpha = 0;
  double lambda = 0;
  double kappa = 0;
  double delta = 0;
  double rate = 0;
  int par_count = 0;
  int npars = options(11);
  if((unsigned int)(npars)==theta.size())
  {
    pars_from_theta = false;
  }
  if(OU_mod)
  {
    arma::ucolvec anc = arma::conv_to<arma::uvec>::from(edge.col(0));
    arma::ucolvec des = arma::conv_to<arma::uvec>::from(edge.col(1));
    
    if(pars_from_theta)
    {
      alpha = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      alpha = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    if(!arma::is_finite(alpha)) return 1e128;
    par_count++;
    //if(alpha>1e-7/Tmax)
    //{
    if(!OUfixedRoot)
    {
      arma::vec distFromRoot = exp(-2*alpha*times);
      arma::vec d1 = distFromRoot(anc-nspecies);
      arma::vec d2 = arma::zeros(nedge);
      d2(externalEdge) = exp(-2*alpha*OU_D(des(externalEdge)));
      d2(not_externalEdge) = distFromRoot(des(not_externalEdge)-nspecies);
      edgelength = d2 - d1;
      rootedge = arma::min(distFromRoot);
    } else
    {
      arma::vec distFromRoot = exp(-2*alpha*times);//%(1 - exp(-2*alpha*(Tmax-times)));
      for(unsigned int looper=0;looper<distFromRoot.size();looper++)
      {
        distFromRoot(looper) *= (1-exp(-2*alpha*(Tmax-times(looper))));
      }
      arma::vec d1 = distFromRoot(anc-nspecies);
      arma::vec d2 = arma::zeros(nedge);
      d2(externalEdge) = exp(-2*alpha*OU_D(des(externalEdge)));// % (1-exp(-2*alpha*(Tmax-OU_D(des(externalEdge)))));
      for(unsigned int looper=0;looper<externalEdge.size();looper++)
      {
        d2(externalEdge(looper)) *= (1-exp(-2*alpha*Tmax));
      }
      d2(not_externalEdge) = distFromRoot(des(not_externalEdge)-nspecies);
      edgelength = d2 - d1;
      rootedge = arma::min(distFromRoot);
    }
    //}
  }
  if(lambda_mod)
  {
    if(pars_from_theta)
    {
      lambda = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      lambda = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    //if((lambda<1) & (lambda>1e-7))
    //{
    edgelength = edgelength * lambda;
    edgelength(externalEdge) = edgelength(externalEdge) + (1-lambda)*dist_des(externalEdge);
    //}
    par_count++;
  }
  if(kappa_mod)
  {
    if(pars_from_theta)
    {
      kappa = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      kappa = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    //if((kappa<1) & (kappa>1e-6))
    //{
    edgelength = exp(kappa*log(edgelength));
    //}
    par_count++;
  }
  if(delta_mod)
  {
    if(pars_from_theta)
    {
      delta = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      delta = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    edgelength = (exp(delta*log(dist_des)) - exp(delta*log(dist_anc)))*exp((1-delta)*log(Tmax));  
    par_count++;
  }
  if(EB_mod)
  {
    if(pars_from_theta)
    {
      rate = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      rate = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    if(rate!=0)
    {
      edgelength = (exp(rate*dist_des)-exp(rate*dist_anc))/rate;
    }
    par_count++;
  }
  if(par_count>0)
  {
    edgevec = edgelength(externalEdge);
  }
  
  int nvar = options(2);
  int nn = options(8);
  int verbose = options(4);
  int nind = options(9);
  int REML = options(10);
  int ret_YY = options(13);
  arma::mat AB = convert_pars(theta,options,Tmin);
  if(AB(0,0)==1.0123456789)
  {
    //Rcout << "pars conversion failed";
    return 1e128;
  }
  arma::mat phylocov = AB(0,0,arma::size(nvar,nvar));
  arma::mat itpainv = arma::zeros(nvar,nvar);
  arma::mat inverse = arma::zeros(nvar,nvar);
  
  List sub_mats = inv_subset(phylocov,subset_list);
  
  arma::mat XU = arma::zeros(nvar*nn,nvar);
  arma::mat XY = arma::zeros(nvar*nn,1);
  arma::mat YU = arma::zeros(nn,nvar);
  arma::mat XX = arma::zeros(nvar*nn,nvar);
  arma::vec YY = arma::zeros(nn);
  arma::vec lnW = arma::zeros(nn);
  arma::mat YU_XY = arma::zeros(1,nvar);
  int u = 0;
  int e = 0;
  int uc_length = 0;
  int uci = 0;
  double len = 0;
  double determinant = 0;
  arma::uvec zero_vec(1);
  zero_vec(0)=0;
  for(u=0;u<nind;u++)
  {
    len = edgevec(u);
    arma::uvec Ka = species_subset(u);
    
    arma::mat Bainv = sub_mats(tip_combn(u));
    /*if(Bainv(0,0)==1.0123456789)
    {
    Rcout << "Ba inverse failed";
    return 1e128;
    }*/
    Bainv = Bainv / len;
    //int ma = y.size();
    arma::mat Ma = Bainv*phylocov.rows(Ka);
    
    XU.rows((u*nvar)+Ka) = Ma;
    XY.submat((u*nvar)+Ka,zero_vec) = Bainv * y.subvec(ymin(u),ymax(u));
    XX.submat((u*nvar)+Ka,Ka) = Bainv;
    YU.row(u) = trans(y.subvec(ymin(u),ymax(u))) * Ma;
    arma::mat temp_mat = (trans(y.subvec(ymin(u),ymax(u))) * Bainv * y.subvec(ymin(u),ymax(u)));
    YY(u) = temp_mat(0,0);
    lnW(u) = log(det(phylocov.submat(Ka,Ka)*len));
  }
  int i = 0;
  for(e=0;e<(nedge+1);e++) // e=nedge+1 is to add the root edge
  {
    u = nind + e; // index of edge after observations were added as polytomies
    if(e==nedge)
    {
      i = nspecies;
    } else
    {
      i = edge(e,1) + 1;
    }
    arma::uvec uchildren = uchildren_list(e);
    uc_length = uchildren(0);
    for(uci=1;uci<uc_length;uci++)
    {
      XU((u*nvar),0,arma::size(nvar,nvar)) = XU((u*nvar),0,arma::size(nvar,nvar)) + XU((uchildren[uci]*nvar),0,arma::size(nvar,nvar));
      XY((u*nvar),0,arma::size(nvar,1)) = XY((u*nvar),0,arma::size(nvar,1)) + XY(uchildren[uci]*nvar,0,arma::size(nvar,1));
      XX((u*nvar),0,arma::size(nvar,nvar)) = XX((u*nvar),0,arma::size(nvar,nvar)) + XX((uchildren[uci]*nvar),0,arma::size(nvar,nvar));
      YU((u),0,arma::size(1,nvar)) = YU((u),0,arma::size(1,nvar)) + YU(uchildren[uci],0,arma::size(1,nvar));    
      YY(u) = YY(u) + YY(uchildren(uci));
      lnW(u) = lnW(u) + lnW(uchildren(uci));
    }
    if(e==nedge)
    {
      len = rootedge;
    } else if(i<=nspecies)
    {
      len = 0;
    } else
    {
      len = edgelength(e);
    }
    itpainv = arma::eye(nvar,nvar) + (XU((u*nvar),0,arma::size(nvar,nvar))*len);
    inverse = try_inv(itpainv,nvar);
    if(inverse(0,0)==1.0123456789)
    {
      return 1e128;
    }
    determinant = det(itpainv);
    lnW(u) = lnW(u) + log(determinant);
    XY((u*nvar),0,arma::size(nvar,1)) = inverse * XY((u*nvar),0,arma::size(nvar,1));
    XX((u*nvar),0,arma::size(nvar,nvar)) = inverse * XX((u*nvar),0,arma::size(nvar,nvar));
    arma::mat YU_XY = (YU((u),0,arma::size(1,nvar))*XY((u*nvar),0,arma::size(nvar,1)));
    YY(u) = YY(u) - (len * YU_XY(0,0));
    YU((u),0,arma::size(1,nvar)) = YU((u),0,arma::size(1,nvar)) * inverse;
    XU((u*nvar),0,arma::size(nvar,nvar)) = XU((u*nvar),0,arma::size(nvar,nvar)) * inverse;
  }
  if(ret_YY==1)
  {
    return YY(u);
  }
  
  arma::mat anc = try_solve(XX((u*nvar),0,arma::size(nvar,nvar)),XY((u*nvar),0,arma::size(nvar,1)));
  /*for(i=0;i<nvar;i++)
  {
  if(anc(i,0)>1e3 || anc(i,0)<-1e3)
  {
  //Rcout << "anc out of range" << std::endl;
  return 1e128;
  }
  }*/
  arma::mat anc_XY = anc.t() * XY((u*nvar),0,arma::size(nvar,1));
  arma::mat anc_XX = anc.t() * XX((u*nvar),0,arma::size(nvar,nvar)) * anc;
  double minus2ll = 0;
  if(REML==1)
  {
    minus2ll = (log(det(XX((u*nvar),0,arma::size(nvar,nvar)))) + lnW(u) + YY(u) - (2 * anc_XY(0,0)) + anc_XX(0,0)) / 2;
  } else
  {
    minus2ll = (lnW(u) + YY(u) - (2 * anc_XY(0,0)) + anc_XX(0,0)) / 2;
  }
  
  /*if(minus2ll<-1e3)
  {
  return 1e128;
  }*/
  if(verbose==1)
  {
    Rcout << minus2ll << std::endl;
  }
  return minus2ll;
  //return List::create(XU,XX,YU,YY,XY,lnW);
  }

// [[Rcpp::export]]
double threepoint2(arma::vec theta,arma::vec options,arma::vec y,
                   arma::mat edge,arma::vec edgelength,List subset_list,
                   List species_subset,arma::vec tip_combn,arma::vec ymin,
                   arma::vec ymax,arma::vec ind_edge,arma::vec anc_edge,
                   arma::mat fixed_phylocov,arma::mat fixed_phenocov,
                   List inv_phenocovs)
{
  // edgevec is a vector of tip lengths
  // subset_list is a list of unique matrices to invert once
  // species_subset is a list of length nspecies
  // tip_combn is a vector signifyinc which matrix subset corresponds to each species
  // ymin, ymax
  
  // options
  // 0 = correlated
  // 1 = error
  // 2 = nvar
  // 3 = clip
  // 4 = verbose
  // 5 = nspecies
  // 6 = nedge
  // 7 = nob
  // 8 = nedge + nob + 1
  
  int nspecies = options(5);
  int nedge = options(6);
  
  double rootedge = 0;
  
  int nvar = options(2);
  int nn = options(8);
  int verbose = options(4);
  int nind = options(9);
  int REML = options(10);
  arma::mat AB = convert_pars2(theta,options,fixed_phylocov,fixed_phenocov,options(16));
  if(AB(0,0)==1.0123456789)
  {
    //Rcout << "pars conversion failed";
    return 1e128;
  }
  arma::mat phylocov = AB(0,0,arma::size(nvar,nvar));
  arma::mat phenocov = AB(nvar,0,arma::size(nvar,nvar));
  arma::mat itpa = arma::zeros(nvar,nvar);
  arma::mat itpainv = arma::zeros(nvar,nvar);
  
  List sub_mats(subset_list.size());
  if((options(15)==0) & (options(1)>0))
  {
    sub_mats = inv_subset(phenocov,subset_list);
  } else if(options(1)==0)
  {
    sub_mats = inv_subset(phylocov,subset_list);
  }
  arma::mat XU = arma::zeros(nvar*(nedge+1),nvar); // this is p
  arma::mat XY = arma::zeros(nvar*(nedge+1),1); // this is Vr (XW-1Y)
  //arma::mat YU = arma::zeros(1,nvar*(nedge+1)); // this is Ul (YW-1U)
  arma::mat XX = arma::zeros(nvar*(nedge+1),nvar); // this is Vr (XW-1X)
  arma::vec YY = arma::zeros(nedge+1); // this is Q (YW-1Y)
  arma::vec lnW = arma::zeros(nedge+1);
  //arma::mat YU_XY = arma::zeros(1,nvar);
  arma::mat pA = arma::zeros(nvar,nvar);
  arma::mat I = arma::eye(nvar,nvar);
  int i = 0;
  int u = 0;
  int anc = 0;
  int e = 0;
  int uc_length = 0;
  int uci = 0;
  double len = 0;
  double determinant = 0;
  arma::uvec zero_vec(1);
  zero_vec(0)=0;
  unsigned int tip_edge_num = 0;
  // initialization
  for(u=0;u<nind;u++)
  {
    tip_edge_num = ind_edge(u);
    arma::uvec Ka = species_subset(u);
    arma::mat Bainv = arma::zeros(Ka.size(),Ka.size());
    if(options(15)==0)
    {
      Bainv = Rcpp::as<arma::mat>(sub_mats(tip_combn(u)));
      if(Bainv(0,0)==1.0123456789)
      {
        //Rcout << "Ba inverse failed";
        return 1e128;
      }
    } else
    {
      Bainv = Rcpp::as<arma::mat>(inv_phenocovs(u));
    }
    if(options(1)==0) Bainv = Bainv/edgelength(tip_edge_num);
    //int ma = y.size();
    //XU.rows((tip_edge_num*nvar)+Ka) += Ma;
    XX.submat((tip_edge_num*nvar)+Ka,Ka) += Bainv; // this is p = X'W-1X
    XY.submat((tip_edge_num*nvar)+Ka,zero_vec) += Bainv * y.subvec(ymin(u),ymax(u)); // this is V = X'W-1Y
    //YU.submat(zero_vec,(tip_edge_num*nvar)+Ka) += trans(y.subvec(ymin(u),ymax(u))) * Bainv; // this is U = Y'W-1X
    arma::mat temp_mat = trans(y.subvec(ymin(u),ymax(u))) * Bainv * y.subvec(ymin(u),ymax(u)); // this is Q = Y'W-1Y
    YY(tip_edge_num) += temp_mat(0,0);
    lnW(tip_edge_num) -= logdet(Bainv);
  }
  
  // recursion
  for(e=0;e<nedge;e++)
  {
    u = e; // index of des edge after observations were added as polytomies
    anc = anc_edge(e);
    //anc = edge(e,0); // index of anc edge
    
    
    if(e==nedge)
    {
      i = nspecies;
    } else
    {
      i = edge(e,1) + 1;
    }
    if(e==nedge)
    {
      len = rootedge;
    } else if(i<=nspecies)
    {
      len = 0;
    } else
    {
      len = edgelength(e);
    }
    
    pA = XX((u*nvar),0,arma::size(nvar,nvar)); // this is sum_s (X'W-1X)_s
    itpa = I + (len*phylocov)*pA; // this is 
    itpainv = try_inv(itpa,nvar);
    
    if(itpainv(0,0)==1.0123456789)
    {
      return 1e128;
    }
    
    // logdet(I+T*pA)
    lnW(u) += logdet(itpa);
    lnW(anc) += lnW(u);
    // XX is p = X'W-1X
    XX((u*nvar),0,arma::size(nvar,nvar)) = pA * itpainv; // this is p
    XX((anc*nvar),0,arma::size(nvar,nvar)) += XX((u*nvar),0,arma::size(nvar,nvar)); // this adds to pA of the anc
    
    // YY is Q = Y'W-1Y
    arma::mat YU_XY = (trans(XY((u*nvar),0,arma::size(nvar,1)))*itpainv*(len*phylocov)*XY((u*nvar),0,arma::size(nvar,1)));
    YY(u) = YY(u) - (YU_XY(0,0));
    YY(anc) += YY(u);
    
    // XY is V = X'W-1Y
    XY((u*nvar),0,arma::size(nvar,1)) = trans(trans(XY((u*nvar),0,arma::size(nvar,1)))*itpainv); // this is V = X'W-1Y
    XY((anc*nvar),0,arma::size(nvar,1)) += XY((u*nvar),0,arma::size(nvar,1)); // this adds V to anc
    
    // YU is U = Y'W-1X
    //YU(0,(u*nvar),arma::size(1,nvar)) = YU(0,(u*nvar),arma::size(1,nvar)) * itpainv;
    //YU(0,(anc*nvar),arma::size(1,nvar)) += YU(0,(u*nvar),arma::size(1,nvar));
  }
  u = nedge;
  len = rootedge;
  pA = XX((u*nvar),0,arma::size(nvar,nvar)); // this is sum_s (X'W-1X)_s
  itpa = I + (len*phylocov)*pA; // this is 
  itpainv = try_inv(itpa,nvar);
  lnW(u) += logdet(itpa);
  //Rcout << XX((u*nvar),0,arma::size(nvar,nvar)) << std::endl;
  //Rcout << XY((u*nvar),0,arma::size(nvar,1));
  
  XX((u*nvar),0,arma::size(nvar,nvar)) = XX((u*nvar),0,arma::size(nvar,nvar)) * itpainv; // this is p
  
  arma::mat YU_XY = (trans(XY((u*nvar),0,arma::size(nvar,1)))*itpainv*(len*phylocov)*XY((u*nvar),0,arma::size(nvar,1)));
  YY(u) = (YY(u) - (YU_XY(0,0)));
  
  XY((u*nvar),0,arma::size(nvar,1)) = (trans(trans(XY((u*nvar),0,arma::size(nvar,1)))*itpainv)); // this is V = X'W-1Y
  
  //YU(0,(u*nvar),arma::size(1,nvar)) = YU(0,(u*nvar),arma::size(1,nvar)) * itpainv;
  arma::mat root_anc = try_solve(XX((u*nvar),0,arma::size(nvar,nvar)),XY((u*nvar),0,arma::size(nvar,1)));
  arma::mat anc_XY = root_anc.t() * XY((u*nvar),0,arma::size(nvar,1));
  arma::mat anc_XX = root_anc.t() * XX((u*nvar),0,arma::size(nvar,nvar)) * root_anc;
  double minus2ll = 0;
  if(REML==1)
  {
    minus2ll = (logdet(XX((u*nvar),0,arma::size(nvar,nvar))) + lnW(u) + YY(u) - (2 * anc_XY(0,0)) + anc_XX(0,0)) / 2;
  } else
  {
    minus2ll = (lnW(u) + YY(u) - (2 * anc_XY(0,0)) + anc_XX(0,0)) / 2;
  }
  
  if(verbose==1)
  {
    Rcout << minus2ll << std::endl;
  }
  /*if(minus2ll<-1e3)
  {
  return 1e128;
  }*/
  return minus2ll;
  //return List::create(XU,XX,YU,YY,XY,lnW);
}

// [[Rcpp::export]]
List threepoint3(arma::vec theta,arma::vec options,arma::vec y,
                 arma::mat edge,arma::vec edgelength,List subset_list,
                 List species_subset,arma::vec tip_combn,arma::vec ymin,
                 arma::vec ymax,arma::vec ind_edge,arma::vec anc_edge,
                 arma::mat fixed_phylocov,arma::mat fixed_phenocov,
                 List inv_phenocovs)
{
  // edgevec is a vector of tip lengths
  // subset_list is a list of unique matrices to invert once
  // species_subset is a list of length nspecies
  // tip_combn is a vector signifyinc which matrix subset corresponds to each species
  // ymin, ymax
  
  // options
  // 0 = correlated
  // 1 = error
  // 2 = nvar
  // 3 = clip
  // 4 = verbose
  // 5 = nspecies
  // 6 = nedge
  // 7 = nob
  // 8 = nedge + nob + 1
  
  int nspecies = options(5);
  int nedge = options(6);
  
  double rootedge = 0;
  
  int nvar = options(2);
  int nn = options(8);
  int verbose = options(4);
  int nind = options(9);
  int REML = options(10);
  arma::mat AB = convert_pars2(theta,options,fixed_phylocov,fixed_phenocov,options(16));
  if(AB(0,0)==1.0123456789)
  {
    //Rcout << "pars conversion failed";
    return 1e128;
  }
  arma::mat phylocov = AB(0,0,arma::size(nvar,nvar));
  arma::mat phenocov = AB(nvar,0,arma::size(nvar,nvar));
  arma::mat itpa = arma::zeros(nvar,nvar);
  arma::mat itpainv = arma::zeros(nvar,nvar);
  Rcout << phylocov << std::endl << std::endl << phenocov;
  List sub_mats(subset_list.size());
  if((options(15)==0) & (options(1)>0))
  {
    sub_mats = inv_subset(phenocov,subset_list);
  } else if(options(1)==0)
  {
    sub_mats = inv_subset(phylocov,subset_list);
  }
  arma::mat XU = arma::zeros(nvar*(nedge+1),nvar); // this is p
  arma::mat XY = arma::zeros(nvar*(nedge+1),1); // this is Vr (XW-1Y)
  //arma::mat YU = arma::zeros(1,nvar*(nedge+1)); // this is Ul (YW-1U)
  arma::mat XX = arma::zeros(nvar*(nedge+1),nvar); // this is Vr (XW-1X)
  arma::vec YY = arma::zeros(nedge+1); // this is Q (YW-1Y)
  arma::vec lnW = arma::zeros(nedge+1);
  //arma::mat YU_XY = arma::zeros(1,nvar);
  arma::mat pA = arma::zeros(nvar,nvar);
  arma::mat I = arma::eye(nvar,nvar);
  int i = 0;
  int u = 0;
  int anc = 0;
  int e = 0;
  int uc_length = 0;
  int uci = 0;
  double len = 0;
  double determinant = 0;
  arma::uvec zero_vec(1);
  zero_vec(0)=0;
  unsigned int tip_edge_num = 0;
  // initialization
  for(u=0;u<nind;u++)
  {
    tip_edge_num = ind_edge(u);
    arma::uvec Ka = species_subset(u);
    arma::mat Bainv = arma::zeros(Ka.size(),Ka.size());
    if(options(15)==0)
    {
      Bainv = Rcpp::as<arma::mat>(sub_mats(tip_combn(u)));
      if(Bainv(0,0)==1.0123456789)
      {
        //Rcout << "Ba inverse failed";
        return 1e128;
      }
    } else
    {
      Bainv = Rcpp::as<arma::mat>(inv_phenocovs(u));
    }
    if(options(1)==0) Bainv = Bainv/edgelength(tip_edge_num);
    //int ma = y.size();
    //XU.rows((tip_edge_num*nvar)+Ka) += Ma;
    XX.submat((tip_edge_num*nvar)+Ka,Ka) += Bainv; // this is p = X'W-1X
    XY.submat((tip_edge_num*nvar)+Ka,zero_vec) += Bainv * y.subvec(ymin(u),ymax(u)); // this is V = X'W-1Y
    //YU.submat(zero_vec,(tip_edge_num*nvar)+Ka) += trans(y.subvec(ymin(u),ymax(u))) * Bainv; // this is U = Y'W-1X
    arma::mat temp_mat = trans(y.subvec(ymin(u),ymax(u))) * Bainv * y.subvec(ymin(u),ymax(u)); // this is Q = Y'W-1Y
    YY(tip_edge_num) += temp_mat(0,0);
    lnW(tip_edge_num) -= logdet(Bainv);
  }
  
  // recursion
  for(e=0;e<nedge;e++)
  {
    u = e; // index of des edge after observations were added as polytomies
    anc = anc_edge(e);
    //anc = edge(e,0); // index of anc edge
    
    
    if(e==nedge)
    {
      i = nspecies;
    } else
    {
      i = edge(e,1) + 1;
    }
    if(e==nedge)
    {
      len = rootedge;
    } else if(i<=nspecies)
    {
      len = 0;
    } else
    {
      len = edgelength(e);
    }
    
    pA = XX((u*nvar),0,arma::size(nvar,nvar)); // this is sum_s (X'W-1X)_s
    itpa = I + (len*phylocov)*pA; // this is 
    itpainv = try_inv(itpa,nvar);
    
    if(itpainv(0,0)==1.0123456789)
    {
      return 1e128;
    }
    
    // logdet(I+T*pA)
    lnW(u) += logdet(itpa);
    lnW(anc) += lnW(u);
    // XX is p = X'W-1X
    XX((u*nvar),0,arma::size(nvar,nvar)) = pA * itpainv; // this is p
    XX((anc*nvar),0,arma::size(nvar,nvar)) += XX((u*nvar),0,arma::size(nvar,nvar)); // this adds to pA of the anc
    
    // YY is Q = Y'W-1Y
    arma::mat YU_XY = (trans(XY((u*nvar),0,arma::size(nvar,1)))*itpainv*(len*phylocov)*XY((u*nvar),0,arma::size(nvar,1)));
    YY(u) = YY(u) - (YU_XY(0,0));
    YY(anc) += YY(u);
    
    // XY is V = X'W-1Y
    XY((u*nvar),0,arma::size(nvar,1)) = trans(trans(XY((u*nvar),0,arma::size(nvar,1)))*itpainv); // this is V = X'W-1Y
    XY((anc*nvar),0,arma::size(nvar,1)) += XY((u*nvar),0,arma::size(nvar,1)); // this adds V to anc
    
    // YU is U = Y'W-1X
    //YU(0,(u*nvar),arma::size(1,nvar)) = YU(0,(u*nvar),arma::size(1,nvar)) * itpainv;
    //YU(0,(anc*nvar),arma::size(1,nvar)) += YU(0,(u*nvar),arma::size(1,nvar));
  }
  u = nedge;
  len = rootedge;
  pA = XX((u*nvar),0,arma::size(nvar,nvar)); // this is sum_s (X'W-1X)_s
  itpa = I + (len*phylocov)*pA; // this is 
  itpainv = try_inv(itpa,nvar);
  lnW(u) += logdet(itpa);
  //Rcout << XX((u*nvar),0,arma::size(nvar,nvar)) << std::endl;
  //Rcout << XY((u*nvar),0,arma::size(nvar,1));
  
  XX((u*nvar),0,arma::size(nvar,nvar)) = XX((u*nvar),0,arma::size(nvar,nvar)) * itpainv; // this is p
  
  arma::mat YU_XY = (trans(XY((u*nvar),0,arma::size(nvar,1)))*itpainv*(len*phylocov)*XY((u*nvar),0,arma::size(nvar,1)));
  YY(u) = (YY(u) - (YU_XY(0,0)));
  
  XY((u*nvar),0,arma::size(nvar,1)) = (trans(trans(XY((u*nvar),0,arma::size(nvar,1)))*itpainv)); // this is V = X'W-1Y
  
  //YU(0,(u*nvar),arma::size(1,nvar)) = YU(0,(u*nvar),arma::size(1,nvar)) * itpainv;
  arma::mat root_anc = try_solve(XX((u*nvar),0,arma::size(nvar,nvar)),XY((u*nvar),0,arma::size(nvar,1)));
  arma::mat anc_XY = root_anc.t() * XY((u*nvar),0,arma::size(nvar,1));
  arma::mat anc_XX = root_anc.t() * XX((u*nvar),0,arma::size(nvar,nvar)) * root_anc;
  double minus2ll = 0;
  if(REML==1)
  {
    minus2ll = (logdet(XX((u*nvar),0,arma::size(nvar,nvar))) + lnW(u) + YY(u) - (2 * anc_XY(0,0)) + anc_XX(0,0)) / 2;
  } else
  {
    minus2ll = (lnW(u) + YY(u) - (2 * anc_XY(0,0)) + anc_XX(0,0)) / 2;
  }
  
  if(verbose==1)
  {
    Rcout << minus2ll << std::endl;
  }
  /*if(minus2ll<-1e3)
  {
  return 1e128;
  }*/
  return List::create(XU,XX,YY,XY,lnW);
}

// [[Rcpp::export]]
double threepoint_phenocorr(arma::vec theta,arma::vec options,arma::vec y,arma::vec ku,arma::vec iu,arma::mat edge,arma::vec edgelength,Rcpp::List uchildren_list,List subset_list,List species_subset,arma::vec tip_combn,arma::vec ymin,arma::vec ymax,arma::vec models,arma::uvec externalEdge,arma::uvec not_externalEdge,arma::vec dist_anc,arma::vec dist_des,double Tmax,double Tmin,int nmodels,arma::vec lower_bounds,arma::vec upper_bounds,arma::vec OU_D,arma::vec times)
{
  // edgevec is a vector of tip lengths
  // subset_list is a list of unique matrices to invert once
  // species_subset is a list of length nspecies
  // tip_combn is a vector signifyinc which matrix subset corresponds to each species
  // ymin, ymax
  
  // options
  // 0 = correlated
  // 1 = error
  // 2 = nvar
  // 3 = clip
  // 4 = verbose
  // 5 = nspecies
  // 6 = nedge
  // 7 = nob
  // 8 = nedge + nob + 1
  
  int nspecies = options(5);
  int nedge = options(6);
  
  double rootedge = 0;
  bool OU_mod = ((models(0) == 1) | (models(1) == 1));
  bool OUfixedRoot = models(0) == 1;
  bool lambda_mod = models(2) == 1;
  bool kappa_mod = models(3) == 1;
  bool delta_mod = models(4) == 1;
  bool EB_mod = models(5) == 1;
  bool pars_from_theta = true;
  double alpha = 0;
  double lambda = 0;
  double kappa = 0;
  double delta = 0;
  double rate = 0;
  int par_count = 0;
  int npars = options(11);
  if((unsigned int)(npars)==theta.size())
  {
    pars_from_theta = false;
  }
  if(OU_mod)
  {
    arma::ucolvec anc = arma::conv_to<arma::uvec>::from(edge.col(0));
    arma::ucolvec des = arma::conv_to<arma::uvec>::from(edge.col(1));
    if(pars_from_theta)
    {
      alpha = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      alpha = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    par_count++;
    if(alpha>1e-7/Tmax)
    {
      if(!OUfixedRoot)
      {
        arma::vec distFromRoot = exp(-2*alpha*times);
        arma::vec d1 = distFromRoot(anc-nspecies);
        arma::vec d2 = arma::zeros(nedge);
        d2(externalEdge) = exp(-2*alpha*OU_D(des(externalEdge)));
        d2(not_externalEdge) = distFromRoot(des(not_externalEdge)-nspecies);
        edgelength = d2 - d1;
        rootedge = arma::min(distFromRoot);
      } else
      {
        arma::vec distFromRoot = exp(-2*alpha*times)%(1 - exp(-2*alpha*(Tmax-times)));
        arma::vec d1 = distFromRoot(anc-nspecies);
        arma::vec d2 = arma::zeros(nedge);
        d2(externalEdge) = exp(-2*alpha*OU_D(des(externalEdge))) % (1-exp(-2*alpha*(Tmax-OU_D(des(externalEdge)))));
        d2(not_externalEdge) = distFromRoot(des(not_externalEdge)-nspecies);
        edgelength = d2 - d1;
        rootedge = arma::min(distFromRoot);
      }
    }
  }
  if(lambda_mod)
  {
    if(pars_from_theta)
    {
      lambda = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      lambda = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    if((lambda<1) & (lambda>1e-7))
    {
      edgelength = edgelength * lambda;
      edgelength(externalEdge) = edgelength(externalEdge) + (1-lambda)*dist_des(externalEdge);
    }
    par_count++;
  }
  if(kappa_mod)
  {
    if(pars_from_theta)
    {
      kappa = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      kappa = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    if((kappa<1) & (kappa>1e-6))
    {
      edgelength = exp(kappa*log(edgelength));
    }
    par_count++;
  }
  if(delta_mod)
  {
    if(pars_from_theta)
    {
      delta = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      delta = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    if(delta>1e-5)
    {
      edgelength = (exp(delta*log(dist_des)) - exp(delta*log(dist_anc)))*exp((1-delta)*log(Tmax));
    }
    par_count++;
  }
  if(EB_mod)
  {
    if(pars_from_theta)
    {
      rate = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      rate = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    if(rate!=0)
    {
      edgelength = (exp(rate*dist_des)-exp(rate*dist_anc))/rate;
    }
    par_count++;
  }
  
  int nvar = options(2);
  int nn = options(8);
  int verbose = options(4);
  int nind = options(9);
  int REML = options(10);
  int ret_YY = options(13);
  arma::mat AB = convert_pars(theta,options,Tmin);
  if(AB(0,0)==1.0123456789)
  {
    //Rcout << "pars conversion failed";
    return 1e128;
  }
  arma::mat phylocov = AB(0,0,arma::size(nvar,nvar));
  arma::mat phenocov = AB(nvar,0,arma::size(nvar,nvar));
  arma::mat itpainv = arma::zeros(nvar,nvar);
  arma::mat inverse = arma::zeros(nvar,nvar);
  
  List sub_mats = inv_subset(phenocov,subset_list);
  
  arma::mat XU = arma::zeros(nvar*nn,nvar);
  arma::mat XY = arma::zeros(nvar*nn,1);
  arma::mat YU = arma::zeros(nn,nvar);
  arma::mat XX = arma::zeros(nvar*nn,nvar);
  arma::vec YY = arma::zeros(nn);
  arma::vec lnW = arma::zeros(nn);
  arma::mat YU_XY = arma::zeros(1,nvar);
  int u = 0;
  int e = 0;
  int uc_length = 0;
  int uci = 0;
  double len = 0;
  double determinant = 0;
  arma::uvec zero_vec(1);
  zero_vec(0)=0;
  for(u=0;u<nind;u++)
  {
    arma::uvec Ka = species_subset(u);
    arma::mat Bainv = sub_mats(tip_combn(u));
    if(Bainv(0,0)==1.0123456789)
    {
      //Rcout << "Ba inverse failed";
      return 1e128;
    }
    //int ma = y.size();
    arma::mat Ma = Bainv*phylocov.rows(Ka);    
    XU.rows((u*nvar)+Ka) = Ma;
    XY.submat((u*nvar)+Ka,zero_vec) = Bainv * y.subvec(ymin(u),ymax(u));
    XX.submat((u*nvar)+Ka,Ka) = Bainv;
    YU.row(u) = trans(y.subvec(ymin(u),ymax(u))) * Ma;
    arma::mat temp_mat = (trans(y.subvec(ymin(u),ymax(u))) * Bainv * y.subvec(ymin(u),ymax(u)));
    YY(u) = temp_mat(0,0);
    lnW(u) = logdet(phenocov.submat(Ka,Ka));
  }
  //int i = 0;
  for(e=0;e<(nedge+1);e++) // e=nedge+1 is to add the root edge
  {
    u = nind + e; // index of edge after observations were added as polytomies
    /*if(e==nedge)
    {
    i = nspecies;
    } else
    {
    i = edge(e,1) + 1;
    }*/
    arma::uvec uchildren = uchildren_list(e);
    uc_length = uchildren(0);
    for(uci=1;uci<uc_length;uci++)
    {
      XU((u*nvar),0,arma::size(nvar,nvar)) = XU((u*nvar),0,arma::size(nvar,nvar)) + XU((uchildren[uci]*nvar),0,arma::size(nvar,nvar));
      XY((u*nvar),0,arma::size(nvar,1)) = XY((u*nvar),0,arma::size(nvar,1)) + XY(uchildren[uci]*nvar,0,arma::size(nvar,1));
      XX((u*nvar),0,arma::size(nvar,nvar)) = XX((u*nvar),0,arma::size(nvar,nvar)) + XX((uchildren[uci]*nvar),0,arma::size(nvar,nvar));
      YU((u),0,arma::size(1,nvar)) = YU((u),0,arma::size(1,nvar)) + YU(uchildren[uci],0,arma::size(1,nvar));    
      YY(u) = YY(u) + YY(uchildren(uci));
      lnW(u) = lnW(u) + lnW(uchildren(uci));
    }
    if(e==nedge)
    {
      len = rootedge;
    } else
    {
      len = edgelength(e);
    }
    itpainv = arma::eye(nvar,nvar) + (XU((u*nvar),0,arma::size(nvar,nvar))*len);
    inverse = try_inv(itpainv,nvar);
    if(inverse(0,0)==1.0123456789)
    {
      return 1e128;
    }
    determinant = det(itpainv);
    lnW(u) = lnW(u) + log(determinant);
    
    XY((u*nvar),0,arma::size(nvar,1)) = inverse * XY((u*nvar),0,arma::size(nvar,1));
    XX((u*nvar),0,arma::size(nvar,nvar)) = inverse * XX((u*nvar),0,arma::size(nvar,nvar));
    arma::mat YU_XY = (YU((u),0,arma::size(1,nvar))*XY((u*nvar),0,arma::size(nvar,1)));
    YY(u) = YY(u) - (len * YU_XY(0,0));
    YU((u),0,arma::size(1,nvar)) = YU((u),0,arma::size(1,nvar)) * inverse;
    XU((u*nvar),0,arma::size(nvar,nvar)) = XU((u*nvar),0,arma::size(nvar,nvar)) * inverse;
  }
  if(ret_YY==1)
  {
    return YY(u);
  }
  
  arma::mat anc = try_solve(XX((u*nvar),0,arma::size(nvar,nvar)),XY((u*nvar),0,arma::size(nvar,1)));
  /*for(i=0;i<nvar;i++)
  {
  if(anc(i,0)>1e3 || anc(i,0)<-1e3 || anc(i,0)==1.0123456789)
  {
  //Rcout << "anc out of range" << std::endl;
  return 1e128;
  }
  }*/
  
  arma::mat anc_XY = anc.t() * XY((u*nvar),0,arma::size(nvar,1));
  arma::mat anc_XX = anc.t() * XX((u*nvar),0,arma::size(nvar,nvar)) * anc;
  double minus2ll = 0;
  if(REML==1)
  {
    minus2ll = (log(det(XX((u*nvar),0,arma::size(nvar,nvar)))) + lnW(u) + YY(u) - (2 * anc_XY(0,0)) + anc_XX(0,0)) / 2;
  } else
  {
    minus2ll = (lnW(u) + YY(u) - (2 * anc_XY(0,0)) + anc_XX(0,0)) / 2;
  }
  
  if(verbose==1)
  {
    Rcout << minus2ll << std::endl;
  }
  /*if(minus2ll<-1e3)
  {
  return 1e128;
  }*/
  return minus2ll;
  //return List::create(XU,XX,YU,YY,XY,lnW);
  }

// [[Rcpp::export]]
double threepoint_calc_pheno(arma::vec theta,arma::vec options,arma::vec y,arma::vec ku,arma::vec iu,arma::mat edge,arma::vec edgelength,Rcpp::List uchildren_list,List subset_list,List species_subset,arma::vec tip_combn,arma::vec ymin,arma::vec ymax,List phenocovs,List inv_phenocovs,arma::vec models,arma::uvec externalEdge,arma::uvec not_externalEdge,arma::vec dist_anc,arma::vec dist_des,double Tmax,double Tmin,int nmodels,arma::vec lower_bounds,arma::vec upper_bounds,arma::vec OU_D,arma::vec times)
{
  // edgevec is a vector of tip lengths
  // subset_list is a list of unique matrices to invert once
  // species_subset is a list of length nspecies
  // tip_combn is a vector signifyinc which matrix subset corresponds to each species
  // ymin, ymax
  
  // options
  // 0 = correlated
  // 1 = error
  // 2 = nvar
  // 3 = clip
  // 4 = verbose
  // 5 = nspecies
  // 6 = nedge
  // 7 = nob
  // 8 = nedge + nob + 1
  
  int nspecies = options(5);
  int nedge = options(6);
  
  double rootedge = 0;
  bool OU_mod = ((models(0) == 1) | (models(1) == 1));
  bool OUfixedRoot = models(0) == 1;
  bool lambda_mod = models(2) == 1;
  bool kappa_mod = models(3) == 1;
  bool delta_mod = models(4) == 1;
  bool EB_mod = models(5) == 1;
  bool pars_from_theta = true;
  double alpha = 0;
  double lambda = 0;
  double kappa = 0;
  double delta = 0;
  double rate = 0;
  int par_count = 0;
  int npars = options(11);
  if((unsigned int)(npars)==theta.size())
  {
    pars_from_theta = false;
  }
  if(OU_mod)
  {
    arma::ucolvec anc = arma::conv_to<arma::uvec>::from(edge.col(0));
    arma::ucolvec des = arma::conv_to<arma::uvec>::from(edge.col(1));
    if(pars_from_theta)
    {
      alpha = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      alpha = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    par_count++;
    if(alpha>1e-7/Tmax)
    {
      if(!OUfixedRoot)
      {
        arma::vec distFromRoot = exp(-2*alpha*times);
        arma::vec d1 = distFromRoot(anc-nspecies);
        arma::vec d2 = arma::zeros(nedge);
        d2(externalEdge) = exp(-2*alpha*OU_D(des(externalEdge)));
        d2(not_externalEdge) = distFromRoot(des(not_externalEdge)-nspecies);
        edgelength = d2 - d1;
        rootedge = arma::min(distFromRoot);
      } else
      {
        arma::vec distFromRoot = exp(-2*alpha*times)%(1 - exp(-2*alpha*(Tmax-times)));
        arma::vec d1 = distFromRoot(anc-nspecies);
        arma::vec d2 = arma::zeros(nedge);
        d2(externalEdge) = exp(-2*alpha*OU_D(des(externalEdge))) % (1-exp(-2*alpha*(Tmax-OU_D(des(externalEdge)))));
        d2(not_externalEdge) = distFromRoot(des(not_externalEdge)-nspecies);
        edgelength = d2 - d1;
        rootedge = arma::min(distFromRoot);
      }
    }
  }
  if(lambda_mod)
  {
    if(pars_from_theta)
    {
      lambda = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      lambda = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    if((lambda<1) & (lambda>1e-7))
    {
      edgelength = edgelength * lambda;
      edgelength(externalEdge) = edgelength(externalEdge) + (1-lambda)*dist_des(externalEdge);
    }
    par_count++;
  }
  if(kappa_mod)
  {
    if(pars_from_theta)
    {
      kappa = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      kappa = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    if((kappa<1) & (kappa>1e-6))
    {
      edgelength = exp(kappa*log(edgelength));
    }
    par_count++;
  }
  if(delta_mod)
  {
    if(pars_from_theta)
    {
      delta = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      delta = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    if(delta>1e-5)
    {
      edgelength = (exp(delta*log(dist_des)) - exp(delta*log(dist_anc)))*exp((1-delta)*log(Tmax));
    }
    par_count++;
  }
  if(EB_mod)
  {
    if(pars_from_theta)
    {
      rate = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      rate = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    if(rate!=0)
    {
      edgelength = (exp(rate*dist_des)-exp(rate*dist_anc))/rate;
    }
    par_count++;
  }
  int nvar = options(2);
  int nn = options(8);
  int verbose = options(4);
  int nind = options(9);
  int REML = options(10);
  
  int ret_YY = options(13);
  
  arma::mat AB = convert_pars(theta,options,Tmin);
  if(AB(0,0)==1.0123456789)
  {
    //Rcout << "pars conversion failed";
    return 1e128;
  }
  arma::mat phylocov = AB(0,0,arma::size(nvar,nvar));
  arma::mat itpainv = arma::zeros(nvar,nvar);
  arma::mat inverse = arma::zeros(nvar,nvar);
  
  arma::mat XU = arma::zeros(nvar*nn,nvar);
  arma::mat XY = arma::zeros(nvar*nn,1);
  arma::mat YU = arma::zeros(nn,nvar);
  arma::mat XX = arma::zeros(nvar*nn,nvar);
  arma::vec YY = arma::zeros(nn);
  arma::vec lnW = arma::zeros(nn);
  arma::mat YU_XY = arma::zeros(1,nvar);
  int u = 0;
  int e = 0;
  int uc_length = 0;
  int uci = 0;
  double len = 0;
  double determinant = 0;
  arma::uvec zero_vec(1);
  zero_vec(0)=0;
  for(u=0;u<nind;u++)
  {
    arma::mat Bainv = inv_phenocovs(u);
    arma::mat phenocov = phenocovs(u);
    arma::uvec Ka = species_subset(u);
    //int ma = y.size();
    arma::mat Ma = Bainv*phylocov.rows(Ka);    
    XU.rows((u*nvar)+Ka) = Ma;
    XY.submat((u*nvar)+Ka,zero_vec) = Bainv * y.subvec(ymin(u),ymax(u));
    XX.submat((u*nvar)+Ka,Ka) = Bainv;
    YU.row(u) = trans(y.subvec(ymin(u),ymax(u))) * Ma;
    arma::mat temp_mat = (trans(y.subvec(ymin(u),ymax(u))) * Bainv * y.subvec(ymin(u),ymax(u)));
    YY(u) = temp_mat(0,0);
    lnW(u) = log(det(phenocov));
  }
  //int i = 0;
  for(e=0;e<(nedge+1);e++) // e=nedge+1 is to add the root edge
  {
    u = nind + e; // index of edge after observations were added as polytomies
    /*if(e==nedge)
    {
    i = nspecies;
    } else
    {
    i = edge(e,1) + 1;
    }*/
    arma::uvec uchildren = uchildren_list(e);
    uc_length = uchildren(0);
    for(uci=1;uci<uc_length;uci++)
    {
      XU((u*nvar),0,arma::size(nvar,nvar)) = XU((u*nvar),0,arma::size(nvar,nvar)) + XU((uchildren[uci]*nvar),0,arma::size(nvar,nvar));
      XY((u*nvar),0,arma::size(nvar,1)) = XY((u*nvar),0,arma::size(nvar,1)) + XY(uchildren[uci]*nvar,0,arma::size(nvar,1));
      XX((u*nvar),0,arma::size(nvar,nvar)) = XX((u*nvar),0,arma::size(nvar,nvar)) + XX((uchildren[uci]*nvar),0,arma::size(nvar,nvar));
      YU((u),0,arma::size(1,nvar)) = YU((u),0,arma::size(1,nvar)) + YU(uchildren[uci],0,arma::size(1,nvar));    
      YY(u) = YY(u) + YY(uchildren(uci));
      lnW(u) = lnW(u) + lnW(uchildren(uci));
    }
    if(e==nedge)
    {
      len = rootedge;
    } else
    {
      len = edgelength(e);
    }
    itpainv = arma::eye(nvar,nvar) + (XU((u*nvar),0,arma::size(nvar,nvar))*len);
    inverse = try_inv(itpainv,nvar);
    if(inverse(0,0)==1.0123456789)
    {
      return 1e128;
    }
    determinant = det(itpainv);
    lnW(u) = lnW(u) + log(determinant);
    XY((u*nvar),0,arma::size(nvar,1)) = inverse * XY((u*nvar),0,arma::size(nvar,1));
    XX((u*nvar),0,arma::size(nvar,nvar)) = inverse * XX((u*nvar),0,arma::size(nvar,nvar));
    arma::mat YU_XY = (YU((u),0,arma::size(1,nvar))*XY((u*nvar),0,arma::size(nvar,1)));
    YY(u) = YY(u) - (len * YU_XY(0,0));
    YU((u),0,arma::size(1,nvar)) = YU((u),0,arma::size(1,nvar)) * inverse;
    XU((u*nvar),0,arma::size(nvar,nvar)) = XU((u*nvar),0,arma::size(nvar,nvar)) * inverse;
  }
  if(ret_YY==1)
  {
    return YY(u);
  }
  
  arma::mat anc = try_solve(XX((u*nvar),0,arma::size(nvar,nvar)),XY((u*nvar),0,arma::size(nvar,1)));
  /*for(i=0;i<nvar;i++)
  {
  if(anc(i,0)>1e3 || anc(i,0)<-1e3 || anc(i,0)==1.0123456789)
  {
  //Rcout << "anc out of range" << std::endl;
  return 1e128;
  }
  }*/
  arma::mat anc_XY = anc.t() * XY((u*nvar),0,arma::size(nvar,1));
  arma::mat anc_XX = anc.t() * XX((u*nvar),0,arma::size(nvar,nvar)) * anc;
  double minus2ll = 0;
  if(REML==1)
  {
    minus2ll = (log(det(XX((u*nvar),0,arma::size(nvar,nvar)))) + lnW(u) + YY(u) - (2 * anc_XY(0,0)) + anc_XX(0,0)) / 2;
  } else
  {
    minus2ll = (lnW(u) + YY(u) - (2 * anc_XY(0,0)) + anc_XX(0,0)) / 2;
  }
  
  if(verbose==1)
  {
    Rcout << minus2ll << std::endl;
  }
  /*if(minus2ll<-1e3)
  {
  return 1e128;
  }*/
  return minus2ll;
  //return List::create(XU,XX,YU,YY,XY,lnW);
  }

// [[Rcpp::export]]
arma::mat threepoint_predict(arma::vec theta,arma::vec options,arma::vec y,arma::vec ku,arma::vec iu,arma::mat edge,arma::vec edgelength,Rcpp::List uchildren_list,arma::vec models,arma::uvec externalEdge,arma::uvec not_externalEdge,arma::vec dist_anc,arma::vec dist_des,double Tmax,double Tmin,int nmodels,arma::vec lower_bounds,arma::vec upper_bounds,arma::vec OU_D,arma::vec times)
{
  // options
  // 0 = correlated
  // 1 = error
  // 2 = nvar
  // 3 = clip
  // 4 = verbose
  // 5 = nspecies
  // 6 = nedge
  // 7 = nob
  // 8 = nedge + nob + 1
  
  int nspecies = options(5);
  int nedge = options(6);
  
  double rootedge = 0;
  bool OU_mod = ((models(0) == 1) | (models(1) == 1));
  bool OUfixedRoot = models(0) == 1;
  bool lambda_mod = models(2) == 1;
  bool kappa_mod = models(3) == 1;
  bool delta_mod = models(4) == 1;
  bool EB_mod = models(5) == 1;
  bool pars_from_theta = true;
  double alpha = 0;
  double lambda = 0;
  double kappa = 0;
  double delta = 0;
  double rate = 0;
  int par_count = 0;
  int npars = options(11);
  if((unsigned int)(npars)==theta.size())
  {
    pars_from_theta = false;
  }
  if(OU_mod)
  {
    arma::ucolvec anc = arma::conv_to<arma::uvec>::from(edge.col(0));
    arma::ucolvec des = arma::conv_to<arma::uvec>::from(edge.col(1));
    if(pars_from_theta)
    {
      alpha = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      alpha = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    par_count++;
    if(alpha>1e-7/Tmax)
    {
      if(!OUfixedRoot)
      {
        arma::vec distFromRoot = exp(-2*alpha*times);
        arma::vec d1 = distFromRoot(anc-nspecies);
        arma::vec d2 = arma::zeros(nedge);
        d2(externalEdge) = exp(-2*alpha*OU_D(des(externalEdge)));
        d2(not_externalEdge) = distFromRoot(des(not_externalEdge)-nspecies);
        edgelength = d2 - d1;
        rootedge = arma::min(distFromRoot);
      } else
      {
        arma::vec distFromRoot = exp(-2*alpha*times)%(1 - exp(-2*alpha*(Tmax-times)));
        arma::vec d1 = distFromRoot(anc-nspecies);
        arma::vec d2 = arma::zeros(nedge);
        d2(externalEdge) = exp(-2*alpha*OU_D(des(externalEdge))) % (1-exp(-2*alpha*(Tmax-OU_D(des(externalEdge)))));
        d2(not_externalEdge) = distFromRoot(des(not_externalEdge)-nspecies);
        edgelength = d2 - d1;
        rootedge = arma::min(distFromRoot);
      }
    }
  }
  if(lambda_mod)
  {
    if(pars_from_theta)
    {
      lambda = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      lambda = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    if((lambda<1) & (lambda>1e-7))
    {
      edgelength = edgelength * lambda;
      edgelength(externalEdge) = edgelength(externalEdge) + (1-lambda)*dist_des(externalEdge);
    }
    par_count++;
  }
  if(kappa_mod)
  {
    if(pars_from_theta)
    {
      kappa = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      kappa = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    if((kappa<1) & (kappa>1e-6))
    {
      edgelength = exp(kappa*log(edgelength));
    }
    par_count++;
  }
  if(delta_mod)
  {
    if(pars_from_theta)
    {
      delta = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      delta = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    if(delta>1e-5)
    {
      edgelength = (exp(delta*log(dist_des)) - exp(delta*log(dist_anc)))*exp((1-delta)*log(Tmax));
    }
    par_count++;
  }
  if(EB_mod)
  {
    if(pars_from_theta)
    {
      rate = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      rate = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    if(rate!=0)
    {
      edgelength = (exp(rate*dist_des)-exp(rate*dist_anc))/rate;
    }
    par_count++;
  }
  
  int nvar = options(2);
  int nob = options(7);
  int nn = options(8);
  
  arma::mat AB = convert_pars(theta,options,Tmin);
  arma::mat phylocov = AB(0,0,arma::size(nvar,nvar));
  arma::mat phenocov = AB(nvar,0,arma::size(nvar,nvar));
  arma::mat itpainv = arma::zeros(nvar,nvar);
  arma::mat inverse = arma::zeros(nvar,nvar);
  
  arma::mat XU = arma::zeros(nvar*nn,nvar);
  arma::mat XY = arma::zeros(nvar*nn,1);
  arma::mat XX = arma::zeros(nvar*nn,nvar);
  int u = 0;
  int k = 0;
  int e = 0;
  int uc_length = 0;
  int uci = 0;
  double len = 0;
  
  for(u=0;u<nob;u++)
  {
    k = ku(u);
    XU.row((u*nvar)+k) = phylocov.row(k) / phenocov(k,k);
    XY((u*nvar)+k,0) = y(u) / phenocov(k,k);
    XX((u*nvar)+k,k) = 1 / phenocov(k,k);
  }
  
  for(e=0;e<(nedge+1);e++) // e=nedge+1 is to add the root edge
  {
    u = nob + e; // index of edge after observations were added as polytomies
    /*if(e==nedge)
    {
    i = nspecies;
    } else
    {
    i = edge(e,1);
    }*/
    arma::uvec uchildren = uchildren_list(e);
    uc_length = uchildren(0);
    
    for(uci=1;uci<uc_length;uci++)
    {
      XU((u*nvar),0,arma::size(nvar,nvar)) = XU((u*nvar),0,arma::size(nvar,nvar)) + XU((uchildren[uci]*nvar),0,arma::size(nvar,nvar));
      XY((u*nvar),0,arma::size(nvar,1)) = XY((u*nvar),0,arma::size(nvar,1)) + XY(uchildren[uci]*nvar,0,arma::size(nvar,1));
      XX((u*nvar),0,arma::size(nvar,nvar)) = XX((u*nvar),0,arma::size(nvar,nvar)) + XX((uchildren[uci]*nvar),0,arma::size(nvar,nvar));
    }
    if(e==nedge)
    {
      len = rootedge;
    } else
    {
      len = edgelength(e);
    }
    itpainv = arma::eye(nvar,nvar) + (XU((u*nvar),0,arma::size(nvar,nvar))*len);
    inverse = inv(itpainv);
    XY((u*nvar),0,arma::size(nvar,1)) = inverse * XY((u*nvar),0,arma::size(nvar,1));
    XX((u*nvar),0,arma::size(nvar,nvar)) = inverse * XX((u*nvar),0,arma::size(nvar,nvar));
    XU((u*nvar),0,arma::size(nvar,nvar)) = XU((u*nvar),0,arma::size(nvar,nvar)) * inverse;
  }  
  arma::vec anc = try_solve(XX((u*nvar),0,arma::size(nvar,nvar)),XY((u*nvar),0,arma::size(nvar,1)));
  arma::mat XXinv = inv(XX((u*nvar),0,arma::size(nvar,nvar)));
  arma::vec anc_sd = arma::sqrt(XXinv.diag());
  arma::mat ret = arma::zeros(nvar,2);
  ret.col(0) = anc;
  ret.col(1) = anc_sd;
  return ret;
  }

// [[Rcpp::export]]
arma::mat threepoint_nopheno_predict(arma::vec theta,arma::vec options,arma::vec y,arma::vec ku,arma::vec iu,arma::mat edge,arma::vec edgelength,Rcpp::List uchildren_list,arma::vec edgevec,List subset_list,List species_subset,arma::vec tip_combn,arma::vec ymin,arma::vec ymax,arma::vec models,arma::uvec externalEdge,arma::uvec not_externalEdge,arma::vec dist_anc,arma::vec dist_des,double Tmax,double Tmin,int nmodels,arma::vec lower_bounds,arma::vec upper_bounds,arma::vec OU_D,arma::vec times)
{
  // edgevec is a vector of tip lengths
  // subset_list is a list of unique matrices to invert once
  // species_subset is a list of length nspecies
  // tip_combn is a vector signifyinc which matrix subset corresponds to each species
  // ymin, ymax
  
  // options
  // 0 = correlated
  // 1 = error
  // 2 = nvar
  // 3 = clip
  // 4 = verbose
  // 5 = nspecies
  // 6 = nedge
  // 7 = nob
  // 8 = nedge + nob + 1
  
  int nspecies = options(5);
  int nedge = options(6);
  
  double rootedge = 0;
  bool OU_mod = ((models(0) == 1) | (models(1) == 1));
  bool OUfixedRoot = models(0) == 1;
  bool lambda_mod = models(2) == 1;
  bool kappa_mod = models(3) == 1;
  bool delta_mod = models(4) == 1;
  bool EB_mod = models(5) == 1;
  bool pars_from_theta = true;
  double alpha = 0;
  double lambda = 0;
  double kappa = 0;
  double delta = 0;
  double rate = 0;
  int par_count = 0;
  int npars = options(11);
  if((unsigned int)(npars)==theta.size())
  {
    pars_from_theta = false;
  }
  if(OU_mod)
  {
    arma::ucolvec anc = arma::conv_to<arma::uvec>::from(edge.col(0));
    arma::ucolvec des = arma::conv_to<arma::uvec>::from(edge.col(1));
    if(pars_from_theta)
    {
      alpha = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      alpha = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    par_count++;
    if(alpha>1e-7/Tmax)
    {
      if(!OUfixedRoot)
      {
        arma::vec distFromRoot = exp(-2*alpha*times);
        arma::vec d1 = distFromRoot(anc-nspecies);
        arma::vec d2 = arma::zeros(nedge);
        d2(externalEdge) = exp(-2*alpha*OU_D(des(externalEdge)));
        d2(not_externalEdge) = distFromRoot(des(not_externalEdge)-nspecies);
        edgelength = d2 - d1;
        rootedge = arma::min(distFromRoot);
      } else
      {
        arma::vec distFromRoot = exp(-2*alpha*times)%(1 - exp(-2*alpha*(Tmax-times)));
        arma::vec d1 = distFromRoot(anc-nspecies);
        arma::vec d2 = arma::zeros(nedge);
        d2(externalEdge) = exp(-2*alpha*OU_D(des(externalEdge))) % (1-exp(-2*alpha*(Tmax-OU_D(des(externalEdge)))));
        d2(not_externalEdge) = distFromRoot(des(not_externalEdge)-nspecies);
        edgelength = d2 - d1;
        rootedge = arma::min(distFromRoot);
      }
    }
  }
  if(lambda_mod)
  {
    if(pars_from_theta)
    {
      lambda = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      lambda = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    if((lambda<1) & (lambda>1e-7))
    {
      edgelength = edgelength * lambda;
      edgelength(externalEdge) = edgelength(externalEdge) + (1-lambda)*dist_des(externalEdge);
    }
    par_count++;
  }
  if(kappa_mod)
  {
    if(pars_from_theta)
    {
      kappa = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      kappa = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    if((kappa<1) & (kappa>1e-6))
    {
      edgelength = exp(kappa*log(edgelength));
    }
    par_count++;
  }
  if(delta_mod)
  {
    if(pars_from_theta)
    {
      delta = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      delta = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    if(delta>1e-5)
    {
      edgelength = (exp(delta*log(dist_des)) - exp(delta*log(dist_anc)))*exp((1-delta)*log(Tmax));
    }
    par_count++;
  }
  if(EB_mod)
  {
    if(pars_from_theta)
    {
      rate = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      rate = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    if(rate!=0)
    {
      edgelength = (exp(rate*dist_des)-exp(rate*dist_anc))/rate;
    }
    par_count++;
  }
  
  if(par_count>0)
  {
    edgevec = edgelength(externalEdge);
  }
  
  int nvar = options(2);
  int nn = options(8);
  int nind = options(9);
  arma::mat AB = convert_pars(theta,options,Tmin);
  arma::mat phylocov = AB(0,0,arma::size(nvar,nvar));
  arma::mat itpainv = arma::zeros(nvar,nvar);
  arma::mat inverse = arma::zeros(nvar,nvar);
  
  List sub_mats = inv_subset(phylocov,subset_list);
  
  arma::mat XU = arma::zeros(nvar*nn,nvar);
  arma::mat XY = arma::zeros(nvar*nn,1);
  arma::mat YU = arma::zeros(nn,nvar);
  arma::mat XX = arma::zeros(nvar*nn,nvar);
  arma::vec YY = arma::zeros(nn);
  arma::vec lnW = arma::zeros(nn);
  arma::mat YU_XY = arma::zeros(1,nvar);
  int u = 0;
  int e = 0;
  int uc_length = 0;
  int uci = 0;
  double len = 0;
  double determinant = 0;
  arma::uvec zero_vec(1);
  zero_vec(0)=0;
  for(u=0;u<nind;u++)
  {
    len = edgevec(u);
    arma::uvec Ka = species_subset(u);
    
    arma::mat Bainv = sub_mats(tip_combn(u));
    /*if(Bainv(0,0)==1.0123456789)
    {
    Rcout << "Ba inverse failed";
    return 1e128;
    }*/
    Bainv = Bainv / len;
    //int ma = y.size();
    arma::mat Ma = Bainv*phylocov.rows(Ka);
    
    XU.rows((u*nvar)+Ka) = Ma;
    XY.submat((u*nvar)+Ka,zero_vec) = Bainv * y.subvec(ymin(u),ymax(u));
    XX.submat((u*nvar)+Ka,Ka) = Bainv;
    YU.row(u) = trans(y.subvec(ymin(u),ymax(u))) * Ma;
    arma::mat temp_mat = (trans(y.subvec(ymin(u),ymax(u))) * Bainv * y.subvec(ymin(u),ymax(u)));
    YY(u) = temp_mat(0,0);
    lnW(u) = log(det(phylocov.submat(Ka,Ka)*len));
  }
  int i = 0;
  for(e=0;e<(nedge+1);e++) // e=nedge+1 is to add the root edge
  {
    u = nind + e; // index of edge after observations were added as polytomies
    if(e==nedge)
    {
      i = nspecies;
    } else
    {
      i = edge(e,1) + 1;
    }
    arma::uvec uchildren = uchildren_list(e);
    uc_length = uchildren(0);
    for(uci=1;uci<uc_length;uci++)
    {
      XU((u*nvar),0,arma::size(nvar,nvar)) = XU((u*nvar),0,arma::size(nvar,nvar)) + XU((uchildren[uci]*nvar),0,arma::size(nvar,nvar));
      XY((u*nvar),0,arma::size(nvar,1)) = XY((u*nvar),0,arma::size(nvar,1)) + XY(uchildren[uci]*nvar,0,arma::size(nvar,1));
      XX((u*nvar),0,arma::size(nvar,nvar)) = XX((u*nvar),0,arma::size(nvar,nvar)) + XX((uchildren[uci]*nvar),0,arma::size(nvar,nvar));
      YU((u),0,arma::size(1,nvar)) = YU((u),0,arma::size(1,nvar)) + YU(uchildren[uci],0,arma::size(1,nvar));    
      YY(u) = YY(u) + YY(uchildren(uci));
      lnW(u) = lnW(u) + lnW(uchildren(uci));
    }
    if(e==nedge)
    {
      len = rootedge;
    } else if(i<=nspecies)
    {
      len = 0;
    } else
    {
      len = edgelength(e);
    }
    itpainv = arma::eye(nvar,nvar) + (XU((u*nvar),0,arma::size(nvar,nvar))*len);
    inverse = try_inv(itpainv,nvar);
    /*if(inverse(0,0)==1.0123456789)
    {
    Rcout << "itpa inverse failed";
    return 1e128;
    }*/
    determinant = det(itpainv);
    lnW(u) = lnW(u) + log(determinant);
    XY((u*nvar),0,arma::size(nvar,1)) = inverse * XY((u*nvar),0,arma::size(nvar,1));
    XX((u*nvar),0,arma::size(nvar,nvar)) = inverse * XX((u*nvar),0,arma::size(nvar,nvar));
    arma::mat YU_XY = (YU((u),0,arma::size(1,nvar))*XY((u*nvar),0,arma::size(nvar,1)));
    YY(u) = YY(u) - (len * YU_XY(0,0));
    YU((u),0,arma::size(1,nvar)) = YU((u),0,arma::size(1,nvar)) * inverse;
    XU((u*nvar),0,arma::size(nvar,nvar)) = XU((u*nvar),0,arma::size(nvar,nvar)) * inverse;
  }
  arma::vec anc = try_solve(XX((u*nvar),0,arma::size(nvar,nvar)),XY((u*nvar),0,arma::size(nvar,1)));
  arma::mat XXinv = inv(XX((u*nvar),0,arma::size(nvar,nvar)));
  arma::vec anc_sd = arma::sqrt(XXinv.diag());
  arma::mat ret = arma::zeros(nvar,2);
  ret.col(0) = anc;
  ret.col(1) = anc_sd;
  return ret;
  }

// [[Rcpp::export]]
arma::mat threepoint_phenocorr_predict(arma::vec theta,arma::vec options,arma::vec y,arma::vec ku,arma::vec iu,arma::mat edge,arma::vec edgelength,Rcpp::List uchildren_list,List subset_list,List species_subset,arma::vec tip_combn,arma::vec ymin,arma::vec ymax,arma::vec models,arma::uvec externalEdge,arma::uvec not_externalEdge,arma::vec dist_anc,arma::vec dist_des,double Tmax,double Tmin,int nmodels,arma::vec lower_bounds,arma::vec upper_bounds,arma::vec OU_D,arma::vec times)
{
  // edgevec is a vector of tip lengths
  // subset_list is a list of unique matrices to invert once
  // species_subset is a list of length nspecies
  // tip_combn is a vector signifyinc which matrix subset corresponds to each species
  // ymin, ymax
  
  // options
  // 0 = correlated
  // 1 = error
  // 2 = nvar
  // 3 = clip
  // 4 = verbose
  // 5 = nspecies
  // 6 = nedge
  // 7 = nob
  // 8 = nedge + nob + 1
  
  int nspecies = options(5);
  int nedge = options(6);
  
  double rootedge = 0;
  bool OU_mod = ((models(0) == 1) | (models(1) == 1));
  bool OUfixedRoot = models(0) == 1;
  bool lambda_mod = models(2) == 1;
  bool kappa_mod = models(3) == 1;
  bool delta_mod = models(4) == 1;
  bool EB_mod = models(5) == 1;
  bool pars_from_theta = true;
  double alpha = 0;
  double lambda = 0;
  double kappa = 0;
  double delta = 0;
  double rate = 0;
  int par_count = 0;
  int npars = options(11);
  if((unsigned int)(npars)==theta.size())
  {
    pars_from_theta = false;
  }
  if(OU_mod)
  {
    arma::ucolvec anc = arma::conv_to<arma::uvec>::from(edge.col(0));
    arma::ucolvec des = arma::conv_to<arma::uvec>::from(edge.col(1));
    if(pars_from_theta)
    {
      alpha = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      alpha = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    par_count++;
    if(alpha>1e-7/Tmax)
    {
      if(!OUfixedRoot)
      {
        arma::vec distFromRoot = exp(-2*alpha*times);
        arma::vec d1 = distFromRoot(anc-nspecies);
        arma::vec d2 = arma::zeros(nedge);
        d2(externalEdge) = exp(-2*alpha*OU_D(des(externalEdge)));
        d2(not_externalEdge) = distFromRoot(des(not_externalEdge)-nspecies);
        edgelength = d2 - d1;
        rootedge = arma::min(distFromRoot);
      } else
      {
        arma::vec distFromRoot = exp(-2*alpha*times)%(1 - exp(-2*alpha*(Tmax-times)));
        arma::vec d1 = distFromRoot(anc-nspecies);
        arma::vec d2 = arma::zeros(nedge);
        d2(externalEdge) = exp(-2*alpha*OU_D(des(externalEdge))) % (1-exp(-2*alpha*(Tmax-OU_D(des(externalEdge)))));
        d2(not_externalEdge) = distFromRoot(des(not_externalEdge)-nspecies);
        edgelength = d2 - d1;
        rootedge = arma::min(distFromRoot);
      }
    }
  }
  if(lambda_mod)
  {
    if(pars_from_theta)
    {
      lambda = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      lambda = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    if((lambda<1) & (lambda>1e-7))
    {
      edgelength = edgelength * lambda;
      edgelength(externalEdge) = edgelength(externalEdge) + (1-lambda)*dist_des(externalEdge);
    }
    par_count++;
  }
  if(kappa_mod)
  {
    if(pars_from_theta)
    {
      kappa = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      kappa = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    if((kappa<1) & (kappa>1e-6))
    {
      edgelength = exp(kappa*log(edgelength));
    }
    par_count++;
  }
  if(delta_mod)
  {
    if(pars_from_theta)
    {
      delta = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      delta = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    if(delta>1e-5)
    {
      edgelength = (exp(delta*log(dist_des)) - exp(delta*log(dist_anc)))*exp((1-delta)*log(Tmax));
    }
    par_count++;
  }
  if(EB_mod)
  {
    if(pars_from_theta)
    {
      rate = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      rate = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    if(rate!=0)
    {
      edgelength = (exp(rate*dist_des)-exp(rate*dist_anc))/rate;
    }
    par_count++;
  }
  
  int nvar = options(2);
  int nn = options(8);
  int nind = options(9);
  arma::mat AB = convert_pars(theta,options,Tmin);
  arma::mat phylocov = AB(0,0,arma::size(nvar,nvar));
  arma::mat phenocov = AB(nvar,0,arma::size(nvar,nvar));
  arma::mat itpainv = arma::zeros(nvar,nvar);
  arma::mat inverse = arma::zeros(nvar,nvar);
  
  List sub_mats = inv_subset(phenocov,subset_list);
  
  arma::mat XU = arma::zeros(nvar*nn,nvar);
  arma::mat XY = arma::zeros(nvar*nn,1);
  arma::mat YU = arma::zeros(nn,nvar);
  arma::mat XX = arma::zeros(nvar*nn,nvar);
  arma::vec YY = arma::zeros(nn);
  arma::vec lnW = arma::zeros(nn);
  arma::mat YU_XY = arma::zeros(1,nvar);
  int u = 0;
  int e = 0;
  int uc_length = 0;
  int uci = 0;
  double len = 0;
  double determinant = 0;
  arma::uvec zero_vec(1);
  zero_vec(0)=0;
  for(u=0;u<nind;u++)
  {
    arma::uvec Ka = species_subset(u);
    arma::mat Bainv = sub_mats(tip_combn(u));
    arma::mat Ma = Bainv*phylocov.rows(Ka);    
    XU.rows((u*nvar)+Ka) = Ma;
    XY.submat((u*nvar)+Ka,zero_vec) = Bainv * y.subvec(ymin(u),ymax(u));
    XX.submat((u*nvar)+Ka,Ka) = Bainv;
    YU.row(u) = trans(y.subvec(ymin(u),ymax(u))) * Ma;
    arma::mat temp_mat = (trans(y.subvec(ymin(u),ymax(u))) * Bainv * y.subvec(ymin(u),ymax(u)));
    YY(u) = temp_mat(0,0);
    lnW(u) = log(det(phenocov.submat(Ka,Ka)));
  }
  for(e=0;e<(nedge+1);e++) // e=nedge+1 is to add the root edge
  {
    u = nind + e; // index of edge after observations were added as polytomies
    arma::uvec uchildren = uchildren_list(e);
    uc_length = uchildren(0);
    for(uci=1;uci<uc_length;uci++)
    {
      XU((u*nvar),0,arma::size(nvar,nvar)) = XU((u*nvar),0,arma::size(nvar,nvar)) + XU((uchildren[uci]*nvar),0,arma::size(nvar,nvar));
      XY((u*nvar),0,arma::size(nvar,1)) = XY((u*nvar),0,arma::size(nvar,1)) + XY(uchildren[uci]*nvar,0,arma::size(nvar,1));
      XX((u*nvar),0,arma::size(nvar,nvar)) = XX((u*nvar),0,arma::size(nvar,nvar)) + XX((uchildren[uci]*nvar),0,arma::size(nvar,nvar));
      YU((u),0,arma::size(1,nvar)) = YU((u),0,arma::size(1,nvar)) + YU(uchildren[uci],0,arma::size(1,nvar));    
      YY(u) = YY(u) + YY(uchildren(uci));
      lnW(u) = lnW(u) + lnW(uchildren(uci));
    }
    if(e==nedge)
    {
      len = rootedge;
    } else
    {
      len = edgelength(e);
    }
    itpainv = arma::eye(nvar,nvar) + (XU((u*nvar),0,arma::size(nvar,nvar))*len);
    inverse = try_inv(itpainv,nvar);
    /*if(inverse(0,0)==1.0123456789)
    {
    Rcout << "itpa inverse failed";
    return 1e128;
    }*/
    determinant = det(itpainv);
    lnW(u) = lnW(u) + log(determinant);
    XY((u*nvar),0,arma::size(nvar,1)) = inverse * XY((u*nvar),0,arma::size(nvar,1));
    XX((u*nvar),0,arma::size(nvar,nvar)) = inverse * XX((u*nvar),0,arma::size(nvar,nvar));
    arma::mat YU_XY = (YU((u),0,arma::size(1,nvar))*XY((u*nvar),0,arma::size(nvar,1)));
    YY(u) = YY(u) - (len * YU_XY(0,0));
    YU((u),0,arma::size(1,nvar)) = YU((u),0,arma::size(1,nvar)) * inverse;
    XU((u*nvar),0,arma::size(nvar,nvar)) = XU((u*nvar),0,arma::size(nvar,nvar)) * inverse;
  }
  arma::vec anc = try_solve(XX((u*nvar),0,arma::size(nvar,nvar)),XY((u*nvar),0,arma::size(nvar,1)));
  arma::mat XXinv = inv(XX((u*nvar),0,arma::size(nvar,nvar)));
  arma::vec anc_sd = arma::sqrt(XXinv.diag());
  arma::mat ret = arma::zeros(nvar,2);
  ret.col(0) = anc;
  ret.col(1) = anc_sd;
  return ret;
  }

// [[Rcpp::export]]
arma::mat threepoint_calc_pheno_predict(arma::vec theta,arma::vec options,arma::vec y,arma::vec ku,arma::vec iu,arma::mat edge,arma::vec edgelength,Rcpp::List uchildren_list,List subset_list,List species_subset,arma::vec tip_combn,arma::vec ymin,arma::vec ymax,List phenocovs,List inv_phenocovs,arma::vec models,arma::uvec externalEdge,arma::uvec not_externalEdge,arma::vec dist_anc,arma::vec dist_des,double Tmax,double Tmin,int nmodels,arma::vec lower_bounds,arma::vec upper_bounds,arma::vec OU_D,arma::vec times)
{
  // edgevec is a vector of tip lengths
  // subset_list is a list of unique matrices to invert once
  // species_subset is a list of length nspecies
  // tip_combn is a vector signifyinc which matrix subset corresponds to each species
  // ymin, ymax
  
  // options
  // 0 = correlated
  // 1 = error
  // 2 = nvar
  // 3 = clip
  // 4 = verbose
  // 5 = nspecies
  // 6 = nedge
  // 7 = nob
  // 8 = nedge + nob + 1
  
  int nspecies = options(5);
  int nedge = options(6);
  
  double rootedge = 0;
  bool OU_mod = ((models(0) == 1) | (models(1) == 1));
  bool OUfixedRoot = models(0) == 1;
  bool lambda_mod = models(2) == 1;
  bool kappa_mod = models(3) == 1;
  bool delta_mod = models(4) == 1;
  bool EB_mod = models(5) == 1;
  bool pars_from_theta = true;
  double alpha = 0;
  double lambda = 0;
  double kappa = 0;
  double delta = 0;
  double rate = 0;
  int par_count = 0;
  int npars = options(11);
  if((unsigned int)(npars)==theta.size())
  {
    pars_from_theta = false;
  }
  if(OU_mod)
  {
    arma::ucolvec anc = arma::conv_to<arma::uvec>::from(edge.col(0));
    arma::ucolvec des = arma::conv_to<arma::uvec>::from(edge.col(1));
    if(pars_from_theta)
    {
      alpha = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      alpha = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    par_count++;
    if(alpha>1e-7/Tmax)
    {
      if(!OUfixedRoot)
      {
        arma::vec distFromRoot = exp(-2*alpha*times);
        arma::vec d1 = distFromRoot(anc-nspecies);
        arma::vec d2 = arma::zeros(nedge);
        d2(externalEdge) = exp(-2*alpha*OU_D(des(externalEdge)));
        d2(not_externalEdge) = distFromRoot(des(not_externalEdge)-nspecies);
        edgelength = d2 - d1;
        rootedge = arma::min(distFromRoot);
      } else
      {
        arma::vec distFromRoot = exp(-2*alpha*times)%(1 - exp(-2*alpha*(Tmax-times)));
        arma::vec d1 = distFromRoot(anc-nspecies);
        arma::vec d2 = arma::zeros(nedge);
        d2(externalEdge) = exp(-2*alpha*OU_D(des(externalEdge))) % (1-exp(-2*alpha*(Tmax-OU_D(des(externalEdge)))));
        d2(not_externalEdge) = distFromRoot(des(not_externalEdge)-nspecies);
        edgelength = d2 - d1;
        rootedge = arma::min(distFromRoot);
      }
    }
  }
  if(lambda_mod)
  {
    if(pars_from_theta)
    {
      lambda = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      lambda = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    if((lambda<1) & (lambda>1e-7))
    {
      edgelength = edgelength * lambda;
      edgelength(externalEdge) = edgelength(externalEdge) + (1-lambda)*dist_des(externalEdge);
    }
    par_count++;
  }
  if(kappa_mod)
  {
    if(pars_from_theta)
    {
      kappa = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      kappa = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    if((kappa<1) & (kappa>1e-6))
    {
      edgelength = exp(kappa*log(edgelength));
    }
    par_count++;
  }
  if(delta_mod)
  {
    if(pars_from_theta)
    {
      delta = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      delta = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    if(delta>1e-5)
    {
      edgelength = (exp(delta*log(dist_des)) - exp(delta*log(dist_anc)))*exp((1-delta)*log(Tmax));
    }
    par_count++;
  }
  if(EB_mod)
  {
    if(pars_from_theta)
    {
      rate = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(npars+par_count)))+lower_bounds(par_count);
    } else
    {
      rate = (upper_bounds(par_count)-lower_bounds(par_count))/(1+exp(-theta(par_count)))+lower_bounds(par_count);
    }
    if(rate!=0)
    {
      edgelength = (exp(rate*dist_des)-exp(rate*dist_anc))/rate;
    }
    par_count++;
  }
  
  int nvar = options(2);
  int nn = options(8);
  int nind = options(9);
  arma::mat AB = convert_pars(theta,options,Tmin);
  
  arma::mat phylocov = AB(0,0,arma::size(nvar,nvar));
  arma::mat itpainv = arma::zeros(nvar,nvar);
  arma::mat inverse = arma::zeros(nvar,nvar);
  
  arma::mat XU = arma::zeros(nvar*nn,nvar);
  arma::mat XY = arma::zeros(nvar*nn,1);
  arma::mat YU = arma::zeros(nn,nvar);
  arma::mat XX = arma::zeros(nvar*nn,nvar);
  arma::vec YY = arma::zeros(nn);
  arma::vec lnW = arma::zeros(nn);
  arma::mat YU_XY = arma::zeros(1,nvar);
  int u = 0;
  int e = 0;
  int uc_length = 0;
  int uci = 0;
  double len = 0;
  double determinant = 0;
  arma::uvec zero_vec(1);
  zero_vec(0)=0;
  for(u=0;u<nind;u++)
  {
    arma::mat Bainv = inv_phenocovs(u);
    arma::mat phenocov = phenocovs(u);
    
    arma::uvec Ka = species_subset(u);
    arma::mat Ma = Bainv*phylocov.rows(Ka);    
    XU.rows((u*nvar)+Ka) = Ma;
    XY.submat((u*nvar)+Ka,zero_vec) = Bainv * y.subvec(ymin(u),ymax(u));
    XX.submat((u*nvar)+Ka,Ka) = Bainv;
    YU.row(u) = trans(y.subvec(ymin(u),ymax(u))) * Ma;
    arma::mat temp_mat = (trans(y.subvec(ymin(u),ymax(u))) * Bainv * y.subvec(ymin(u),ymax(u)));
    YY(u) = temp_mat(0,0);
    lnW(u) = log(det(phenocov));
  }
  for(e=0;e<(nedge+1);e++) // e=nedge+1 is to add the root edge
  {
    u = nind + e; // index of edge after observations were added as polytomies
    arma::uvec uchildren = uchildren_list(e);
    uc_length = uchildren(0);
    for(uci=1;uci<uc_length;uci++)
    {
      XU((u*nvar),0,arma::size(nvar,nvar)) = XU((u*nvar),0,arma::size(nvar,nvar)) + XU((uchildren[uci]*nvar),0,arma::size(nvar,nvar));
      XY((u*nvar),0,arma::size(nvar,1)) = XY((u*nvar),0,arma::size(nvar,1)) + XY(uchildren[uci]*nvar,0,arma::size(nvar,1));
      XX((u*nvar),0,arma::size(nvar,nvar)) = XX((u*nvar),0,arma::size(nvar,nvar)) + XX((uchildren[uci]*nvar),0,arma::size(nvar,nvar));
      YU((u),0,arma::size(1,nvar)) = YU((u),0,arma::size(1,nvar)) + YU(uchildren[uci],0,arma::size(1,nvar));    
      YY(u) = YY(u) + YY(uchildren(uci));
      lnW(u) = lnW(u) + lnW(uchildren(uci));
    }
    if(e==nedge)
    {
      len = rootedge;
    } else
    {
      len = edgelength(e);
    }
    itpainv = arma::eye(nvar,nvar) + (XU((u*nvar),0,arma::size(nvar,nvar))*len);
    inverse = try_inv(itpainv,nvar);
    /*if(inverse(0,0)==1.0123456789)
    {
    Rcout << "itpa inverse failed";
    return 1e128;
    }*/
    determinant = det(itpainv);
    lnW(u) = lnW(u) + log(determinant);
    XY((u*nvar),0,arma::size(nvar,1)) = inverse * XY((u*nvar),0,arma::size(nvar,1));
    XX((u*nvar),0,arma::size(nvar,nvar)) = inverse * XX((u*nvar),0,arma::size(nvar,nvar));
    arma::mat YU_XY = (YU((u),0,arma::size(1,nvar))*XY((u*nvar),0,arma::size(nvar,1)));
    YY(u) = YY(u) - (len * YU_XY(0,0));
    YU((u),0,arma::size(1,nvar)) = YU((u),0,arma::size(1,nvar)) * inverse;
    XU((u*nvar),0,arma::size(nvar,nvar)) = XU((u*nvar),0,arma::size(nvar,nvar)) * inverse;
  }
  arma::vec anc = try_solve(XX((u*nvar),0,arma::size(nvar,nvar)),XY((u*nvar),0,arma::size(nvar,1)));
  arma::mat XXinv = inv(XX((u*nvar),0,arma::size(nvar,nvar)));
  arma::vec anc_sd = arma::sqrt(XXinv.diag());
  arma::mat ret = arma::zeros(nvar,2);
  ret.col(0) = anc;
  ret.col(1) = anc_sd;
  return ret;
  }

// [[Rcpp::export]]
List recursion(arma::mat edge,int any_species_no_data,arma::vec which_no_species_data,arma::vec iu,arma::vec options,arma::vec species,arma::vec tiplabel)
{
  int e=0,i=0,isnull_length=0,temp_size=0;
  bool pheno = options(1)!=1;
  int nedge = options(6);
  int nspecies = options(5);
  int nob = options(7);
  int nind = options(9);
  List uchildren(nedge+1);
  for(e=0;e<=nedge;e++)
  {
    if (e==nedge)
    {
      i = nspecies;
    } else
    {
      i = edge(e,1);
    }
    if(pheno)
    {
      if((i+1)<=nspecies)
      {
        if(any_species_no_data==1)
        {
          arma::uvec isnull = find(which_no_species_data==(i+1));
          isnull_length = isnull.size();
          if(isnull_length>0)
          {
            arma::uvec temp(1);
            temp(0)=1;
            uchildren(e) = temp;      
          } else
          {
            arma::uvec temp = find(species==tiplabel(i));
            temp_size = temp.size()+1;
            arma::uvec temp2(temp_size);
            temp2(0) = temp_size;
            temp2.subvec(1,temp_size-1) = temp;
            uchildren(e) = temp2;          
          }
        } else
        {
          arma::uvec temp = find(species==tiplabel(i));
          temp_size = temp.size()+1;
          arma::uvec temp2(temp_size);
          temp2(0) = temp_size;
          temp2.subvec(1,temp_size-1) = temp;
          uchildren(e) = temp2;
        }
      } else
      {
        uchildren(e) = condition_f(iu,i,nspecies,edge,nind);
      }
    } else
    {
      if(any_species_no_data==1)
      {
        arma::uvec isnull = find(which_no_species_data==(i+1));
        isnull_length = isnull.size();
        if(isnull_length>0)
        {
          arma::uvec temp(1);
          temp(0)=1;
          uchildren(e) = temp;      
        } else
        {
          uchildren(e) = condition_f(iu,i,nspecies,edge,nob);
        }
      } else
      {
        uchildren(e) = condition_f(iu,i,nspecies,edge,nob);
      }
    }
  }
  return uchildren;
}