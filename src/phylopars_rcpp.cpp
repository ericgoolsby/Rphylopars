// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

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
arma::mat cDot(arma::mat A,arma::mat B)
{
  return A*B;
}

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
double threepoint(arma::vec theta,arma::vec options,arma::vec y,arma::vec ku,arma::vec iu,arma::mat edge,arma::vec edgelength,Rcpp::List uchildren_list,arma::vec non_optim_transform,arma::vec models,arma::uvec externalEdge,arma::uvec not_externalEdge,arma::vec dist_anc,arma::vec dist_des,double Tmax,double Tmin,int nmodels,arma::vec lower_bounds,arma::vec upper_bounds,arma::vec OU_D,arma::vec times)
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
    return 1e32;
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
        return 1e32;
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
      return 1e32;
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
    return 1e32;
  }*/

  if(verbose==1)
  {
    Rcout << minus2ll << std::endl;
  }
  return minus2ll;
}

// [[Rcpp::export]]
double threepoint_nopheno(arma::vec theta,arma::vec options,arma::vec y,arma::vec ku,arma::vec iu,arma::mat edge,arma::vec edgelength,Rcpp::List uchildren_list,arma::vec edgevec,List subset_list,List species_subset,arma::vec tip_combn,arma::vec ymin,arma::vec ymax,arma::vec non_optim_transform,arma::vec models,arma::uvec externalEdge,arma::uvec not_externalEdge,arma::vec dist_anc,arma::vec dist_des,double Tmax,double Tmin,int nmodels,arma::vec lower_bounds,arma::vec upper_bounds,arma::vec OU_D,arma::vec times)
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
    if(!arma::is_finite(alpha)) return 1e32;
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
    return 1e32;
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
      return 1e32;
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
        return 1e32;
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
      return 1e32;
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
    return 1e32;
  }*/
  if(verbose==1)
  {
    Rcout << minus2ll << std::endl;
  }
  return minus2ll;
  //return List::create(XU,XX,YU,YY,XY,lnW);
}

// [[Rcpp::export]]
double threepoint_phenocorr(arma::vec theta,arma::vec options,arma::vec y,arma::vec ku,arma::vec iu,arma::mat edge,arma::vec edgelength,Rcpp::List uchildren_list,List subset_list,List species_subset,arma::vec tip_combn,arma::vec ymin,arma::vec ymax,arma::vec non_optim_transform,arma::vec models,arma::uvec externalEdge,arma::uvec not_externalEdge,arma::vec dist_anc,arma::vec dist_des,double Tmax,double Tmin,int nmodels,arma::vec lower_bounds,arma::vec upper_bounds,arma::vec OU_D,arma::vec times)
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
    return 1e32;
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
      return 1e32;
    }
    //int ma = y.size();
    arma::mat Ma = Bainv*phylocov.rows(Ka);    
    XU.rows((u*nvar)+Ka) = Ma;
    XY.submat((u*nvar)+Ka,zero_vec) = Bainv * y.subvec(ymin(u),ymax(u));
    XX.submat((u*nvar)+Ka,Ka) = Bainv;
    YU.row(u) = trans(y.subvec(ymin(u),ymax(u))) * Ma;
    arma::mat temp_mat = (trans(y.subvec(ymin(u),ymax(u))) * Bainv * y.subvec(ymin(u),ymax(u)));
    YY(u) = temp_mat(0,0);
    lnW(u) = log(det(phenocov.submat(Ka,Ka)));
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
        return 1e32;
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
      return 1e32;
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
    return 1e32;
  }*/
  return minus2ll;
  //return List::create(XU,XX,YU,YY,XY,lnW);
}

// [[Rcpp::export]]
double threepoint_calc_pheno(arma::vec theta,arma::vec options,arma::vec y,arma::vec ku,arma::vec iu,arma::mat edge,arma::vec edgelength,Rcpp::List uchildren_list,List subset_list,List species_subset,arma::vec tip_combn,arma::vec ymin,arma::vec ymax,List phenocovs,List inv_phenocovs,arma::vec non_optim_transform,arma::vec models,arma::uvec externalEdge,arma::uvec not_externalEdge,arma::vec dist_anc,arma::vec dist_des,double Tmax,double Tmin,int nmodels,arma::vec lower_bounds,arma::vec upper_bounds,arma::vec OU_D,arma::vec times)
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
    return 1e32;
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
        return 1e32;
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
      return 1e32;
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
    return 1e32;
  }*/
  return minus2ll;
  //return List::create(XU,XX,YU,YY,XY,lnW);
}

// [[Rcpp::export]]
arma::mat threepoint_predict(arma::vec theta,arma::vec options,arma::vec y,arma::vec ku,arma::vec iu,arma::mat edge,arma::vec edgelength,Rcpp::List uchildren_list,arma::vec non_optim_transform,arma::vec models,arma::uvec externalEdge,arma::uvec not_externalEdge,arma::vec dist_anc,arma::vec dist_des,double Tmax,double Tmin,int nmodels,arma::vec lower_bounds,arma::vec upper_bounds,arma::vec OU_D,arma::vec times)
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
arma::mat threepoint_nopheno_predict(arma::vec theta,arma::vec options,arma::vec y,arma::vec ku,arma::vec iu,arma::mat edge,arma::vec edgelength,Rcpp::List uchildren_list,arma::vec edgevec,List subset_list,List species_subset,arma::vec tip_combn,arma::vec ymin,arma::vec ymax,arma::vec non_optim_transform,arma::vec models,arma::uvec externalEdge,arma::uvec not_externalEdge,arma::vec dist_anc,arma::vec dist_des,double Tmax,double Tmin,int nmodels,arma::vec lower_bounds,arma::vec upper_bounds,arma::vec OU_D,arma::vec times)
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
      return 1e32;
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
      return 1e32;
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
arma::mat threepoint_phenocorr_predict(arma::vec theta,arma::vec options,arma::vec y,arma::vec ku,arma::vec iu,arma::mat edge,arma::vec edgelength,Rcpp::List uchildren_list,List subset_list,List species_subset,arma::vec tip_combn,arma::vec ymin,arma::vec ymax,arma::vec non_optim_transform,arma::vec models,arma::uvec externalEdge,arma::uvec not_externalEdge,arma::vec dist_anc,arma::vec dist_des,double Tmax,double Tmin,int nmodels,arma::vec lower_bounds,arma::vec upper_bounds,arma::vec OU_D,arma::vec times)
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
      return 1e32;
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
arma::mat threepoint_calc_pheno_predict(arma::vec theta,arma::vec options,arma::vec y,arma::vec ku,arma::vec iu,arma::mat edge,arma::vec edgelength,Rcpp::List uchildren_list,List subset_list,List species_subset,arma::vec tip_combn,arma::vec ymin,arma::vec ymax,List phenocovs,List inv_phenocovs,arma::vec non_optim_transform,arma::vec models,arma::uvec externalEdge,arma::uvec not_externalEdge,arma::vec dist_anc,arma::vec dist_des,double Tmax,double Tmin,int nmodels,arma::vec lower_bounds,arma::vec upper_bounds,arma::vec OU_D,arma::vec times)
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
      return 1e32;
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