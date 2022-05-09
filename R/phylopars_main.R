convert_to_means <- function(df,index_col = 1,sort_vec,FUN = function(X) mean(X,na.rm=TRUE))
{
  ret <- matrix(NA,ncol = ncol(df)-1,nrow = length(unique(df[,index_col])))
  index <- df[,index_col]
  col_counter <- 1
  cols <- colnames(df)[-index_col]
  for(i in (1:ncol(df))[-index_col])
  {
    ret[,col_counter] <- sapply(split(df[,i],index),FUN = FUN)
    col_counter <- col_counter + 1
  }
  rownames(ret) <- names(split(df[,i],index))
  colnames(ret) <- cols
  if(missing(sort_vec)) return(ret) else return(ret[sort_vec,,drop=F])
}

simtraits <- function(ntaxa=15,ntraits=4,nreps=1,nmissing=0,tree,v,anc,intraspecific,model="BM",parameters,nsim=1,return.type="data.frame")
{
  if(model=="OU") model <- "OUfixedRoot"
  if(nmissing>(ntaxa*ntraits*nreps)) nmissing <- round(runif(1,ntaxa*ntraits*nreps-1))
  if(missing(tree))
  {
    tree <- pbtree(n=ntaxa)
  } else
  {
    tree <- tree[c("edge","tip.label","edge.length","Nnode")]
    class(tree) <- "phylo"
    ntaxa <- length(tree$tip.label)
  }
  tree <- reorder(tree,"postorder")
  perm_tree <- tree
  if(model!="BM") tree <- transf.branch.lengths(phy = tree,model = model,parameters = parameters)$tree
  if(missing(v))
  {
    v <- matrix(0,ntraits,ntraits)
    npars <- length(v[upper.tri(v,TRUE)])
    pars <- rnorm(npars)
    v[upper.tri(v,TRUE)] <- pars
    v <- t(v)%*%v
  } else ntraits = length(diag(v))
  if(missing(anc))
  {
    anc <- rep(0,ntraits)
  } else if(length(anc)==1) anc <- rep(anc,ntraits)
  anc <- as.double(anc)
  
  if(missing(intraspecific)) if(nreps>1) intraspecific <- 0.1 else intraspecific <- 0
  if(length(intraspecific)==1)
  {
    opt <- 1
    #intraspecific <- matrix(rep(intraspecific,ntraits*ntaxa),nrow = ntaxa,ncol = ntraits)
    intraspecific <- diag(rep(intraspecific,ntraits),nrow = ntraits)
  } else if(length(intraspecific)==ntraits)
  {
    opt <- 2
    #intraspecific <- t(matrix(rep(intraspecific,ntaxa),nrow=ntraits))
    intraspecific <- diag(intraspecific,nrow = ntraits)
  } else
  {
    opt <- 3
  }
  #anc_mat <- matrix(1,ntaxa) %*% anc
  
  # copied from geiger package, which is set to be archived on May 8, 2022
  
  phe <- array(data = 0,dim=c(length(tree$tip.label)+tree$Nnode,ntraits,nsim))

  for(i in nrow(tree$edge):1)
  {
    anc_node <- tree$edge[i,1]
    des_node <- tree$edge[i,2]
    phe[des_node,,] <- t(mvrnorm(n = nsim,mu = rep(0,ntraits),Sigma = v*tree$edge.length[i])) + phe[anc_node,,]
  }
  #Xall <- sim.char(phy = tree,par = v,nsim = nsim,root = 0)
  Xall <- phe[1:length(tree$tip.label),,,drop=FALSE]
  dimnames(Xall) <- list(tree$tip.label,paste("V",1:ntraits,sep=""),NULL)
  #colnames(Xall) <- paste("V",1:ntraits,sep="")
  #if(nreps==1 & nmissing==0 & nsim==1)
  #{
  #  if(return.type=="matrix") return(list(trait_data=Xall[,,1,drop=FALSE],tree=perm_tree,sim_tree=tree,original_X=Xall[,,1,drop=FALSE])) else
  #    return(list(trait_data=data.frame(species=rownames(Xall[,,1,drop=FALSE]),Xall[,,1]),tree=perm_tree,sim_tree=tree,original_X=Xall[,,1,drop=FALSE]))
  #} else if(nreps==1 & nmissing==0) 
  #{
  #  if(return.type=="matrix")
  #  {
  #    return(list(trait_data=lapply(apply(Xall,3,function(X) list(X)),function(X) X[[1]]),tree=perm_tree,sim_tree=tree,original_X=lapply(apply(Xall,3,function(X) list(X)),function(X) X[[1]])))
  #  } else
  #    return(list(trait_data=lapply(apply(Xall,3,function(X) list(X)),function(X) data.frame(species=rownames(X[[1]]),X[[1]])),tree=perm_tree,sim_tree=tree,original_X=lapply(apply(Xall,3,function(X) list(X)),function(X) X[[1]])))
  #}
  
  X <- original_X <- rep(list(matrix(0,ntaxa*nreps,ntraits)),nsim)
  for(j in 1:nsim)
  {
    #for(ii in 1:nrow(Xall[,,j,drop=FALSE]))
    #{
    #  Xall[,,j] <- Xall[,,j] + anc_mat
    #}
    #Xall[,,j] <- as.matrix(Xall[,,j,drop=FALSE]) + anc_mat
    Xall[,,j] <- scale(Xall[,,j],-anc,FALSE)
    original_X[[j]] <- Xall[,,j,drop=FALSE]
    simdat_j <- matrix(NA,nrow = ntaxa*nreps,ncol = ntraits)
    
    #if(nreps==1)
    #{
    #  X[[j]][1:(ntraits*ntaxa)] <- original_X[[j]][1:(ntraits*ntaxa)] <- Xall[,,j]
    #} else
    #{
    
      for(ii in 1:nrow(Xall[,,j,drop=FALSE]))
      {
        X[[j]][(1:nreps-1)*ntaxa+ii,] <- mvrnorm(n = nreps,mu = Xall[ii,,j,drop=FALSE],Sigma = intraspecific)
      }
    
      #simdat_j <- t(apply(as.matrix(Xall[,,j,drop=FALSE]),1,function(X) mvrnorm(n = nreps,mu = X,Sigma = intraspecific)))
      #for(jj in 1:nreps)
      #{
        #X[[j]][1:ntaxa + (jj-1)*(ntaxa),] <- rnorm(n = ntraits*ntaxa,mean = Xall[,,j],sd = intraspecific)
        #X[[j]][1:ntaxa + (jj-1)*(ntaxa),] <- simdat_j[,1:ntraits*jj]
      #  X[[j]][1:ntaxa + (jj-1)*(ntaxa),] <- simdat_j[,(1:ntraits-1)*nreps+jj]
        
        #return(X)
      #}
    #}
    X[[j]][sample(1:length(X[[j]]),nmissing)] <- NA
    colnames(X[[j]]) <- paste("V",1:ncol(X[[j]]),sep = "")
    species <- rep(rownames(Xall[,,j,drop=FALSE]),nreps)
    rownames(X[[j]]) <- 1:nrow(X[[j]])
    X[[j]] <- data.frame(species=species,X[[j]])
    if(nreps==1) rownames(X[[j]]) <- species
  }
  if(nsim==1) list(trait_data=X[[1]],tree=perm_tree,sim_tree=tree,original_X=original_X[[1]]) else
    list(trait_data=X,tree=perm_tree,sim_tree=tree,original_X=original_X)
}

write.phylopars <- function(trait_data,tree,data_file,tree_file,species_identifier="species")
{ 
  trait_data <- as.data.frame(trait_data)
  if(!any(tolower(colnames(trait_data))=="species")) 
  {
    if(species_identifier %in% colnames(trait_data))
    {
      colnames(trait_data)[which(colnames(trait_data)==species_identifier)] <- "species"
    } else stop("trait_data must be a data.frame with a column named species")
  }
  species_col <- match("species",tolower(colnames(trait_data)))[1]
  cols <- rep(FALSE,ncol(trait_data))
  for(i in 1:ncol(trait_data))
  {
    cols[i] <- is.numeric(trait_data[,i])
  }
  trait_data <- data.frame(species=trait_data[,species_col],trait_data[,which(cols)])
  colnames(trait_data)[1] <- ""
  species <- unique(trait_data[,1])
  if(!all(species %in% tree$tip.label))
  {
    trait_data <- trait_data[-which(!(species %in% tree$tip.label)),]
  }
  write.table(x = trait_data,file = trait_data,sep="\t",row.names = FALSE,quote = FALSE,na = "")
  write.tree(tree,tree_file)
}