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

optim2 <- function(args)
{
  names(args)[1] <- "theta"
  f <- function(theta,args)
  {
    AB <- Rphylopars:::convert_pars2(theta,args$options,args$fixed_phylocov,args$fixed_phenocov,check1=0)
    A <- AB[1:args$options["nvar"],1:args$options["nvar"]]
    B <- AB[1:args$options["nvar"]+args$options["nvar"],1:args$options["nvar"]]
    if(args$options["pheno_error"] > 0)
    {
      tr <- try(as.matrix(nearPD(B,corr = FALSE,keepDiag = FALSE)$mat),silent=TRUE)
      B <- if(class(tr)=="try-error") B else tr
    }
    tr <- try(as.matrix(nearPD(A,corr = FALSE,keepDiag = FALSE)$mat),silent=TRUE)
    A <- if(class(tr)=="try-error") A else tr
    args[[1]] <- Rphylopars:::depost(A,B,options=args$options)
    do.call(args[[2]],args[-2])
  }
  optim_args <- args
  names(optim_args)[[1]] <- "par"
  do.call(optim,optim_args)
}

optim_BFGS_NM <- function(args,optim_limit=50)
{
  abs_tol <- 1e-1
  args_BFGS <- c(args,list(method="BFGS"))#,control=list(ndeps=rep(sqrt(.Machine$double.eps),length(args$par)),maxit=length(args$par)*200,abstol=1e-5)))
  #o <- try(suppressWarnings(do.call(optim,args_BFGS)),silent=TRUE)
  o <- try(suppressWarnings(optim2(args_BFGS)),silent=TRUE)
  if(class(o)!="try-error")
  {
    args$par <- o$par
  } else
  {
    #o <- suppressWarnings(do.call(optim,args))
    o <- try(suppressWarnings(optim2(args)),silent=TRUE)
    args$par <- o$par
  }
  #o2 <- suppressWarnings(do.call(optim,args))
  o2 <- try(suppressWarnings(optim2(args)),silent=TRUE)
  counter <- 1
  while(abs(o$value-o2$value)>abs_tol & counter<optim_limit)
  {
    counter <- counter+1
    args_BFGS$par <- o2$par
    #o <- try(suppressWarnings(do.call(optim,args_BFGS)),silent=TRUE)
    o <- try(suppressWarnings(optim2(args_BFGS)),silent=TRUE)
    if(class(o)!="try-error")
    {
      args$par <- o$par
    } else
    {
      #o <- suppressWarnings(do.call(optim,args))
      o <- try(suppressWarnings(optim2(args)),silent=TRUE)
      args$par <- o$par        
    }
    if(abs(o$value-o2$value)<=abs_tol) return(o)
    #o2 <- suppressWarnings(do.call(optim,args))
    o2 <- try(suppressWarnings(optim2(args)),silent=TRUE)
  }
  if(abs(o$value-o2$value)>abs_tol)
  {
    warning(paste("Difference between BFGS and Nelder-Mead:",abs(o$value-o2$value),"\nMay not have reached optimal solution. Consider running phylopars.rerun()."))
  }
  return(o2)
}

phylopars.rerun <- function(PPE)
{
  o2 <- optim_BFGS_NM(PPE$optim_args,PPE$optim_limit)
  PPE$optim_args$par <- o2$par
  if(PPE$options["REML"]==1) minus2ll <- o2$value*2+((PPE$options["nob"]-PPE$options["nvar"])*log(2*pi)) else minus2ll <- o2$value*2+((PPE$options["nob"])*log(2*pi))
  PPE$minus2ll <- minus2ll
  PPE$phylopars_logl <- o2$value
  PPE$theta <- o2$par
  ret <- list(minus2ll=minus2ll,phylopars_logl=o2$value,theta=o2$par,options=PPE$options,norms=PPE$norms,offsets=PPE$offsets)
  PPE$pars <- get_final_pars(ret)
  class(PPE) <- "phylopars"
  PPE
}

phylopars <- function(trait_data,tree,model="BM",pheno_error=TRUE,phylo_correlated=TRUE,pheno_correlated=FALSE,calc_pheno=FALSE,calc_pheno_auto_n=20,calc_pheno_phenocovs,use_means=FALSE,species_identifier="species",verbose=FALSE,phylocov_start,phenocov_start,theta_start,model_start,skip_optim=FALSE,REML=TRUE,optim_limit=50,BM_first=TRUE,usezscores=TRUE)
{
  if(use_means & !calc_pheno) calc_pheno <- TRUE
  if(!missing(calc_pheno_auto_n) | !missing(calc_pheno_phenocovs)) calc_pheno <- TRUE
  f_args <- as.list(match.call())
  ret_YY <- FALSE
  tree <- reorder(tree,"postorder")
  n <- nspecies <- length(tree$tip.label)
  OU_D <- numeric(nspecies)
  names(OU_D) <- tree$tip.label
  phy <- reorder(tree,"pruningwise")
  times <- pruningwise.branching.times(phy)
  Tmax <- max(times)
  
  if(BM_first & (model[1]!="BM" | length(model)>1))
  {
    temp_args <- f_args
    temp_args[["model"]] <- "BM"
    temp_args <- temp_args[2:length(temp_args)]
    BM_p <- do.call(phylopars,temp_args)
    raw_BM <- BM_p
    temp_args[["phylocov_start"]] <- BM_p$pars[[1]]
    temp_args[["phenocov_start"]] <- BM_p$pars[[2]]
    temp_args[["BM_first"]] <- FALSE
    strt <- c(alpha=1e-7/Tmax,alpha=1e-7/Tmax,lambda=1e-7,kappa=1e-6,delta=1e-6,rate=0)
    strt1 <- double()
    if("OU" %in% model) strt1 <- c(strt1,alpha=strt[1])
    if("lambda" %in% model) strt1 <- c(strt1,lambda=strt[3])
    if("kappa" %in% model) strt1 <- c(strt1,kappa=strt[4])
    if("delta" %in% model) strt1 <- c(strt1,delta=strt[5])
    if("EB" %in% model) strt1 <- c(strt1,rate=strt[6])
    temp_args[["model"]] <- model
    temp_args[["model_start"]] <- strt1
    BM_p <- do.call(phylopars,temp_args)
  }
  models <- c(OUfixedRoot=as.integer(("OU" %in% model) | ("OUfixedRoot" %in% model)),OUrandomRoot=as.integer("OUrandomRoot" %in% model),
              lambda=as.integer("lambda" %in% model),
              kappa=as.integer("kappa" %in% model),delta=as.integer("delta" %in% model),EB=as.integer("EB" %in% model))
  if(models[1]==1 & models[2]==2) stop("Cannot have both OUfixedRoot and OUrandomRoot.")
  
  if(any(grepl("OU",model)))
  {
    if("OU" %in% model)
    {
      warning("OU should be set to either OUfixedRoot or OUrandomRoot. Assigning OU to OUfixedRoot.")
    } else if(!(("OUfixedRoot" %in% model) | ("OUrandomRoot" %in% model))) stop("Unknown OU model specification.")
    if (!is.ultrametric(phy))# & missing(OU_D)) # adapated from phylolm
    {
      externalEdge <- which(phy$edge[,2] <= nspecies)
      OU_D <- numeric(nspecies) # adjustments to external branck lengths
      stop("Tree is not ultrametric (required for OU).") # Currently doesn't work
      flag <- 1
      dis <- pruningwise.distFromRoot(phy) # has all nodes
      OU_D <- max(dis[1:nspecies]) - dis[1:nspecies]
      OU_D <- OU_D - mean(OU_D)
      des <- phy$edge[, 2]
      externalEdge = (des <= n)
      
      phy$edge.length[externalEdge] <- phy$edge.length[externalEdge] + OU_D[des[externalEdge]]
      tree <- reorder(phy,"postorder")
    } else
      #  if(missing(OU_D))
      #{
      #  OU_D <- numeric(nspecies)
      #} else
    {
      if (length(OU_D)!=nspecies) stop("OU_D should be a vector with one term for each tip in the tree.")
      if (is.null(names(OU_D))) stop("OU_D is lacking names (tip labels).")
      ordr <- match(phy$tip.label, names(OU_D))
      if (sum(is.na(ordr))>0) stop("names of OU_D do not match the tree tip labels.")
      OU_D <- OU_D[ordr,drop=F]
    }
  }
  
  #if(missing(OU_D)) OU_D <- numeric(nspecies)
  
  anc <- tree$edge[, 1]
  des <- tree$edge[, 2]
  #externalEdge <- (tree$edge[,2] <= nspecies)
  
  #externalEdge <- match(1:nspecies,tree$edge[,2])-1
  #not_externalEdge <- match((nspecies+2):dim(tree$edge)[1],tree$edge[,2])-1
  externalEdge <- tree$edge[,2]<=nspecies
  not_externalEdge <- which(!externalEdge)-1
  externalEdge <- which(externalEdge)-1
  #not_externalEdge <- which(!externalEdge)-1
  #externalEdge <- which(externalEdge)-1
  #if(!pheno_error)
  #{
  externalEdge <- match(1:nspecies,tree$edge[,2])-1
  #}
  #not_externalEdge <- match((nspecies+2):(tree$Nnode+nspecies),tree$edge[,2])-1
  distFromRoot <- numeric(phy$Nnode + nspecies)
  for (i in length(phy$edge.length):1) distFromRoot[phy$edge[i, 2]] <- distFromRoot[phy$edge[i,1]] + phy$edge.length[i]
  Tmean <- mean(distFromRoot[1:nspecies])
  Tmin <- min(distFromRoot[1:nspecies])
  dist_anc <- distFromRoot[anc]
  dist_des <- distFromRoot[des]
  starting.values.default <- c(alpha=0.5/Tmax,alpha=0.5/Tmax,lambda=0.5,kappa=0.5,delta=0.5,rate=-1/Tmax)
  bounds.default = matrix(c(1e-7/Tmax,50/Tmax,1e-7/Tmax,50/Tmax,1e-7,1,1e-6,1,1e-5,3,-3/Tmax,0), ncol=2, byrow=TRUE)
  #starting.values.default = c(0,0,1,1,1,0)
  #bounds.default = matrix(c(-1e-6,50/Tmax,-1e-6,50/Tmax,-1e-6,1+1e-6,-1e-6,1+1e-6,-1e-6,3,-3/Tmax,1+1e-6), ncol=2, byrow=TRUE)
  rownames(bounds.default) = c("alpha","alpha","lambda","kappa","delta","rate")
  if(!missing(model_start)) 
  {
    model_index <- which(!is.na(match(names(starting.values.default),names(model_start))))
    starting.values.default[model_index] <- model_start
  }
  colnames(bounds.default) = c("min","max")
  starting.values.default <- -log(-(-bounds.default[,2]+starting.values.default)/(-bounds.default[,1]+starting.values.default))
  lower_bounds <- bounds.default[as.logical(models),1]
  upper_bounds <- bounds.default[as.logical(models),2]
  nmodels <- length(which(models==1))
  model_args <- list(lower_bounds=lower_bounds,upper_bounds=upper_bounds,models=models,externalEdge=externalEdge,not_externalEdge=not_externalEdge,dist_anc=dist_anc,dist_des=dist_des,Tmax=Tmax,Tmin=Tmin,nmodels=nmodels,OU_D=OU_D,times=times)
  
  clip <- FALSE
  colnames(trait_data)[which(colnames(trait_data)==species_identifier)] <- "species"
  data_not_tree <- name.check(tree,data.names = unique(trait_data$species))
  if(length(data_not_tree)>1)
  {
    warning(paste("\nDropping trait_data for",data_not_tree[[2]],"(species not found in tree)"),immediate. = TRUE)
    data_not_tree <- data_not_tree[[2]]
    trait_data <- trait_data[-which(!is.na(match(trait_data$species,data_not_tree))),]
  }
  nvar <- ncol(trait_data)-1
  species_col <- which(colnames(trait_data)=="species")
  blank_rows <- which(apply(trait_data[,(1:ncol(trait_data))[-species_col],drop=FALSE],1,function(X) all(is.na(X))))
  if(length(blank_rows>0))
  {
    trait_data <- trait_data[-blank_rows,]
  }
  #missing_taxa <- name.check(tree,data.names = unique(trait_data$species))
  #if(length(missing_taxa) > 1) tree <- drop.tip(tree,missing_taxa[[1]])
  
  trait_data$species <- factor(trait_data$species, levels=tree$tip.label)
  taxa <- tree$tip.label
  ns <- convert_to_means(trait_data,index_col = 1,sort_vec = tree$tip.label,FUN = function(X) length(which(!is.na(X))))
  if(max(ns,na.rm=TRUE)==1)
  {
    if(missing(pheno_error)) pheno_error <- FALSE
    #cat("Single observation per species. Setting pheno_error to FALSE.")
    if(calc_pheno & missing(calc_pheno_phenocovs))
    {
      warning("Single observation per species cannot be used to estimate within-species variance. Setting calc_pheno to FALSE.",immediate. = TRUE)
      calc_pheno <- FALSE
    }
    if(pheno_error) warning("Single observation per species. Consider setting pheno_error to FALSE.",immediate. = TRUE)
  }
  
  if(calc_pheno==TRUE) pheno_correlated <- TRUE
  if(pheno_correlated) pheno_error <- TRUE
  if(!missing(calc_pheno_phenocovs))
    if(!is.list(calc_pheno_phenocovs)) calc_pheno_phenocovs <- rep(list(calc_pheno_phenocovs),length(tree$tip.label))
  if(phylo_correlated) npars <- (nvar^2-nvar)/2+nvar else npars <- nvar
  if(pheno_correlated & !calc_pheno) npars <- npars + ((nvar^2-nvar)/2+nvar) else if(pheno_error & !calc_pheno) npars <- npars + nvar
  trait_data <- trait_data[,c(which(colnames(trait_data)=="species"),which(colnames(trait_data)!="species"))]
  
  trait_data[,1:nvar+1] <- apply(trait_data[,1:nvar+1,drop=F],2,as.double)
  if(!pheno_error)
  {
    rawX <- convert_to_means(trait_data,1,sort_vec = tree$tip.label,FUN = function(X) mean(X,na.rm=TRUE))
    trait_data <- data.frame(species=rownames(rawX),rawX)
  }
  
  means <- apply(trait_data[,2:ncol(trait_data),drop=FALSE],2,function(X) mean(X,na.rm=TRUE))
  sds <- apply(trait_data[,2:ncol(trait_data),drop=FALSE],2,function(X) sd(X,na.rm=TRUE))
  sds[sds==0] <- 1
  ztrait_data <- trait_data
  if(usezscores & !ret_YY)
  {
    for(i in 2:(nvar+1))
    {
      ztrait_data[,i] <- (trait_data[,i] - means[i-1]) / sds[i-1]
    }
    offsets <- means
    norms <- sds
  } else
  {
    offsets <- rep(0.0,nvar)
    norms <- rep(1.0,nvar)
  }
  X <- phylocurve:::convert_to_means(ztrait_data,1,tree$tip.label)
  Xtrait_data <- phylocurve:::convert_to_means(trait_data,1,tree$tip.label)
  vars <- phylocurve:::convert_to_means(trait_data,1,tree$tip.label,function(X) var(X,na.rm=TRUE))
  vars[vars==0] <- 1e-6
  pooled <- colSums((ns-1)*vars/apply(ns,2,function(X) sum(X[X>1 & !is.na(X)])-length(X[X>1 & !is.na(X)])),na.rm=TRUE)
  if(missing(phenocov_start) & all(!is.na(pooled))) 
  {
    if(all(pooled>0))
    {
      phenocov_start <- matrix(0,nvar,nvar)
      diag(phenocov_start) <- pooled
    }
  }
  
  if(missing(calc_pheno_phenocovs) & calc_pheno)
  {
    vars <- vars[tree$tip.label,]
    calc_pheno_phenocovs <- rep(list(matrix(0,nvar,nvar)),nspecies)
    for(i in 1:nspecies)
    {
      diag(calc_pheno_phenocovs[[i]]) <- vars[i,]
    }
    counts <- ns[tree$tip.label,]
    counts[counts<1 | is.na(counts)] <- 1
    for(i in 1:length(calc_pheno_phenocovs))
    {
      is_na <- diag(calc_pheno_phenocovs[[i]])==0 | which(counts[i,]<calc_pheno_auto_n)
      diag(calc_pheno_phenocovs[[i]])[is_na] <- pooled[is_na]
      if(use_means)
      {
        diag(calc_pheno_phenocovs[[i]]) <- diag(calc_pheno_phenocovs[[i]]) / counts[i,]
      }
    }
  }
  if(calc_pheno) 
  {
    if(is.null(names(calc_pheno_phenocovs))) names(calc_pheno_phenocovs) <- tree$tip.label
  }
  
  species_no_data <- apply(X,1,function(X) all(is.na(X)))
  any_species_no_data <- any(species_no_data)
  which_species_no_data <- which(species_no_data)
  null_species <- factor(names(which_species_no_data),levels=tree$tip.label)
  which_species_no_data <- match(null_species,tree$tip.label)
  null_species_index <- match(null_species,trait_data$species)
  
  if(all(complete.cases(X)) & missing(phylocov_start))
  {
    Y <- apply(X,2,pic,phy=multi2di(drop.tip(tree,tip = names(X[,1])[is.na(X[,1])])))
    phylocov_start <- t(Y)%*%Y/nrow(Y) * norms%*%t(norms)
  } else
  {
    pstart <- try(suppressWarnings(mlest(X)$sigmahat),silent=TRUE)
    if(class(pstart)!="try-error")
    {
      if(all(!is.na(pstart))) phylocov_start <- pstart * norms %*% t(norms)
    }
  }
  A <- matrix(0,nvar,nvar)
  for (i in 1:ncol(X)){
    tmp <- pic(X[,i][!is.na(X[,i])],multi2di(drop.tip(tree,tip = names(X[,i])[is.na(X[,i])])))
    A[i,i] <- sum(tmp^2)/length(tmp)
  }
  
  B <- apply(X,2,function(X) sample_sd(X)^2)
  greater_than_1 <- apply(ns,2,function(X) max(X)>1)
  B[!greater_than_1] <- 1e-6
  if(nvar>1)
  {
    B <- diag(B)
  } else
  {
    B <- as.matrix(B)
  }
  
  nedge <- dim(tree$edge)[1]
  nob <-  sum(!is.na(ztrait_data[,2:(nvar+1)]))
  uind <- which(!is.na(ztrait_data[,2:(nvar+1)]))
  nspecies <- length(tree$tip.label)
  y <- as.matrix(ztrait_data[,2:(nvar+1)])[uind]
  ku <- iu <- double()
  for(i in 1:nvar)
  {
    ku <- c(ku,rep(i,sum(!is.na(ztrait_data[,i+1]))))
    iu <- c(iu,ztrait_data$species[which(!is.na(ztrait_data[,i+1]))])
  }
  if(!pheno_error | pheno_correlated)
  {
    nind <- nrow(trait_data)
  } else
  {
    nind <- length(y)
  }
  ns <- c(nspecies,nedge,nvar,nob)
  options <- c(phylo_correlated=as.integer(phylo_correlated),pheno_error=as.integer(pheno_error)+as.integer(pheno_correlated)+as.integer(calc_pheno),nvar=as.integer(nvar)
               ,clip=as.integer(clip),verbose=as.integer(verbose),nspecies=as.integer(nspecies),
               nedge=as.integer(nedge),nob=as.integer(nob),nn=as.integer(nob+nedge+1),nind=nind,REML=as.integer(REML),npars=npars)
  ku <- ku-1
  iu <- iu-1
  edge <- tree$edge-1
  edgelength <- tree$edge.length
  uchildren <- vector("list",nedge+1)
  for(e in 0:nedge)
  {
    if (e==nedge)
    {  
      i = nspecies
    } else { i = edge[e+1,2] } # i=index of child node
    if(!pheno_error | pheno_correlated)
    {
      if((i+1)<=nspecies)
      {
        if(any(which_species_no_data==(i+1)))
        {
          uchildren[[e+1]] <- rbind(1)
        } else
        {
          temp <- cbind(which(trait_data$species==tree$tip.label[i+1]))
          uchildren[[e+1]] <- rbind(length(temp)+1,temp - 1)          
        }
      } else
      {
        uchildren[[e+1]] = condition_f(iu,i,nspecies,edge,nind)
      }
    } else
    {
      if(any(which_species_no_data==(i+1)))
      {
        uchildren[[e+1]] <- rbind(1)
      } else
      {
        uchildren[[e+1]] <- condition_f(iu,i,nspecies,edge,nob)
      }
    }
  }
  
  if(!pheno_error | pheno_correlated)
  {
    options["nn"] <- options["nind"]+options["nedge"]+1
    al <- apply(trait_data[,1:nvar+1,drop=F],1,function(X) which(!is.na(X))-1)
    if(class(al)=="matrix")
    {
      un <- list(al[,1])
      al <- rep(un,nrow(trait_data))
    } else
    {
      un <- unique(al)
    }
    inds <- lapply(un,function(Y) as.double(which(sapply(al,function(X) identical(X,Y)))))
    tip_combn <- vector("double",length=nrow(trait_data))
    for(i in 1:length(inds))
    {
      tip_combn[inds[[i]]] <- i-1
    }
    ylen <- as.double(sapply(al,length))
    ymin <- Reduce('+',ylen,1,accumulate=TRUE)[1:length(ylen)]
    ymax <- ymin+ylen
    ymax[ymin>ymax] <- ymin[ymin>ymax]
    ymin <- ymin-1
    ymax <- ymax-2
    ymin[ymin>(sum(ylen)-1)] <- (sum(ylen)-1)
    ymax[ymax>(sum(ylen)-1)] <- (sum(ylen)-1)
    ymax[ymin>ymax] <- ymin[ymin>ymax]
    
    if(!pheno_error)
    {
      edgevec <- tree$edge.length[match(match(trait_data$species,tree$tip.label),tree$edge[,2])]
      if(any(species_no_data) & FALSE)
      {
        edgevec <- edgevec[names(edgevec)!=null_species]
      }
    }
    y <- as.double(t(ztrait_data[,1:nvar+1,drop=FALSE]))
    y <- y[!is.na(y)]
    if(any_species_no_data)
    {
      zero_un <- match(0,sapply(un,length))
      if(length(zero_un)>0 & !is.na(zero_un))
      {
        un[[zero_un]] <- 0
        zero_al <- which(sapply(al,length)==0)
        al[zero_al] <- un[zero_un]
        tip_combn[which(tip_combn==(zero_un-1))] <- zero_un-1
      }
    }
  }
  if(use_means)
  {
    mean_data <- convert_to_means(trait_data,1,tree$tip.label)
    mean_data <- data.frame(species=rownames(mean_data),mean_data)
    f_args[["use_means"]] <- FALSE
    f_args[["trait_data"]] <- mean_data
    f_args[["calc_pheno"]] <- TRUE
    calc_pheno_phenocovs <- lapply(calc_pheno_phenocovs,function(X) X)
    f_args[["calc_pheno_phenocovs"]] <- calc_pheno_phenocovs
    return(do.call(phylopars,f_args[2:length(f_args)]))
  }
  start0 <- get_starting_pars(A,B,options)
  if(!missing(phylocov_start) | !missing(phenocov_start))
  {
    if(!missing(phylocov_start)) A <- phylocov_start else A <- convert_pars(start0,options)[1:nvar,]
    if(!missing(phenocov_start)) B <- phenocov_start else B <- convert_pars(start0,options)[1:nvar+nvar,]
    start0 <- depost(A,B,list(options=options,norms=norms,theta=start0))
  } else if(!missing(theta_start))
  {
    start0 <- theta_start
  }
  
  if(calc_pheno)
  {
    phenocovs <- calc_pheno_phenocovs[match(match(trait_data$species,tree$tip.label),match(names(calc_pheno_phenocovs),tree$tip.label))]
    inv_phenocovs <- vector("list",length(phenocovs))
    phenocovs <- lapply(phenocovs,function(X) X / (norms%*%t(norms)))
    lens <- sapply(al,length)
    for(i in 1:length(phenocovs))
    {
      phenocovs[[i]] <- phenocovs[[i]][al[[i]]+1,al[[i]]+1,drop=FALSE]
      inv_phenocovs[[i]] <- try_inv(phenocovs[[i]],lens[i])      
    }
  }
  constant <- ifelse(REML,((nob-nvar)*log(2*pi)),((nob)*log(2*pi)))
  options <- c(options,constant = constant,ret_YY=0)
  start0 <- c(start0,starting.values.default[models==1])
  if(calc_pheno)
  {
    args <- list(par=start0,fn=threepoint_calc_pheno,options = options,y = y,ku = ku,iu = iu,edge = edge,edgelength = edgelength,uchildren_list = uchildren,subset_list = un,species_subset = al,tip_combn = tip_combn,ymin = ymin,ymax = ymax,phenocovs = phenocovs, inv_phenocovs = inv_phenocovs)
  } else if(!pheno_error)
  {
    args <- list(par=start0,fn=threepoint_nopheno,options = options,y = y,ku = ku,iu = iu,edge = edge,edgelength = edgelength,uchildren_list = uchildren,edgevec = edgevec,subset_list = un,species_subset = al,tip_combn = tip_combn,ymin = ymin,ymax = ymax)
  } else if(pheno_correlated)
  {
    args <- list(par=start0,fn=threepoint_phenocorr,options = options,y = y,ku = ku,iu = iu,edge = edge,edgelength = edgelength,uchildren_list = uchildren,subset_list = un,species_subset = al,tip_combn = tip_combn,ymin = ymin,ymax = ymax)
  } else
  {
    args <- list(par=start0,fn=threepoint,options = options,y = y,ku =ku,iu = iu,edge = edge,edgelength = edgelength,uchildren_list = uchildren)
  }
  args <- c(args,model_args)
  if(ret_YY)
  {
    args$options["ret_YY"] <- 1
    return(do.call(args$fn,c(list(theta=start0),args[-1:-2])))
  }
  if(skip_optim)
  {
    o2 <- list(par=start0,value=do.call(args$fn,c(list(theta=start0),args[-1:-2])))
  } else
  {
    o2 <- optim_BFGS_NM(args,optim_limit)
    if(BM_first & (model[1]!="BM" | length(model)>1))
    {
      if(BM_p$phylopars_logl < o2$value)
      {
        o2$par <- BM_p$theta
        o2$value <- BM_p$phylopars_logl
      }
    }
  }
  
  if(REML) minus2ll <- o2$value*2+((nob-nvar)*log(2*pi)) else minus2ll <- o2$value*2+((nob)*log(2*pi))
  ret <- list(minus2ll=minus2ll,phylopars_logl=o2$value,theta=o2$par,options=options,norms=norms,offsets=offsets)
  PPE <- c(list(pars=get_final_pars(ret),tree=tree,trait_data=trait_data,ztrait_data=ztrait_data,optim_limit=optim_limit),ret)
  PPE <- c(PPE,list(optim_args=args))
  
  if(calc_pheno)
  {
    calc_pheno_phenocovs <- calc_pheno_phenocovs[sapply(calc_pheno_phenocovs,function(X) if(any(is.na(X))) FALSE else TRUE)]
    
    PPE$pars$Phenotypic <- Reduce('+',calc_pheno_phenocovs) / length(calc_pheno_phenocovs)
    PPE$calc_pheno_phenocovs <- calc_pheno_phenocovs
  }
  if(any(models>0))
  {
    PPE$model <- (PPE$optim_args$upper_bounds - PPE$optim_args$lower_bounds)/(1+exp(-PPE$theta[(PPE$options["npars"]+1):length(PPE$theta)])) + PPE$optim_args$lower_bounds
    if(any(abs(PPE$model-PPE$optim_args$upper_bounds)<1e-6 | abs(PPE$model-PPE$optim_args$lower_bounds)<1e-6)) warning("Parameter estimates appear at bounds.")
  }
  class(PPE) <- "phylopars"
  return(PPE)
}

sample_sd <- function(x)
{
  x <- x[!is.na(x)]
  (sqrt(var(x)*(length(x)-1)/length(x)))
}

##### converts covariance matrices to a vector for unconstrained optimization
get_starting_pars <- function(A,B,options)
{
  pheno_error <- as.logical(options["pheno_error"])
  phylo_correlated <- as.logical(options["phylo_correlated"])
  pheno_correlated <- options["pheno_error"] == 2
  calc_pheno <- options["pheno_error"] == 3
  nvar <- options["nvar"]
  if(pheno_error & FALSE)
  {
    phylogeneticfraction <- .5
  } else
  {
    phylogeneticfraction <- 1
  }
  A <- A*phylogeneticfraction
  if(pheno_error) B <- B#*(1-phylogeneticfraction)
  diag(A) <- abs(diag(A))
  if(pheno_error) diag(B) <- abs(diag(B))
  A <- as.matrix(nearPD(A,keepDiag = FALSE)$mat)
  B <- as.matrix(nearPD(B,keepDiag = FALSE)$mat)
  if(calc_pheno)
  {
    pheno_error <- pheno_correlated <- FALSE
  }
  if(length(A)>1)
  {
    if(phylo_correlated)
    {
      A <- diag(log(diag(chol(A))))
    } else
    {
      A <- diag(log(diag(A)))
    }
    if(pheno_error) B <- diag(log(diag(chol(B))))*2
  } else
  {
    if(phylo_correlated)
    {
      A <- log(chol(A))
    }
    else
    {
      A <- log(A)
    }
    if(pheno_error) B <- log(chol(B))*2
  }
  
  if(pheno_error & phylo_correlated)
  {
    if(nvar>1)
    {
      if(pheno_correlated)
      {
        start0 <- c(A[upper.tri(A,TRUE)],B[upper.tri(B,TRUE)])
      } else
      {
        start0 <- c(A[upper.tri(A,TRUE)],diag(B))        
      }
    } else
    {
      start0 <- c(A,B)
    }
  } else if(pheno_error)
  {
    if(nvar>1)
    {
      if(pheno_correlated)
      {
        start0 <- c(diag(A),B[upper.tri(B,TRUE)])
      } else
      {
        start0 <- c(diag(A),diag(B))
      }
    } else
    {
      start0 <- c(A,B)
    }
  } else if(phylo_correlated)
  {
    if(nvar>1)
    {
      start0 <- A[upper.tri(A,TRUE)]
    } else
    {
      start0 <- c(A)
    }
  } else
  {
    if(nvar>1)
    {
      start0 <- diag(A)
    } else
    {
      start0 <- A
    }
  }
  return(start0)
}

##### convert theta back to covariance matrices
get_final_pars <- function(PPE,theta)
{
  if(!is.null(PPE$optim_args$fixed_phylocov))
  {
    AB <- Rphylopars:::convert_pars2(PPE$theta,PPE$options,
          if(is.null(PPE$optim_args$fixed_phylocov)) diag(PPE$options["nvar"]) else PPE$optim_args$fixed_phylocov,
          if(is.null(PPE$optim_args$fixed_phenocov)) diag(PPE$options["nvar"]) else PPE$optim_args$fixed_phenocov,
          PPE$optim_args$options["tree_height"])
    A <- AB[1:PPE$options["nvar"],1:PPE$options["nvar"]]
    B <- AB[1:PPE$options["nvar"]+PPE$options["nvar"],1:PPE$options["nvar"]]
  
    #if(PPE$optim_args$options["estim_phylocov"]==0) 
      A <- A * PPE$norms %*% t(PPE$norms)
    #if(PPE$optim_args$options["estim_phenocov"]==0) 
      B <- B * PPE$norms %*% t(PPE$norms)
    
    rownames(A) <- rownames(B) <- colnames(A) <- colnames(B) <- colnames(PPE$pars[[1]])
    return(list(Phylogenetic=A,Phenotypic=B))
  }
  
  if(missing(theta))
  {
    theta <- PPE$theta
  }
  A <- B <- matrix(0,PPE$options["nvar"],PPE$options["nvar"])
  phylo_correlated <- as.logical(PPE$options["phylo_correlated"])
  pheno_error <- as.logical(PPE$options["pheno_error"])
  pheno_correlated <- PPE$options["pheno_error"] == 2
  calc_pheno <- PPE$options["pheno_error"] == 3
  if(calc_pheno) pheno_error <- pheno_corelated <- FALSE
  norms <- PPE$norms
  if(phylo_correlated)
  {
    A[upper.tri(A,TRUE)] <- theta[1:(length(A[upper.tri(A,TRUE)]))]
  } else
  {
    diag(A) <- theta[1:(length(diag(A)))]
  }
  A[lower.tri(A,TRUE)] <- t(A)[lower.tri(A,TRUE)]
  A[upper.tri(A,FALSE)] <- 0
  
  if(dim(A)[1]>1 & phylo_correlated)
  {
    if(phylo_correlated)
    {
      A[2,1] <- A[2,1]*exp(A[1,1])
    }
  }
  diag(A) <- exp(diag(A))
  if(phylo_correlated) A <- A%*%t(A)
  
  if(pheno_error)
  {
    if(phylo_correlated & !pheno_correlated)
    {
      from <- (1+(length(A[upper.tri(A,TRUE)])))
      to <- (length(A[upper.tri(A,TRUE)])+length(diag(B)))
    } else if(phylo_correlated & pheno_correlated)
    {
      from <- (1+(length(A[upper.tri(A,TRUE)])))
      to <- (length(A[upper.tri(A,TRUE)])*2)
    } else
    {
      from <- (1+length(diag(A)))
      if(pheno_correlated)
      {
        to <- (length(diag(A))+length(B[upper.tri(B,TRUE)]))
      } else to <- (length(diag(A))+length(diag(B)))
    }
    if(!pheno_correlated)
    {
      diag(B) <- theta[from:to]
    } else
    {
      B[upper.tri(B,TRUE)] <- theta[from:to]
    }
    B[lower.tri(B,TRUE)] <- t(B)[lower.tri(B,TRUE)]
    B[upper.tri(B,FALSE)] <- 0
    
    if(dim(B)[1]>1 & pheno_correlated)
    {
      if(pheno_correlated)
      {
        B[2,1] <- B[2,1]*exp(B[1,1])
      }
    }
    diag(B) <- exp(diag(B))
    if(pheno_correlated) B <- B%*%t(B)
  } else
  {
    diag(B) <- 0
  }
  if(as.logical(PPE$options["clip"]))
  {
    e <- try(eigen(A),silent=TRUE)
    if(class(e)!="try-error")
    {
      w <- e$values
      v <- e$vectors
      minw <- min(w)
      if(minw < 1e-12)
      {
        w[w<1e-12] <- 1e-12
        A <- v%*%diag(w)%*%t(v)
      }
    }
    e <- try(eigen(B),silent=TRUE)
    if(class(e)!="try-error")
    {
      w <- e$values
      v <- e$vectors
      minw <- min(w)
      if(minw < 1e-12)
      {
        w[w<1e-12] <- 1e-12
        B <- v%*%diag(w)%*%t(v)
      }
    }
  }
  
  A <- A*norms%*%t(norms)
  B <- B*norms%*%t(norms)
  rownames(A) <- rownames(B) <- colnames(A)
  return(list(Phylogenetic=A,Phenotypic=B))
}

simtraits <- function(ntaxa=15,ntraits=4,nreps=1,nmissing=0,tree,v,anc,intraspecific,model="BM",parameters,nsim=1,return.type="data.frame")
{
  if(nmissing>(ntaxa*ntraits*nreps)) nmissing <- round(runif(1,ntaxa*ntraits*nreps-1))
  if(missing(tree))
  {
    tree <- pbtree(n=ntaxa)
  } else ntaxa <- length(tree$tip.label)
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
  
  if(missing(intraspecific)) intraspecific <- 0.1
  if(length(intraspecific)==1)
  {
    opt <- 1
    intraspecific <- matrix(rep(intraspecific,ntraits*ntaxa),nrow = ntaxa,ncol = ntraits)
  } else if(length(intraspecific)==ntraits)
  {
    opt <- 2
    intraspecific <- t(matrix(rep(intraspecific,ntaxa),nrow=ntraits))
  } else opt <- 3
  anc_mat <- matrix(1,ntaxa) %*% anc
  Xall <- sim.char(phy = tree,par = v,nsim = nsim)
  colnames(Xall) <- paste("V",1:ntraits,sep="")
  if(nreps==1 & nmissing==0 & nsim==1)
  {
    if(return.type=="matrix") return(list(trait_data=Xall[,,1],tree=perm_tree,sim_tree=tree)) else
      return(list(trait_data=data.frame(species=rownames(Xall[,,1]),Xall[,,1]),tree=perm_tree,sim_tree=tree))
  } else if(nreps==1 & nmissing==0) 
  {
    if(return.type=="matrix")
    {
      return(list(trait_data=lapply(apply(Xall,3,function(X) list(X)),function(X) X[[1]]),tree=perm_tree,sim_tree=tree))
    } else
      return(list(trait_data=lapply(apply(Xall,3,function(X) list(X)),function(X) data.frame(species=rownames(X[[1]]),X[[1]])),tree=perm_tree,sim_tree=tree))
  }
  
  X <- original_X <- rep(list(matrix(0,ntaxa*nreps,ntraits)),nsim)
  for(j in 1:nsim)
  {
    Xall[,,j] <- Xall[,,j] + anc_mat
    if(nreps==1)
    {
      X[[j]][1:(ntraits*ntaxa)] <- original_X[[j]][1:(ntraits*ntaxa)] <- Xall[,,j]
    } else
    {
      for(jj in 1:nreps)
      {
        original_X[[j]] <- Xall[,,j]
        X[[j]][1:ntaxa + (jj-1)*(ntaxa),] <- rnorm(n = ntraits*ntaxa,mean = Xall[,,j],sd = intraspecific)
      }
    }
    X[[j]][sample(1:length(X[[j]]),nmissing)] <- NA
    colnames(X[[j]]) <- paste("V",1:ncol(X[[j]]),sep = "")
    species <- rep(rownames(Xall[,,j]),nreps)
    rownames(X[[j]]) <- 1:nrow(X[[j]])
    X[[j]] <- data.frame(species=species,X[[j]])
    if(nreps==1) rownames(X[[j]]) <- species
  }
  if(nsim==1) list(trait_data=X[[1]],tree=perm_tree,sim_tree=tree,original_X=original_X[[1]]) else
    list(trait_data=X,tree=perm_tree,sim_tree=tree,original_X=original_X)
}

write.phylopars <- function(trait_data,tree,data_file,tree_file,species_identifier="species")
{ 
  rnms <- rownames(trait_data)
  rownames(trait_data) <- NULL
  trait_data <- as.data.frame(trait_data)
  if(is.null(colnames(trait_data))) colnames(trait_data) <- paste("V",1:ncol(trait_data),sep="")
  featurenames <- colnames(trait_data)
  if(length(featurenames[which(featurenames==species_identifier)])==0)
  {
    if(is.null(rnms))
    {
      warning("No tip labels and no row names. Assuming trait_data are in order of tips with no within-species replicates.")
      rnms <- tree$tip.label
    }
    trait_data[[species_identifier]] <- rnms
    featurenames <- colnames(trait_data)
  }
  featurenames[which(featurenames==species_identifier)] <- "species"
  colnames(trait_data) <- featurenames
  featurenames <- featurenames[featurenames!="species"]
  namecheck <- name.check(tree,data.names=trait_data$species)
  if((namecheck!="OK")[[1]])
  {
    if(length(namecheck$data_not_tree)>0)
    {
      warning(paste("trait_data for",paste(unique(namecheck$data_not_tree),collapse=", "),"not in tree. Pruning from trait_dataset."))
      drop <- which(!is.na(match(trait_data$species,namecheck$data_not_tree)))
      trait_data <- trait_data[-drop,]
    }
    if(length(namecheck$tree_not_data)>0)
    {
      missing_species <- namecheck$tree_not_data
      for(i in 1:length(missing_species))
      {
        trait_data <- rbind(trait_data,NA)
        trait_data$species[nrow(trait_data)] <- missing_species[i]
      }
    }
  }
  if(ncol(trait_data)==2)
  {
    trait_data[,-match("species",colnames(trait_data))] <- as.double(trait_data[,-match("species",colnames(trait_data))])
  } else
  {
    trait_data[,-match("species",colnames(trait_data))] <- apply(trait_data[,-match("species",colnames(trait_data)),drop=F],2,as.double)
  }
  species <- trait_data$species
  spcs <- unique(species)
  trait_data <- trait_data[,-which(colnames(trait_data)=="species"),drop=FALSE]
  out <- matrix(nrow=length(unique(species)),ncol=ncol(trait_data))
  for(i in 1:ncol(trait_data))
  {
    for(j in 1:length(unique(species)))
    {
      out[j,i] <- paste(trait_data[,i][species==spcs[j]],collapse = ";")
      out[j,i] <- gsub("NA","",out[j,i])
      semicolons <- which(is.finite(match(strsplit(out[j,i],"")[[1]],";")))
      if(length(semicolons)>0)
      {
        if(semicolons[length(semicolons)]==length(strsplit(out[j,i],"")[[1]]))
        {
          strng <- strsplit(out[j,i],"")[[1]]
          out[j,i] <- paste(strng[1:(length(strng)-1)],collapse = "")
        }
        if(semicolons[1]==1) out[j,i] <- sub(";","",out[j,i])
      }
      while(grepl(";;",out[j,i]))
      {
        out[j,i] <- sub(";;",";",out[j,i])
      }
      semicolons <- which(is.finite(match(strsplit(out[j,i],"")[[1]],";")))
      if(length(semicolons)>0)
      {
        if(semicolons[length(semicolons)]==length(strsplit(out[j,i],"")[[1]]))
        {
          strng <- strsplit(out[j,i],"")[[1]]
          out[j,i] <- paste(strng[1:(length(strng)-1)],collapse = "")
        }
        if(semicolons[1]==1) out[j,i] <- sub(";","",out[j,i])
      }
    }
  }
  rownames(out) <- unique(species)
  write.table(out,data_file,sep="\t",quote=FALSE,na = "",col.names=NA)
  write.tree(tree,tree_file)
}

phylopars.pca <- function(PPE,mode="cov",pgls_means=FALSE)
{
  # adapted from phyl.pca (phytools package)
  if(PPE$options[2]!=0) warning("Phylogeneic PCA does not take measurement error into account when calculating trait loadings.",immediate. = TRUE)
  VP <- PPE$pars[[1]]
  tree <- PPE$tree
  trait_data <- PPE$trait_data
  if(pgls_means) 
  {
    Y <- phylopars.predict(PPE,nodes = NA,verbose = FALSE)[[1]]
  } else
  {
    Y <- convert_to_means(trait_data,1,tree$tip.label)
  }
  Y <- Y[tree$tip.label,]
  n <- nrow(Y)
  m <- ncol(Y)
  a <- phylopars.predict(PPE,nodes = length(tree$tip.label)+1,tips = NA,verbose = FALSE)[[1]]
  if (mode == "corr") {
    Y = Y/matrix(rep(sqrt(diag(VP)), n), n, m, byrow = T)
    VP = VP/(sqrt(diag(VP)) %*% t(sqrt(diag(VP))))
    threepoint_calc <- three.point.compute(transf.branch.lengths(tree)$tree,Y)
    a <- matrix(threepoint_calc$P1/threepoint_calc$vec11, m, 1)
  }
  
  es = eigen(VP)
  result <- list()
  result$Eval <- diag(es$values)
  result$Evec <- es$vectors
  A <- matrix(rep(a, n), n, m, byrow = T)
  result$S <- (Y - A) %*% result$Evec
  Ccv <- three.point.compute(transf.branch.lengths(tree)$tree,Y-A,result$S)$QP/(n - 1)
  result$L <- matrix(, m, m, dimnames = list(colnames(Y), paste("PC", 
                                                                1:ncol(Y), sep = "")))
  for (i in 1:m) for (j in 1:m) result$L[i, j] <- Ccv[i, j]/sqrt(VP[i, 
                                                                    i] * result$Eval[j, j])
  class(result) <- "phyl.pca"
  result
}

print.phylopars <- function(x, ...)
{
  PPE <- x
  
  cat("Phylogenetic trait variance-covariance\n")
  print(PPE$pars[[1]])
  cat("\n")
  if(x$options["pheno_error"]>0)
  {
    cat("Phenotypic trait variance-covariance\n")
    print(PPE$pars[[2]])
    
    nvar <- PPE$options["nvar"]
    xs <- as.list(PPE$trait_data[,1:nvar+1,drop=FALSE])
    finite <- lapply(xs,is.finite)
    xs <- lapply(xs,function(X) X[is.finite(X)])
    
    featurecount <- unlist(lapply(xs,length))
    featuresum <- unlist(lapply(xs,function(X) mean(sum(X),na.rm=TRUE)))
    featuremean <- featuresum/featurecount
    feature2sum <- unlist(lapply(xs,function(X) mean(sum(X^2),na.rm=TRUE)))
    featurevar <- feature2sum/featurecount - featuremean^2
    cat("\n% variance explained by phlogeny\n")
    percent_var <- (1-diag(PPE$pars[[2]])/featurevar)*100
    names(percent_var) <- colnames(PPE$pars[[1]])
    print(percent_var)
    
  }
  
  if(length(PPE$model)>0)
  {
    cat(paste("",names(PPE$optim_args$models[PPE$optim_args$models==1])," model: ",names(PPE$model)," = ",format(PPE$model),"\n",sep=""))
    cat("\n")
  } else cat("Brownian motion model\n\n")
  
}

summary.phylopars <- function(object, ...)
{
  PPE <- object
  
  
  if(length(PPE$model)>0)
  {
    cat(paste(names(PPE$optim_args$models[PPE$optim_args$models==1])," model: ",names(PPE$model)," = ",format(PPE$model),"\n",sep=""))
    cat("\n")
  } else cat("Brownian motion model\n\n")
  
  ret <- matrix(0,nrow(PPE$pars[[1]]),4,dimnames = list(colnames(PPE$pars[[1]]),c("phylogenetic mean","phylogenetic sd","phenotypic sd",
                                                                                  "% variance explained by phylogeny")))
  ret[,1] <- round(phylopars.predict(PPE,nodes = length(PPE$tree$tip.label)+1,NA,FALSE)[[1]],4)
  ret[,2] <- round(sqrt(diag(PPE$pars[[1]])),4)
  ret[,3] <- round(sqrt(diag(PPE$pars[[2]])),4)
  nvar <- PPE$options["nvar"]
  xs <- as.list(PPE$trait_data[,1:nvar+1,drop=FALSE])
  finite <- lapply(xs,is.finite)
  xs <- lapply(xs,function(X) X[is.finite(X)])  
  featurecount <- unlist(lapply(xs,length))
  featuresum <- unlist(lapply(xs,function(X) mean(sum(X),na.rm=TRUE)))
  featuremean <- featuresum/featurecount
  feature2sum <- unlist(lapply(xs,function(X) mean(sum(X^2),na.rm=TRUE)))
  featurevar <- feature2sum/featurecount - featuremean^2
  percent_var <- (1-diag(PPE$pars[[2]])/featurevar)*100
  ret[,4] <- round(percent_var,2)
  ret
}

depost <- function(A,B,PPE,convert=TRUE,options)
{
  if(missing(PPE))
  {
    PPE <- list(options=options,norms=rep(1,options["nvar"]),theta=numeric(options["npars"]))
  }
  if(missing(options)) options <- PPE$options
  phylo_correlated <- as.logical(PPE$options["phylo_correlated"])
  pheno_error <- as.logical(PPE$options["pheno_error"])
  pheno_correlated <- PPE$options["pheno_error"] == 2
  calc_pheno <- PPE$options["pheno_error"] == 3
  estim_phylocov <- as.logical(options["estim_phylocov"])
  estim_phenocov <- as.logical(options["estim_phenocov"])
  skip_phenocov <- if(!is.na(estim_phenocov)) !estim_phenocov else FALSE
  skip_phylocov <- if(!is.na(estim_phylocov)) !estim_phylocov else FALSE
  if(calc_pheno) pheno_error <- pheno_corelated <- FALSE
  theta <- PPE$theta
  if(convert) A <- as.matrix(A) / (PPE$norms%*%t(PPE$norms))
  if(!skip_phenocov)
  {
    if(pheno_error & !pheno_correlated)
    {
      if(convert) B <- as.matrix(B) / (PPE$norms%*%t(PPE$norms))
      if(length(B)>1)
      {
        diag(B) <- log(diag(B))
        theta[(length(theta)-length(diag(B))+1):length(theta)] <- diag(B)
      } else
      {
        theta[length(theta)] <- log(B)
      }
    } else if(pheno_correlated)
    {
      B <- B / (PPE$norms%*%t(PPE$norms))
      if(length(B)>1)
      {
        B <- t(chol(B))
        diag(B) <- log(diag(B))
        B[2,1] <- B[2,1]/exp(B[1,1])
        B[upper.tri(B,FALSE)] <- t(B)[upper.tri(B,FALSE)]
        theta[(length(theta)-length(which(upper.tri(B,TRUE)))+1):length(theta)] <- B[upper.tri(B,TRUE)]
      } else
      {
        B <- log(B)
        theta[length(theta)] <- B
      }
    }
  }
  if(!skip_phylocov)
  {
    if(phylo_correlated)
    {
      A <- as.matrix(nearPD(A)$mat)
      A <- t(chol(A))
    }
    if(length(A)>1)
    {
      diag(A) <- log(diag(A))
    } else
    {
      A <- log(A)
    }
    if(phylo_correlated & length(A)>1)
    {
      A[2,1] <- A[2,1]/exp(A[1,1])
    }
    A[upper.tri(A,FALSE)] <- t(A)[upper.tri(A,FALSE)]
    if(phylo_correlated)
    {
      theta[1:length(which(upper.tri(A,TRUE)))] <- A[upper.tri(A,TRUE)]
    } else
    {
      if(length(A)>1)
      {
        theta[1:length(diag(A))] <- diag(A)
      } else
      {
        theta[1] <- A
      }
    }
  }
  theta
}

phylopars.predict <- function(PPE,nodes,tips,verbose=FALSE)
{
  usezscores <- TRUE
  theta <- PPE$theta
  trait_data <- PPE$trait_data
  tree <- PPE$tree
  taxa <- tree$tip.label
  options <- PPE$options
  nvar <- options["nvar"]
  nind <- options["nind"]
  means <- PPE$offsets
  sds <- PPE$norms
  pheno_error <- as.logical(options["pheno_error"])
  pheno_correlated <- options["pheno_error"]==2
  calc_pheno <- options["pheno_error"]==3
  if(calc_pheno) pheno_error <- pheno_correlated <- TRUE
  ztrait_data <- PPE$ztrait_data
  X <- convert_to_means(ztrait_data,1,tree$tip.label)
  #X <- apply(ztrait_data[,2:ncol(ztrait_data),drop=FALSE],2,function(X) tapply(X,ztrait_data[,1],function(Y) mean(Y,na.rm=TRUE)))
  
  species_no_data <- apply(X,1,function(X) all(is.na(X)))
  any_species_no_data <- any(species_no_data)
  which_species_no_data <- which(species_no_data)
  null_species <- factor(names(which_species_no_data),levels=tree$tip.label)
  which_species_no_data <- match(null_species,tree$tip.label)
  null_species_index <- match(null_species,trait_data$species)
  
  data_species_factor <- as.integer(factor(trait_data$species,levels=tree$tip.label))
  tree_species_factor <- as.integer(factor(tree$tip.label,levels=tree$tip.label))
  
  nedge <- dim(tree$edge)[1]
  nob <-  sum(!is.na(ztrait_data[,2:(nvar+1)]))
  uind <- which(!is.na(ztrait_data[,2:(nvar+1)]))
  nspecies <- length(tree$tip.label)
  
  if(missing(tips))
  {
    tips <- 1:nspecies
  } else
  {
    if(is.character(tips))
    {
      tips <- match(unique(tips),tree$tip.label)
    }
  }
  if(missing(nodes))
  {
    nodes <- nspecies:PPE$options["nedge"]+1
  }
  reconstruct <- unique(c(tips,nodes))
  reconstruct <- reconstruct[reconstruct>0 & reconstruct<=(PPE$options["nedge"]+1)]
  reconstruct <- reconstruct[is.finite(reconstruct) & is.numeric(reconstruct)]
  reconstruct.length <- length(reconstruct)
  y <- as.matrix(ztrait_data[,2:(nvar+1)])[uind]
  ku <- iu <- double()
  for(i in 1:nvar)
  {
    ku <- c(ku,rep(i,sum(!is.na(ztrait_data[,i+1]))))
    iu <- c(iu,ztrait_data$species[which(!is.na(ztrait_data[,i+1]))])
  }
  if(!pheno_error | pheno_correlated)
  {
    options["nn"] <- options["nind"]+options["nedge"]+1
    al <- apply(trait_data[,1:nvar+1,drop=F],1,function(X) which(!is.na(X))-1)
    if(class(al)=="matrix")
    {
      un <- list(al[,1])
      al <- rep(un,nrow(trait_data))
    } else
    {
      un <- unique(al)
    }
    inds <- lapply(un,function(Y) as.double(which(sapply(al,function(X) identical(X,Y)))))
    tip_combn <- vector("double",length=nrow(trait_data))
    for(i in 1:length(inds))
    {
      tip_combn[inds[[i]]] <- i-1
    }
    ylen <- as.double(sapply(al,length))
    ymin <- Reduce('+',ylen,1,accumulate=TRUE)[1:length(ylen)]
    ymax <- ymin+ylen
    ymax[ymin>ymax] <- ymin[ymin>ymax]
    ymin <- ymin-1
    ymax <- ymax-2
    ymin[ymin>(sum(ylen)-1)] <- (sum(ylen)-1)
    ymax[ymax>(sum(ylen)-1)] <- (sum(ylen)-1)
    ymax[ymin>ymax] <- ymin[ymin>ymax]
    
    if(!pheno_error)
    {
      edgevec <- tree$edge.length[match(match(trait_data$species,tree$tip.label),tree$edge[,2])]
      if(any(species_no_data) & FALSE)
      {
        edgevec <- edgevec[names(edgevec)!=null_species]
      }
    }
    y <- as.double(t(ztrait_data[,1:nvar+1]))
    y <- y[!is.na(y)]
    if(any_species_no_data)
    {
      zero_un <- match(0,sapply(un,length))
      if(length(zero_un)>0 & !is.na(zero_un))
      {
        un[[zero_un]] <- 0
        zero_al <- which(sapply(al,length)==0)
        al[zero_al] <- un[zero_un]
        tip_combn[which(tip_combn==(zero_un-1))] <- zero_un-1
      }
    }
  } else
  {
    nind <- length(y)
  }
  ns <- c(nspecies,nedge,nvar,nob)
  ku <- ku-1
  iu <- iu-1
  
  predict_mat <- predict_sd_mat <- matrix(NA,nedge+1,nvar)
  counter <- 0
  for(temp_root in 1:(nedge+1))
  {
    if(temp_root %in% reconstruct)
    {
      counter <- counter + 1
      if(verbose) cat(paste(counter/reconstruct.length*100," %\n"))
      if(temp_root==(nspecies+1)) 
      {
        temp_tree <- tree
      } else if(temp_root<=nspecies)
      {
        temp_tree <- reroot2(tree,temp_root)
      } else
      {
        temp_tree <- reorder(root(tree,node = temp_root,resolve.root = TRUE),"postorder")
      }
      edgelength <- temp_tree$edge.length
      edge <- temp_tree$edge
      
      edge <- edge-1
      uchildren <- recursion(edge,as.integer(any_species_no_data),which_species_no_data,iu,options,data_species_factor,tree_species_factor)
      
      if(calc_pheno)
      {
        temp <- threepoint_calc_pheno_predict(theta,
                                              options = options,y = y,ku = ku,iu = iu,edge = edge,edgelength = edgelength,
                                              uchildren_list = uchildren,subset_list = un,species_subset = al,tip_combn = tip_combn,
                                              ymin = ymin,ymax = ymax,phenocovs = PPE$optim_args$phenocovs,inv_phenocovs = PPE$optim_args$inv_phenocovs,
                                              lower_bounds=PPE$optim_args$lower_bounds, upper_bounds=PPE$optim_args$upper_bounds, models=PPE$optim_args$models,
                                              externalEdge = PPE$optim_args$externalEdge, not_externalEdge=PPE$optim_args$not_externalEdge,
                                              dist_anc=PPE$optim_args$dist_anc,dist_des=PPE$optim_args$dist_des,Tmax=PPE$optim_args$Tmax,
                                              Tmin=PPE$optim_args$Tmin,
                                              nmodels=PPE$optim_args$nmodels,OU_D=PPE$optim_args$OU_D,times=PPE$optim_args$times)
        
        
      } else if(pheno_error & !pheno_correlated)
      {
        temp <- threepoint_predict(theta,options,y,ku,iu,edge,edgelength,uchildren,
                                   lower_bounds=PPE$optim_args$lower_bounds, upper_bounds=PPE$optim_args$upper_bounds, models=PPE$optim_args$models,
                                   externalEdge = PPE$optim_args$externalEdge, not_externalEdge=PPE$optim_args$not_externalEdge,
                                   dist_anc=PPE$optim_args$dist_anc,dist_des=PPE$optim_args$dist_des,Tmax=PPE$optim_args$Tmax,
                                   Tmin=PPE$optim_args$Tmin,
                                   nmodels=PPE$optim_args$nmodels,OU_D=PPE$optim_args$OU_D,times=PPE$optim_args$times)
      } else if(pheno_correlated)
      {
        temp <- threepoint_phenocorr_predict(theta,
                                             options = options,y = y,ku = ku,iu = iu,edge = edge,edgelength = edgelength,
                                             uchildren_list = uchildren,subset_list = un,species_subset = al,tip_combn = tip_combn,
                                             ymin = ymin,ymax = ymax,
                                             lower_bounds=PPE$optim_args$lower_bounds, upper_bounds=PPE$optim_args$upper_bounds, models=PPE$optim_args$models,
                                             externalEdge = PPE$optim_args$externalEdge, not_externalEdge=PPE$optim_args$not_externalEdge,
                                             dist_anc=PPE$optim_args$dist_anc,dist_des=PPE$optim_args$dist_des,Tmax=PPE$optim_args$Tmax,
                                             Tmin=PPE$optim_args$Tmin,
                                             nmodels=PPE$optim_args$nmodels,OU_D=PPE$optim_args$OU_D,times=PPE$optim_args$times)
      } else
      {
        edgevec <- temp_tree$edge.length[match(match(trait_data$species,temp_tree$tip.label),temp_tree$edge[,2])]
        temp <- threepoint_nopheno_predict(theta,
                                           options = options,y = y,ku = ku,iu = iu,edge = edge,edgelength = edgelength,
                                           uchildren_list = uchildren,edgevec = edgevec, subset_list = un,species_subset = al,tip_combn = tip_combn,
                                           ymin = ymin,ymax = ymax,
                                           lower_bounds=PPE$optim_args$lower_bounds, upper_bounds=PPE$optim_args$upper_bounds, models=PPE$optim_args$models,
                                           externalEdge = PPE$optim_args$externalEdge, not_externalEdge=PPE$optim_args$not_externalEdge,
                                           dist_anc=PPE$optim_args$dist_anc,dist_des=PPE$optim_args$dist_des,Tmax=PPE$optim_args$Tmax,
                                           Tmin=PPE$optim_args$Tmin,
                                           nmodels=PPE$optim_args$nmodels,OU_D=PPE$optim_args$OU_D,times=PPE$optim_args$times)
      }
      predict_mat[temp_root,] <- temp[,1]
      predict_sd_mat[temp_root,] <- temp[,2]
    }
  }
  offsets <- PPE$offsets
  norms <- PPE$norms
  predict_mat <- t(apply(predict_mat,1,function(X) X*norms+offsets))
  predict_sd_mat <- t(apply(predict_sd_mat,1,function(X) X*norms))
  if(nrow(predict_mat)==nvar)
  {
    predict_mat <- t(predict_mat)
    predict_sd_mat <- t(predict_sd_mat)
  }
  rownames(predict_mat) <- paste("n",1:(nedge+1),sep="")
  rownames(predict_mat)[1:nspecies] <- tree$tip.label
  rownames(predict_sd_mat) <- rownames(predict_mat)
  colnames(predict_mat) <- colnames(predict_sd_mat) <- colnames(trait_data)[which(colnames(trait_data)!="species")]
  return(list(predicted=predict_mat[complete.cases(predict_mat),],predicted_sd=predict_sd_mat[complete.cases(predict_sd_mat),]))
}

phylopars.likelihood <- function(PPE,phylocov,phenocov,theta)
{
  PPE$options["verbose"] <- 0
  options <- PPE$options
  nvar <- options["nvar"]
  nob <- options["nob"]
  if(missing(theta))
  {
    if(missing(phenocov)) phenocov <- PPE$pars[[2]]
    if(missing(phylocov)) phylocov <- PPE$pars[[1]]
    theta <- depost(phylocov,phenocov,PPE)
  }
  temp <- do.call(PPE$optim_args$fn,c(list(theta=theta),PPE$optim_args[c(3:length(PPE$optim_args))]))  
  if(options["REML"]==1)  minus2ll <- as.double(temp*2+((nob-nvar)*log(2*pi))) else minus2ll <- as.double(temp*2+((nob)*log(2*pi)))
  return(list(minus2ll=minus2ll,phylopars_logl=temp))
}

phylopars.crossvalidate <- function(PPE,plot=FALSE,verbose=FALSE)
{
  if(class(PPE)=="phylopars.lm") PPE <- PPE$PPE
  cophen <- cophenetic.phylo(PPE$tree)
  usezscores <- TRUE
  theta <- PPE$theta
  trait_data <- PPE$trait_data
  tree <- PPE$tree
  taxa <- tree$tip.label
  options <- PPE$options
  nvar <- options["nvar"]
  nind <- options["nind"]
  means <- PPE$offsets
  sds <- PPE$norms
  pheno_error <- as.logical(options["pheno_error"])
  pheno_correlated <- options["pheno_error"]==2
  calc_pheno <- options["pheno_error"]==3
  if(calc_pheno) pheno_error <- pheno_correlated <- TRUE
  ztrait_data <- PPE$ztrait_data  
  var_names <- colnames(ztrait_data)[1:nvar+1]
  nedge <- dim(tree$edge)[1]
  nob <-  sum(!is.na(ztrait_data[,2:(nvar+1),drop=FALSE]))
  uind <- which(!is.na(ztrait_data[,2:(nvar+1),drop=FALSE]))
  nspecies <- length(tree$tip.label)
  y <- as.matrix(ztrait_data[,2:(nvar+1),drop=FALSE])[uind]
  ku <- iu <- KU <- IU <- double()
  for(i in 1:nvar)
  {
    KU <- c(KU,rep(i,sum(!is.na(ztrait_data[,i+1,drop=FALSE]))))
    IU <- c(IU,ztrait_data$species[which(!is.na(ztrait_data[,i+1,drop=FALSE]))])
  }
  
  cvfeaturematrix <- matrix(NA,nedge+1,nvar)
  counter <- 0
  for(temp_root in 1:nspecies)
  {
    counter <- counter + 1
    if(verbose) cat(paste("\n",counter,"...",sep=""))
    for(j in 1:nvar)
    {
      ztrait_data <- PPE$ztrait_data
      trait_data <- PPE$trait_data
      ztrait_data[ztrait_data$species==tree$tip.label[temp_root],j+1] <- NA
      trait_data[trait_data$species==tree$tip.label[temp_root],j+1] <- NA
      X <- convert_to_means(ztrait_data,1,tree$tip.label)
      #X <- apply(ztrait_data[,2:ncol(ztrait_data),drop=FALSE],2,function(X) tapply(X,ztrait_data[,1],function(Y) mean(Y,na.rm=TRUE)))
      X <- X[PPE$tree$tip.label,,drop=FALSE]
      species_no_data <- apply(X,1,function(X) all(is.na(X)))
      any_species_no_data <- any(species_no_data)
      which_species_no_data <- which(species_no_data)
      null_species <- factor(names(which_species_no_data),levels=tree$tip.label)
      which_species_no_data <- match(null_species,tree$tip.label)
      null_species_index <- match(null_species,trait_data$species)
      
      data_species_factor <- as.integer(factor(ztrait_data$species,levels=tree$tip.label))
      tree_species_factor <- as.integer(factor(tree$tip.label,levels=tree$tip.label))
      nob <-  options["nob"] <- sum(!is.na(ztrait_data[,2:(nvar+1),drop=FALSE]))
      uind <- which(!is.na(ztrait_data[,2:(nvar+1),drop=FALSE]))
      y <- as.matrix(ztrait_data[,2:(nvar+1),drop=FALSE])[uind]    
      iu <- IU[IU!=temp_root | KU!=j]
      ku <- KU[IU!=temp_root | KU!=j]
      if(!pheno_error | pheno_correlated)
      {
        options["nind"] <- nind <- length(which(apply(ztrait_data[,1:nvar+1,drop=FALSE],1,function(X) any(is.finite(X)))))
        y <- as.double(t(ztrait_data[,1:nvar+1,drop=FALSE]))
        y <- y[!is.na(y)]
        
        options["nn"] <- options["nind"]+options["nedge"]+1
        al <- apply(trait_data[,1:nvar+1,drop=F],1,function(X) which(!is.na(X))-1)
        if(class(al)=="matrix")
        {
          un <- list(al[,1])
          al <- rep(un,nrow(trait_data))
        } else
        {
          un <- unique(al)
        }
        inds <- lapply(un,function(Y) as.double(which(sapply(al,function(X) identical(X,Y)))))
        tip_combn <- vector("double",length=nrow(trait_data))
        for(i in 1:length(inds))
        {
          tip_combn[inds[[i]]] <- i-1
        }
        ylen <- as.double(sapply(al,length))
        ymin <- Reduce('+',ylen,1,accumulate=TRUE)[1:length(ylen)]
        ymax <- ymin+ylen
        ymax[ymin>ymax] <- ymin[ymin>ymax]
        ymin <- ymin-1
        ymax <- ymax-2
        ymin[ymin>(sum(ylen)-1)] <- (sum(ylen)-1)
        ymax[ymax>(sum(ylen)-1)] <- (sum(ylen)-1)
        ymax[ymin>ymax] <- ymin[ymin>ymax]
        
        if(any_species_no_data)
        {
          zero_un <- match(0,sapply(un,length))
          if(length(zero_un)>0 & !is.na(zero_un))
          {
            un[[zero_un]] <- 0
            zero_al <- which(sapply(al,length)==0)
            al[zero_al] <- un[zero_un]
            tip_combn[which(tip_combn==(zero_un-1))] <- zero_un-1
          }
        }
        
      } else
      {
        nind <- options["nind"] <- length(y)
        options["nn"] <- options["nind"]+options["nedge"]+1
      }
      ku <- ku-1
      iu <- iu-1
      if(temp_root==(nspecies+1)) 
      {
        temp_tree <- tree
      } else if(temp_root<=nspecies)
      {
        temp_tree <- reroot2(tree,temp_root)
      } else
      {
        temp_tree <- reorder(root(tree,node = temp_root,resolve.root = TRUE),"postorder")
      }
      
      
      edgelength <- temp_tree$edge.length
      edge <- temp_tree$edge
      
      if(temp_root<=nspecies) temp_tree$edge.length[(nedge-1):nedge] <- temp_tree$edge.length[nedge:(nedge-1)]
      
      
      nspecies <- length(temp_tree$tip.label)
      phy <- reorder(temp_tree,"pruningwise")
      times <- pruningwise.branching.times(phy)
      Tmax <- max(times)
      
      OU_D <- PPE$optim_args$OU_D[match(temp_tree$tip.label,tree$tip.label)]
      anc <- temp_tree$edge[, 1]
      des <- temp_tree$edge[, 2]
      externalEdge <- temp_tree$edge[,2]<=nspecies
      not_externalEdge <- which(!externalEdge)-1
      externalEdge <- which(externalEdge)-1
      externalEdge <- match(1:nspecies,temp_tree$edge[,2])-1
      distFromRoot <- numeric(phy$Nnode + nspecies)
      for (i in length(phy$edge.length):1) distFromRoot[phy$edge[i, 2]] <- distFromRoot[phy$edge[i,1]] + phy$edge.length[i]
      Tmean <- mean(distFromRoot[1:nspecies])
      Tmin <- min(distFromRoot[1:nspecies])
      dist_anc <- distFromRoot[anc]
      dist_des <- distFromRoot[des]
      
      model_args <- list(lower_bounds=PPE$optim_args$lower_bounds,upper_bounds=PPE$optim_args$upper_bounds,models=PPE$optim_args$models,
                         externalEdge=externalEdge,not_externalEdge=not_externalEdge,dist_anc=dist_anc,dist_des=dist_des,
                         Tmax=Tmax,Tmin=Tmin,nmodels=PPE$optim_args$nmodels,
                         OU_D=OU_D,times=times)
      
      
      edge <- edge-1
      uchildren <- recursion(edge,as.integer(any_species_no_data),which_species_no_data,iu,options,data_species_factor,tree_species_factor)
      
      if(calc_pheno)
      {
        temp_phenocovs <- PPE$optim_args$phenocovs
        temp_inv_phenocovs <- PPE$optim_args$inv_phenocovs
        temp_mat_nums <- which(ztrait_data$species==tree$tip.label[temp_root])
        for(mat_i in 1:length(temp_mat_nums))
        {
          rows_cols <- (1:ncol(temp_phenocovs[[temp_mat_nums[mat_i]]]))
          drop_rows_cols <- which(colnames(temp_phenocovs[[temp_mat_nums[mat_i]]])==var_names[j])
          if(length(drop_rows_cols)>0) rows_cols <- rows_cols[-drop_rows_cols]
          temp_phenocovs[[temp_mat_nums[mat_i]]] <- temp_phenocovs[[temp_mat_nums[mat_i]]][rows_cols,rows_cols]
          temp_inv_phenocovs[[temp_mat_nums[mat_i]]] <- temp_inv_phenocovs[[temp_mat_nums[mat_i]]][rows_cols,rows_cols]    
        }
        temp_args <- c(list(theta=theta,
                            options = options,y = y,ku = ku,iu = iu,edge = edge,edgelength = edgelength,
                            uchildren_list = uchildren,subset_list = un,species_subset = al,tip_combn = tip_combn,
                            ymin = ymin,ymax = ymax,phenocovs = temp_phenocovs,inv_phenocovs = temp_inv_phenocovs),
                       model_args)
        temp <- do.call(threepoint_calc_pheno_predict,temp_args)    
      } else if(pheno_error & !pheno_correlated)
      {
        temp_args <- c(list(theta=theta,options=options,y=y,ku=ku,iu=iu,edge=edge,edgelength=edgelength,uchildren=uchildren),
                       model_args)
        temp <- do.call(threepoint_predict,temp_args)
      } else if(pheno_correlated)
      {
        temp_args <- c(list(theta=theta,
                            options = options,y = y,ku = ku,iu = iu,edge = edge,edgelength = edgelength,
                            uchildren_list = uchildren,subset_list = un,species_subset = al,tip_combn = tip_combn,
                            ymin = ymin,ymax = ymax),
                       model_args)
        temp <- do.call(threepoint_phenocorr_predict,temp_args)
      } else
      {
        edgevec <- temp_tree$edge.length[match(match(ztrait_data$species,temp_tree$tip.label),temp_tree$edge[,2])]
        if(any(species_no_data) & FALSE)
        {
          edgevec <- edgevec[names(edgevec)!=null_species]
        }
        
        temp_args <- c(list(theta=theta,
                            options = options,y = y,ku = ku,iu = iu,edge = edge,edgelength = edgelength,
                            uchildren_list = uchildren,edgevec = edgevec, subset_list = un,species_subset = al,tip_combn = tip_combn,
                            ymin = ymin,ymax = ymax),
                       model_args)
        temp <- do.call(threepoint_nopheno_predict,temp_args)
      }      
      cvfeaturematrix[temp_root,j] <- temp[j,1]
    }
  }
  if(verbose) cat("\n")
  offsets <- PPE$offsets
  norms <- PPE$norms
  cvfeaturematrix <- t(apply(cvfeaturematrix,1,function(X) X*norms+offsets))
  if(nrow(cvfeaturematrix)==nvar)
  {
    cvfeaturematrix <- t(cvfeaturematrix)
  }
  rownames(cvfeaturematrix) <- paste("n",1:(nedge+1),sep="")
  rownames(cvfeaturematrix)[1:nspecies] <- tree$tip.label
  colnames(cvfeaturematrix) <- colnames(ztrait_data)[which(colnames(ztrait_data)!="species")]
  cvfeaturematrix <- cvfeaturematrix[1:nspecies,,drop=FALSE]
  featurematrix <- convert_to_means(PPE$trait_data,1,tree$tip.label)
  delta <- cvfeaturematrix - featurematrix
  absdelta <- abs(delta)
  absbias <- apply(delta,2,function(X) mean(X,na.rm=TRUE))
  relbias <- apply(delta/featurematrix,2,function(X) mean(X,na.rm=TRUE))
  
  abspheno_errors <- apply(absdelta,2,function(X) mean(X,na.rm=TRUE))
  relpheno_errors <- apply(absdelta/featurematrix,2,function(X) mean(X,na.rm=TRUE))
  xs <- as.list(PPE$trait_data[,1:nvar+1,drop=FALSE])
  finite <- lapply(xs,is.finite)
  xs <- lapply(xs,function(X) X[is.finite(X)])
  
  featurecount <- unlist(lapply(xs,length))
  featuresum <- unlist(lapply(xs,function(X) mean(sum(X),na.rm=TRUE)))
  featuremean <- featuresum/featurecount
  feature2sum <- unlist(lapply(xs,function(X) mean(sum(X^2),na.rm=TRUE)))
  featurevar <- feature2sum/featurecount - featuremean^2
  reldeltas <- delta/featurematrix
  
  r <- apply(featurematrix,2,function(X) length(which(is.finite(X))))
  pheno_error_plots <- alldeltas <- plot_factors <- values <- nnest <- w <- meanest <- vcompr <- deltas_h1 <- deltas_h2 <- reldeltas_h2 <- reldeltas_h1 <- vector("list",nvar)
  
  tab <- matrix(0,2,3,dimnames = list(c("mean bias","mean error"),c("evolutionary model","mean model","nearest neighbor model")))
  ret <- list()
  
  for(i in 1:nvar)
  {
    ret[[i]] <- tab
    vcompr[[i]] <- t(featurematrix[is.finite(featurematrix[,i]),i])
    w[[i]] <- matrix(1/(r[i]-1),r[i],r[i])
    diag(w[[i]]) <- 0
    meanest[[i]] <- cDot(vcompr[[i]],w[[i]])
    deltas_h1[[i]] <- meanest[[i]]-vcompr[[i]]
    reldeltas_h1[[i]] <- deltas_h1[[i]]/vcompr[[i]]
    for(j in 1:length(vcompr[[i]]))
    {
      taxa <- colnames(vcompr[[i]])
      taxa_num <- taxa[j]
      taxa_j <- match(taxa_num,taxa)
      temp_cophen <- names(sort(cophen[colnames(vcompr[[i]]),colnames(vcompr[[i]])][taxa_j,-taxa_j]))
      inc <- 1
      nnest[[i]][j] <- featurematrix[temp_cophen[inc],i]
      while(is.na(nnest[[i]][j]))
      {
        inc <- inc + 1
        nnest[[i]][j] <- featurematrix[temp_cophen[inc],i]
        if(inc>nspecies) break
      }
    }
    deltas_h2[[i]] <- nnest[[i]] - vcompr[[i]]
    reldeltas_h2[[i]] <- deltas_h2[[i]]/vcompr[[i]]
    values[[i]] <- delta[,i][is.finite(delta[,i])]
    alldeltas[[i]] <- c(values[[i]],deltas_h1[[i]],deltas_h2[[i]])
    #if(plot)
    #{
    #  plot_factors[[i]] <- factor(c(rep("Phylopars",length(values[[i]])),rep("Mean",length(values[[i]])),
    #                                rep("NearestNeighbor",length(values[[i]]))))
    #  x_range <- range(alldeltas[[i]])
    #  x_range[1] <- x_range[1]-x_range[1]*.2
    #  x_range[1] <- x_range[1]+x_range[1]*.2
    #  suppressWarnings(print(qplot(alldeltas[[i]], colour=plot_factors[[i]],geom="density",xlim=x_range,ylim=c(0,1),xlab = "error",ylab="Kernel density",main=colnames(ztrait_data)[i+1])+theme_bw()+theme(legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    #                                                                                                                                                                                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))+geom_vline(xintercept=0, linetype="dotted")))
    #}
    ret[[i]][1,] <- c(mean(values[[i]]),mean(deltas_h1[[i]]),mean(deltas_h2[[i]]))
    ret[[i]][2,] <- c(mean(abs(values[[i]])),mean(abs(deltas_h1[[i]])),mean(abs(deltas_h2[[i]])))
  }
  names(ret) <- colnames(ztrait_data)[1:nvar+1]
  mean_bias <- sapply(values,mean)
  mae <- sapply(values,function(X) mean(abs(X)))
  rmse <- sapply(values,function(X) sqrt(mean(X^2))) 
  ret <- list(cross_validation=ret,plot_values=values,plot_deltas=alldeltas,cvfeaturematrix=cvfeaturematrix,featurematrix=featurematrix)
  corrs <- matrix(data = NA,nrow = nvar,ncol = 3,dimnames = list(colnames(cvfeaturematrix),c("R2","N","p")))
  for(i in 1:nvar)
  {
    temp_lm <- lm(cvfeaturematrix[,i]~featurematrix[,i])
    corrs[i,1] <- summary(temp_lm)$r.squared
    corrs[i,2] <- length(which(complete.cases(cbind(featurematrix[,i],cvfeaturematrix[,i]))))
    corrs[i,3] <- summary(temp_lm)$coefficients[2,4]
    if(plot) for(i in 1:nvar)
    {
      plot(featurematrix[,i],cvfeaturematrix[,i],xlab = paste("Observed mean:",colnames(cvfeaturematrix)[i]),ylab=paste("Predicted mean:",colnames(cvfeaturematrix)[i]))
      abline(temp_lm)
    }
  }
  ret$corrs <- corrs
  class(ret) <- "crossvalidate"
  ret
}

plot.crossvalidate <- function(x, ...)
{
  crossvalidate <- x
  for(i in 1:ncol(x$featurematrix))
  {
    plot(x$featurematrix[,i],x$cvfeaturematrix[,i],xlab = paste("Observed mean:",colnames(x$cvfeaturematrix)[i]),ylab=paste("Predicted mean:",colnames(x$cvfeaturematrix)[i]))
    abline(lm(x$cvfeaturematrix[,i]~x$featurematrix[,i]))
  }
  #values <- crossvalidate$plot_values
  #alldeltas <- crossvalidate$plot_deltas
  #plot_factors <- vector("list",length(values))
  #for(i in 1:length(values))
  #{
  #  plot_factors[[i]] <- factor(c(rep("Phylopars",length(values[[i]])),rep("Mean",length(values[[i]])),
  #                                rep("NearestNeighbor",length(values[[i]]))))
  #  x_range <- range(alldeltas[[i]])
  #  x_range[1] <- x_range[1]-x_range[1]*.2
  #  x_range[1] <- x_range[1]+x_range[1]*.2
  #  suppressWarnings(print(qplot(alldeltas[[i]], colour=plot_factors[[i]],geom="density",xlim=x_range,ylim=c(0,1),xlab = "error",ylab="Kernel density",main=names(crossvalidate$cross_validation)[i])+theme_bw()+theme(legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
  #                                                                                                                                                                                                                     panel.background = element_blank(), axis.line = element_line(colour = "black"))+geom_vline(xintercept=0, linetype="dotted")))
  #}
}

print.crossvalidate <- function(x, ...)
{
  #print(x$cross_validation)
  print(x$corrs)
}

pval <- function(r,tree)
{
  N <- length(tree$tip.label)
  r <- abs(r)
  t <- r*sqrt((N-2)/(1-r^2))
  p <- 2*pt(t,N-2-1,lower.tail = FALSE)
  p[p>1] <- 1
  p
}

reroot2 <- function (tree, temp_root) 
{
  tree <- multi2di(tree,random = FALSE)
  tt <- splitTree(tree, list(node = temp_root, bp = 0))
  p <- tt[[1]]
  d <- tt[[2]]
  p <- root(p, outgroup = "NA", resolve.root = T)
  bb <- which(p$tip.label == "NA")
  p$edge.length[which(p$edge[, 2] == bb)] <- 0
  cc <- p$edge[which(p$edge[, 2] == bb), 1]
  dd <- setdiff(p$edge[which(p$edge[, 1] == cc), 2], bb)
  temp_tree <- paste.tree(p, d)
  nedge <- nrow(temp_tree$edge)
  des <- match(temp_tree$edge[(temp_tree$edge[nedge,1]==temp_tree$edge[,1]),2],temp_tree$edge[,2])
  root <- match(temp_root,temp_tree$edge[des,2])
  others <- (1:length(des))[-root]
  temp_tree$edge.length[des[others]] <- temp_tree$edge.length[des[others]] + temp_tree$edge.length[des[root]]
  temp_tree$edge.length[des[root]] <- 1e-16
  temp_tree <- reorder(di2multi(temp_tree),"postorder")
  return(temp_tree)
}

phylopars.lm <- function()
{
  args <- as.list(match.call())
  args <- args[3:length(args)]
  colnames(trait_data)[which(colnames(trait_data)==species_identifier)] <- "species"
  trait_data$species <- factor(trait_data$species, levels=tree$tip.label)
  trait_data <- trait_data[,c(which(colnames(trait_data)=="species"),which(colnames(trait_data)!="species"))]
  original_data <- trait_data
  original_option <- getOption("na.action")
  options(na.action="na.pass")
  mod.mat <- model.matrix(object = formula,data = trait_data)
  intercept <- attr(terms(formula,data = trait_data),"intercept")==1
  if(!intercept) stop("Intercept-free PGLS not currently supported.")
  y_var <- model.frame(formula,data=trait_data)
  var_name <- colnames(y_var)[1]
  y_var <- y_var[,1,drop=FALSE]
  mod.mat <- cbind(mod.mat,y_var)
  colnames(mod.mat)[ncol(mod.mat)] <- var_name
  if(intercept)
  {
    trait_data <- data.frame(species=trait_data$species,mod.mat[,2:ncol(mod.mat)])
    colnames(trait_data) <- c("species",colnames(mod.mat)[2:ncol(mod.mat)])
  } else
  {
    trait_data <- data.frame(species=trait_data$species,mod.mat[,1:ncol(mod.mat)])
    colnames(trait_data) <- c("species",colnames(mod.mat)[1:ncol(mod.mat)])
  }
  options(na.action = original_option)
  args$trait_data <- trait_data
  PPE <- do.call(phylopars,args)
  n <- nspecies <- as.integer(PPE$options["nspecies"])
  means <- phylopars.predict(PPE,nodes = nspecies + 1,NA)[[1]]
  trait_data <- PPE$trait_data
  df.int <- as.integer(intercept)
  k <- ncol(PPE$pars[[1]])
  rdf <- n - k
  covX <- PPE$pars[[1]]
  npred <- ncol(covX)-1
  y_pos <- ncol(covX)
  
  if(ncol(covX)==1 & intercept)
  {
    R2 <- 0
    ts <- ps <- SEs <- NA
    
  } else
  {
    coefs <- solve(covX[1:npred,1:npred,drop=FALSE])%*%covX[1:npred,y_pos,drop=FALSE]
    R2 <- as.double(sum(covX[1:npred,y_pos,drop=FALSE] * coefs) / covX[y_pos,y_pos,drop=FALSE])
  }
  R2adj <- 1-(1-R2)*(n-df.int)/(rdf)
  SST <- as.double(covX[y_pos,y_pos]) * (n-1)
  
  SSreg <- SST * R2
  SSres <- SST - SSreg
  MSres <- SSres / ((rdf))
  if(any(names(PPE$model)=="alpha")) MSres <- 2 * PPE$model["alpha"] * MSres
  sigma <- sqrt(MSres)
  if(!(ncol(covX)==1 & intercept))
  {
    SEs <- sqrt(diag(solve((covX)[1:npred,1:npred,drop=FALSE]) * MSres / (n-1) ))
    ts <- coefs / SEs
    ps <- 2*(1-pt(abs(ts),rdf))
  }
  if(intercept==1)
  {
    if(ncol(covX)==1 & intercept) coefs <- setNames(means[y_pos],"(Intercept)") else
    {
      coefs <- as.double(c(Intercept=means[y_pos] - means[1:npred] %*% solve(covX[1:npred,1:npred,drop=FALSE])%*%covX[1:npred,y_pos,drop=FALSE],coefs))
      SEs <- c(NA,as.double(SEs))
      ts <- c(NA,as.double(ts))
      ps <- c(NA,as.double(ps))
      names(coefs) <- c("(Intercept)",colnames(covX)[1:npred])
    }
  } else names(coefs) <- colnames(covX)[1:npred]
  Fstat <- rdf / (k-df.int)*R2 / (1-R2)
  pval <- as.double(pf(Fstat,k-df.int,rdf,lower.tail = FALSE))
  logdet <- three.point.compute(tree,cbind(setNames(rep(1,n),tree$tip.label)))$logd
  ll <- -n/2 * log(2*pi) - n/2 * log((n-k) * MSres/n) - logdet/2 - n/2
  if(any(is.na(trait_data))) ll <- NA
  ret <- list(coefficients=coefs,SEs=SEs,ts=ts,ps=ps,R2=R2,R2adj=R2adj,sigma=sigma,Fstat=Fstat,pval=pval,df1=k,df2=rdf,dims=list(N=n,p=npred,REML=PPE$options["REML"],df.int=df.int),model=formula,SST=SST,SSres=SSres,SSreg=SSreg,logLik=ll,PPE=PPE,original_data=original_data,covX=covX)
  class(ret) <- "phylopars.lm"
  ret
}

anova.phylopars.lm <- function(object,...)
{
  trait_data <- object$original_data
  covX <- object$covX
  SST <- object$SST
  tlabels <- attr(terms(object$model),"term.labels")
  k <- length(tlabels)
  n <- object$dims$N
  y_pos <- ncol(covX)
  NR <- length(tlabels) + 1
  rss <- resdf <- rep(NA, NR)
  rss[1] <- object$SST
  resdf[1] <- n - 1
  coefs <- coef(object)
  vars <- double()
  accum <- 0
  for (i in 1:(k)) {
    fmla <- as.formula(paste(colnames(object$PPE$pars[[1]])[ncol(object$PPE$pars[[1]])],"~",paste(tlabels[1:i],collapse="+")))
    rdf <-  length(colnames(model.matrix(fmla,data=trait_data))[-1]) - length(vars)
    accum <- accum + rdf
    vars <- colnames(model.matrix(fmla,data=trait_data))[-1]
    coefs <- solve(covX[1:accum,1:accum,drop=FALSE])%*%covX[1:accum,y_pos,drop=FALSE]
    R2 <- as.double(sum(covX[1:accum,y_pos,drop=FALSE] * coefs) / covX[y_pos,y_pos,drop=FALSE])
    SSreg <- SST * R2
    SSres <- SST - SSreg
    rss[i + 1] <- SSres
    resdf[i + 1] <- n - rdf
    resdf[i+1] <- resdf[i]-rdf
  }
  ss <- c(abs(diff(rss)), object$SSres)
  df <- c(abs(diff(resdf)), n - object$df1)
  ms <- ss/df
  fval <- ms/ms[NR]
  P <- pf(fval, df, df[NR], lower.tail = FALSE)
  table <- data.frame(df, ss, ms, f = fval, P)
  table[length(P), 4:5] <- NA
  dimnames(table) <- list(c(tlabels, "Residuals"), c("Df", 
                                                     "Sum Sq", "Mean Sq", "F value", "Pr(>F)"))
  structure(table, heading = c("Analysis of Variance Table\nSequential SS",
                               paste("Response:", deparse(object$model[[2L]]))), 
            class = c("anova", "data.frame"))
}

formals(phylopars.lm) <- c(alist(formula = ),formals(phylopars))

coef.phylopars.lm <- function(object,...)
{
  object$coefficients
}

print.phylopars.lm <- function(x,...)
{
  dd <- x$dims
  mCall <- x$call
  cat("Generalized least squares fit by ")
  cat(ifelse(dd$REML == 1, "REML\n", "maximum likelihood\n"))
  cat("  Model:", deparse(x$model), "\n")
  if(length(x$PPE$model)>0)
  {
    cat(paste("\n  ",names(x$PPE$optim_args$models[x$PPE$optim_args$models==1])," model: ",names(x$PPE$model)," = ",format(x$PPE$model),sep=""))
    cat("\n\n")
  } else cat("\nBrownian motion model\n")
  cat("  Log-", ifelse(dd$REML == 1, "restricted-", 
                       ""), "likelihood: ", format(x$logLik), "\n", sep = "")
  cat("\nCoefficients:\n")
  print(coef(x))
  cat("\n")
  
  cat("Degrees of freedom:", dd[["N"]] - dd[["p"]]-dd[["df.int"]], "\n")
  cat("Residual standard error:", format(x$sigma), "\n")
  invisible(x)
}

logLik.phylopars <- function(object,...)
{
  val <- -object$minus2ll/2
  attr(val,"nall") <- object$options["nob"]
  attr(val,"nobs") <- object$options["nob"]
  attr(val,"df") <- length(object$theta)
  class(val) <- "logLik"
  val
}

logLik.phylopars.lm <- function(object,...)
{
  val <- object$logLik
  attr(val,"nall") <- object$dims$N
  attr(val,"nobs") <- object$dims$N
  attr(val,"df") <- object$df1 + 1
  class(val) <- "logLik"
  val
}

summary.phylopars.lm <- function(object,...)
{
  tTable <- data.frame(coef(object), object$SEs, object$ts, object$ps)
  dimnames(tTable) <- list(names(coef(object)), c("Value", "Std.Error", 
                                                  "t-value", "p-value"))
  dd <- object$dims
  mCall <- object$call
  cat("Generalized least squares fit by ")
  cat(ifelse(dd$REML == 1, "REML\n", "maximum likelihood\n"))
  cat("  Model:", deparse(object$model), "\n")
  if(length(object$PPE$model)>0)
  {
    cat(paste("\n",names(object$PPE$optim_args$models[object$PPE$optim_args$models==1])," model: ",names(object$PPE$model)," = ",format(object$PPE$model),sep=""))
    cat("\n")
  } else cat("\nBrownian motion model\n")
  aux <- logLik(object)
  object$BIC <- BIC(aux)
  object$AIC <- AIC(aux)
  
  print(data.frame(AIC = object$AIC, BIC = object$BIC, logLik = as.vector(object$logLik), 
                   row.names = " "))
  cat("\nCoefficients:\n")
  colnames(tTable) <- c("Estimate", "Std.Err", "Z value", "Pr(>z)")
  printCoefmat(tTable,P.values = TRUE,has.Pvalue = TRUE,na.print = "")
  
  cat("\n")
  cat("Residual standard error:", format(object$sigma), "\n")
  cat("Degrees of freedom:", dd[["N"]] - 
        dd[["p"]] - dd[["df.int"]], "\n")
  cat("Multiple R-squared:", format(object$R2), "Adjusted R-squared:", format(object$R2adj), "\n")
  cat("F-statistic:", format(object$Fstat), "on", format(object$df1-1), "and", format(object$df2), "DF, p-value:", format(object$pval), "\n")
  invisible(object)
}