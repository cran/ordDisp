specification  <- function(formula,data,n_bs){
  
  formula <- paste(formula)
  name_y <- formula[2]
  
  ### y ###  
  if(!(name_y %in% names(data))){
    stop("Response y undefined", call.=FALSE)
  }
  y <- data[,name_y]
  
  ### X und Z ### 
  name_xz <- gsub(" ","",formula[3])
  
  h0 <- "^([^\\|]+?)\\|?([^\\|]*)$"
  test_formula <- grepl(h0,name_xz)
  if(!test_formula){
    stop("Formula is not correctly specified!")
  }
  
  h1 <- "^([^\\|]+?)\\|?$"
  h2 <- "(s[(])(.+)([)])"
  nodisp <- grepl(h1,name_xz)
  if(nodisp){
    warning("No dispersion effects. Proportional odds model will be fitted.")
    names_x <- unlist(strsplit(name_xz,"\\+"))
    
    # smooth 
    ind2x <- grepl(h2,names_x)
    names_sx <- sub(h2,"\\2",names_x[ind2x])
    names_x[ind2x] <- names_sx
    
    if(any(!(names_x %in% names(data)))){
      wrong <- paste(names_x[which(!(names_x %in% names(data)))],collapse="and")
      stop(paste("Covariate",wrong,"undefined!"),call.=FALSE)
    }
    X <- X0 <- data[,names_x]
    Z <- Z0 <- NULL 
    
    # factors 
    for(j in names_x){
      if(is.factor(X[,j])){
        desj <- model.matrix(formula(paste0("~", j)), data=X)
        X <- as.data.frame(cbind(X[,-which(names(X)==j), drop=F], desj[,-1, drop=F]))
      }
    }
    
    # smooth
    for(j in names_x){
      if(ind2x[names_x==j]){
        desj  <- bs(X[,j], df=n_bs)
        desj <- apply(desj, 2, function(x){x-max(x)/2})
        colnames(desj) <- paste0(j, seq(1:n_bs))
        X <- as.data.frame(cbind(X[,-which(names(X)==j), drop=F], desj))
      }
    }
    
  } else{
    xz <- unlist(strsplit(name_xz,"\\|"))
    names_x <- unlist(strsplit(xz[1],"\\+"))
    names_z <- unlist(strsplit(xz[2],"\\+"))
    
    # smooth 
    ind2x <- grepl(h2,names_x)
    names_sx <- sub(h2,"\\2",names_x[ind2x])
    names_x[ind2x] <- names_sx
    ind2z <- grepl(h2,names_z)
    names_sz <- sub(h2,"\\2",names_z[ind2z])
    names_z[ind2z] <- names_sz
    
    if(any(!(names_x %in% names(data)))){
      wrong <- paste(names_x[which(!(names_x %in% names(data)))],collapse="and")
      stop(paste("Covariate",wrong,"undefined!"),call.=FALSE)
    }
    if(any(!(names_z %in% names(data)))){
      wrong <- paste(names_z[which(!(names_z %in% names(data)))],collapse="and")
      stop(paste("Covariate",wrong,"undefined!"),call.=FALSE)
    }
    X <- X0 <- data[,names_x, drop=F]
    Z <- Z0 <- data[,names_z, drop=F]
    
    # factors 
    for(j in names_x){
      if(is.factor(X[,j])){
        desj <- model.matrix(formula(paste0("~", j)), data=X)
        X <- as.data.frame(cbind(X[,-which(names(X)==j), drop=F], desj[,-1, drop=F]))
      }
    }
    for(j in names_z){
      if(is.factor(Z[,j])){
        desj <- model.matrix(formula(paste0("~", j)), data=Z)
        Z <- as.data.frame(cbind(Z[,-which(names(Z)==j), drop=F], desj[,-1, drop=F]))
      }
    }
    
    # smooth 
    for(j in names_x){
      if(ind2x[names_x==j]){
        desj  <- bs(X[,j], df=n_bs)
        desj <- apply(desj, 2, function(x){x-max(x)/2})
        colnames(desj) <- paste0(j, seq(1:n_bs))
        X <- as.data.frame(cbind(X[,-which(names(X)==j), drop=F], desj))
      }
    }
    for(j in names_z){
      if(ind2z[names_z==j]){
        desj  <- bs(Z[,j], df=n_bs)
        desj <- apply(desj, 2, function(x){x-max(x)/2})
        colnames(desj) <- paste0(j, seq(1:n_bs))
        Z <- as.data.frame(cbind(Z[,-which(names(Z)==j), drop=F], desj))
      }
    }
  }
  return(list("y"=y,
              "X"=X,
              "Z"=Z,
              "X0"=X0,
              "Z0"=Z0))
}
