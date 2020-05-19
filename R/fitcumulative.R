fitcumulative <- function(y,X,Z,scaling,middle,m,reverse,...){ 

  n  <- nrow(X) # observations
  k  <- length(unique(y)) # categories 
  px <- ncol(X) # covariates in X
  
  if(is.null(Z)){
    
    # generate formula 
    colnames(X)    <- namesX    <- paste0(names(X),"x")
    f01 <- paste(namesX,collapse="+")
    formula0 <- formula(paste("y~",f01))
    
    # estimation with vglm 
    mod <- vglm(formula0,
                family=cumulative(parallel=FALSE~1, reverse=reverse),
                data=X,
                control=vglm.control(checkwz=FALSE), ...)
    
    return(list("mod"=mod,"m"=m))
    
    
  } else{
    
    
    pz <- ncol(Z) # covariates in Z  
    if(reverse){
      fac <- -1 
    } else{
      fac <- 1 
    }
    
    # generate design-matrix  
    Zbig <- Z[,rep(1:pz,each=k-1)]
    
    if(k%%2!=0){ # odd
      if(middle){
        m <- floor(k/2)+1
      }
      for(i in 0:(pz-1)){
        if(!scaling){
          if(m>1){
            Zbig[,(1:(m-1))+i*(k-1)] <- fac*-Zbig[,(1:(m-1))+i*(k-1)]
            Zbig[,(m:(k-1))+i*(k-1)] <- fac*Zbig[,(m:(k-1))+i*(k-1)]
          }
        } else{
          if(m>1){
            for(r in 1:(m-1)){
              Zbig[,r+i*(k-1)] <- fac*-((m-r-1)+0.5)*Zbig[,r+i*(k-1)]
            }
          }
          for(r in m:(k-1)){
            Zbig[,r+i*(k-1)] <- fac*+((r-m)+0.5)*Zbig[,r+i*(k-1)]
          }
        }
      }
    }
    
    if(k%%2==0){ # even
      if(middle){
        m <- k/2
      }
      for(i in 0:(pz-1)){
        if(!scaling){
          if(m>1){
            Zbig[,(1:(m-1))+i*(k-1)]     <- fac*-Zbig[,(1:(m-1))+i*(k-1)]
            Zbig[,((m+1):(k-1))+i*(k-1)] <- fac*Zbig[,((m+1):(k-1))+i*(k-1)]
          }
          Zbig[,m+i*(k-1)] <- 0
        } else{
          for(r in 1:m){
            Zbig[,r+i*(k-1)] <- fac*-(m-r)*Zbig[,r+i*(k-1)]
          }
          if(m < (k-1)){
            for(r in (m+1):(k-1)){
              Zbig[,r+i*(k-1)] <- fac*+(r-m)*Zbig[,r+i*(k-1)]
            }
          }
        }
      }
    }
    
    # labels 
    h1 <- paste0(rep(names(Z),each=k-1),"z")
    h2 <- rep(1:(k-1),times=pz)
    colnames(Zbig) <- namesZbig <- paste0(h1,h2)
    colnames(Z)    <- namesZ    <- paste0(names(Z),"z")
    colnames(X)    <- namesX    <- paste0(names(X),"x")
    
    DM <- data.frame(X,Zbig,Z)
    
    # generate formulas 
    f11 <- paste(namesX,collapse="+")
    f12 <- paste(namesZ,collapse="+")
    formula1 <- formula(paste("y~",f11,"+",f12))
    
    f21 <- paste(namesZbig,collapse="+")
    formula2 <- formula(paste("~",f11,"+",f12,"+",f21))
    
    formula3 <- c()
    for(i in 0:(pz-1)){
      f31 <- paste(namesZbig[(1:(k-1))+(k-1)*i],collapse="+")
      f32 <- formula(paste(namesZ[i+1],"~",f31))
      formula3 <- c(formula3,f32)
    }
    
    # estimation with vglm 
    mod <- vglm(formula1,
                family=cumulative(parallel=FALSE~1, reverse=reverse),
                xij=formula3,
                form2=formula2,
                data=DM,
                checkwz=FALSE, ...)
    
    return(list("mod"=mod,"m"=m))
    
  }
  
}