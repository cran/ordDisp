plotordDispsmooth <- function(x, names_s, cex, title){
  
  par(mfrow=c(length(names_s),2))
  for(i in 1:length(names_s)){
    names_X <- names(x@X)[grepl(names_s[i], names(x@X))]
    names_Z <- names(x@Z)[grepl(names_s[i], names(x@Z))]
    
    newnames <- substr(names(coef(x)), start = 1, stop = nchar(names(coef(x))) - 1)
    
    fgamma <- fbeta <- numeric(100)
    KIgamma <- KIbeta <- matrix(0, nrow=100, ncol=2)
    covmat <- summary(x)@cov.unscaled 
    
    if(length(names_X)==length(names_Z)){
      des <- bs(seq(min(x@X0[,names_s[i]]),max(x@X0[,names_s[i]]),length=100), df=length(names_X))
      des <- apply(des, 2, function(x){x-max(x)/2})
      beta <- gamma <- numeric(length(names_X))
      for(j in 1:length(names_X)){
        h1 <- which(newnames==names_X[j])
        gamma[j] <- coef(x)[h1[2]]
        beta[j]  <- coef(x)[h1[1]]
      }
      fgamma <- des%*%gamma
      fbeta  <- des%*%beta
      
      covbeta <- covmat[which(newnames%in%names_X)[1:length(names_X)], 
                        which(newnames%in%names_X)[1:length(names_X)]]
      varbeta <- des%*%covbeta%*%t(des)
      KIbeta[,1] <- fbeta-1.96*sqrt(diag(varbeta))
      KIbeta[,2] <- fbeta+1.96*sqrt(diag(varbeta))
      
      covgamma <- covmat[which(newnames%in%names_Z)[(length(names_X)+1):(length(names_X)*2)], 
                         which(newnames%in%names_Z)[(length(names_X)+1):(length(names_X)*2)]]
      vargamma <- des%*%covgamma%*%t(des)
      KIgamma[,1] <- fgamma-1.96*sqrt(diag(vargamma))
      KIgamma[,2] <- fgamma+1.96*sqrt(diag(vargamma))
      
      
    } else{
      if(length(names_X)<=1){
        if(length(names_X)==0){
          fbeta <- rep(0,100) 
        } else{
          h1 <- which(newnames==names_X)
          des<- seq(min(x@X0[,names_s[i]]),max(x@X0[,names_s[i]]),length=100)
          beta  <- coef(x)[h1]
          fbeta <- des*beta
          covbeta <- covmat[h1, h1]
          varbeta <- des*covbeta*des
          KIbeta[,1] <- fbeta-1.96*sqrt(varbeta)
          KIbeta[,2] <- fbeta+1.96*sqrt(varbeta)
        }
        des <- bs(seq(min(x@Z0[,names_s[i]]),max(x@Z0[,names_s[i]]),length=100), df=length(names_Z))
        des <- apply(des, 2, function(x){x-max(x)/2})
        gamma <- numeric(length(names_Z))
        for(j in 1:length(names_Z)){
          h1 <- which(newnames==names_Z[j])
          gamma[j] <- coef(x)[h1]
        }
        fgamma   <- des%*%gamma
        covgamma <- covmat[which(newnames%in%names_Z), which(newnames%in%names_Z)]
        vargamma <- des%*%covgamma%*%t(des)
        KIgamma[,1] <- fgamma-1.96*sqrt(diag(vargamma))
        KIgamma[,2] <- fgamma+1.96*sqrt(diag(vargamma))
      } else{
        if(length(names_Z)==0){
          fgamma <- rep(0,100)
        } else{
          h1 <- which(newnames==names_Z)
          des<- seq(min(x@Z0[,names_s[i]]),max(x@Z0[,names_s[i]]),length=100)
          gamma  <- coef(x)[h1]
          fgamma <- des*gamma
          covgamma <- covmat[h1, h1]
          vargamma <- des*covgamma*des
          KIgamma[,1] <- fgamma-1.96*sqrt(vargamma)
          KIgamma[,2] <- fgamma+1.96*sqrt(vargamma)
        }
        des  <- bs(seq(min(x@X0[,names_s[i]]),max(x@X0[,names_s[i]]),length=100), df=length(names_X))
        des  <- apply(des, 2, function(x){x-max(x)/2})
        beta <- numeric(length(names_X))
        for(j in 1:length(names_X)){
          h1 <- which(newnames==names_X[j])
          beta[j] <- coef(x)[h1]
        }
        fbeta   <- des%*%beta
        covbeta <- covmat[which(newnames%in%names_X), which(newnames%in%names_X)]
        varbeta <- des%*%covbeta%*%t(des)
        KIbeta[,1] <- fbeta-1.96*sqrt(diag(varbeta))
        KIbeta[,2] <- fbeta+1.96*sqrt(diag(varbeta))
      }
    }
      
    ymax <- max(c(KIgamma,KIbeta))
    ymin <- min(c(KIgamma,KIbeta))
      
    plot(seq(min(x@X0[,names_s[i]]),max(x@X0[,names_s[i]]),length=100), fbeta, type="l", lwd=cex, 
        xlab=names_s[i], ylab=bquote("f("*beta*")"), cex.axis=cex, cex.lab=cex, ylim=c(ymin, ymax))
    lines(seq(min(x@X0[,names_s[i]]),max(x@X0[,names_s[i]]),length=100), KIbeta[,1], type="l", lwd=cex, lty="dashed")
    lines(seq(min(x@X0[,names_s[i]]),max(x@X0[,names_s[i]]),length=100), KIbeta[,2], type="l", lwd=cex, lty="dashed")
    plot(seq(min(x@X0[,names_s[i]]),max(x@X0[,names_s[i]]),length=100), fgamma, type="l", lwd=cex, 
        xlab=names_s[i], ylab=bquote("f("*alpha*")"), cex.axis=cex, cex.lab=cex, ylim=c(ymin, ymax))
    lines(seq(min(x@X0[,names_s[i]]),max(x@X0[,names_s[i]]),length=100), KIgamma[,1], type="l", lwd=cex, lty="dashed")
    lines(seq(min(x@X0[,names_s[i]]),max(x@X0[,names_s[i]]),length=100), KIgamma[,2], type="l", lwd=cex, lty="dashed")
      
    if(!is.null(title)){
      mtext(title,side=3,line=2,cex=cex)
    }
  }
  par(mfrow=c(1,1))
}