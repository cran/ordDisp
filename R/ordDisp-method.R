#' Visualization of Estimated Effects 
#'
#'@description 
#'A function to visualize the estimated effects of the location-shift model or rating-scale model accounting for response styles
#'(RSRS) obtained by \code{ordDisp}. The function returns a two-dimensional plot of the tupel \eqn{(exp{\alpha},exp{\beta}}. 
#'It is optional to include pointwise 95\% confidence intervals represented by stars, where the horizontal and vertical 
#'length correspond to the confidence intervals of \eqn{exp{\alpha}} (dispersion or response-style effect) and \eqn{exp{\beta}} 
#'(location or content-related effect).  
#'
#'@param x Object of class \code{\link[ordDisp]{ordDisp}}
#'@param names Names of the variables that shall be plotted 
#'@param colorvec Vector of colors that are used for plotting (same length as names)
#'@param reference Optional name of reference with estimate \eqn{(\alpha,\beta)=(0,0)} (for categorical covariates)
#'@param labels Optional names that are used as labels in the plot (same length as names)
#'@param cex Global argument to set the size of all the labels in the plot 
#'@param KI If true, pointwise 95\% confidence intervals are included in the plot 
#'@param KIfactor Ratio that is used to plot the stars that represent confidence intervals (only if \code{KI=TRUE})
#'@param title Optional title that is added to the plot 
#'@param ... Further arguments passed to or from other methods
#'
#'@author Moritz Berger <moritz.berger@stat.uni-muenchen.de> \cr \url{http://www.statistik.lmu.de/~mberger/}
#'
#'@references
#'Tutz, Gerhard and Berger, Moritz (2016a): Response Styles in Rating Scales - Simultaneous Modelling of 
#'Content-Related Effects and the Tendency to Middle or Extreme Categories, 
#'Journal of Educational and Behavioral Statistics 41(3), 239-268.
#' 
#'Tutz, Gerhard and Berger, Moritz (2016b): Seperating Location and Dispersion in Ordinal Regression Models, 
#'LMU Muenchen, Institut fuer Statistik, Technical Report 190 
#'
#'@seealso \code{\link[ordDisp]{ordDisp}}
#'
#'@examples 
#'data(reti)
#'
#'mod <- ordDisp(RET~SM+DIAB+GH+BP|SM+DIAB,data=reti,family="cumulative")
#'plot(mod,names=c("SM","DIAB"),colorvec=c(1,2))
#'plotvglm(mod)
#'
#'@docType methods
#'@export
#'@rdname plotordDisp
#'@importFrom grDevices grey 
#'@importFrom graphics abline axis box mtext par plot.new plot.window points polygon text 
#'

plotordDisp <- function(x,names,colorvec,reference=NULL,labels=NULL,cex=2,KI=FALSE,KIfactor=10/11,title=NULL,...){
  
  # check input 
  if(missing(names)){
    stop("Names of variables to be plotted are missing.")
  }
  if(missing(colorvec)){
    colorvec <- rep(1,length(names)+length(reference))
  }
  if(!all(names%in%names(x@X))){
    wrong <- paste(names[which(!(names%in%names(x@X)))],collapse="and")
    stop(paste("Covariate",wrong,"undefined!"),call.=FALSE)
  }
  
  newnames <- substr(names(coef(x)), start = 1, stop = nchar(names(coef(x))) - 1)
  
  beta <- gamma <- numeric(length(names))
  
  for(i in 1:length(names)){
    h1 <- which(newnames==names[i])
    gamma[i] <- exp(coef(x)[h1[2]])
    beta[i]  <- exp(coef(x)[h1[1]])
  }
  
  if(KI){
    betaKI <- gammaKI <- matrix(NA,nrow=length(names),ncol=2)
    stats  <- summaryvglm(x)@coef3
    
    for(i in 1:length(names)){
      h1 <- which(newnames==names[i])
      gammaKI[i,1] <- exp(stats[h1[2],1]-1.96*stats[h1[2],2])
      gammaKI[i,2] <- exp(stats[h1[2],1]+1.96*stats[h1[2],2])
      
      betaKI[i,1] <- exp(stats[h1[1],1]-1.96*stats[h1[1],2])
      betaKI[i,2] <- exp(stats[h1[1],1]+1.96*stats[h1[1],2])
    }
  }
  
  
  plot.new()
  if(!is.null(title)){
    par(mar=c(6,6,5,2))
  } else{
    par(mar=c(6,6,2,2))
  }
  
  if(KI){
    plot.window(ylim=c(min(1,min(beta),min(betaKI)),max(max(beta),max(betaKI),1)),xlim=c(min(1,min(gamma),min(gammaKI)),max(max(gamma),max(gammaKI),1)))
  } else{
    plot.window(ylim=c(min(1,min(beta))-0.01,max(max(beta),1)+0.01),xlim=c(min(1,min(gamma))-0.15,max(max(gamma),1)+0.15))
  }
  box()
  axis(1,cex.axis=cex)
  axis(2,cex.axis=cex)
  abline(h=1,col=grey(0.8),lwd=2.5)
  abline(v=1,col=grey(0.8),lwd=2.5)
  
  if(KI){
    for(i in 1:length(names)){
      
      x <- c(gammaKI[i,1],gammaKI[i,1]+(gamma[i]-gammaKI[i,1])*(KIfactor),gamma[i],gamma[i]+(gamma[i]-gammaKI[i,1])*(1-KIfactor),
             gammaKI[i,2],gamma[i]+(gamma[i]-gammaKI[i,1])*(1-KIfactor),gamma[i],gammaKI[i,1]+(gamma[i]-gammaKI[i,1])*(KIfactor),
             gammaKI[i,1])
      y <- c(beta[i],betaKI[i,1]+(beta[i]-betaKI[i,1])*(KIfactor),betaKI[i,1],betaKI[i,1]+(beta[i]-betaKI[i,1])*(KIfactor),beta[i],
             beta[i]+(beta[i]-betaKI[i,1])*(1-KIfactor),betaKI[i,2],beta[i]+(beta[i]-betaKI[i,1])*(1-KIfactor),beta[i])
      polygon(x,y,col=grey(0.9))
    }
  }
  
  for(i in 1:length(names)){
    points(gamma[i],beta[i],pch=19,cex=cex,col=colorvec[i])
    if(!is.null(labels)){
      text(gamma[i]+0.01,beta[i],labels[i],adj=c(0,-0.2),cex=cex,col=colorvec[i])
    } else{
      text(gamma[i]+0.01,beta[i],names[i],adj=c(0,-0.2),cex=cex,col=colorvec[i])
    }
  }
  if(!is.null(title)){
    mtext(title,side=3,line=2,cex=cex)
  }
  mtext(expression(exp(alpha)),side=1,line=4,cex=cex)
  mtext(expression(exp(beta)),side=2,line=3,cex=cex)
  
  if(!is.null(reference)){
    points(1,1,pch=19,cex=cex,col=colorvec[length(colorvec)])
    text(1+0.01,1,reference,adj=c(0,0),cex=cex,col=colorvec[length(colorvec)])
  }
}

setMethod("plot",signature(x="ordDisp"),plotordDisp)



