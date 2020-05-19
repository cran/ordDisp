#' Visualization of Estimated Effects 
#'
#'@description 
#'A function to visualize the estimated effects of the location-shift model or rating-scale model accounting for response styles
#'(RSRS) obtained by \code{ordDisp}. In case of linear effects, the function returns a two-dimensional plot of the tupel \eqn{(exp{\alpha},exp{\beta}}. 
#'It is optional to include pointwise 95\% confidence intervals represented by stars, where the horizontal and vertical 
#'length correspond to the confidence intervals of \eqn{exp{\alpha}} (dispersion or response-style effect) and \eqn{exp{\beta}} 
#'(location or content-related effect). In case of smooth effects, the function returns two plots of the fitted (non-linear) functions \eqn{f(\beta)} 
#'and \eqn{f(\alpha)}. 
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
#'@author Moritz Berger <moritz.berger@imbie.uni-bonn.de> \cr \url{https://www.imbie.uni-bonn.de/personen/dr-moritz-berger/}
#'
#'@references
#'Tutz, Gerhard and Berger, Moritz (2016): Response Styles in Rating Scales - Simultaneous Modelling of 
#'Content-Related Effects and the Tendency to Middle or Extreme Categories, 
#'Journal of Educational and Behavioral Statistics 41(3), 239-268.
#' 
#'Tutz, Gerhard and Berger, Moritz (2017): Seperating Location and Dispersion in Ordinal Regression Models,
#'Econometrics and Statistics 2, 131-148.
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
#'mod2 <- ordDisp(RET~SM+s(DIAB)+GH+BP|SM+DIAB+GH+BP, data=reti, 
#'family="cumulative", n_bs=4, scaling=FALSE)
#'plot(mod2, names=c("DIAB"))
#'
#'@docType methods
#'@export
#'@rdname plotordDisp
#'@importFrom grDevices grey 
#'@importFrom graphics abline axis box mtext par plot.new plot.window points polygon text lines 
#'@importFrom splines bs
#'

plotordDisp <- function(x,names,colorvec,reference=NULL,labels=NULL,cex=1,KI=FALSE,KIfactor=10/11,title=NULL,...){
  
  # check input 
  if(missing(names)){
    stop("Names of variables to be plotted are missing.")
  }

  terms <- paste(x@form)[3]
  terms <- gsub(" ","", terms)
  terms <- unique(strsplit(terms, "\\+|\\~|\\|")[[1]])
  # smooth 
  h2 <- "(s[(])(.+)([)])"
  vars <- sub(h2, "\\2", terms)
  is_smooth <- grepl(h2, terms)
  vars_s <- sub(h2,"\\2", terms[is_smooth])
  if(!all(names%in%vars)){
    wrong <- paste(names[which(!(names%in%vars))],collapse="and")
    stop(paste("Covariate",wrong,"undefined!"),call.=FALSE)
  }
  
  names_s <- names[names%in%vars_s]
  if(length(names_s)>0){
    for(v in names_s){
      plotordDispsmooth(x, names_s, cex, title)
    } 
    names_plot <- names[!names%in%vars_s]
  } else{
    names_plot <- names 
  }
  
  if(length(names_plot)>0){
    
    names_plot <- names(x@X)[grepl(paste0(names_plot, collapse="|"), names(x@X))]
    names_Z <- names(x@Z)[grepl(paste0(names_plot, collapse="|"), names(x@Z))]
    
    if(missing(colorvec)){
      colorvec <- rep(1,length(names_plot)+length(reference))
    }
  
    newnames <- substr(names(coef(x)), start = 1, stop = nchar(names(coef(x))) - 1)
    
    beta <- gamma <- numeric(length(names_plot))
    
    for(i in 1:length(names_plot)){
      h1 <- which(newnames==names_plot[i])
      if(length(h1)==1){
        if(!names_plot[i]%in%names_Z){
          gamma[i] <- 1
          beta[i]  <- exp(coef(x)[h1[1]])
        } else{
          beta[i]  <- 1
          gamma[i] <- exp(coef(x)[h1[1]])
        }
      } else{
        gamma[i] <- exp(coef(x)[h1[2]])
        beta[i]  <- exp(coef(x)[h1[1]])
      }
    }
    
    if(KI){
      betaKI <- gammaKI <- matrix(NA,nrow=length(names_plot),ncol=2)
      stats  <- summaryvglm(x)@coef3
      
      for(i in 1:length(names_plot)){
        h1 <- which(newnames==names_plot[i])
        if(length(h1)==1){
          if(!names_plot[i]%in%names_Z){
            gammaKI[i,1] <- gammaKI[i,2] <- 1
            betaKI[i,1] <- exp(stats[h1[1],1]-1.96*stats[h1[1],2])
            betaKI[i,2] <- exp(stats[h1[1],1]+1.96*stats[h1[1],2])
          } else{
            betaKI[i,1] <- betaKI[i,2] <- 1
            gammaKI[i,1] <- exp(stats[h1[1],1]-1.96*stats[h1[1],2])
            gammaKI[i,2] <- exp(stats[h1[1],1]+1.96*stats[h1[1],2])
          }
        } else{
          gammaKI[i,1] <- exp(stats[h1[2],1]-1.96*stats[h1[2],2])
          gammaKI[i,2] <- exp(stats[h1[2],1]+1.96*stats[h1[2],2])
          
          betaKI[i,1] <- exp(stats[h1[1],1]-1.96*stats[h1[1],2])
          betaKI[i,2] <- exp(stats[h1[1],1]+1.96*stats[h1[1],2])
        }
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
      for(i in 1:length(names_plot)){
        
        x <- c(gammaKI[i,1],gammaKI[i,1]+(gamma[i]-gammaKI[i,1])*(KIfactor),gamma[i],gamma[i]+(gamma[i]-gammaKI[i,1])*(1-KIfactor),
               gammaKI[i,2],gamma[i]+(gamma[i]-gammaKI[i,1])*(1-KIfactor),gamma[i],gammaKI[i,1]+(gamma[i]-gammaKI[i,1])*(KIfactor),
               gammaKI[i,1])
        y <- c(beta[i],betaKI[i,1]+(beta[i]-betaKI[i,1])*(KIfactor),betaKI[i,1],betaKI[i,1]+(beta[i]-betaKI[i,1])*(KIfactor),beta[i],
               beta[i]+(beta[i]-betaKI[i,1])*(1-KIfactor),betaKI[i,2],beta[i]+(beta[i]-betaKI[i,1])*(1-KIfactor),beta[i])
        polygon(x,y,col=grey(0.9))
      }
    }
    
    for(i in 1:length(names_plot)){
      points(gamma[i],beta[i],pch=19,cex=cex,col=colorvec[i])
      if(gamma[i]<1){
        if(!is.null(labels)){
          text(gamma[i]+0.01,beta[i],labels[i],adj=c(0,-0.2),cex=cex,col=colorvec[i])
        } else{
          text(gamma[i]+0.01,beta[i],names_plot[i],adj=c(0,-0.2),cex=cex,col=colorvec[i])
        }
      } else{
        if(!is.null(labels)){
          text(gamma[i]-0.01,beta[i],labels[i],adj=c(1,-0.2),cex=cex,col=colorvec[i])
        } else{
          text(gamma[i]-0.01,beta[i],names_plot[i],adj=c(1,-0.2),cex=cex,col=colorvec[i])
        }
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
}

setMethod("plot",signature(x="ordDisp"),plotordDisp)







