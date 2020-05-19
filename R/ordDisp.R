#' Separating Location and Dispersion in Ordinal Regression Models 
#'
#'@description
#'A function to estimate the location-shift model or rating-scale model accounting for response styles (RSRS) 
#'for the regression analysis of ordinal responses. The model allows to account for differing variability in 
#'subgroups of the population. The model explicitely links varying disperion (or response behaviour) to explanatory 
#'variables (metric, binary, ordinal and/or nominal). The basic models are described in Tutz and Berger (2016) and Tutz and Berger (2017). 
#'
#'@param formula Object of class \code{\link{formula}}: a symbolic description of the model to be fitted. See details.
#'@param data Data.frame of class \code{\link{data.frame}} containing the variables of the model.
#'@param family Type of link function that is used to link the mean responses to the linear predictors of the model; 
#'ordDisp currently allows only one out of \code{"cumulative"} and \code{"acat"}. See details.
#'@param scaling If true, the thresholds of the location-shift model are shifting by using scale values for the widening of 
#'the intervals between two thresholds. 
#'@param middle If true, the model expects a symmetric response of the form 'strongly disagree','moderatly disagree',...,
#''moderatly agree','strongly agree'. 
#'@param m Middle category of the (non-symmetric) response, chosen for the model. Only relevant, if \code{middle=FALSE}.
#'@param n_bs Number of inner B-spline basis functions for smooth components (see details).  
#'@param reverse Argument of the family function passed to \code{\link{vglm}}. 
#'@param ... Further arguments passed to or from other methods
#'
#'@details 
#'The \link{formula} has to have the form \code{response ~ x-variables|z-variables}, where \code{response} is the name of 
#'the ordinal response variable, \code{x-variables} are the terms that specify the location (or content-related) effects 
#'of the model and \code{z-variables} are the terms that specify the dispersion (or response-style) effects.
#'
#'If all the variables are entered in both parts of the model, the right hand side of the formula can, for example, 
#'have the form \code{x1+...+xp|x1+...+xp}. If the second part is omitted, a simple model without dispersion 
#'(or response-style) effects is fitted.
#'
#'The function allows for smooth (non-linear) effects in the x-variables and/or the z-variables. Smooth effects are specified 
#'by entering s(x) and/or s(z) into the formula. The functions are fitted using \code{n_bs} B-spline basis functions.
#'
#'Function \code{ordDisp} internally calls \code{\link{vglm}} from package \code{\link{VGAM}}. Argument \code{family} is passed to \code{vglm}. 
#'Currently two link functions are implemented
#'
#'\itemize{
#'\item \code{"cumulative"} to estimate a cumulative model of the form 
#'\deqn{P(y\leq r)/P(y>r)=eta_r}
#'\item \code{"acat"} to estimate a adjacent-categories model of the form 
#'\deqn{P(y=r+1)/P(y=r)=eta_r}
#'}
#'
#'
#'@return 
#'Object of class \code{ordDisp} which inherits from \code{\link{vglm}}. The object comprises all the slots of an 
#'\code{"vglm"}-object and in addition the following components: 
#'
#'\item{outercall}{The matched call of \code{ordDisp}.}
#'\item{X}{Design matrix of x-variables.}
#'\item{Z}{Design matrix of z-variables.}
#'
#'All the methods implemented for objects of class \code{vglm}, like \code{print}, \code{summary}, \code{predict} 
#'and \code{plot} can be applied. 
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
#'@seealso \code{\link[VGAM]{summaryvglm}}, \code{\link[VGAM]{predictvglm}}, \code{\link[ordDisp]{plotordDisp}}
#'
#'@examples 
#'data(reti)
#'
#'mod <- ordDisp(RET~SM+DIAB+GH+BP|SM+DIAB,data=reti,family="cumulative")
#'summary(mod)
#'
#'mod2 <- ordDisp(RET~SM+s(DIAB)+GH+BP|SM+DIAB+GH+BP, data=reti, 
#'family="cumulative", n_bs=4, scaling=FALSE)
#'summary(mod2)
#'
#'@import VGAM
#'@importFrom methods new 
#'@importFrom splines bs 
#'@export
#'  

ordDisp <- function(formula,
                    data, 
                    family=c("cumulative","acat"),
                    scaling=TRUE,
                    middle=TRUE,
                    m=NULL,
                    n_bs=6, 
                    reverse=FALSE,
                    ...){
  
  # check input 
  if(missing(formula)){
    stop("Argument formula is missing with no default.")
  }
  if(missing(data)){
    stop("Argument data is missing with no default.")
  }
  if(!middle && is.null(m)){
    stop("If middle is FALSE, m has to be specified.")
  }
  family <- match.arg(family)
  
  #predefinition
  comp <- specification(formula,data,n_bs)  
  y    <- comp$y
  X    <- comp$X
  Z    <- comp$Z 
  
  # estimation 
  if(family=="cumulative"){
    output <- fitcumulative(y,X,Z,scaling,middle,m,reverse,...)
  } 
  if(family=="acat"){
    output <- fitadjacent(y,X,Z,scaling,middle,m,reverse,...)
  }
  outercall <- match.call()
  form      <- formula

  # define class 
  if(is.null(Z)){
    to_return  <- new("ordDisp",output$mod,outercall=outercall,form=form,X=X,X0=comp$X0)
  } else{
    to_return  <- new("ordDisp",output$mod,outercall=outercall,form=form,X=X,Z=Z,X0=comp$X0,Z0=comp$Z0,m=output$m)
  }
  
  return(to_return)
  
}


