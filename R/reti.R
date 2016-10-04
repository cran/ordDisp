#' @name reti
#' @title Example Retinopathy  
#' @description The data set contains information about persons with retinopathy. In the 6-year followup study on diabetes 
#' and retinopathy status the interesting question is how the retinopathy status is associated wie several risk factors. 
#' @docType data
#' @usage data(reti)
#' @format A data frame containing 613 observations on 5 variables:
#' 
#' \code{RET } retinopathy status (1:no retinopathy, 2:nonproliferative retinopathy, 3:advanced retinopathy or
#' blind) 
#' 
#' \code{SM } smoker (1:yes, 0:no)
#' 
#' \code{DIAB } diabetes duration in years
#' 
#' \code{GH } glycosylated hemoglobin measured in percent
#' 
#' \code{BP } diastolic blood pressure in mmHg
#' 
#' 
#' @references 
#' 
#' Bender and Grouven (1998): Using binary logistic regression models for ordinal data with nonproportional
#' odds, J. Clin. Epidemiol., 51, 809-816.
#' 
#' 
#' @examples 
#' data(reti)
#' table(reti$RET)
#'   
#'  
NULL

