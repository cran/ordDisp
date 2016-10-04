specification  <- function(formula,data){
  
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
  nodisp <- grepl(h1,name_xz)
  if(nodisp){
    warning("No dispersion effects. Proportional odds model will be fitted.")
    names_x <- unlist(strsplit(name_xz,"\\+"))
    if(any(!(names_x %in% names(data)))){
      wrong <- paste(names_x[which(!(names_x %in% names(data)))],collapse="and")
      stop(paste("Covariate",wrong,"undefined!"),call.=FALSE)
    }
    X <- data[,names_x]
    Z <- NULL 
  } else{
    xz <- unlist(strsplit(name_xz,"\\|"))
    names_x <- unlist(strsplit(xz[1],"\\+"))
    names_z <- unlist(strsplit(xz[2],"\\+"))
    if(any(!(names_x %in% names(data)))){
      wrong <- paste(names_x[which(!(names_x %in% names(data)))],collapse="and")
      stop(paste("Covariate",wrong,"undefined!"),call.=FALSE)
    }
    if(any(!(names_z %in% names(data)))){
      wrong <- paste(names_z[which(!(names_z %in% names(data)))],collapse="and")
      stop(paste("Covariate",wrong,"undefined!"),call.=FALSE)
    }
    X <- data[,names_x]
    Z <- data[,names_z]
  }
  return(list("y"=y,
              "X"=X,
              "Z"=Z))
}