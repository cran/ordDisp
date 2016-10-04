#'@rdname ordDisp
#'@export

setClass("ordDisp",representation(outercall="call",X="data.frame",Z="data.frame",m="numeric"),contains="vglm",
         prototype=prototype(new("vglm")))