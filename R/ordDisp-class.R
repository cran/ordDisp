#'@rdname ordDisp
#'@export

setClass("ordDisp",representation(outercall="call",form="formula",X="data.frame",Z="data.frame",X0="data.frame",Z0="data.frame",m="numeric"),contains="vglm",
         prototype=prototype(new("vglm")))