

trenchDet<- function(c){
  
 if(!is.numeric(c)){
   stop("Please provide a numeric vector")
 }else if(length(c) < 2){
   stop("Please provide a vector of at least length 2")
 }else{
   logdet <- trenchDetcpp(c)$logdet
 }
  
 return(logdet)
  
}