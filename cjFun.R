#' Clear working session
#' 
#' Function to transform P-values into symbols of significance (***)
#' @param s A vector of P-values
#' @return TRUE for clear
#' 
#' @author C.J.
#' @export
clearSession <- function(){
  lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
}


#' Write table
#' 
#' Function to write.table with regular parameter
#' @param data Data.frame to be written
#' @param path Directory to write
#' @return NULL
#' 
#' @author C.J.
#' @export
write.tsv<-function(x = "", file = ""){
  write.table(x, file = file, append = F, quote = F, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T, qmethod = 'escape')
}



