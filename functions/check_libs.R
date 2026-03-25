check_libs <- function(libs){
  
  # check if libraries already installed
  check <- unlist(lapply(libs, require, character.only=TRUE))
  
  # install missing libraries
  if(sum(!check)){ 
    install.packages(libs[!check])
    tmp<-lapply(libs[!check], require, character.only=TRUE)
  }
}