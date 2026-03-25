w_eCDF <- function(x,scale=0) {
  #Compute empirical weighted cummulative distribution function
  x <- sort(x)
  n <- length(x)
  vals <- unique(x)
  
  # calculate weighted eCDF
  weight <- 2/(1+exp(scale*vals))
  w_pdf <- tabulate(match(x, vals)) * weight
  res <- cumsum(w_pdf)/sum(w_pdf)
  # plot(vals,res)
  
  # create function
  res2 <- approxfun(vals, res, method="constant", yleft=0, yright=1, f=0, ties = "ordered")
  class(res2) <- c("ecdf", "stepfun", class(res2))
  assign("nobs", n, envir = environment(res2))
  attr(res2, "call") <- sys.call()
  res2
}
