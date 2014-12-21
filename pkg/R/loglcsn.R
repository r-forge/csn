loglcsn <- function(x, mu, sigma, gamma, nu, delta){
  if (is.matrix(x))
    result <- rep(0,nrow(x))
  if (length(mu)==1){
    result <- rep(0,length(x))
    x <-as.matrix(x)
  }
  if (!is.matrix(x)){
    x <-matrix(x,1,length(x))
    result <-rep(0,nrow(x))
  }
  if (length(mu) != ncol(x)) 
    stop("mismatch of dimensions 'x' and 'mu'")
  if (is.matrix(gamma) && nrow(gamma) != length(nu))
    stop("incorrect 'gamma'")
  if (is.matrix(gamma) && ncol(gamma) != length(mu))
    stop("incorrect 'gamma'")
  if (is.matrix(sigma) && nrow(sigma) != length(mu))
    stop("incorrect 'sigma'")
  if (is.matrix(sigma) && ncol(sigma) != length(mu))
    stop("incorrect 'sigma'")
  if (is.matrix(delta) && nrow(delta) != length(nu))
    stop("incorrect 'delta'")
  if (is.matrix(delta) && ncol(delta) != length(nu))
    stop("incorrect 'delta'")
  for (i in 1:nrow(x)){
    f1 <- log(pmnorm(0,nu,delta+gamma%*%sigma%*%t(gamma)))
    f2 <- log(pmnorm(as.vector(gamma%*%(x[i,]-mu)), nu, delta))
    f3 <- dmnorm(x[i,], mu, sigma,log=TRUE)
    result[i] <- -f1+f2+f3
  }
    result <- sum(result)
  return(result)
}