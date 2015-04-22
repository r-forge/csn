\name{rcsn}
\alias{rcsn}
\title{
Random number generation
}
\description{
  Random number generation of the closed-skew normal distribution
}
\usage{
rcsn(k, mu = rep(0, n), sigma, gamma, nu = rep(0, q), delta)
}
\arguments{
  \item{k}{the number of random  numbers to be generated
}
  \item{mu}{a numeric vector representing the location parameter of the
        distribution; it must be of length \code{n}, as defined above
}
  \item{sigma}{a positive definite matrix representing the scale parameter 
        of the distribution; a vector of length 1 is also allowed
}
  \item{gamma}{a matrix representing the skewness parameter
        of the distribution; a vector of length 1 is also allowed 
}
  \item{nu}{a numeric vector allows for closure with conditional densities;
        it must be of length \code{q}, as defined above
}
  \item{delta}{a positive definite matrix allows for closure with the
        marginal densities; a vector of length 1 is also allowed 
}
}
\details{
Function rcsn makes use of rmvnorm from package mvtnorm;
}
\value{
\code{rcsn} returns a matrix of \code{k} rows of random vectors
}
\seealso{
\code{\link[mvtnorm]{rmvnorm}}
}
\examples{
mu <- c(1,2,3)
sigma <- matrix(c(2,-1,0,-1,2,-1,0,-1,2),3)
gamma <- matrix(c(0,1,0,2,2,3),2,3)
nu <- c(1,3)
delta <- matrix(c(1,1,1,2),2)
x1 <- rcsn(100, mu, sigma, gamma, nu, delta)
x2 <- rcsn(100,5,9,1,0,0.05)
}

