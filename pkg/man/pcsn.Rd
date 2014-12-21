\name{pcsn}
\alias{pcsn}
\title{
The cumulative distribution function
}
\description{
The cumulative distribution function of the closed-skew normal distribution
}
\usage{
pcsn(x, mu, sigma, gamma, nu, delta)
}
\arguments{
  \item{x}{this is either a vector of length \code{n} or 
           a matrix with \code{n} columns, where \code{n=ncol(sigma)},
           giving the coordinates of the point(s) where the cdf must
           be evaluated
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
Function pcsn makes use of pmnorm from package mnormt
}
\value{
\code{pcsn} returns a vector of cdf values
}
\seealso{
\code{\link[mnormt]{pmnorm}}
}
\examples{
x1 <- seq(4,6,by = 0.1)
x2 <- x1+sin(x1)
x3 <- x1-cos(x1)
x <- cbind(x1,x2,x3)
mu <- c(1,2,3)
sigma <- matrix(c(2,-1,0,-1,2,-1,0,-1,2),3)
gamma <- matrix(c(0,1,0,2,2,3),2,3)
nu <- c(1,3)
delta <- matrix(c(1,0,1,2),2)
pcsn(6,5,9,1,0,0.05)
pcsn(c(3,4,5),mu,sigma,gamma,nu,delta)
pcsn(x,mu,sigma,gamma,nu,delta)
}
