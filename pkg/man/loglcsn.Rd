\name{loglcsn}
\alias{loglcsn}
\title{
The log-likelihood function
}
\description{
  The log-likelihood function of the closed-skew normal distribution
  }
\usage{
loglcsn(x, mu, sigma, gamma, nu, delta)
}
\arguments{
  \item{x}{this is either a vector of length \code{n} or 
           a matrix with \code{n} columns, where \code{n=ncol(sigma)},
           giving the coordinates of the point(s) where the density must
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
Function loglcsn makes use of pmvnorm and dmvnorm from package mvtnorm
}
\value{
  \code{loglcsn} returns a sum of log-transformed density values
}
\seealso{
\code{\link[mvtnorm]{pmvnorm}}, \code{\link[mvtnorm]{dmvnorm}}
}
\examples{
x <- cbind(seq(3,9,length=100),seq(7,13,length=100))
mu <- c(5,7)
sigma <- matrix(c(1,0.2,0.2,4),2)
gamma <- matrix(c(4,0,0,5),2)
nu <- c(-2,6)
delta <- matrix(c(1,0,0,1),2)
L <- loglcsn(x, mu, sigma, gamma, nu, delta)
}

