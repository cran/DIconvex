\name{DIconvex}
\alias{DIconvex}
\title{Finding patterns of monotonicity and convexity in two-dimensional data}
\description{
This package takes as input x values \eqn{x_1,\ldots,x_n}, as well as lower \eqn{L_1,\ldots,L_n}, and upper bounds \eqn{U_1,\ldots,U_n}. It maximizes \eqn{\sum _{i=1}^{n}f_i, \, f_i\in \{0,1\}} such that there exists at least one convex increasing (decreasing) set of values \eqn{L_j\le y_j\le U_j, j\in C}, where \eqn{C} is the set of indices \eqn{i=1,\ldots,n} for which \eqn{f_i=1}.
}
\usage{
DIconvex(x, lower, upper, increasing = FALSE, epsim = 0, epsic = 0,visual=TRUE)}

\arguments{
  \item{x}{a numeric vector containing a set of points. The elements of \code{x} have to be positive and ranked in ascending order. The vector \code{x} can not contain duplicate data.  }
  \item{lower}{a numeric vector of the same length as \code{x} containing the lower limit points. The elements of the vector \code{lower} have to be non-negative and finite.}
  \item{upper}{a numeric vector of the same length as \code{x} containing the upper limit points. The elements of the vector \code{upper} have to be non-negative and finite. Furthermore, \eqn{L_i\le U_i, i=1,\ldots,n.}}
  \item{increasing}{a boolean value determining whether to look for an increasing or decreasing pattern. The default value is FALSE.  }
  \item{epsim}{a non-negative value controlling the monotonicity conditions, \eqn{y_{i+1}-y_{i}\le (\ge)epsim, \, i=1,\ldots,n-1}. The default value is 0.}
  \item{epsic}{a positive value controlling the convexity condition. For \eqn{\alpha_i:=(x_i-x_{i+1})/(x_{i-1}-x_{i+1})} the condition imposed is \eqn{y_i- \alpha _i y_{i+1}-(1-\alpha_i)y_{i-1}\le epsic, \, i=2,\ldots,n-1}. The default value is 0.}
 \item{visual}{a boolean value indicating whether a visual representation of the solution is desired. Here a solution is depicted for all values of x, with linearly interpolated y if \eqn{i \notin C}. The default value is TRUE. }
}
\details{
The package \code{DIconvex} is solved as a linear program facilitating \code{\link[lpSolveAPI:solve.lpExtPtr]{lpSolveAPI}}.
It lends itself to applications with financial options data. Given a dataset of call or put options, the function maximizes the number of data points such that there exists at least one set of arbitrage-free fundamental option prices within bid and ask spreads.

For this particular application, \code{x} is the vector of strike prices, \code{lower} represents the vector of bid prices and \code{upper} represents the vector of ask prices.
}

\value{
a list containing:

a vector containing \eqn{f_1,\ldots,f_n}.

a vector containing \eqn{y_j, \, j \in C}.

a single integer value containing the status code of the underlying linear program. For the interpretation of status codes please see \code{\link[lpSolveAPI:solve.lpExtPtr]{lpSolveAPI}} R documentation. The value 0 signifies success.
}

\author{

{Liudmila Karagyaur <liudmila.karagyaur@usi.ch>}{}

 {Paul Schneider <paul.schneider@usi.ch>}{}


}



\examples{
x = c(315, 320, 325, 330, 335, 340, 345, 350)
upper = c(0.5029714, 0.5633280, 0.6840411, 0.8751702, 3.0000000, 1.5692708, 2.3237279, 3.5207998)
lower = c(0.2514857, 0.4325554, 0.4325554, 0.6236845, 2.5000000, 1.1870125, 1.9414696, 3.1385415)

DIconvex(x, lower, upper, increasing = TRUE)

x = c(340, 345, 350, 355, 360, 365)
lower = c(2.7661994, 1.3177168, 1.5029454, 0.1207069, 0.1207069, 0.1207069)
upper = c(3.1383790, 1.5088361, 1.6236522, 0.3721796, 0.1810603, 0.2514727)

DIconvex(x, lower, upper, increasing = FALSE)
}
