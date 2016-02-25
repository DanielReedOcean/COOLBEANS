\name{Exchange.sensitivity}
\alias{Exchange.sensitivity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bottom water oxygen sensitivity to vertical exchange between surface and bottom waters
%%  ~~function to do ... ~~
}
\description{
Calculates the rate of change of bottom water oxygen with respect to changes in vertical exchange between surface and bottom waters.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Exchange.sensitivity(BOD, Exchange, L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BOD}{
  Benthic oxygen demand (mol O2/[m2 y])
%%     ~~Describe \code{BOD} here~~
}
  \item{Exchange}{
  Vertical exchange coefficient between surface and bottom waters (/y)
%%     ~~Describe \code{Exchange} here~~
}
  \item{L}{
  Thickness of bottom water layer (m)
%%     ~~Describe \code{L} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 Returns the rate of change of bottom water oxygen with respect to changes in vertical exchange ([mol O2 * y]/m3)

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (BOD, Exchange, L) 
{
    return(BOD/(L * Exchange^2))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
