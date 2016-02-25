\name{Loading.sensitivity}
\alias{Loading.sensitivity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bottom water oxygen sensitivity to nutrient loading
%%  ~~function to do ... ~~
}
\description{
Calculates the rate of change of bottom water oxygen with respect to changes in nitrogen loading.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Loading.sensitivity(Depth, Exchange, L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Depth}{
  Water depth (m)
%%     ~~Describe \code{Depth} here~~
}
  \item{Exchange}{
  Vertical exchange between surface and bottom waters (y)
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
Returns the rate of change of bottom water oxygen with respect to changes in nitrogen loading ([mol O2 * y]/[mol N * m2])
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
function (Depth, Exchange, L) 
{
    gamma <- p.remin(Depth)
    return(-10 * ((77/60) * gamma - (17/60))/(Exchange * L * 
        ((733/1500) * gamma + (1007/1500))))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
