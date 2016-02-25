\name{C.burial}
\alias{C.burial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Carbon burial in sediments
%%  ~~function to do ... ~~
}
\description{
Calculates the amount of carbon buried in sediments for a given depth according to Middelburg et al. (1997).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
C.burial(z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{
  Water depth (m)
%%     ~~Describe \code{z} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns carbon burial rate (mol/[m2 y])
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Jack J. Middelburg, Karline Soetaert, Peter M.J. Herman (1997) Empirical relationships for use in global diagenetic models, Deep Sea Research Part I: Oceanographic Research Papers, 44(2), 327–344
%% ~put references to the literature/web site here ~
}
\author{
Daniel C. Reed
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
function (z) 
{
    return(10 * 4.4 * 10^{
        -0.84672973 - 0.00061506 * z
    })
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
