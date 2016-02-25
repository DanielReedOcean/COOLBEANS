\name{C.flux}
\alias{C.flux}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Organic carbon flux to underlying sediments
%%  ~~function to do ... ~~
}
\description{
Calculates the depositional flux of organic carbon to underlying sediments for given conditions.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
C.flux(NPP, z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{NPP}{
  Net primary production (mol C/[m2 y])
%%     ~~Describe \code{NPP} here~~
}
  \item{z}{
  Water depth (m)
%%     ~~Describe \code{z} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Return the depositional flux of organic carbon to sediments (mol C/[m2 y]).
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
function (NPP, z) 
{
    return(NPP * ((1 - 0.17) * exp(-0.018 * z) + 0.17 * exp(-0.00046 * 
        z)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
