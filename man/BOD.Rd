\name{BOD}
\alias{BOD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Benthic oxygen demand
%%  ~~function to do ... ~~
}
\description{
Calculate the benthic oxygen demand of a coastal sediment for given conditions.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BOD(z, NPP, C.flux.extra = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{
  Water depth (m)
%%     ~~Describe \code{z} here~~
}
  \item{NPP}{
  Net primary production (mol C/[m2 y])
%%     ~~Describe \code{NPP} here~~
}
  \item{C.flux.extra}{
  Additional organic carbon fluxes to sediments (mol C/[m2 y])
%%     ~~Describe \code{C.flux.extra} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns benthic oxygen demand of sediment (mol O2/[m2 y])
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
function (z, NPP, T, C.flux.extra = 0) 
{
    C.flux.sed <- C.flux(NPP, z) + C.flux.extra
    Oxygen.demand <- (((77/60) * p.remin(z) - (17/60)) * C.flux.sed)
    if (any(NPP <= 0)) 
        warning(paste("Net primary production is ", NPP[NPP <= 
            0], " - must be positive.\n", sep = ""))
    if (any(z <= 0)) 
        warning(paste("Water depth is ", z[z <= 0], " - must be positive.\n", 
            sep = ""))
    return(Oxygen.demand)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
