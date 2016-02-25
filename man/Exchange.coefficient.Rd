\name{Exchange.coefficient}
\alias{Exchange.coefficient}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Vertical exchange coefficient between surface and bottom waters
%%  ~~function to do ... ~~
}
\description{
Calculates the rate constant of exchange between surface and bottom waters.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Exchange.coefficient(BOD, O2.surf, O2.bottom, L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BOD}{
  Benthic oxygen demand (mol O2/[m2 y])
%%     ~~Describe \code{BOD} here~~
}
  \item{O2.surf}{
  Surface oxygen concentration (mol/m3)
%%     ~~Describe \code{O2.surf} here~~
}
  \item{O2.bottom}{
  Bottom oxygen concentration (mol/m3)
%%     ~~Describe \code{O2.bottom} here~~
}
  \item{L}{
  Thickness of bottom layer (m)
%%     ~~Describe \code{L} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns vertical exchange coefficient between surface and bottom waters (/y)
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
function (BOD, O2.surf, O2.bottom, L) 
{
    Exchange <- BOD/((O2.surf - O2.bottom) * L)
    mask <- (0.05 > (O2.surf - O2.bottom)/O2.surf | O2.bottom > 
        O2.surf)
    Exchange[mask] <- BOD[mask]/(0.1 * O2.surf[mask] * L)
    return(Exchange)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
