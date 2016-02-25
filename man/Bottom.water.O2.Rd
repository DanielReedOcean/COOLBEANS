\name{Bottom.water.O2}
\alias{Bottom.water.O2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bottom water oxygen
%%  ~~function to do ... ~~
}
\description{
Calculate the bottom water oxygen concentration for a given set of conditions. This is the chief function of the COOLBEANS package with others aiding in calculating benthic oxygen demand and the vertical exchange coefficient required by this function.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Bottom.water.O2(O2.surf, BOD, Exchange, L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{O2.surf}{
  Oxygen concentration in surface waters (mol/m3)
%%     ~~Describe \code{O2.surf} here~~
}
  \item{BOD}{
  Benthic oxygen demand (mol O2/[m2 y])
%%     ~~Describe \code{BOD} here~~
}
  \item{Exchange}{
  Vertical exchange coefficient (/y)
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
Return bottom water oxygen concentration (mol/m3)
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
function (O2.surf, BOD, Exchange, L) 
{
    O2.bottom <- O2.surf - BOD/(Exchange * L)
    return(O2.bottom)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
