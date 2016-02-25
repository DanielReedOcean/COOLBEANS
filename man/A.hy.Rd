\name{A.hy}
\alias{A.hy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximum hypoxic area
}
\description{
This function returns the largest hypoxic area that can be support for the given conditions.}
\usage{
A.hy(NPP, Exchange, O2.surf, Depth, Area, L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{NPP}{
Net primary production (mol C/[m2 y])
%%     ~~Describe \code{NPP} here~~
}
  \item{Exchange}{
Exchange coefficient between surface and bottom waters (/y)
%%     ~~Describe \code{Exchange} here~~
}
  \item{O2.surf}{
Surface oxygen concentration (mol/m3)
%%     ~~Describe \code{O2.surf} here~~
}
  \item{Depth}{
Water depth (m)
%%     ~~Describe \code{Depth} here~~
}
  \item{Area}{
Area of primary production (m2)
%%     ~~Describe \code{Area} here~~
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
Returns the largest hypoxic area supported by given conditions (m2)
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
function (NPP, Exchange, O2.surf, Depth, Area, L) 
{
    M.sed <- C.flux(NPP, Depth) * Area
    gamma <- p.remin(Depth)
    return(M.sed * (gamma * 77/60 - 17/60)/(L * Exchange * (O2.surf - 
        0.063)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
