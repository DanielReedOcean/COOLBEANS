\name{p.remin}
\alias{p.remin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Proportion of organic matter remineralised
%%  ~~function to do ... ~~
}
\description{
Calculates the proportion of organic matter deposited in sediments that is ultimately remineralised using empirical functions from Middelburg et al. (1997).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
p.remin(z)
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
Returns the proportion of organic matter that is remineralised in a sediment of a given depth (unitless).
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
\code{\link{C.remin}},
\code{\link{C.burial}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (z) 
{
    return(C.remin(z)/(C.remin(z) + C.burial(z)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
