\name{signTest}
\alias{signTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Exact Sign Test with Confidence Intervals
}
\description{
Uses \code{binom.exact} to create sign test with confidence intervals.
}
\usage{
signTest(x,  stat = c("pos-neg", "prop pos"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
numeric vector
}
  \item{stat}{
statistic for estimates and confidence intervals,
"pos-neq"= proportion positive - proportion negative, "prop pos"= proportion positive 
(all proportions are out of non-zero values)
}
  \item{\dots}{
arguments passed to \code{binom.exact}
}
}
\details{
Exact sign test, using exact binomial test on the proportion positive out of the non-zero values.
}
\value{
A list of class 'htest'
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{%% ~put references to the literature/web site here }
\author{
Michael P. Fay
}
\note{
The sign test can be interpreted as a test that the median is zero assuming continuous data. If you want to test on the median without making continuity assumptions use \code{\link{medianTest}}.
}
%\seealso{%\code{\link{help}}, ~~~}
\examples{
x<-c(rep(-1,10),rep(0,60),rep(1,30))
signTest(x)
# sample median is zero, 
# and not surprisingly the median test 
# properly gives a large p-value
medianTest(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}