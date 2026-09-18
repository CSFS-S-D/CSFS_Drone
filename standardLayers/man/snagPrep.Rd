\name{snagPrep}
\alias{snagPrep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
snagPrep
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Prepares a lidar point cloud for snag detection. The function thins dense point clouds 
and makes an attempt to scale Intensity data. Most importantly, it tries to find intensity
thresholds between bole and branch returns and foliage returns.
}
\usage{
add_LasSnagAttributes(pc, ttops)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pc}{
  %%     ~~Describe \code{x} here~~
  An LAS point cloud with return intensity information
  }
  \item{ttops}{
  An simple feature collection of tree data (usually from cloud2trees).
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
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
##--	or standard data sets, see data().

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
