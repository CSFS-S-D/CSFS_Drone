\name{snagMapper}
\alias{snagMapper}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
snagMapper
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A function to find snags in a lidar point cloud.
}
\usage{
snagMapper(pc, Lint, Uint)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pc}{
  %%     ~~Describe \code{x} here~~
  An LAS point cloud with return intensity information
  }
  \item(Lint){
  The lower intensity threshold for the point cloud. The intensity threshold between 
  bole and branch returns and foliage returns.
  }
  \item(Unit){
  The upper intensity threshold for the point cloud. This threshold attempts to
  separate older snags that have lost all their bark and are very reflective from
  younger snags that still have bark.
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
