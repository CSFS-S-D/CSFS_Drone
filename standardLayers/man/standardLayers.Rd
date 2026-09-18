\name{standardLayers}
\alias{standardLayers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
standardLayers
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Creates and saves raster and vector files that are useful to foresters from a set of point clouds.
}
\usage{
standardLayers(las_dir, output="./", bounds=NULL, ortho=NULL, dtm_res=1,
chm_res=1, ws=cloud2trees::itd_ws_functions()$log_fn)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{las_dir}{
  %%     ~~Describe \code{x} here~~
  Path to a directory containing a set of las files to be processed.
  }
  \item{output}{
  Path to a directory where files will be written.
  }
  \item{bounds}{
  Path to a polygon file with project boundaries. Products will be
  spatially unconstrained if no bounds are provided.
  }
  \item{ortho}{
  Path to an ortomosaic. If provided, this function will crop an
  orthomosaic to the project area. It's for convenience really.
  }
  \item{chm_res}{
  numeric. canopy height model resolution.
  }
  \item{dtm_res}{
  numeric. Digital terrain model resolution.
  }
  \item{ws}{
  numeric or function. Window size for individual tree detection via lidR::lmf().
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
