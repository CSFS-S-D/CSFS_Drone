\name{smoothRetile}
\alias{smoothRetile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
smoothRetile
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Homogenizes point cloud density in a set of overlapping .las files and retiles 
the data into files that do not overlap.
}
\usage{
smoothRetile(las_dir, output_dir, chunk_size = 200, res = 1, pct = 1, density = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{las_dir}{
%%     ~~Describe \code{x} here~~
  Path to a directory with a set of overlapping .las files.
  }
  \item{output_dir}{
  Path to a directory where output .las tiles will be written.
  }
  \item{chunk_size}{
  Length in the units of the data of the sides of square output tiles (or something like that).
  }
  \item{res}{
  numeric. The resolution of the grid used to filter the point cloud. See \link{decimate_points} for details.
  }
  \item{pct}{
  numeric. The proportion of points in the original point cloud to keep.
  }
  \item{density}{
  numeric. Optionally, a target final point density. Set/leave pct = 1 if using.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Point clouds derived from aerial sources often have non-uniform point density 
due to flight line and mission overlap, which can cause erroneous patterns in 
data products. This function is meant to produce a point cloud with uniform point
density. If a target density is not provided, it sets the target point density to  
the proportion specified with "pct" of the lowest point density of the input files.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
An LAScatalog object
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
