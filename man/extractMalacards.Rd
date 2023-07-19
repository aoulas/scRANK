\name{extractMalacards}
\alias{extractMalacards}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract prior knowledge from MalaCards
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
extractMalacards(disease,files,path)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
%%     ~~Describe \code{x} here~~
  \item{disease}{
  \code{The name of the disease under study}
  }
  \item{files}{
  \code{A MalaCards html file name to parse}
  }
  \item{path}{
  \code{Path where html files are stored and by defult will be used to also write output files}
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{disease }{A Character string with the name of disease}
\item{files }{A vector containing characters strings of file name(s) of downloaded html disease page(s) from MalaCards}
\item{path }{A character strings with the path to the directory where html files are stored. By defult the same path will also be used to write output files}
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
