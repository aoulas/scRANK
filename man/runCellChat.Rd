\name{runCellChat}
\alias{runCellChat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Run CellChat for diease and control samples for each cell type in the experiment
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
runCellChat(seuratObject,labels,cellIDs,LablesUniq=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seuratObject}{
  \code{The SEURAT object to analyse}
  }
  \item{labels}{
  \code{Define the name of the SEURAT object column containing disease-control labels (if data is not annotated this needs to be defined, otherwise if will take the default value 'label')}
  }
  \item{cellIDs}{
  \code{Define the name of the SEURAT object column containing cell cluster IDs (if data is not annotated this needs to be defined, otherwise if will take the default value 'celltype')}
  }
  \item{LablesUniq}{
  \code{Define the labels (2) to perform pairwise comparisons}
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{seuratObject }{A SEURAT object obtained from teh basic analysis funtion}
\item{labels }{A  Character String with the name of the column in the SEURAT object defining the sample labels to be used for DEG identification.}
\item{cellIDs }{A  Character String with the name of the column in the SEURAT object defining the annotated cell IDs.}
\item{LablesUniq }{A  Character vector with the names of sample labels to be used for DEG identification (default=NULL). If NULL the labels parameter will be used to define the sample classes.}  
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
