\name{runBasicAnalysis}
\alias{runBasicAnalysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Run basic analysis for a scRNA-Seq dataset
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
runBasicAnalysis(
  disease,
  path,
  annotate = TRUE,
  userlabel,
  usercelltype,
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{disease}{
  %%     ~~Describe \code{x} here~~
  \code{The name of the disease under study}
  }
  \item{path}{
  \code{The path to the directory containing the data}
  }
  \item{annotate}{
  \code{Define whether the data is annotated or if annoation is required (default=TRUE). If FALSE then annotation will be performed automatically using the 100 most significant marker genes for each cell type across three marker gene databases ("PanglaoDB_Augmented_2021","CellMarker_Augmented_2021","Tabula_Sapiens").}
  }
  \item{userlabel}{
  \code{Define the name of the SEURAT object column containing disease-control labels (if data is not annotated this needs to be defined, otherwise if will take the default value 'label')}
  }
  \item{usercelltype}{
  \code{Define the name of the SEURAT object column containing cell cluster IDs (if data is not annotated this needs to be defined, otherwise if will take the default value 'celltype')}
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The path paremter should point to a directory containing either of the following: 1) The data supplied in seperate subdirectories for all the samples in the analysis. The names given to the sample subdirectories will be used later for labeling the samples and therefore have to be structured correctly in the form of: unique label followed by a number (e.g., Disease1, Disease2, Disease3, Control1, Control2, Control3). 2) The data in an already inetgrated format in a single subdirectory within the directory provided by the path parameter.

Data can be provided in any of the following supported formats:
1) barcodes.tsv, genes.tsv, matrix.mtx, 2) H5 format, 3) SEURAT objects. If a metadata file exists with further annontations it may be supplied in the individual subdirectory/ies as a file with the name "meta.txt".

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{disease }{A Character string with the name of disease}
\item{path }{A Character string with a valid path to the directory containing the scRNA-Seq data.}
\item{annotate }{A Boolean (TRUE/FALSE) denoting whether to annotate the data (default=TRUE).}
\item{userlabel }{A  Character string with the name of the SEURAT object column which contains the information labels with respect to the conditions of the experiment (disease-control labels).}
\item{usercelltype }{A  Character string with the name of the SEURAT object column which contains the annotated cell cluster IDs.}
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
