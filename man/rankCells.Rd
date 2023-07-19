\name{rankCells}
\alias{rankCells}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Rank all cells in the anlaysis by mapping to prior knowledge
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rankCells(
  seuratObject,
  path,
  scan="Cell",
  priorknowledgePathsKEGG,
  priorknowledgePathsGO,
  priorknowledgePathsMSIG,
  priorknowledgePathsWiki,
  priorknowledgePathsReact,
  priorknowledgeMOA,
  labels,
  cellIDs,
  checkdrug=TRUE
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{seuratObject}{
  \code{The SEURAT object to analyse}
  }
  \item{path}{
  \code{Full path to where the data data was extracted}
  }
  \item{scan}{
  \code{Define whether to scan single-cell or simulate bulk RNA (takes two values defualt="Cell" or "Bulk")}
  }
  \item{priorknowledgePathsKEGG}{
  \code{Define the prior knowdledge from KEGG pathways}
  }
  \item{priorknowledgePathsGO}{
  \code{Define the prior knowdledge from GO Biological processes}
  }
  \item{priorknowledgePathsMSIG}{
  \code{Define the prior knowdledge from MSIG pathways}
  }
  \item{priorknowledgePathsWiki}{
  \code{Define the prior knowdledge from WIKI pathways}
  }
  \item{priorknowledgePathsReact}{
  \code{Define the prior knowdledge from REACTOME pathways}
  }
  \item{priorknowledgeMOA}{
  \code{Define the prior knowdledge drugs or drug MOAs from CMAP}
  }
  \item{labels}{
  \code{Define the name of the SEURAT object column containing disease-control labels (if data is not annotated this needs to be defined, otherwise if will take the default value 'label')}
  }
  \item{cellIDs}{
  \code{Define the name of the SEURAT object column containing cell cluster IDs (if data is not annotated this needs to be defined, otherwise if will take the default value 'celltype')}
  }
  \item{checkdrug}{
  \code{Define whether to use drug names (MalaCards approach=default), or dug mode of actions (MOAs) (Hypothesis-driven)}
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Works with the output of the basic analysis function which is in the form of a SEURAT object. Sample labels have to be supplied in order to perform the analysis for differentially expressed genes (DEGs) identification. Genes are used to predict pathways using five different databases supported by enrichR (KEGG, GO, REACTOME, MSIG and WIKI). Drug repurposing is also performed in order to obtain the drugs (and their MOAs) that reverse the expression profiles obtained by the DEGs. Finally this function performs a mapping step and returns the rankings per cell-type across pathwways and drugs.
}
\value{
%%  ~Describe the value returned
\item{seuratObject }{A SEURAT object obtained from the basic analysis funtion}
\item{scan }{A Character string the type of analysis to perform (Cell of Bulk).}
\item{priorknowledgePathsKEGG }{A  Character vector with the prior knowledge information from KEGG.}
\item{priorknowledgePathsGO }{A  Character vector with the prior knowledge information from GO .}
\item{priorknowledgePathsMSIG }{A  Character vector with the prior knowledge information from MSIG.}
\item{priorknowledgePathsWiki }{A  Character vector with the prior knowledge information from WIKI.}
\item{priorknowledgePathsReact }{A  Character vector with the prior knowledge information from REACTOME.}
\item{priorknowledgeMOA }{A  Character vector with the names of drugs of drug MOAs from CMAP to be used as prior knowledge.}
\item{labels }{A  Character String with the name of the column in the SEURAT object defining the sample labels to be used for DEG identification.}
\item{cellIDs }{A  Character String with the name of the column in the SEURAT object defining the annotated cell IDs.}
\item{checkdrug }{A Boolean (TRUE/FALSE) denoting whether to use drug names (default=TRUE) or alternatively MOAs (FALSE).}
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
