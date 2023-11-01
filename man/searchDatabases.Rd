\name{searchDatabases}
\alias{searchDatabases}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Search Individual Databases for Malacards Terms
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
searchDatabases(
  disease,
  path,
  scenario="Malacards",
  checkdrug = TRUE,
  keywordsWikiUser,
  keywordsKEGGUser,
  keywordsGOUser,
  keywordsMSIGUser,
  keywordsReactUser,
  keywordsMOAUser
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
  \item{scenario}{
  \code{Define which scenario/approcah to use (takes two values: default="Malacards", or "Hypothesis")}
  }
  \item{checkdrug}{
  \code{Define whether to use drug names for the Malacards scenario (default=TRUE), or dug mode of actions (MOAs) (Hypothesis scenario)}
  }
  \item{keywordsWikiUser}{
  \code{Define the prior knowdledge from WIKI pathways (only for Hypothesis-driven scenario)}
  }
  \item{keywordsKEGGUser}{
  \code{Define the prior knowdledge from KEGG pathways (only for Hypothesis-driven scenario)}
  }
  \item{keywordsGOUser}{
  \code{Define the prior knowdledge from GO Biological processes (only for Hypothesis-driven scenario)}
  }
  \item{keywordsMSIGUser}{
  \code{Define the prior knowdledge from MSIG pathways (only for Hypothesis-driven scenario)}
  }
  \item{keywordsReactUser}{
  \code{Define the prior knowdledge from REACTOME pathways (only for Hypothesis-driven scenario)}
  }
  \item{keywordsMOAUser}{
  \code{Define the prior knowdledge from CMAP drug MOAs (only for Hypothesis-driven scenario)}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Two options for this function to work 1) For the Malacards scenario, this function works with the output of the extractMalacards() function. This extracts prior knowledge output files which should be saved in the same path/diretory as the one containing the data. 2) For the Hypothesis scenario the user has to supply the prior knowledge. The function will then scan the relevant databases and return the results for inspection.
}
}
\value{
\item{disease }{A Character string with the name of disease}
\item{path }{A Character string with a valid path to the directory containing the scRNA-Seq data.}
\item{scenario }{A  Character string with the name of scenarion (default="Malacards") alternative option="Hypothesis".}
\item{checkdrug }{A Boolean (TRUE/FALSE) denoting whether use drug names (default=TRUE) alternative option for MOAs=FALSE.}
\item{keywordsWikiUser }{A  Character vector with the names of WIKI pathways to be used as prior knowledge.}
\item{keywordsKEGGUser }{A  Character vector with the names of KEGG pathways to be used as prior knowledge.}
\item{keywordsGOUser }{A  Character vector with the names of GO Biological processes to be used as prior knowledge.}
\item{keywordsMSIGUser }{A  Character vector with the names of MSIG pathways to be used as prior knowledge.}
\item{keywordsReactUser }{A  Character vector with the names of REACTOME pathways to be used as prior knowledge.}
\item{keywordsMOAUser }{A  Character vector with the names of drug MOAs from CMAP to be used as prior knowledge.}
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
