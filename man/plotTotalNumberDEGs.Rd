\name{plotTotalNumberDEGs}
\alias{plotTotalNumberDEGs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the total DEGs output
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotTotalNumberDEGs(filename, title = "Total DEGs Rankings")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
  %%     ~~Describe \code{filename} here~~
  \code{Name of the file with Total Number of DEGs output}
  }
  \item{title}{
%%     ~~Describe \code{title} here~~
  \code{Title for the plot (Default = "Total DEGs Rankings")}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{filename }{A character string of the file name with the total DEGs output}
\item{title }{A character string to use as title for the plot (Default = 'Total DEGs Rankings')}
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
function (filename, title = "Total DEGs Rankings")
{
    TotalNumberDEGs <- read.table(filename, sep = "\t", header = T,
        as.is = T, row.names = 1)
    TotalNumberDEGs$CellID <- rownames(TotalNumberDEGs)
    p <- ggplot(data = TotalNumberDEGs, aes(x = reorder(CellID,
        TotalNumberDEGs, FUN = max, decreasing = FALSE), y = TotalNumberDEGs)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        theme(plot.title = element_text(size = 11, )) + ggtitle(title) +
        xlab("")
    p <- p + theme(plot.title = element_text(size = 15, face = "bold"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15,
            face = "bold"), legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))
    p <- p + coord_flip()
    print(p)
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
