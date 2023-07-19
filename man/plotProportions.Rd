\name{plotProportions}
\alias{plotProportions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the cell proportions output
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotProportions(filename, title = "Cell Proportion Rankings")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
%%     ~~Describe \code{filename} here~~
  \code{Name of the file with the cell proportions output}
}
  \item{title}{
%%     ~~Describe \code{title} here~~
  \code{Title for the plot (Default = 'Cell Proportion Rankings')}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{filename }{A character string of the file name with the cell proportions output}
\item{title }{A character string to use as title for the plot (Default = 'Cell Proportion Rankings')}
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
function (filename, title = "Cell Proportion Rankings")
{
    proportions <- read.table(filename, sep = "\t", header = T,
        as.is = T, row.names = 1)
    colnames(proportions) <- c("CellIDs", "Condition", "Proportion")
    proportions2 <- proportions \%>\% arrange(CellIDs) \%>\% group_by(CellIDs) \%>\%
        mutate(rel_inc = abs(Proportion - lag(Proportion, default = first(Proportion))))
    p <- ggplot(data = proportions2, aes(x = reorder(CellIDs,
        rel_inc, FUN = max, decreasing = TRUE), y = Proportion,
        color = Condition, fill = Condition)) + geom_bar(stat = "identity",
        position = position_dodge()) + geom_segment(aes(x = 8,
        y = max(Proportion), xend = 10, yend = max(Proportion)),
        arrow = arrow(length = unit(0.5, "cm"), ends = "last")) +
        geom_label(aes(x = 9, y = max(Proportion) - 0.04, label = "Dec. Diff. in Cell Prop."),
            colour = "white", fontface = "bold", vjust = -0.25,
            show.legend = FALSE) + theme(plot.title = element_text(size = 11,
        )) + ggtitle(title) + xlab("") + ylab("Proportion of cells")
    p <- p + theme(plot.title = element_text(size = 15, face = "bold"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15,
            face = "bold"), legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))
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
