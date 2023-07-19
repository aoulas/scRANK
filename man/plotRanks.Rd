\name{plotRanks}
\alias{plotRanks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the cell rankings output
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotRanks(filename, title = "Rankings")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
%%     ~~Describe \code{filename} here~~
  \code{Name of the file with Rankings output}
}
  \item{title}{
%%     ~~Describe \code{title} here~~
  \code{Title for the plot (Default = 'Rankings')}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{filename }{A character string of the file name with the cell rankings output}
\item{title }{A character string to use as title for the plot (Default = 'Rankings')}
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
function (filename, title = "Rankings")
{
    Ranks <- read.table(filename, sep = "\t", header = T, as.is = T,
        row.names = 1)
    Ranks <- round(Ranks, 2)
    Ranks <- Ranks[vapply(Ranks, function(x) length(unique(x)) >
        1, logical(1L))]
    CellIDS <- rownames(Ranks)
    Ranks <- as.data.frame(sapply(Ranks, rank))
    MeanRanks <- as.data.frame(apply(Ranks, 1, mean))
    rownames(MeanRanks) <- CellIDS
    MeanRanks <- cbind(CellIDS, MeanRanks)
    MeanRanks <- MeanRanks[order(MeanRanks$`apply(Ranks, 1, mean)`,
        decreasing = T), ]
    indexBulk <- which(MeanRanks$CellIDS == "Bulk")
    colsbulk <- c(rep("black", length(CellIDS)))
    colsbulk[indexBulk] <- "red"
    Ranks <- Ranks \%>\% gather(Parameter, Rank)
    Ranks <- cbind(CellIDS, Ranks)
    Ranks$CellIDS <- gsub("Lymph.Node.cd1c.positive.Myeloid.Dendritic.Cell",
        "Myeloid.Dendritic.Cell", Ranks$CellIDS)
    p <- ggplot(data = Ranks, aes(x = reorder(CellIDS, Rank,
        FUN = mean, decreasing = TRUE), y = Rank)) + geom_boxplot(outlier.colour = NA,
        color = colsbulk) + guides(color = "none") + geom_dotplot(aes(color = Parameter,
        fill = Parameter), binaxis = "y", stackdir = "center",
        dotsize = 0.3, position = position_dodge(), binwidth = 0.5) +
        stat_summary(fun.y = mean, aes(shape = "average Rank"),
            geom = "point", size = 2, color = "red", fill = "red") +
        scale_shape_manual("", values = c(`average Rank` = 24)) +
        scale_y_continuous(breaks = seq(0, max(Ranks$Rank), 2)) +
        theme_bw() + theme(axis.text = element_text(size = 15),
        axis.text.y = element_text(color = colsbulk), axis.title = element_text(size = 15,
            face = "bold"), legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) + ggtitle(title) +
        xlab("") + ylab("Rank")
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
