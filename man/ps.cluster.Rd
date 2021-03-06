\name{ps.cluster}
\alias{ps.cluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the prediction strength of a clustering model
}
\description{
This function computes the prediction strength of a clustering model as published in R. Tibshirani and G. Walther 2005.
}
\usage{
ps.cluster(cl.tr, cl.ts, na.rm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cl.tr}{
Clusters membership as defined by the original clustering model, i.e. the one that was not fitted on the dataset of interest.
}
  \item{cl.ts}{
Clusters membership as defined by the clustering model fitted on the dataset of interest.
}
  \item{na.rm}{
\code{TRUE} if missing values should be removed, \code{FALSE} otherwise.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{ps }{the overall prediction strength (minimum of the prediction strengths at cluster level).}
  \item{ps.cluster }{Prediction strength for each cluster}
  \item{ps.individual }{Prediction strength for each sample.}
}
\references{
R. Tibshirani and G. Walther (2005) "Cluster Validation by Prediction Strength", \emph{Journal of Computational and Graphical Statistics}, \bold{14}(3):511--528.
}
\author{
Benjamin Haibe-Kains
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
%%\examples{
%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ clustering }
