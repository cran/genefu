\name{gene70}
\alias{gene70}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the 70 genes prognosis profile (GENE70) as published by van't Veer et al. 2002
}
\description{
This function computes signature scores and risk classifications from gene expression values following the algorithm used for the 70 genes prognosis profile (GENE70) as published by van't Veer et al. 2002.
}
\usage{
gene70(data, logged2 = TRUE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Matrix of gene expressions with samples in rows and probes in columns, dimnames being properly defined.
}
  \item{logged2}{
\code{FALSE} if gene expressions must be transformed in log2, \code{TRUE} otherwise.
}
  \item{verbose}{
\code{TRUE} to print informative messages, \code{FALSE} otherwise.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
	\item{score }{Continuous signature scores}
	\item{risk }{Binary risk classification, \code{1} being high risk and \code{0} being low risk.}
}
\references{
L. J. van't Veer and H. Dai and M. J. van de Vijver and Y. D. He and A. A. Hart and M. Mao and H. L. Peterse and K. van der Kooy and M. J. Marton and A. T. Witteveen and G. J. Schreiber and R. M. Kerkhiven and C. Roberts and P. S. Linsley and R. Bernards and S. H. Friend (2002) "Gene Expression Profiling Predicts Clinical Outcome of Breast Cancer", \emph{Nature}, \bold{415}:530--536.
}
\author{
Benjamin Haibe-Kains
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[genefu]{nki}}
}
\examples{
## load GENE70 signature
data(sig.gene70)
## load NKI dataset
data(nki)
## compute relapse score
rs.nki <- gene70(data=data.nki)
table(rs.nki$risk)
## note that the discrepancies compared to the original publication are closed to the official cutoff, raising doubts on its exact value
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ prognosis }
%%\keyword{ breast cancer }
