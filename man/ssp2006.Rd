\name{ssp2006}
\alias{ssp2006}
\docType{data}
\title{
SSP2006 classifier for identification of breast cancer molecular subtypes (Hu et al 2006)
}
\description{
List of parameters defining the SSP2006 classifier for identification of breast cancer molecular subtypes (Hu et al 2006). 
}
\usage{data(ssp2006)}
\format{
List of parameters for SSP2006:
\describe{
\item{\code{centroids} }{Gene expression centroids for each subtype.}
\item{\code{centroids.map} }{Mapping for centroids.}
\item{\code{method.cor} }{Method of correlation used to compute distance to the centroids.}
\item{\code{method.centroids} }{Method used to compute the centroids.}
\item{\code{std} }{Method of standardization for gene expressions.}
\item{\code{mins} }{Minimum number of samples within each cluster allowed during the fitting of the model.}
}
}
%%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%%}
\source{
\url{http://www.biomedcentral.com/1471-2164/7/96}
}
\references{
Hu, Zhiyuan and Fan, Cheng and Oh, Daniel and Marron, JS and He, Xiaping and Qaqish, Bahjat and Livasy, Chad and Carey, Lisa and Reynolds, Evangeline and Dressler, Lynn and Nobel, Andrew and Parker, Joel and Ewend, Matthew and Sawyer, Lynda and Wu, Junyuan and Liu, Yudong and Nanda, Rita and Tretiakova, Maria and Orrico, Alejandra and Dreher, Donna and Palazzo, Juan and Perreard, Laurent and Nelson, Edward and Mone, Mary and Hansen, Heidi and Mullins, Michael and Quackenbush, John and Ellis, Matthew and Olopade, Olufunmilayo and Bernard, Philip and Perou, Charles (2006) "The molecular portraits of breast tumors are conserved across microarray platforms", \emph{BMC Genomics}, \bold{7}(96)
}
\examples{
data(ssp2006)
str(ssp2006)
}
\keyword{data}