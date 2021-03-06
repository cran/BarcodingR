% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FMFtheta12.R
\name{FMFtheta12}
\alias{FMFtheta12}
\title{Calculate Intraspecific and Interspecific Variation}
\usage{
FMFtheta12(ref)
}
\arguments{
\item{ref}{object of class "DNAbin" used as a reference dataset, which contains taxon information.}
}
\value{
a data frame containing intraspecific (sd, theta1) and interspefic variation (mean) of all species,
and their corresponding nearest neighbor (NN).
}
\description{
Calculation intraspecific variation (sd) of the potential species theta1, and mean interspecific
distance (here, the mean distance between the potential species and its nearest neighbor theta2) 
(fuzzy-set based method,slightly modified from Zhang et al. 2012). The calculation was done for all species in the 
reference dataset.
}
\examples{

data(TibetanMoth) 
ref<-as.DNAbin(as.character(TibetanMoth[1:50,]))
FMF.theta12<-FMFtheta12(ref)
FMF.theta12
}
\references{
Zhang, A. B., C. Muster, H.B. Liang, C.D. Zhu, R. Crozier, P. Wan, J. Feng, R. D. Ward.(2012). A fuzzy-set-theory-based approach 
to analyse species membership in DNA barcoding. Molecular Ecology, 21(8):1848-63.
}
\author{
Ai-bing ZHANG, PhD. CNU, Beijing, CHINA, contact at zhangab2008 (at) mail.cnu.edu.cn.
}
\keyword{FMFtheta12}
