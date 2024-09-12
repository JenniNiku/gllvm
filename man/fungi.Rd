\name{fungi}
\docType{data}
\alias{fungi}
\title{Wood-decaying fungi data}
\description{
Dataset of 1666 binary observations for 215 fungal species, in different 53 European Beech forests spread across 8 regions.
}
\usage{
data(fungi)
}
\format{
\describe{
\item{Y}{A data frame with the presence-absences of 215 fungal species measured at 1666 logs.}
\item{X}{A data frame of 8 predictor variables.}
\describe{
\item{DBH.CM}{Diameter at breast height (cm)} 
\item{AVERDP}{Decay stage of logs on a 1-5 scale}
\item{CONNECT10}{Connectivity of the surrounding forest at 10km scale}
\item{TEMPR}{Annual temperature range in degrees Celsius}
\item{PRECIP}{Annual precipitation in milimeters}
\item{log.Area}{ln(area in hectares) of reserves}
\item{REGION}{Site groups identified based on spatial clusters}
\item{RESERVE}{Site name}
}
\item{TR}{A data frame of the traits used in Abrego et al. (2022).}
\item{tree}{The phylogenetic tree.}
\item{C}{The phylogenetic covariance matrix.}
\item{dist}{The phylogenetic distance matrix.}
\describe{
}
}

}
\details{
Observations of fungi species inhabiting European beech logs, in different European countries. The countries have been grouped in eight different regions. Logs were surveyed in 53 different reserves (or sites). Included environment and trait covariates are limited to those analyzed in the original article, though more are available in the published dataset on datadryad.org.
}

\references{
Abrego, N., Bässler, C., Christensen, M., and Heilmann‐Clausen, J. (2022). Traits and phylogenies modulate the environmental responses of wood‐inhabiting fungal communities across spatial scales. Journal of Ecology, 110(4), 784-798.

Abrego, N., Bässler, C., Christensen, M., and Heilmann-Clausen, J. (2022). Data and code from: Traits and phylogenies modulate the environmental responses of wood-inhabiting fungal communities across spatial scales [Dataset]. Dryad. https://doi.org/10.5061/dryad.t76hdr82r
}

\examples{
\dontrun{
data(fungi)
Y <- fungi$Y
X <- fungi$X
TR <- fungi$TR
C <- fungi$C
dist <- fungi$dist

model <- gllvm(y = Y, X = cbind(int = 1, X), TR = TR, 
        formula = ~DBH.CM + AVERDP + I(AVERDP^2) + CONNECT10 + TEMPR + PRECIP + 
        log.AREA + (DBH.CM + AVERDP + I(AVERDP^2) + CONNECT10 + TEMPR + PRECIP + 
        log.AREA):(FB.type + Sp.log.vol.µ3 + Lifestyle), 
        family = "binomial", num.lv = 0, studyDesign = X[,c("REGION", "RESERVE")], 
                     colMat = list(C, dist = dist), colMat.rho.struct = "term", 
                     row.eff = ~(1 | REGION/RESERVE), sd.errors = FALSE, 
                     randomX = ~int + DBH.CM + AVERDP + I(AVERDP^2) + 
                                CONNECT10 + TEMPR + PRECIP + log.AREA, 
                     beta0com = TRUE, nn.colMat = 10, maxit = 20000)
}
}