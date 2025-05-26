\name{Skabbholmen}
\docType{data}
\alias{Skabbholmen}
\title{Skabbholmen island data}
\description{
Dataset of ordinal observations of plants, on the island Skabbholmen in the Stocholm archipelago. Includes 65 unique sites and 70 species, surveyed in two different years.
}
\usage{
data(Skabbholmen)
}
\format{
\describe{
\item{Y}{A data frame with ordinal of 70 plant species measured at 126 plots.}
\item{X}{A matrix of 2 predictor variables at 126 plots.}
\item{species}{A matrix of full species names and abbreviations used in the community data (Y).}
}

}
\details{
Observations of vascular plant cover in 126 one-square-meter plots divided over four transects. The ordinal responses are on a five-degree Hult-Sernander-Du Rietz scale, and were originally recorded by Wolfgang and Cramer (1987) and additionally analyzed by ter Braak (1987). There is a total of 64 unique sites, that were surveyed in two different years (1978 and 1984), but two plots were only surveyed in one year (thus bringing the total number of rows in the data to 126). The plots were located on an elevation gradient, running from the shoreline to the edge of old-growth forest. Elevation to the shoreline was recorded in centimeters during the sampling in 1978.

This dataset was published with permission from the CANOCO FORTRAN package example datasets.
}

\references{
ter Braak, C.J.F. and Smilauer, P.  (1998).  CANOCO reference manual and user's guide to CANOCO for Windows: software for canonical community ordination (version 4). Microcomputer Power, New York, New York, USA.

ter Braak, C.J.F. (1987). The analysis of vegetation-environment relationships by canonical correspondence analysis. Vegetatio, 69(1), 69-77.

Cramer, W. & Hytteborn, H. (1987). The separation of fluctuation and long-term change in vegetation dynamics of a rising seashore. Vegetatio, 69, 157–167.
}
\examples{
\dontrun{
data(Skabholmen)
Y <- Skabbholmen$Y
X <- Skabbholmen$X
model <- gllvm(y = Y, X = X, studyDesign = X,
    lv.formula = ~scale(Elevation) + I(Year-1978),
    num.RR = 2, 
    family = "ordinal",
    zeta.struc="common",
    row.eff=~(1|transectID))
}
}