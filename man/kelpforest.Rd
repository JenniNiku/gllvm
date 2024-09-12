\name{kelpforest}
\docType{data}
\alias{kelpforest}
\title{Kelp Forest community data}
\description{
Dataset of 836 percent cover records for 132 algae and invertebrates, from 2000 to 2020 in the Santa Barbara Channel, California, USA.
}
\usage{
data(fungi)
}
\format{
\describe{
\item{Y}{A data frame with the percent cover records of 132 sessile invertebrates and understory macroalgae measured at 836 permanent transects.}
\item{X}{A data frame of study design variables and predictors.}
\describe{
\item{SITE}{Kelp forest site} 
\item{YEAR}{Sampling year}
\item{TRANSECT}{Permanent transect identifying number (unique within site), nested within sites}
\item{KELP_FRONDS}{A number of stipes of giant kelp}
\item{PERCENT_ROCKY}{percent rock coverage}
}
\item{SPinfo}{A data frame of the species information including species names, group and taxonomy.}
\describe{
\item{SP_CODE}{Species codes corresponding species names in abundance matrix} 
\item{GROUP}{Species group; algae, invertebrate or plant}
\item{COMMON_NAME, SCIENTIFIC_NAME}{Species' common and scientific names}
\item{TAXON_KINGDOM,TAXON_PHYLUM,TAXON_CLASS,TAXON_ORDER,TAXON_FAMILY, TAXON_GENUS}{Species taxonomic information}
}
}

}
\details{
Data consists of the cover of sessile invertebrates and understory macroalgae starting from year 2000 to 2020. The presence of invertebrates and macroalgae are recorded using a point contact method; 80 uniformly spaced points along permanent 40m x 2m transects, such that a percent cover of a given species on a transect is estimated as the fraction of total points at which that species was present. Data are from SBC LTERs kelp forest monitoring program in the Santa Barbara Channel, California, USA. Detailed description of the data available in the reference and the website  https://sbclter.msi.ucsb.edu/data/catalog/package/?package=knb-lter-sbc.15
}

\references{
Reed, D, R. Miller. 2023. SBC LTER: Reef: Kelp Forest Community Dynamics: Cover of sessile organisms, Uniform Point Contact ver 33. Environmental Data Initiative. https://doi.org/10.6073/pasta/0af1a5b0d9dde5b4e5915c0012ccf99c. Accessed 2024-09-11.
}

\examples{
\dontrun{
data(kelpforest)
Y <- kelpforest$Y
X <- kelpforest$X
SPinfo <- kelpforest$SPinfo
}
}