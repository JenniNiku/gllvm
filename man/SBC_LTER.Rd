\name{SBC_LTER}
\docType{data}
\alias{SBC_LTER}
\title{SBC LTER: Reef: Kelp Forest Community Dynamics: Cover of sessile organisms, Uniform Point Contact}
\description{
These data describe the cover of sessile invertebrates, understory macroalgae, and bottom substrate types as determined by a uniform point contact method. The presence of over 150 taxa of sessile invertebrates and macroalgae are recorded at 80 uniformly spaced points along permanent 40m x 2m transects. Multiple species can be recorded at any given point. Percent cover of a given species on a transect can be estimated from UPC observations as the fraction of total points at which that species was present x 100. The total percent cover of all species combined using this method can exceed 100\%; however, the percent cover of any single species cannot exceed 100\%. This specific version of the data includes 61 species of macroalgae, 69 species of sessile invertebrates, and two species of plants.
}
\usage{
data(SBC_LTER)
}
\format{
\describe{
\item{SPinfo}{Information on the species included.}
\item{Xenv}{A data frame containing the year, site name and transect number, and two environmental covariates; number of kelp fronds and rockiness of the seabed.}
\item{Yabund}{A data frame containing the recorded percent covers for the included macroalgae, invertebrates, and plants.}
}

}
\details{
These data are part of SBC LTERs kelp forest monitoring program, which began in 2000 and was designed to track long-term patterns in species abundance and diversity of reef-associated organisms in the Santa Barbara Channel, California, USA. The sampling locations in this dataset include nine reef sites along the mainland coast of the Santa Barbara Channel and at two sites on the north side of Santa Cruz Island. These sites reflect several oceanographic regimes in the channel and vary in distance from sources of terrestrial runoff. Data collection began in 2000 and this dataset is updated annually.

The time period of data collection varied among the 11 kelp forest sites. Sampling at BULL, CARP, and NAPL began in 2000, sampling at the other 6 mainland sites (AHND, AQUE, IVEE, GOLB, ABUR, MOHK) began in 2001 (transects 3, 5, 6, 7, 8 at IVEE were added in 2011). Data collection at the two Santa Cruz Island sites (SCTW and SCDI) began in 2004.
}

\references{
Reed, D. and R. Miller. 2023. SBC LTER: Reef: Kelp Forest Community Dynamics: Cover of sessile organisms, Uniform Point Contact ver 33. Environmental Data Initiative. https://doi.org/10.6073/pasta/0af1a5b0d9dde5b4e5915c0012ccf99c (Accessed: 2023-12-01).
}