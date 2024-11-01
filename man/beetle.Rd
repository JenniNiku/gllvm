\name{beetle}
\docType{data}
\alias{beetle}
\title{ground beetle assemblages}
\description{
These data describe the abundance of ground beetle assemblages from the Scottish landscape alongside with the environmental data and species traits. The data includes abundances of 68 species of ground beetle species.
}
\usage{
data(beetle)
}
\format{
\describe{
\item{Y}{A data frame of species composition of the ground-beetle assemblages.}
\item{X}{A data frame of study design variables and environmental data.}
\describe{
\item{SiteCode}{Unique id for sample.} 
\item{Landuse}{Land use type.}
\item{Grid}{Grid where sample were collected.}
\item{Area}{sampling area in Scotland.}
\item{Samplingyear}{Sampling year}
\item{Texture}{1, peat; 2, peaty loam; 3, loamy sand; 4, sandy loam; 5, sandy clay loam; 6, sandy silt loam; 7, silty clay}
\item{Org}{organic content (\% loss of organic content on ignition), log10 transformed}
\item{pH}{soil pH}
\item{AvailP}{available P (mg/L), log10 transformed}
\item{AvailK}{available K (mg/L)}
\item{Moist}{ percentage moisture content}
\item{Bare}{ percentage cover estimate of bare ground in 11 1-m2 quadrats, arcsine transformed}
\item{Litter}{ percentage cover estimate of litter cover in 11 1-m2 quadrats, log10 transformed}
\item{Bryophyte}{ percentage cover estimate of bryophytes in 11 1-m2 quadrats, arcsine transformed}
\item{Plants.m2}{ number of reproductive stems (ﬂowering or fruiting) in 11 1-m2 quadrats}
\item{Canopy}{ height canopy height (cm) in 11 1-m2 quadrats}
\item{Stemdensity}{ number of stems (ramets) in 100 cm2}
\item{Biom_l5}{ dry mass (g) of biomass 0–5 cm from soil surface in 400 cm2}
\item{Biom_m5}{ dry mass (g) of biomass .5 cm from soil surface in 400 cm2, log10 transformed}
\item{Reprobiom}{ biomass of reproductive parts (ﬂowers and fruits) in 100 cm2, log10 transformed}
\item{Elevation}{ elevation in m a.s.l.}
\item{Management}{ management intensity index (see Materials and Methods: Environmental data)}
}
\item{TR}{A data frame of the species names and species morphological and life trait characteristics.}
\describe{
\item{SPECIES}{Species names.} 
\item{CODE}{Species codes corresponding species names in abundance matrix}
\item{LYW}{ diameter of the eye, measured from above}
\item{LAL}{ length of the antenna}
\item{LPW}{ maximum width of the pronotum}
\item{LPH}{ maximum depth (“vaulting”) of the pronotum}
\item{LEW}{ maximum width of the elytra}
\item{LFL}{ length of the metafemur (with the articulation segments), from the coxa to the apex}
\item{LTR}{ length of the metatrochanter}
\item{LRL}{ length of the metatarsi}
\item{LFW}{ maximum width of the metafemur}
\item{LTL}{ total length (length of the pronotum in the medial line plus length of the elytra, from the medial ridge of the scutellum to the apex)}
\item{CLG}{ color of the legs (1, pale; 2, black; 3, metallic)}
\item{CLB}{ color of the body (1, pale; 2, black; 3, metallic)}
\item{WIN}{ wing development (1, apterous or brachypterous; 2, dimorphic; 3, macropterous)}
\item{PRS}{ shape of the pronotum (1, oval; 2, cordiform; 3, trapezoidal)}
\item{OVE}{ overwintering (1, only adults; 2, adults and larvae or only larvae)}
\item{FOA}{ food of the adult (1, mostly Collembola; 2, generalist predator; 3, mostly plant material)}
\item{DAY}{ daily activity (1, only diurnal; 2, nocturnal)}
\item{BRE}{ breeding season (1, spring; 2, summer; 3, autumn or winter)}
\item{EME}{ main period of emergence of the adults (1, spring; 2, summer; 3, autumn)}
\item{ACT}{ main period of adult activity (1, autumn; 2, summer only)}
}
}

}
\details{

Beetles were sampled with two parallel rows of nine pitfall traps (diameter 7.5 cm, 2 m apart) at each site, starting in early May.

Detailed description of the data available in the reference and in the Ecological Archives E082-012.
}

\references{
Ignacio Ribera, Sylvain Dolédec, Iain S. Downie, and Garth N. Foster. 2001. Effect of land disturbance and stress on species traits of ground beetle assemblages. Ecology 82:1112-1129.
}

\examples{
\dontrun{
data(beetle)
# Abundance matrix
Y <- beetle$Y
# Environmental data
X <- beetle$X
# Species traits
TR <- beetle$TR
}
}