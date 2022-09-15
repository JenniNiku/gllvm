\name{eSpider}
\docType{data}
\alias{eSpider}
\title{Hunting spider data}
\description{
Extended dataset of counts of hunting spiders in a dune area in the Netherlands, measured at 100 pitfall traps.
}
\usage{
data(eSpider)
}
\format{
\describe{
\item{abund}{A data frame with abundances of 12 hunting spider species measured at 100 sites.}
\item{X}{A matrix of 26 predictor variables at 28 of the 100 sites.}
\item{nonNA}{An vector of indices indicating at which sites the predictor variables were measured.}
}

}
\details{
Counts of hunting spiders in a dune area in the Netherlands, measured with 100 different pitfall traps. This dataset was published with permission from the CANOCO FORTRAN package (version 4 or higher) example datasets.

Species names have been abbreviated, corresponding to: 
Alopacce = \emph{Alopecosa accentuata}, Alopcune = \emph{Alopecosa cuneata}, Alopfabr = \emph{Alopecosa fabrilis}, Arctlute = \emph{Arctosa lutetiana}, Arctperi = \emph{Arctosa perita}, Auloalbi = \emph{Alonia albimana}, Pardlugu = \emph{Pardosa lugubris}, Pardmont = \emph{Pardosa monticola}, Pardnigr = \emph{Pardosa nigriceps}, Pardpull = \emph{Pardosa pullata}, Trocterr = \emph{Trochosa terricola}, Zoraspin = \emph{Zora spinimana}.

Enviromental measurements were taken at 28 of the 100 pitfall traps measuring soil properties(Water content: "conWate", Humus content: "conHumu", Acidity (pH-KCl)), vegetation ("BareSand": percentage bare sand, "FallTwig": cover on the ground by leaves and twigs, "CovMoss": cover by mosses and lichens, "CovHerb": cover by the herb and grass layer (including maximum height, minimum height, "CovCala": cover by \emph{Calamagrostis epigejos}, cover by \emph{Carex arenaria}, "CovFest": cover by \emph{Festuca ovina}, "CovCory": cover by \emph{Corynephorus canescens}, "CovUrti": cover by \emph{Urtica dioica}, "CovMoeh": cover by \emph{Moehringia trinervia}), "CovShru": cover by the shrub layer (minimum and maximum height, and "CovLigu": cover by \emph{Ligustrum vulgare}), "CovTree": cover by the tree layer (including maximum height, cover by \emph{Populus tremula} or \emph{Crataegus monogyna})), and light properties ("LuxGrey": lux at equal grey sky, "LuxSun": lux at cloudless sky, "LuxRef": lux by reflection of the soil surface)

The original publication of Canonical Correspondence Analysis used standardized versions of the \code{log} and \code{log1p} transformed predictors "ConWate", "BareSand", "FallTwig",  "CovMoss", "CovHerb", "LuxRef".
}

\references{
ter Braak, C.J.F. and Smilauer, P.  (1998).  CANOCO reference manual and user's guide to CANOCO for Windows: software for canonical community ordination (version 4). Microcomputer Power, New York, New York, USA.
ter Braak, C.J.F. (1986). Canonical correspondence analysis: a new eigenvector technique for multivariate direct gradient analysis. Ecology, 67(5), 1167-1179.
Van der Aart, P. J. M. and Smeenk-Enserink, N. (1975). Correlations between distributions of hunting spiders (\emph{Lycosidae}, \emph{Ctenidae}) and environmental characteristics in a dune area. Netherlands Journal of Zoology, 25(1), 1-45.
}
\examples{
data(eSpider)
Y <- eSpider$abund[eSpider$nonNA, ]
X <- eSpider$X[eSpider$nonNA, ]
model <- gllvm(y = Y, X = X, 
    lv.formula = ~log(ConWate) + log1p(BareSand) + log1p(FallTwig) + 
    log1p(CovMoss) + log1p(CovHerb) + log(LuxRef),
    num.RR = 2, 
    family = "poisson")
}