\name{eSpider}
\docType{data}
\alias{spider}
\title{Hunting spider data}
\description{
Extended dataset of counts of hunting spiders in a dune area in the Netherlands, measured at 100 pitfall traps.
}
\usage{
data(spider)
}
\format{
\describe{
\item{abund}{A data frame with abundances of 12 hunting spider species measured at 100 sites.
}
\item{X}{A matrix of 26 predictor variables at 28 of the 100 sites.
}
\item{nonNA}{An vector of indices indicating at which sites the predictor variables were measured.
}
}

}}
\description{
Counts of hunting spiders in a dune area in the Netherlands, measured with 100 different pitfall traps. This dataset was published with permission from the CANOCO FORTRAN package (version 4 or higher) example datasets.

Species names have been abbreviated, corresponding to: 
Alopacce = \emph{Alopecosa accentuata}, Alopcune = \emph{Alopecosa cuneata}, Alopfabr = \emph{Alopecosa fabrilis}, Arctlute = \emph{Arctosa lutetiana}, Arctperi = \emph{Arctosa perita}, Auloalbi = \emph{Alonia albimana}, Pardlugu = \emph{Pardosa lugubris}, Pardmont = \emph{Pardosa monticola}, Pardnigr = \emph{Pardosa nigriceps}, Pardpull = \emph{Pardosa pullata}, Trocterr = \emph{Trochosa terricola}, Zoraspin = \emph{Zora spinimana}.

Enviromental measurements were taken at 28 of the 100 pitfall traps measuring soil properties(Water content: "conWate", Humus content: "conHumu", Acidity (pH-KCl)), vegetation (percentage bare sand, cover on the ground by leaves and twigs, cover by mosses and lichens, cover by the herb and grass layer (including maximum height, minimum height, cover by \emph{Calamagrostis epigejos}, cover by \emph{Carex arenaria}, cover by \emph{Festuca ovina}, cover by \emph{Corynephorus canescens}, cover by \emph{Urtica dioica}, cover by \emph{Moehringia trinervia}), cover by the shrub layer (minimum and maximum height, and cover by \emph{Ligustrum vulgare}), cover by the tree layer (including maximum height, cover by \emph{Populus tremula} or \emph{Crataegus monogyna})), and light properties (lux at equal grey sky, lux at cloudless sky, lux by reflection of the soil surface)
}
\references{
Smeenk-Enserink, N., & Van der Aart, P. J. M. (1974). Correlations between distributions of hunting spiders (\emph{Lycosidae}, \emph{Ctenidae}) and environmental characteristics in a dune area. Netherlands Journal of Zoology, 25(1), 1-45.
}