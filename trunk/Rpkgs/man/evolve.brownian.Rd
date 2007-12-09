\name{evolve.brownian}
\alias{evolve.brownian}
\title{ Trait evolution under Brownian motion }
\description{
  Trait evolution under Brownian motion model
}
\usage{
evolve.brownian(phy, value = 0, var = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{ phylo object }
  \item{value}{ mean }
  \item{var}{ variance }
}
\value{
  Data frame of trait values with row.names corresponding to phylo$tip.label
}
\author{ Steve Kembel <skembel@berkeley.edu }
