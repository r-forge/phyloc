\name{writetraits}
\alias{writetraits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Write a Phylocom traits formatted file}
\description{
  Write a Phylocom traits formatted file
}
\usage{
writetraits(trt, file = "", bin = NULL, sigd = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trt}{ Data frame containing trait data }
  \item{file}{ Filename path }
  \item{bin}{ Vector index of trait columns to be treated as binary  }
  \item{sigd}{ Significant digits for output }
}

\references{ Webb et al. 2007. Phylocom. }
\author{ David Ackerly <dackerly@berkeley.edu> and Steve Kembel <skembel@berkeley.edu> }
