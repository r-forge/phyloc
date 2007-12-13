
## not bothering to check for zero branch lengths:
##   consensus is that this isn't very important,
##  and that it's simple enough to do
##   any(EdgeLength(x)==0) if necessary
hasPoly <- function(object) {
  degree <- tabulate(edges(object)[, 1])
  struc <- any(degree > 2)
  return(struc)
}

hasSingles <- function(object) {
  degree <- tabulate(edges(object)[, 1])
  any(degree == 1)
}
