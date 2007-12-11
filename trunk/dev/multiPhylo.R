"[.multiPhylo" <- function(x, i)
{
    class(x) <- NULL
    structure(x[i], class = "multiPhylo")
}
