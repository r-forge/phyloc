"all.r2tbl" <-
function (phylogeny) 
{
    r2t.bl <- rep(NA, length(phylogeny$tip.label))
    for (i in 1:length(r2t.bl)) {
        n.edge <- which(phylogeny$edge[, 2] == as.character(i))
        sum.so.far <- phylogeny$edge.length[n.edge]
        while (phylogeny$edge[n.edge, 1] != "-1") {
            n.edge <- which(phylogeny$edge[, 2] == phylogeny$edge[n.edge, 
                1])
            sum.so.far <- sum.so.far + phylogeny$edge.length[n.edge]
        }
        r2t.bl[i] <- sum.so.far
    }
    return(r2t.bl)
}

