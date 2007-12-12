"phylo.overlap" <-
function(assemblage.one, assemblage.two, cm, method=c("shared","unique")){

        # returns the TBL PD score for a subset of species defined
        # either as occuring in assemblages 1 AND 2 or 1 NOT 2.
        # Note that phylo.overlap(1, 2, cm, "shared") is equivalent
        # to phylo.overlap(2, 1, cm, "shared") but that
        # phylo.overlap(1, 2, cm, "unique") is not the same as
        # phylo.overlap(2, 1, cm, "unique").

        method <- match.arg(method)
        switch(method,
                shared = {taxon.subset <-  assemblage.one[assemblage.one %in% assemblage.two]},
                unique = {taxon.subset <-  assemblage.one[! assemblage.one %in% assemblage.two]})
        if(length(taxon.subset) > 0){
		pd <- pd.calc(cm, taxon.subset)$pd
        } else {
                pd <- 0
        }
	return(pd)
}

