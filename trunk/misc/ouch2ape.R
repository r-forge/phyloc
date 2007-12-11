ouch2ape<-function(node, ancestor, time, species)
{
	node[is.na(species)]<--(1:sum(is.na(species)))
	node[!is.na(species)]<-1:sum(!is.na(species))
	new.ancestor<-node[as.numeric(ancestor)]
	edge<-cbind(new.ancestor[-1], node[-1])
	edge.length<-numeric(length(time)-1)
	for(i in 1:length(edge.length))
		edge.length[i]<-time[i+1]-time[as.numeric(ancestor[i+1])];
	tip.label<-species[!is.na(species)]
	mode(edge) <- "character"
    mode(tip.label) <- "character"
    obj <- list(edge = edge, edge.length = edge.length, tip.label=tip.label)
    class(obj) <- "phylo"
    obj
}
