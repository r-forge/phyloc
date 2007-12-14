plotNodeNumbers.phylo <- function(phy)
{
	phy$node.label <- (length(phy$tip.label)+1):max(phy$edge);
	phy$tip.label <- paste(1:length(phy$tip.label), phy$tip.label, sep='_');
	
	plot.phylo(phy, show.node.label=TRUE);	
}
