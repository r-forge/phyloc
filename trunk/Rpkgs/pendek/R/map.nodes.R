"map.nodes" <-
function(bi.tree, poly.tree){
# Written by Andy Purvis
# Maps nodes in a binary tree produced by force.binary.tree onto the polytomous parent tree.
# Also produces the weights later to be used in the linear regression

number.of.binary.nodes<-length(bi.tree$tip.label)-1
nnodes.in.poly.tree<-length(poly.tree$edge[,1])-length(poly.tree$tip.label)+1

maps.from<-1 #Initialise node index for binary tree
maps.onto<-array(NA,c(number.of.binary.nodes,4))

for (i in 1:nnodes.in.poly.tree)
{
	poly.node.name<-as.character(-i)
	n.desc<-length(which(poly.tree$edge[,1]==poly.node.name))
	for (j in 1:(n.desc-1))
	{
		maps.onto[maps.from,1]<-poly.node.name
		maps.onto[maps.from,3]<-1/(n.desc-1) #weights for use in regression later
		maps.onto[maps.from,4]<-ifelse(j==1,TRUE,FALSE) #whether clade is in both trees
		maps.from<-maps.from+1
	}	
}

maps.onto[,2]<-as.character(-seq(1:number.of.binary.nodes))

maps.onto
}

