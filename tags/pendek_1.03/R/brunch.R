"brunch" <-
function(phylogeny, dataset, y, x, cherries.only=FALSE, arbitrary.bl=0){

#Generates BRUNCH-like contrasts in vars from a phylogeny, which can
#contain polytomies, but only when x shows variation.

variable.list<-c(y,x)
dataset<-strip.spp.not.in.tree(dataset,phylogeny)
dataset<-strip.missing.data(dataset,variable.list)
phylogeny<-strip.phylo(phylogeny,dataset)
y.data<-dataset[,names(dataset)==y]
names(y.data)<-dataset[,1]
x.data<-dataset[,names(dataset)==x]
names(x.data)<-dataset[,1]

plum.list<-find.plums(phylogeny)

#Generate contrasts and add each in turn to list; adjust branch lengths and data if necessary
for (i in 1:length(plum.list))
{
	y.plum<-as.numeric(y.data[is.element(names(y.data),plum.list[[i]]$tips)])
	x.plum<-as.numeric(x.data[is.element(names(x.data),plum.list[[i]]$tips)])
	plum.list[[i]]<-c(plum.list[[i]],list(y.plum),list(x.plum))
	names(plum.list[[i]])[3]<-y
	names(plum.list[[i]])[4]<-x
	contrast<-brunch.contrast(phylogeny, plum.list[[i]],arbitrary.bl=arbitrary.bl)
	plum.list[[i]]<-c(plum.list[[i]],contrast)
	
	#Now drop tips from data set
	y.data<-y.data[-match(plum.list[[i]]$to.drop, names(y.data))]
	x.data<-x.data[-match(plum.list[[i]]$to.drop, names(x.data))]
	
	#Amend branch length if appropriate
	if (!is.null(plum.list[[i]]$bl.adjustment))
	{
		tip.number<-which(phylogeny$tip.label==plum.list[[i]]$to.drop[1])
		edge.number<-which(phylogeny$edge[,2]==as.character(tip.number))
		phylogeny$edge.length[edge.number]<-plum.list[[i]]$bl.adjustment	
	}
	
	#Put new tip name and new data into phylogeny and data set if appropriate
	if (!is.null(plum.list[[i]]$nodal.y))
	{
		new.name<-paste(plum.list[[i]]$to.drop[1],"+",sep="")
		phylogeny$tip.label[which(phylogeny$tip.label==plum.list[[i]]$to.drop[1])]<-new.name
		y.data<-c(y.data,plum.list[[i]]$nodal.y)
		names(y.data)[length(y.data)]<-new.name
		x.data<-c(x.data,plum.list[[i]]$nodal.x)
		names(x.data)[length(x.data)]<-new.name
		plum.list[[i]]$to.drop<-plum.list[[i]]$to.drop[-1]
	}

}

#Cycle through dropping tips if necessary
to.drop<-unlist(lapply(plum.list,function(x){x$to.drop}))
phylogeny<-try(drop.tip(phylogeny,to.drop),silent=TRUE)


#If not cherries only, then iterate until no variance or no data left in x
if (cherries.only==FALSE)
{
	x.var.left<-ifelse(length(unique(x.data))>1,TRUE,FALSE)
	while (x.var.left==TRUE)
	{
		more.plums<-find.plums(phylogeny)
		#Generate contrasts and add each in turn to list; adjust branch lengths and data if necessary
		for (i in 1:length(more.plums))
		{
			y.plum<-as.numeric(y.data[is.element(names(y.data),more.plums[[i]]$tips)])
			x.plum<-as.numeric(x.data[is.element(names(x.data),more.plums[[i]]$tips)])
			more.plums[[i]]<-c(more.plums[[i]],list(y.plum),list(x.plum))
			names(more.plums[[i]])[3]<-y
			names(more.plums[[i]])[4]<-x
			contrast<-brunch.contrast(phylogeny, more.plums[[i]],arbitrary.bl=arbitrary.bl)
			more.plums[[i]]<-c(more.plums[[i]],contrast)
			
			#Now drop tips from data set
			y.data<-y.data[-match(more.plums[[i]]$to.drop, names(y.data))]
			x.data<-x.data[-match(more.plums[[i]]$to.drop, names(x.data))]
	
			#Amend branch length if appropriate
			if (!is.null(more.plums[[i]]$bl.adjustment))
			{
				tip.number<-which(phylogeny$tip.label==more.plums[[i]]$to.drop[1])
				edge.number<-which(phylogeny$edge[,2]==as.character(tip.number))
				phylogeny$edge.length[edge.number]<-more.plums[[i]]$bl.adjustment	
			}
	
			#Put new tip name and new data into phylogeny and data set if appropriate
			if (!is.null(more.plums[[i]]$nodal.y))
			{
				new.name<-paste(more.plums[[i]]$to.drop[1],"+",sep="")
				phylogeny$tip.label[which(phylogeny$tip.label==more.plums[[i]]$to.drop[1])]<-new.name
				y.data<-c(y.data,more.plums[[i]]$nodal.y)
				names(y.data)[length(y.data)]<-new.name
				x.data<-c(x.data,more.plums[[i]]$nodal.x)
				names(x.data)[length(x.data)]<-new.name
				more.plums[[i]]$to.drop<-more.plums[[i]]$to.drop[-1]
			}
			
		}
		
		#Cycle through dropping tips if necessary
		to.drop<-unlist(lapply(more.plums,function(x){x$to.drop}))
		phylogeny<-try(drop.tip(phylogeny,to.drop),silent=TRUE)

						
		plum.list<-c(plum.list,more.plums)
		x.var.left<-ifelse(length(unique(x.data))>1,TRUE,FALSE)

	}
	
}
print(paste(length(plum.list),"contrasts generated in all (N.B. Not all need be informative)"))
plum.list
}

