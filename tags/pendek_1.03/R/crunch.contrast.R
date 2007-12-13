"crunch.contrast" <-
function(phylogeny, plum){

#Computes a CRUNCH contrast.

branch.lengths<-rep(NA,length(plum$tips))
for (i in 1:length(plum$tips)) branch.lengths[i]<-get.branch.length(phylogeny, plum$tips[i])


contrast.weights<-rep(NA,length(plum[[3]])) #construct vector to hold contrast weights
mean.x<-mean(plum[[4]])
contrast.weights<-ifelse((plum[[4]]>=mean(plum[[4]])),(1/sum(plum[[4]]>=mean.x)),(-1/sum(plum[[4]]<mean.x)))

if (length(plum$tips)==2)
{
	#Contrast at a bifurcation
	root.y<-sum((plum[[3]]/branch.lengths))/sum(1/branch.lengths)
	root.x<-sum((plum[[4]]/branch.lengths))/sum(1/branch.lengths)
	root.bl.adjustment<-1/sum(1/branch.lengths)
	y.contrast<-sum(plum[[3]]*contrast.weights)
	x.contrast<-sum(plum[[4]]*contrast.weights)
	variance<-sum(branch.lengths)
	to.drop<-plum$tips
}
	
if (length(plum$tips)>2)
	{
	#Contrast at a polytomy
	#At present, does nothing different from bifurcations - this needs more work!!!
	root.y<-sum((plum[[3]]/branch.lengths))/sum(1/branch.lengths)
	root.x<-sum((plum[[4]]/branch.lengths))/sum(1/branch.lengths)
	root.bl.adjustment<-1/sum(1/branch.lengths)
	y.contrast<-sum(plum[[3]]*contrast.weights)
	x.contrast<-sum(plum[[4]]*contrast.weights)
	variance<-sum(branch.lengths)
	to.drop<-plum$tips
	}


contrast<-list(nodal.y=root.y, nodal.x=root.x, bl.adjustment=root.bl.adjustment, y.contrast=y.contrast, x.contrast=x.contrast,
	variance=variance, to.drop=to.drop)
	
contrast
}

