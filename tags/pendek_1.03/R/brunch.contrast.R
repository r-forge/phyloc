"brunch.contrast" <-
function(phylogeny, plum, arbitrary.bl=0){

#Computes a BRUNCH contrast if there is one.

branch.lengths<-rep(NA,length(plum$tips))
for (i in 1:length(plum$tips)) branch.lengths[i]<-get.branch.length(phylogeny, plum$tips[i])

if (var(plum[[4]])==0)
{
	#No variance in X within node, so compute average
	root.y<-sum((plum[[3]]/branch.lengths))/sum(1/branch.lengths)
	root.x<-sum((plum[[4]]/branch.lengths))/sum(1/branch.lengths)
	root.bl.adjustment<-1/sum(1/branch.lengths)
	y.contrast<-NA
	x.contrast<-NA
	variance<-NA
	to.drop<-plum$tips
}

if (var(plum[[4]])!=0)
{
	#Variance in X within node, so compute contrast
	contrast.weights<-rep(NA,length(plum[[3]])) #construct vector to hold contrast weights
	mean.x<-mean(plum[[4]])
	contrast.weights<-ifelse((plum[[4]]>=mean(plum[[4]])),(1/sum(plum[[4]]>=mean.x)),(-1/sum(plum[[4]]<mean.x)))

	if (length(plum$tips)==2)
	{
	#Contrast at a bifurcation
	root.y<-NULL
	root.x<-NULL
	root.bl.adjustment<-NULL
	y.contrast<-sum(plum[[3]]*contrast.weights)
	x.contrast<-sum(plum[[4]]*contrast.weights)
	variance<-sum(branch.lengths)
	to.drop<-plum$tips
	}
	
	if (length(plum$tips)>2)
	{
	#Contrast at a polytomy
	root.y<-NULL
	root.x<-NULL
	root.bl.adjustment<-NULL
	va<-arbitrary.bl
	vb<-arbitrary.bl
	branch.lengths<-branch.lengths-arbitrary.bl
	a<-plum[[4]]>=mean(plum[[4]])
	ya<-sum((plum[[3]][a]/branch.lengths[a]))/sum(1/branch.lengths[a])
	xa<-sum((plum[[4]][a]/branch.lengths[a]))/sum(1/branch.lengths[a])
	yb<-sum((plum[[3]][!a]/branch.lengths[!a]))/sum(1/branch.lengths[!a])
	xb<-sum((plum[[4]][!a]/branch.lengths[!a]))/sum(1/branch.lengths[!a])
	va<-va+(1/sum(1/branch.lengths[a]))
	vb<-vb+(1/sum(1/branch.lengths[!a]))
	
	y.contrast<-ya-yb
	x.contrast<-xa-xb
	variance<-va+vb
	to.drop<-plum$tips
	}
}

contrast<-list(nodal.y=root.y, nodal.x=root.x, bl.adjustment=root.bl.adjustment, y.contrast=y.contrast, x.contrast=x.contrast,
	variance=variance, to.drop=to.drop)
	
contrast
}

