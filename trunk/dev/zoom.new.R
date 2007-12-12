 zoom.new<-function(x, wait=FALSE)
 {
click<-"Complete tree"
repeat {

 plot(x, sub=paste("Node :", click))
 N.tip<-Ntip(x)
 N.node<-Nnode(x)
 sub<-subtrees(x, wait=wait)
 coor<-plot.phylo.coor(x)
 tips<-x$tip.label
 nodes<-x$node.label
 if (is.null(x$node.label)) nodes<-(N.tip+1):(N.tip+N.node)
 labs<-c(tips, nodes)
  
 
	
  click<-identify(coor[,1], coor[,2], labels=labs, n=1)
  	 if (click > N.tip)
 		 {
 		 for (i in 1:length(sub)) if (sub[[i]]$name==click) break
 		 x11()
 		 x<-sub[[i]]
 		 }
 	 else cat("this is a tip, you have to choose a node\n")


 }
}


