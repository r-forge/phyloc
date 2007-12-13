zoom.new<-function(x, wait=FALSE)

{
sub<-subtrees(x, wait=wait)
y<-NULL

repeat {

split.screen(c(1,2))
screen(2)
if (is.null(y)) plot(x)
    else plot(y,sub=paste("Node :", click))
screen(1)
plot(x,sub="Complete tree")


 N.tip<-Ntip(x)
 N.node<-Nnode(x)

 coor<-plot.phylo.coor(x)
 tips<-x$tip.label
 nodes<-x$node.label
 if (is.null(x$node.label)) nodes<-(N.tip+1):(N.tip+N.node)
 labs<-c(rep("",N.tip), nodes)



  click<-identify(coor[,1], coor[,2], labels=labs, n=1)
  #plot.window(xlim=c(0,10), ylim=c(0,10))
  #loc<-locator(n=1)
  	 if (click > N.tip)
 		 {
        close.screen(c(1,2),all.screens = TRUE)
        split.screen(c(1,2))
        screen(1) #selects the screen to plot in
        plot(x, sub="Complete tree") # plots x in screen 1 (left)
        screen(2)
 		for (i in 1:length(sub)) if (sub[[i]]$name==click) break		
        y<-sub[[i]]
        #plot(y,sub=paste("Node :", click))
        ##here starts the function to create the clickable points for choices
        #coor.y<-plot.phylo.coor(y)
        #N.tip.y<-Ntip(y)    
	   
 		 }
 	 else cat("this is a tip, you have to choose a node\n")


 }
}



