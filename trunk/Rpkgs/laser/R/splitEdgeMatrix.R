splitEdgeMatrix <- function(phy, node)
{
  x <- branching.times(phy); #necessary later in cbind z section
  rootnode <- length(phy$tip.label) + 1 ;
  phy$tag <-rep(1, nrow(phy$edge))  #to tag subtrees

  if (node >= rootnode) #tags all descendants of a particular internal node
  {
    node.desc <- node
    pos<-1
    phy$tag[phy$edge[,1] == node.desc[1]]<-2
    while(pos != (length(node.desc)+1))
    {
      temp <- node.sons(phy, node.desc[pos])
      temp<-temp[temp > rootnode];

      for (k in 1:length(temp)){
        phy$tag[phy$edge[,1]== temp[k]] <- 2
      }
      node.desc<-c(node.desc, temp)
      pos<-pos+1
    }
  }
  else if (node > 0)    #tags terminal nodes only.
    phy$tag[phy$edge[,2]==node] <- 2
    
  #at this point, all branches are tagged by subtree 2 or 1 (defined
  # by removing ONLY branches in 2
  
  z <- cbind(phy$edge, gsr(phy$edge[,1], names(x), x), phy$edge.length, phy$phenotype, phy$tag)
  z<- matrix(as.numeric(z), dim(z))
  
  #elements of z are as follows:
  #z[,1:2]<-edge; z[,3]<-branching_time; z[,4] <- branch.length; z[,5]<-phenotype, z[,6]<-tag.
  z<-as.data.frame(z)
  return(z)
}
