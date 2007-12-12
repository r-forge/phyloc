subtrees<-function(tree)
{
N.tip<-Ntip (tree)
N.node<-Nnode(tree)
#if (is.unrooted(tree))
  limit<-N.tip+N.node
print("limit")
print(limit)
sub<-list(N.node)
u<-0
  for (k in (N.tip+2):limit)
  {
  u<-u+1
  fils<-NULL
  pere<-res <- k 
	repeat
	{
	for (i in 1: length(pere)) fils<-c(fils, tree$edge[,2][tree$edge[,1]==pere[i]])
	res<-c(res, fils)
      pere<-fils
	fils<-NULL
	if (length(pere)==0) break
	}
  len<-res[res>N.tip] 
  len.tip<-res[res<N.tip+1] 
  vec<-1:length(tree$tip.label)
  len.tip.stay<-setdiff(vec, len.tip)
  tree2<-drop.tip(tree, len.tip.stay)
  sub[[u]]<-tree2

  sub[[u]]$name<-k

  }
return(sub)
}
