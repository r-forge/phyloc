"pandit" <-
function (tree, class="phylo", type="s", quiet=FALSE, model=NULL, p=0.3) {
	
	if (class!="treeshape" & class!="phylo") {
		stop("'class' must be one of 'phylo' or 'treeshape'\n")
	}

	if (class(tree)=='list'){
		tree=unlist(tree)
	}
	err=options()$show.error.messages
	warn=options()$warn
	options(warn=-1)
	options(show.error.messages=FALSE)
		
      if(inherits(try(open(url("http://www.google.com")), silent = TRUE),
             "try-error")) {
		cat("No connexion available\n")
		return(NULL)
      }	
	options(warn=warn)
	options(show.error.messages=err)
	
	
	number=0
	res=list()
#Recuperation des arbres :
	if (!quiet) {cat("Connecting to Pandit...")}
		
	
	for (i in 1:length(tree)) {
		URL="http://www.sanger.ac.uk/cgi-bin/Pfam/download_hmm.pl?acc=PF"
		aux=as.character(tree[i])
		if (nchar(aux)==1){ aux=paste("0000", aux, sep="")}
		if (nchar(aux)==2){ aux=paste("000", aux, sep="")}
		if (nchar(aux)==3){ aux=paste("00", aux, sep="")}
		if (nchar(aux)==4){ aux=paste("0", aux, sep="")}
		URL = paste(URL, aux, sep="")
		
		if (type=="s") { URL = paste(URL,"&tree=tree&type=seed", sep="") }
		if (type=="e") { URL = paste(URL,"&tree=tree&type=full", sep="") }
		
		text=""

            err=options()$show.error.messages
		warn=options()$warn
		options(warn=-1)
		options(show.error.messages=FALSE)
		text=try(scan(file=URL, what="", quiet=TRUE, comment.char=""), silent=TRUE)
		options(warn=warn)
		options(show.error.messages=err)
		
 		
		if (length(text)!=1) {
			text[1]="("
			#text=text[-length(text)]
			tmp<-""
			for (i in 1:length(text)) {
				tmp<-paste(tmp, text[i], sep="")
			}
			text <- tmp
			#return(text)
			phy=read.tree2(text=text)
                        if (is.null(phy)) {
				return(NULL)
			}
			if (class=="treeshape") {
				tmp=as.treeshape.phylo(phy, model, p)
				if (identical(tmp, NULL)==FALSE) {
					number=number+1
					res[[number]]=tmp
				}
			} 
			if (class=="phylo") {
				number=number+1
				res[[number]]=phy
			}
		}
	}
	
	#Transformation des arbres
	if (number==0){
		if (!quiet) {
			cat("\nno object of class \"", class)
			cat(" \"\n")
		}
		return(NULL)
	}
	
	if (number==1){
		if (!quiet) {
			cat("\n 1 object of class \"", class)
			cat(" \"\n")
		}
		return(res[[1]])
	}
	

	if (!quiet) {
		cat("\n", number)
		cat(" objects of class \"", class)
		cat(" \"\n")
	}
	return(res)

}

