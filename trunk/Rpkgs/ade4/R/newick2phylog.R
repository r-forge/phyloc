"newick2phylog" <- function (x.tre, add.tools = TRUE, call =match.call()) {
    complete <- function(x.tre) {
        # Si la cha�ne est en plusieurs morceaux elle est rassembl�e
        if (length(x.tre) > 1) {
            w <- ""
            for (i in 1:length(x.tre)) w <- paste(w, x.tre[i], 
                sep = "")
            x.tre <- w
        }
        # Si les parenth�ses gauches et droites ont des effectifs diff�rents -> out
        ndroite <- nchar(gsub("[^)]","",x.tre))
        ngauche <- nchar(gsub("[^(]","",x.tre))
        if (ndroite !=ngauche) stop (paste (ngauche,"( versus",ndroite,")"))
        # on doit trouver un ;
        if (regexpr(";", x.tre) == -1) 
            stop("';' not found")
        # Tous les commentaires entre [] sont supprim�s
        i <- 0
        kint <- 0
        kext <- 0
        arret <- FALSE
        if (regexpr("\\[", x.tre) != -1) {
            x.tre <- gsub("\\[[^\\[]*\\]", "", x.tre, ext = FALSE)
        }
        x.tre <- gsub(" ", "", x.tre, ext = FALSE)
        # On ne peut supprimer les . qui sont dans les distances !
        # x.tre <- gsub("[.]","_", x.tre, ext = FALSE)
        while (!arret) {
            i <- i + 1
            # examen de la cha�ne par couple de charact�res
            if (substr(x.tre, i, i) == ";") 
                arret <- TRUE
            # (, c'est une feuille sans label
            if (substr(x.tre, i, i + 1) == "(,") {
                kext <- kext + 1
                add <- paste("Ext", kext, sep = "")
                x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre, 
                  i + 1), sep = "")
                i <- i + 1
            }
            # ,, c'est une feuille sans label
            else if (substr(x.tre, i, i + 1) == ",,") {
                kext <- kext + 1
                add <- paste("Ext", kext, sep = "")
                x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre, 
                  i + 1), sep = "")
                i <- i + 1
            }
            # ,) c'est une feuille sans label
            else if (substr(x.tre, i, i + 1) == ",)") {
                kext <- kext + 1
                add <- paste("Ext", kext, sep = "")
                x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre, 
                  i + 1), sep = "")
                i <- i + 1
            }
            # (: c'est une feuille sans label avec distance
            else if (substr(x.tre, i, i + 1) == "(:") {
                kext <- kext + 1
                add <- paste("Ext", kext, sep = "")
                x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre, 
                  i + 1), sep = "")
                i <- i + 1
            }
            # ,: c'est une feuille sans label avec distance
            else if (substr(x.tre, i, i + 1) == ",:") {
                kext <- kext + 1
                add <- paste("Ext", kext, sep = "")
                x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre, 
                  i + 1), sep = "")
                i <- i + 1
            }
            # ), c'est un noeud sans label
            else if (substr(x.tre, i, i + 1) == "),") {
                kint <- kint + 1
                add <- paste("I", kint, sep = "")
                x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre, 
                  i + 1), sep = "")
                i <- i + 1
            }
            # )) c'est un noeud sans label
            else if (substr(x.tre, i, i + 1) == "))") {
                kint <- kint + 1
                add <- paste("I", kint, sep = "")
                x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre, 
                  i + 1), sep = "")
                i <- i + 1
            }
            # ): c'est un noeud sans label avec distance
            else if (substr(x.tre, i, i + 1) == "):") {
                kint <- kint + 1
                add <- paste("I", kint, sep = "")
                x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre, 
                  i + 1), sep = "")
                i <- i + 1
            }
            # ); c'est la racine sans label
            else if (substr(x.tre, i, i + 1) == ");") {
                add <- "Root"
                x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre, 
                  i + 1), sep = "")
                i <- i + 1
            }
        }
        # extraction de l'information non structurelle
        lab.points <- strsplit(x.tre, "[(),;]")[[1]]
        lab.points <- lab.points[lab.points != ""]
        # recherche de la pr�sence des longueurs
        no.long <- (regexpr(":", lab.points) == -1)
        # si il n'y avait aucune longueur
        if (all(no.long)) {
            lab.points <- paste(lab.points, ":", c(rep("1", length(no.long) - 
                1), "0.0"), sep = "")
        }
        # si il y en vait partout sauf � la racine
        else if (no.long[length(no.long)]) {
            lab.points[length(lab.points)] <- paste(lab.points[length(lab.points)], 
                ":0.0", sep = "")
        }
        # si il y en a et il n'y en a pas -> out
        else if (any(no.long)) {
            print(x.tre)
            stop("Non convenient data leaves or nodes with and without length")
        }
        w <- strsplit(x.tre, "[(),;]")[[1]]
        w <- w[w != ""]
        leurre <- make.names(w, unique = TRUE)
        leurre <- gsub("[.]","_", leurre, ext = FALSE)
        for (i in 1:length(w)) {
            old <- paste(w[i])
            x.tre <- sub(old, leurre[i], x.tre, ext = FALSE)
        }
        # extraction des labels et des longueurs
        w <- strsplit(lab.points, ":")
        label <- function(x) {
            # ici on peut travailler sur les labels
            lab <- x[1]
            lab <- gsub("[.]","_", lab, ext = FALSE)
            return (lab)
        }
        
        longueur <- function(x) {
            long <- x[2]
            return (long)
        }

        labels <- unlist(lapply(w, label))
        longueurs <- unlist(lapply(w, longueur))
        # ici on peut travailler sur les labels
        labels <- make.names(labels, TRUE)
        labels <- gsub("[.]","_", labels, ext = FALSE)
        w <- labels
        for (i in 1:length(w)) {
            new <- w[i]
            x.tre <- sub(leurre[i], new, x.tre, ext = FALSE)
        }
        # on les a remis � leur place
        cat <- rep("", length(w))
        for (i in 1:length(w)) {
            new <- w[i]
            if (regexpr(paste(")", new, sep = ""), x.tre, ext = FALSE) != 
                -1) 
                cat[i] <- "int"
            else if (regexpr(paste(",", new, sep = ""), x.tre, 
                ext = FALSE) != -1) 
                cat[i] <- "ext"
            else if (regexpr(paste("(", new, sep = ""), x.tre, 
                ext = FALSE) != -1) 
                cat[i] <- "ext"
            else cat[i] <- "unknown"
        }
        return(list(tre = x.tre, noms = labels, poi = as.numeric(longueurs), 
            cat = cat))
    }
    res <- complete(x.tre)
    poi <- res$poi
    nam <- res$noms
    names(poi) <- nam
    cat <- res$cat
    res <- list(tre = res$tre)
    res$leaves <- poi[cat == "ext"]
    names(res$leaves) <- nam[cat == "ext"]
    res$nodes <- poi[cat == "int"]
    names(res$nodes) <- nam[cat == "int"]
    listclass <- list()
    dnext <- c(names(res$leaves), names(res$nodes))
    listpath <- as.list(dnext)
    names(listpath) <- dnext
    x.tre <- res$tre
    while (regexpr("[(]", x.tre) != -1) {
        a <- regexpr("([^()]*)", x.tre, ext = FALSE)
        n1 <- a[1] + 1
        n2 <- n1 - 3 + attr(a, "match.length")
        chasans <- substring(x.tre, n1, n2)
        chaavec <- paste("(", chasans, ")", sep = "")
        nam <- unlist(strsplit(chasans, ","))
        w1 <- strsplit(x.tre, chaavec, ext = FALSE)[[1]][2]
        parent <- unlist(strsplit(w1, "[,);]", ext = FALSE))[1]
        listclass[[parent]] <- nam
        x.tre <- gsub(chaavec, "", x.tre, ext = FALSE)
        w2 <- which(unlist(lapply(listpath, function(x) any(x[1] == 
            nam))))
        for (i in w2) {
            listpath[[i]] <- c(parent, listpath[[i]])
        }
    }
    res$parts <- listclass
    res$paths <- listpath
    dnext <- c(res$leaves, res$nodes)
    names(dnext) <- c(names(res$leaves), names(res$nodes))
    res$droot <- unlist(lapply(res$paths, function(x) sum(dnext[x])))
    res$call <- call
    class(res) <- "phylog"
    if (!add.tools) return(res)
    return(newick2phylog.addtools(res))
}


"hclust2phylog" <- function (hc, add.tools = TRUE) {
    if (!inherits(hc, "hclust")) 
        stop("'hclust' object expected")
    labels.leaves <- make.names(hc$labels, TRUE)
    nnodes <- nrow(hc$merge)
    labels.nodes <- paste("Int", 1:nnodes, sep = "")
    l.bra <- matrix("$", nnodes, 2)
    for (i in nnodes:1) {
        for (j in 1:2) {
            if (hc$merge[i, j] < 0) 
                l.bra[i, j] <- as.character(hc$height[i])
            else l.bra[i, j] <- as.character(hc$height[i] - hc$height[hc$merge[i, 
                j]])
        }
    }
    l.eti <- matrix("$", nnodes, 2)
    for (i in nnodes:1) {
        for (j in 1:2) {
            if (hc$merge[i, j] > 0) 
                l.eti[i, j] <- labels.nodes[hc$merge[i, j]]
            else l.eti[i, j] <- labels.leaves[-hc$merge[i, j]]
        }
    }
    tre <- paste("(", l.eti[nnodes, 1], ":", l.bra[nnodes, 1], 
        ",", l.eti[nnodes, 2], ":", l.bra[nnodes, 2], ")Root:0.0;", 
        sep = "")
    for (j in (nnodes - 1):1) {
        w <- paste("(", l.eti[j, 1], ":", l.bra[j, 1], ",", l.eti[j, 
            2], ":", l.bra[j, 2], ")", labels.nodes[j], ":", 
            sep = "")
        tre <- gsub(paste(labels.nodes[j], ":", sep = ""), w, 
            tre, ext = FALSE)
    }
    res <- newick2phylog(tre, add.tools, call=match.call())
    return(res)
}


taxo2phylog <- function (taxo, add.tools = FALSE, root = "Root", abbrev = TRUE) 
{
    if (!inherits(taxo, "taxo")) 
        stop("Object 'taxo' expected")
    nc <- ncol(taxo)
    for (k in 1:nc) { 
        w <- as.character(k)
        w <- paste("l", w, sep="")
        w1 <- levels(taxo[,k])
        if (abbrev) w1 <- abbreviate(w1)
        levels(taxo[,k]) <- paste(w, w1,sep="")
    }
    leaves.names <- row.names(taxo)
    res <- paste(root,";",sep="")
    x <- taxo[, nc]
    xred <- as.character(levels(x))
    w <- "("
    for (i in xred) w <- paste(w, i, ",", sep = "")
    res <- paste(w, ")", res, sep = "")
    res <- sub(",)", ")", res, ext = FALSE)
    for (j in nc:1) {
        x <- taxo[, j]
        if (j>1) y <- taxo[, j - 1] else y <- as.factor(leaves.names)
        for (k in 1:nlevels(x)) {
            w <- "("
            old <- as.character(levels(x)[k])
            yred <- unique(y[x == levels(x)[k]])
            yred <- as.character(yred)
            for (i in yred) w <- paste(w, i, ",", sep = "")
            w <- paste(w, ")", old, sep = "")
            w <- sub(",)", ")", w, ext = FALSE)
        res <- gsub(old, w, res, ext = FALSE)
        }
    }
    return(newick2phylog(res, add.tools, call = match.call()))
}

   
"newick2phylog.addtools" <- function(res, tol =1e-07) {
    nleaves <- length(res$leaves) # nombre de feuilles
    nnodes <- length(res$nodes)    # nombre de noeuds
    node.names <- names(res$nodes) # noms des feuilles
    leave.names <- names(res$leaves) # noms des noeuds    
    dimnodes<-unlist(lapply(res$parts,length)) # nombres de descendants imm�diats de chaque noeud
    effnodes <- dimnodes # recevra le nombre de descendants total de chaque noeud
    wnodes <- lgamma(dimnodes+1) 
    # recevra le logarithme du nombre de permuations compatibles
    # avec la sous-arborescence associ�e � chaque noeud
    


    # les matrices de proximit� #
    a <- matrix(0, nleaves, nleaves)
    ia <- as.numeric(col(a))
    ja <- as.numeric(row(a))
    a <- cbind(ia, ja)[ia < ja, ]
    # a contient la liste des couples de feuilles
    floc1 <- function(x) {
        # x est un couple de num�ros de deux feuilles i, avec i<j
        # Cette fonction renvoie 
        # resw - la distance � la racine du premier anc�tre commun de deux feuilles
        # resa - l'inverse des produits des nombres de descendants des noeuds 
        # rencontr�s sur le plus court chemin entre les deux feuilles
        c1 <- rev(res$paths[[x[1]]])
        c2 <- rev(res$paths[[x[2]]])
        commonnodes <- c1[c1 %in% c2]
        resw <- res$droot[commonnodes[1]]
        d1 <- c1[! (c1 %in% c2)][-1]
        d2 <- c2[! (c2 %in% c1)][-1]
        pathij <- c(d1,d2,commonnodes[1])
        resa <- 1/prod(unlist(dimnodes[pathij]))
        return(c(resw,resa))
    }
    b <- apply(a, 1, floc1)
    names(b) <- NULL
    w <- matrix(0, nleaves, nleaves)
    w[col(w) < row(w)] <- b[1,]
    w <- w + t(w)
    diag(w) <- res$droot[leave.names]
    dimnames(w) <- list(leave.names,1:nleaves)

    res$Wmat <- w
    
    #############################
    # la composante Wmat contient la matrice W des distances racine-premier anc�tre commun
    #############################
    
    w <- diag(res$Wmat)
    w <- matrix(w, nleaves, nleaves)
    w <- w + t(w) - 2 * res$Wmat
    w <- as.dist(sqrt(w))
    attr(w, "Labels") <- leave.names
    
    res$Wdist <- w
    #############################
    # la composante Wdist contient la matrice des racines des distances nodales
    # qui forment une distance euclidienne
    #############################
    w <- res$Wmat
    w <- w / sum(w)
    w <- bicenter.wt(w)
    w <- eigen(w,sym=TRUE)
    res$Wvalues <- w$values[-nleaves]*nleaves
    w <- as.data.frame(qr.Q(qr(scale(w$vectors, scale = FALSE)))[,-nleaves]*sqrt(nleaves))
    row.names(w) <- leave.names
    names(w) = paste("W",1:(nleaves-1),sep="")
    res$Wscores <- w
    
    
    w <- matrix(0, nleaves, nleaves)
    w[col(w) < row(w)] <- b[2,]
    w <- w + t(w)
    # On rajoute la diagonale pour que A soit bistochastique
    floc2 <- function(x) {
        # cette fonction renvoie pour une feuille la fr�quence des repr�sentations
        # compatibles qui placent cette feuille tout en haut ou tout en bas
        c1 <- rev(res$paths[[x]])
        c1 <- c1[-1] # premier ancetre, second ancetre, ..., racine
        resw <- dimnodes[c1] # ordre des noeuds
        resw <- 1/prod(unlist(resw))
        return(resw)
    }
    diag(w) <- unlist(lapply(leave.names,floc2))
    dimnames(w) <- list(leave.names,1:nleaves)
    res$Amat <- w
    #############################
    # la composante Amat contient la matrice des probabilit�s
    # pour une feuille d'�tre juste au dessus d'une autre
    # dans l'ensemble des permutations compatibles
    #############################
    # double centrage
    w <- bicenter.wt(w)
    # diagonalisation
    eig <- eigen (w, sym = TRUE)
    w0 <- abs(eig$values)/max(abs(eig$values))
    w0 <- which(w0<tol)
    if (length(w0)==0) stop ("abnormal output : no null eigenvalue")
    if (length(w0)==1) w0 <- (1:nleaves)[-w0]
    else if (length(w0)>1) {
        # on ajoute le vecteur d�riv� de 1n 
        w <- cbind(rep(1,nleaves),eig$vectors[,w0])
        # on orthonormalise l'ensemble
        w <- qr.Q(qr(w))
        # on met les valeurs propres � 0
        eig$values[w0] <- 0
        # on remplace les vecteurs du noyau par une base orthonorm�e contenant 
        # en premi�re position le parasite
        eig$vectors[,w0] <- w[,-ncol(w)]
        # on enl�ve la position du parasite
        w0 <- (1:nleaves)[-w0[1]]
    }
    rank <- length(w0)
    res$Avalues <- eig$values[w0]*nleaves
    #############################
    # la composante Avalues contient les valeurs propres de QAQ
    #############################
    res$Adim <- sum(res$Avalues>tol)
    #############################
    # la composante Adim contient le nombre de valeurs propres positives
    # associ�es � la composante positive de la variance
    #############################
    w <- eig$vectors[,w0]*sqrt(nleaves)
    w <- data.frame(w)
    row.names(w) <- leave.names
    names(w) <- paste("A",1:rank,sep="")
    res$Ascores <- w
    #############################
    # la composante Ascores contient une base orthoborm�e de l'orthogonal de n
    # pour la pond�ration uniforme. Elle d�finit un phylogramme
    #############################

    # Compl�ment : la valeur des noeuds #    

    floc3 <- function(k) {
        # k est un num�ro de noeud
        # x est un vecteur comportant un nom de noeud et des noms de descendants 
        # de ce noeud. 
        # A la fin parts wnodes contient le logarithme
        # du nombre de permutations compatibles de chaque sous-arbre
        # et effnodes contient le nombre de descendants de chaque sous-arbre
        y <- res$parts[[k]]
        x <- y[y%in%names(res$nodes)]
        n1 <- names(res$parts)[k]
        if (length(x)<=0) return(NULL)
        effnodes[n1] <<- effnodes[n1] - length(x) + sum(effnodes[x])
        wnodes[n1] <<- wnodes[n1] + sum(wnodes[x])
        return(NULL)
    }
    
    lapply(1:length(res$parts),floc3)
    # typolo.value <- 1-exp(wnodes-lgamma(effnodes+1)) abandon
    
    ####res$Aparam <- data.frame(x1=I(dimnodes), x2=I(effnodes), x3=I(wnodes), x4=I(typolo.value))
    res$Aparam <- data.frame(ndir=dimnodes, nlea=effnodes, lnperm=I(wnodes))
    #############################
    # la composante Aparam est un data.frame de param�tre sur l'ensemble des noeuds
    # x1 = nombre de descendants directs
    # x2 = nombre de feuilles descendantes 
    # x3 = log du nombre de permutations compatibles avec la phylog�nie extraite
    # x4 = 1-rapport du nombre de permutations compatibles sur le nombre de permutations totales
    # pour la phylog�nie extraite dans ce noeud
    # cet indice vaut 0 si le noeud est final et est maximal � la racine
    # attention il ne vaut pas 1 mais 1-epsilon quand il est affich� 1
    #############################

    # Compl�ment : la base B #
    w1 <- matrix(0, nleaves, nnodes)
    ####x1 <- res$Aparam$x2 #le nombre de feuilles descendantes
    x1 <- res$Aparam$lnperm #on trie sur le log des permutations
    # on calcule une matrice auxiliaire pour avoir la liste des feuilles descendantes
    # pour chacun des noeuds
    dimnames(w1) <- list(leave.names, names(x1))
    for (i in leave.names) {
        ancetres <- res$paths[[i]]
        ancetres <- rev(ancetres)[-1]#rev(ancetres[-1])[-1]
        w1[i, ancetres] <- 1
    }
    w1 <- cbind(w1, diag(1, nleaves))
    dimnames(w1)[[2]] <- c(names(x1),leave.names)
    x1 <- c(x1, rep(-1,nleaves))
    names(x1) <-dimnames(w1)[[2]]
    # La matrice w1 contient 1 en i-j si la feuille i descend du noeud j

    ######################################
    # on construit une famille d'indicatrices de classes
    # Une ar�te de l'arborescence est un lien de descendance
    # Chaque noeud et chaque feuille (� l'expection de la racine) a un seul ascendant
    # Il y a n+f-1 ar�tes. Le noeud j a m(j) descendants
    # Les feuilles n'en n'ont pas. Donc m(1)+m(2)+ ... + m(n) = n+f-1
    # Il y a n+f-1 ar�tes r�parties en n blocs.
    # Il y a donc n+f-1-n=f-1 descendants indicateurs DI quand on enl�ve une ar�te descendante par noeud
    # Rien n'est conserv� pour un noeud avec un seul descendant
    # Pour chaque DI on utilise l'indicatrice de la classe des feuilles descendant de cet noeud
    # la composante Bindica contient f-1 indicatrices de classes de feuilles
    # names (w) contient des noms de descendants
    # nomuni contient les noms de DI pour l'�tiquetage final
    ####################################
    funnoe <- function (noeud) {
        # renvoie pour un noeud une liste dont chaque composante est un descendant imm�diat du noeud
        # caract�ris� par la liste des feuilles qui en descendent sous forme de matrice 
        # d'indicatrices. Le dernier descendant imm�diat du noeud est �limin�.
        x <- res$parts[[noeud]] # les descendants imm�diats
        xval <- x1[x] # le nombre de feuilles descendantes des descendants
        xval <- rev(sort(xval)) # tri�e
        x <- names(xval) # on r�cup�re lesquels
        x <- x[-length(x)] # on enl�ve le dernier
        if (length(x) ==0) return(NULL)
        if (length(x) ==1) xmat <- matrix(w1[,x],ncol=1,dimnames=list(leave.names,noeud))
        else {
           xmat <- w1[,x]
           dimnames(xmat)[[2]] <- rep(noeud, ncol(xmat))
        }
        return (list(xmat, x))
        # les noms des colonnes de xmat rep�te le nom du noeud
        # dans y on a le nom des descendants retenus
    }

    nomuni <- NULL
    w <- matrix(1,nleaves,1)
    dimnames(w) <- list(leave.names, "un")
    for (i in names(x1)[1:nnodes]) {
        provi <- funnoe(i)
        if (!is.null(provi)) {
            w <-cbind(w, provi[[1]])
            nomuni <- c(nomuni,provi[[2]])
        }
    }
    w <- w[,-1]
    nomrepet <- dimnames(w)[[2]]
    names(nomrepet) <- nomuni
    dimnames(w)[[2]] <- nomuni
    names(nomuni) <- nomuni
    #############################
    # Les indicatrices sont class�es par ordre d�croissant
    # de xtQWQx la variance phylog�n�tique formelle de l'indicatrice centr�e
    # Bindica n'a qu'une valeur p�dagogique et ne sert pas explicitement
    # mais la proc�dure est simple
    # 1) d�finition des indicatrices, il y en a toujours f-1
    # 2) rangement par valeur d�croissante de la forme quadratique
    # Ce rangement est conserv� dans res$Bindica
    # les valeurs du crit�re de rangement dans Bvalues
    # 3) rajout de 1n devant
    # 4) orthonormalisation
    # on obtient toujours une base orthonorm�e de l'orthogonal de 1n
    #############################

    w.val <- x1[nomuni]
    # trie par ordre descendant
    w.val <- rev(sort(w.val))
    # lesquels
    w <- w[,names(w.val)]
    # nomrepet / w sont tri�s
    nomrepet <- nomrepet[names(w.val)]    
    res$Bindica <- as.data.frame(w)
    w <- cbind(rep(1,nleaves),w)
    w <- qr.Q(qr(w))
    w <- w[, -1] * sqrt(nleaves)
    w <- data.frame(w)
    row.names(w) <- leave.names
    names(w) <- paste("B",1:(nleaves-1),sep="")
    res$Bscores <- w
    ### res$Bvalues <- w.val
    lw <- lapply(node.names, function (x) which(nomrepet==x))
    names(lw) <- node.names
    fun1 <- function (x) {
        if (length(x)==0) return("x")
        if (length(x)==1) return(as.character(x))
        y <- x[1]
        for(k in 2:length(x)) y <- paste(y,x[k],sep="/")
        return(y)
    }       
    lw <- unlist(lapply(lw, fun1))
    res$Blabels <- lw
    return(res)
}
