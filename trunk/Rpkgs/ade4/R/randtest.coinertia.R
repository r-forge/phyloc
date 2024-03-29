"randtest.coinertia" <- function(xtest, nrepet=999, fixed=0, ...) {
  nrepet<-nrepet+1
  if (!inherits(xtest,"dudi"))
    stop("Object of class dudi expected")
  if (!inherits(xtest,"coinertia"))
    stop("Object of class 'coinertia' expected")
  appel<-as.list(xtest$call)
  dudiX<-eval(appel$dudiX,sys.frame(0))
  dudiY<-eval(appel$dudiY,sys.frame(0))
  X<-dudiX$tab
  X.cw<-dudiX$cw
  X.lw<-dudiX$lw
  appelX<-as.list(dudiX$call)
  apx<-appelX$df
  Xinit<-eval(appelX$df,sys.frame(0))
  if (appelX[[1]] == "dudi.pca") {        
    appelX$scale<-eval(appelX$scale,sys.frame(0))
    appelX$center<-eval(appelX$center,sys.frame(0))
    if (is.null(appelX$scale)) appelX$scale<-TRUE
    if (is.null(appelX$center)) appelX$center<-TRUE
    if(!(is.logical(appelX$center))) stop("Not implemented for decentred PCA: read the documentation file.")
    if (appelX$center == FALSE && appelX$scale == FALSE) typX<-"nc"
    if (appelX$center == FALSE && appelX$scale == TRUE) typX<-"cs"
    if (appelX$center == TRUE  && appelX$scale == FALSE) typX<-"cp"
    if (appelX$center == TRUE  && appelX$scale == TRUE) typX<-"cn"
  } else if (appelX[[1]] == "dudi.coa") {
    typX<-"fc"
  } else if (appelX[[1]] == "dudi.fca") {
    typX<-"fc"
  } else if (appelX[[1]] == "dudi.acm") {
    typX<-"cm"
    Xinit <- acm.disjonctif(Xinit)
  } else stop ("Not yet available")
  Y<-dudiY$tab
  Y.cw<-dudiY$cw
  Y.lw<-dudiY$lw
  appelY<-as.list(dudiY$call)
  apy<-appelY$df
  Yinit<-eval(appelY$df,sys.frame(0))
  if (appelY[[1]] == "dudi.pca") {        
    appelY$scale<-eval(appelY$scale,sys.frame(0))
    appelY$center<-eval(appelY$center,sys.frame(0))
    if (is.null(appelY$scale)) appelY$scale<-TRUE
    if (is.null(appelY$center)) appelY$center<-TRUE
    if(!(is.logical(appelY$center))) stop("Not implemented for decentred PCA: read the documentation file.")
    if (appelY$center == FALSE && appelY$scale == FALSE) typY<-"nc"
    if (appelY$center == FALSE && appelY$scale == TRUE) typY<-"cs"
    if (appelY$center == TRUE  && appelY$scale == FALSE) typY<-"cp"
    if (appelY$center == TRUE  && appelY$scale == TRUE) typY<-"cn"
  } else if (appelY[[1]] == "dudi.coa") {
    typY<-"fc"
  } else if (appelY[[1]] == "dudi.fca") {
    typY<-"fc"
  } else if (appelY[[1]] == "dudi.acm") {
    typY<-"cm"
    Yinit <- acm.disjonctif(Yinit)
  } else stop ("Not yet available")
  if (identical(all.equal(X.lw,Y.lw),TRUE)) {
    if (identical(all.equal(X.lw,rep(1/nrow(X), nrow(X))),TRUE)) {
      isim<-testertrace(nrepet, X.cw, Y.cw, X, Y, nrow(X), ncol(X), ncol(Y))
    } else {
      if (fixed==0) {
        cat("Warning: non uniform weight. The results from simulations\n")
        cat("are not valid if weights are computed from analysed data.\n")
        isim<-testertracenu(nrepet, X.cw, Y.cw, X.lw, X, Y, nrow(X), ncol(X), ncol(Y), Xinit, Yinit, typX, typY)
	if(typX=="cm") {isim[-1]=isim[-1]/ncol(eval(appelX$df,sys.frame(0)))}
	if(typY=="cm") {isim[-1]=isim[-1]/ncol(eval(appelY$df,sys.frame(0)))}
      } else if (fixed==1) {
        cat("Warning: non uniform weight. The results from permutations\n")
        cat("are valid only if the row weights come from the fixed table.\n")
        cat("The fixed table is table X : ")
        print(apx)
        isim<-testertracenubis(nrepet, X.cw, Y.cw, X.lw, X, Y, nrow(X), ncol(X), ncol(Y), Xinit, Yinit, typX, typY, fixed)
	if(typY=="cm") {isim[-1]=isim[-1]/ncol(eval(appelY$df,sys.frame(0)))
	}
      } else if (fixed==2) {
        cat("Warning: non uniform weight. The results from permutations\n")
        cat("are valid only if the row weights come from the fixed table.\n")
        cat("The fixed table is table Y : ")
        print(apy)
        isim<-testertracenubis(nrepet, X.cw, Y.cw, X.lw, X, Y, nrow(X), ncol(X), ncol(Y), Xinit, Yinit, typX, typY, fixed)
	if(typX=="cm") {isim[-1]=isim[-1]/ncol(eval(appelX$df,sys.frame(0)))}
	
      }
      else if (fixed>2) {stop ("Error : fixed must be =< 2")}
    }
                                        # On calcule le RV a partir de la coinertie
    isim<-isim/sqrt(sum(dudiX$eig^2))/sqrt(sum(dudiY$eig^2))
    obs<-isim[1]
    return(as.randtest(isim[-1],obs,call=match.call()))
  } else {
    stop ("Equal row weights expected")
  }
}
