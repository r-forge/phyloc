# C.Heibl 29.06.2007

# E.Paradis' DNA mining function optimized for MacOSX

# seq: the sequence file you want to test
# format: do not change anything here
# itree: you might specify an intial tree
# exclude: models to be excluded
# execname: do not change anything here
# path2exec: specify path to executable, but default should work

# DO NOT CLOSE R WHILE RUNNING PHYML!!!

phymltest.mac <- function (seq, format = "sequential", itree = NULL, exclude = NULL, execname="phyml_macOSX", path2exec = FALSE) 
{
	# fix path
	if (!is.character(path2exec)){
		path2exec <- system(paste("locate", execname), intern=TRUE)
		path2exec <- path2exec[1]
		path2exec <- gsub(execname, "", path2exec)
		}
	cat(paste("\npath.paml was set to", path2exec, "\n"))
	setwd(path2exec)
	write.dna.phylip(seq, file="R", path=NULL)
	write.tree(itree, file="INITALTREE")
    outfile <- "R.phylip_phyml_stat.txt"
    inp <- "R.phylip"
    if (file.exists(outfile)) 
        inp <- c(inp, "A")
    if (file.exists("R.phylip_phyml_tree.txt")) 
        inp <- c(inp, "A")
    if (format != "interleaved") 
        inp <- c(inp, "I")
    if (!is.null(itree)) 
        inp <- c(inp, "U", "INITALTREE")
    N <- length(.phymltest.model)
    input.model <- list(c(rep("M", 5), "Y"), c(rep("M", 5), "V", 
        rep("Y", 2)), c(rep("M", 5), "R", "A", rep("Y", 2)), 
        c(rep("M", 5), "R", "A", "Y", "V", rep("Y", 2)), c(rep("M", 
            6), "T", rep("Y", 2)), c(rep("M", 6), "T", "Y", "V", 
            rep("Y", 2)), c(rep("M", 6), "T", "Y", "R", "A", 
            rep("Y", 2)), c(rep("M", 6), "T", "Y", "R", "A", 
            "Y", "V", rep("Y", 2)), c(rep("M", 7), "Y"), c(rep("M", 
            7), "V", rep("Y", 2)), c(rep("M", 7), "R", "A", rep("Y", 
            2)), c(rep("M", 7), "V", "Y", "R", "A", rep("Y", 
            2)), c("M", "T", rep("Y", 2)), c("M", "T", "Y", "V", 
            rep("Y", 2)), c("M", "T", "Y", "R", "A", rep("Y", 
            2)), c("M", "T", "Y", "V", "Y", "R", "A", rep("Y", 
            2)), c("T", rep("Y", 2)), c("T", "Y", "V", rep("Y", 
            2)), c("T", "Y", "R", "A", rep("Y", 2)), c("T", "Y", 
            "V", "Y", "R", "A", rep("Y", 2)), c(rep("M", 2), 
            "T", rep("Y", 2)), c(rep("M", 2), "T", "Y", "V", 
            rep("Y", 2)), c(rep("M", 2), "T", "Y", "R", "A", 
            rep("Y", 2)), c(rep("M", 2), "T", "Y", "R", "A", 
            "Y", "V", rep("Y", 2)), c(rep("M", 3), "Y"), c(rep("M", 
            3), "V", rep("Y", 2)), c(rep("M", 3), "R", "A", rep("Y", 
            2)), c(rep("M", 3), "V", "Y", "R", "A", rep("Y", 
            2)))
    loglik <- numeric(N)
    names(input.model) <- names(loglik) <- .phymltest.model
    if (is.null(path2exec)) exec <- execname
    else exec <- paste(path2exec, execname, sep = "./")
    imod <- if (is.null(exclude)) 1:N
    else {(1:N)[!.phymltest.model %in% exclude]}
    for (i in imod) {
        if (i == 2) {
            if (length(inp) == 1) 
                inp <- c(inp, rep("A", 2))
            else if (inp[2] != "A") 
                inp <- c(inp[1], rep("A", 2), inp[2:length(inp)])
        }
        cat(c(inp, input.model[[i]]), file = "f", sep = "\n")
        system(paste(exec, "f", sep = " < "))
        loglik[i] <- scan("R.phylip_phyml_lk.txt", 
            quiet = TRUE)
    }
    unlink("f")
    loglik <- loglik[imod]
    class(loglik) <- "phymltest"
    loglik
}
