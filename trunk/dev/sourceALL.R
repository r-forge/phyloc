d <- dir() # get all files in dev/
Rsrc <- d[grep("\\.R$", d)] # select only *.R files
Rsrc <- Rsrc[-which(Rsrc == "sourceALL.R")] # remove the current file ;)
cat("start sourcing R files...\n\n")
for (x in Rsrc) try(source(x))
cat("Done\n")
