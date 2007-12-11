d <- dir()
Rsrc <- d[grep("\\.R$", d)]
cat("start sourcing R files...")
for (x in Rsrc) try(source(x))
cat("done\n\n")
