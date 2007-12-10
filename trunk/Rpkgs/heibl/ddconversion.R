###	converting degrees to decimals 
dms2dec <- function(filename)
{
dms <- read.table(file = "filename", header=FALSE)
dec.long <- dms[1] + dms[2]/60 + dms[3]/3600 
dec.lat <- dms[5] + dms[6]/60 +dms[7]/3600
i <- 1
n <- dim(dms)
repeat{	x <- dms[i,4]=="W"
	if(x == FALSE)
	{dec.long[2,] <- -dec.long[2,]} else{}
	i <- i+1 
	if (i<n[1]) next 
	if (i==n[1]+1) break}
i <- 1
n <- dim(dms)
repeat{	x <- dms[i,8]=="S"
	if(x == FALSE)
	{dec.lat[i,] <- -dec.lat[i,]} else{}
	i <- i+1 
	if (i<n[1]) next 
	if (i==n[1]+1) break}
dec <- matrix(c(dec.long$V1, dec.lat$V5), ncol=2)
write.table(dec, file="decimal.txt", col.names=FALSE, row.names=FALSE, sep="\t", dec=",")
}

###	converting decimals to degrees

dec2dms <- function(input.file)

{
dec <- read.table(file=input.file, header=FALSE, sep="\t", dec="," )
deg.long <- as.integer(dec$V1)
min2.long <- (dec$V1 - deg.long)*60
min.long <- as.integer(min2$V1)
sec.long <- (min2.long - min.long)*60
sec.long <- round(sec.long, digits=0)
min.long <- abs(min.long)
sec.long <- abs(sec.long)
dms.long <- matrix(c(deg.long, min.long, sec.long), ncol=3)
deg.lat <- as.integer(dec$V2)
min2.lat <- (dec$V2 - deg.lat)*60
min.lat <- as.integer(min2.lat)
sec.lat <- (min2.lat - min.lat)*60
sec.lat <- round(sec.lat, digits=0)
min.lat <- abs(min.lat)
sec.lat <- abs(sec.lat)
dms.lat <- matrix(c(deg.lat, min.lat, sec.lat), ncol=3)
i <- 1
n <- dim(dms)
a <- dms.long[i,1]<0
if(a == TRUE)
{x <- "W"} else{x <- "E"}
i <- i+1
repeat{	a <- dms.long[i,1]<0
	if(a == TRUE)
	{x <- c(x,"W")} else{x <- c(x, "E")}
	i <- i+1 
	if (i<n[1]) next 
	if (i==n[1]+1) break}
i <- 1
n <- dim(dms)
b <- dms.lat[i,1]<0
if(b == TRUE)
{y <- "S"} else{y <- "N"}
i <- i+1
repeat{	b <- dms.lat[i,1]<0
	if(b == TRUE)
	{y <- c(y,"S")} else{y <- c(y, "N")}
	i <- i+1 
	if (i<n[1]) next 
	if (i==n[1]+1) break}
x <- as.factor(x)
y <- as.factor(y)
dms.long <- matrix(c(abs(deg.long), min.long, sec.long), ncol=3)
dms.lat <- matrix(c(abs(deg.lat), min.lat, sec.lat), ncol=3)
com <- list(dms.long, x, dms.lat, y)
write.table(com, file="degrees_minutes_seconds.txt", col.names=FALSE, row.names=FALSE, sep="\t", dec=",")
}


