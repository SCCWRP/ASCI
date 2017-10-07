"dapply"<-
function(X, MARGINX, Y = X, MARGINY = MARGINX, FUN, ..., half=T)
{
	if(is.character(FUN))
		FUN <- get(FUN, mode = "function")
	else if(mode(FUN) != "function") {
		farg <- substitute(FUN)
		if(mode(farg) == "name")
			FUN <- get(farg, mode = "function")
		else stop(paste("\"", farg, "\" is not a function", sep = ""))
		
	}
	Y <- Y
	MARGINY <- MARGINY
	half <- half && is.all.equal(X,Y) && is.all.equal(MARGINX,MARGINY)
#  SET UP FOR X
	outX <- seq(1, length(dim(X)))[MARGINX]
	if(!is.array(X)) {
		X <- as.matrix(X)
		if(MARGINX[1] == 0)
			MARGINX <- -1
		if(MARGINX == -1)
			outX <- NULL
		else outX <- 1
	}
	else {
		if(MARGINX[1] == 0)
			MARGINX <-  - seq(along = dim(X))
		if(length(MARGINX) == length(dim(X))) {
			X <- array(X, c(dim(X), 1), c(ldimnames(X), list(NULL))
				)
			if(MARGINX[1] < 0)
				outX <- NULL
		}
	}
	dX <- dim(X)
	if(length(class(X))) {
		if(length(dX) == 2)
			X <- as.matrix(X)
		else X <- as.array(X)
	}
	dnX <- ldimnames(X)
	dnX.call <- dnX[ - MARGINX]
	dnX <- dnX[outX]
	dimnames(X) <- NULL
	newX <- aperm(X, c(seq(1, length(dX))[ - MARGINX], seq(1, length(dX))[
		MARGINX]))
	subdimX <- dX[ - MARGINX]
	dim(newX) <- c(prod(subdimX), prod(dX[MARGINX]))
#   SET UP FOR Y
	outY <- seq(length(dim(Y)))[MARGINY]
	if(!is.array(Y)) {
		Y <- as.matrix(Y)
		if(MARGINY[1] == 0)
			MARGINY <- -1
		if(MARGINY == -1)
			outY <- NULL
		else outY <- 1
	}
	else {
		if(MARGINY[1] == 0)
			MARGINY <-  - seq(along = dim(Y))
		if(length(MARGINY) == length(dim(Y))) {
			Y <- array(Y, c(dim(Y), 1), c(ldimnames(Y), list(NULL))
				)
			if(MARGINY[1] < 0)
				outY <- NULL
		}
	}
	dY <- dim(Y)
	if(length(class(Y))) {
		if(length(dY) == 2)
			Y <- as.matrix(Y)
		else Y <- as.array(Y)
	}
	dnY <- ldimnames(Y)
	dnY.call <- dnY[ - MARGINY]
	dnY <- dnY[outY]
	dimnames(Y) <- NULL
	newY <- aperm(Y, c(seq(1, length(dY))[ - MARGINY], seq(1, length(dY))[
		MARGINY]))
	subdimY <- dY[ - MARGINY]
	dim(newY) <- c(prod(subdimY), prod(dY[MARGINY]))	
#    REPLICATE newX ncol(newY) times AND newY ncol(newX) times
#    if half then choose unique pairs only
	nX <- ncol(newX)
	nY <- ncol(newY)
	if(half){
		nXY <- nX * (nX-1)/2
		if(nXY==0) return()
		sel <- row(matrix(0,nX,nX)) > col(matrix(0,nX,nX))
	} else {
		nXY <- nX * nY
		sel <- rep(T,nXY)
	}
	NewX <- newX[, rep(seq(nX), nY)[sel], drop = F]
	NewY <- newY[, rep(seq(nY), rep(nX, nY))[sel], drop = F]
#    COMPUTATIONS
	ans <- vector("list", nXY)
	if(length(subdimX) > 1) {
		if(length(subdimY) > 1) {
			for(i in 1:nXY)
				ans[i] <- list(FUN(array(NewX[, i], subdimX,
					dnX.call), array(NewY[, i], subdimY,
					dnY.call), ...))
		}
		else {
			dimnames(NewY) <- list(dnY.call[[1]], NULL)
			for(i in 1:nXY)
				ans[i] <- list(FUN(array(NewX[, i], subdimX,
					dnX.call), NewY[, i], ...))
		}
	}
	else {
		if(length(subdimY) > 1) {
			dimnames(NewX) <- list(dnX.call[[1]], NULL)
			for(i in 1:nXY)
				ans[i] <- list(FUN(NewX[, i], array(NewY[,
					i], subdimY, dnY.call), ...))
		}
		else {
			dimnames(NewX) <- list(dnX.call[[1]], NULL)
			dimnames(NewY) <- list(dnY.call[[1]], NULL)
			for(i in 1:nXY)
				ans[i] <- list(FUN(NewX[, i], NewY[, i], ...))
			
		}
	}
#  OUTPUT
	list.names <- NULL
	if(half && (length(unlist(dnX))>0|length(unlist(dnY))>0))
		list.names <- paste(
			expandg(dnX,dX[outX])[rep(seq(nX), nY)[sel]],
			expandg(dnY,dY[outY])[rep(seq(nY), rep(nX, nY))[sel]],
			sep="-")
	ans.names <- names(ans[[1]])
	ans.dimnames <- ldimnames(ans[[1]])
	ret.list <- is.recursive(ans[[1]])
	ret.array <- is.array(ans[[1]])
	dimarray <- dim(ans[[1]])
	first.length <- length(ans[[1]])
	if(!ret.list)
		for(i in 1:nXY)
			if(length(ans[[i]]) != first.length) ret.list <- T
	if(ret.array)
		for(i in 1:nXY)
			if(!is.all.equal(dimarray, dim(ans[[i]]))) ret.array <-
					F
	if(ret.array)
		d <- dim(ans[[1]])
	if(!ret.list)
		ans <- unlist(ans, recursive = F)
	if(length(ans) == 0) return()
	if(length(ans) == nXY){
		if(half){ names(ans) <- list.names; return(ans) }
		return(simplify(array(ans, c(dX[outX], dY[outY]),
			dimnames = c(dnX, dnY))))
		}
	if(length(ans) %% nXY == 0) {
		if(half){
		if(ret.array) return(simplify(array(ans, c(d, nXY),
				c(ans.dimnames, list(list.names)))))
		return(simplify(array(ans, c(length(ans)/nXY, nXY),
				list(ans.names, list.names))))
		}
		if(ret.array)
			return(simplify(array(ans, c(d, dX[outX], dY[outY]),
				c(ans.dimnames, dnX, dnY))))
		return(simplify(array(ans, c(length(ans)/nXY, dX[outX],
				dY[outY]), c(list(ans.names), dnX, dnY))))
	}
	return(ans)
}

"usei"<- function(n=2) eval(substitute(i, sys.parent(n)))

"expandg" <- function(args,d,sep="|"){
#  args is a list each of whose entries is either a vector of length d[i]
#  or NULL.  Expand the grid guided by d when the entry is NULL and paste.
#  Simplified from expand.grid()
	ndim <- length(args)
	started <- FALSE
	individual.rep <- 1
	for(i in seq(1,ndim)) {
		if(length(args[[i]])==d[i]){
			group.rep <- prod(d[ - (1:i)])
			this <- rep(rep(args[[i]], rep(individual.rep,d[i])),
				group.rep)
			if(!started){ result <- this; started <- TRUE }
			else result <- paste(result,this,sep=sep)
		}
		individual.rep <- individual.rep * d[i]
	}
	if(!started) return(rep("",prod(d)))
	return(result)
}

"ldimnames"<- function(x) {
	d <- dimnames(x)
	if(is.null(d))
		d <- rep(list(NULL), length(dim(x)))
	d
}
"is.all.equal"<- function(target, current)
if(is.logical(temp <- all.equal(target, current)) && temp) TRUE else FALSE

"makedimnames"<- function(x, wh = seq(length(dim(x)))) {
	if(!is.array(x)){
		names(x) <- paste(substitute(x),".",seq(length(x)),sep="")
		return(x)
	}
	d <- dim(x)
	n <- length(dim(x))
	if(is.null(dimnames(x))) dimnames(x) <- rep(list(NULL), n)
	for(i in wh)
		dimnames(x)[[i]] <- paste(substitute(x), i, ".",
			seq(d[i]), sep = "")
	x
}
"simplify"<- function(a) {
	if(length(unlist(dimnames(a))) == 0)
		dimnames(a) <- NULL
	if(length(dim(a)) == 1) {
		aa <- as.vector(a)
		names(aa) <- unlist(dimnames(a))
		return(aa)
	}
	else a
}
