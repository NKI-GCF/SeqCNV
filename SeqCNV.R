findCovFiles <- function(pattern,...) {
	list.files(pattern=pattern, ...)
}

loadCovData <- function(files, gc=NULL, mapa=NULL, black=NULL, excludechr=NULL, datacol=5) {
	d <- lapply(files, read.delim, header=F)
	cnames <- removeCommonFix(files)
	names(d) <- cnames

	#checks
	if(length(unique(sapply(d,nrow))) != 1) {
		stop("Files have variable number of lines and are not compatible")
	}

	#naturally sort data 
	require(gtools)
	d <- lapply(d, function(x) {
		x[,1] <- factor(x[,1], levels=mixedsort(levels(x[,1])))
		x[order(x[,1], x[,2]),]
		})

	#check row order in data files
	for (i in 2:length(files)) {
		if(!isTRUE(all.equal(d[[1]][,1:3], d[[i]][,1:3], check.attributes=FALSE))) stop(paste("Datafiles 1 and", i,"are not in the same order"))
	}

	cov <- sapply(d, function(x) x[!(x[,1] %in% excludechr),datacol])
	anno <- data.frame(droplevels(d[[1]][!(d[[1]][,1] %in% excludechr),1:3]))
	colnames(anno) <- c("chr","start","end")
	annoid <- paste(anno$chr, anno$start, anno$end, sep=":")
	if(!is.null(gc)) {
		anno$gc <- gc[match(annoid, paste(gc$chr, gc$start, gc$end, sep=":")),"gc"]
	}
	if(!is.null(mapa)) {
		anno$mapa <- mapa[match(annoid, paste(mapa$chr, mapa$start, mapa$end, sep=":")),"mapa"]
	}
	if(!is.null(black)) {
		suppressPackageStartupMessages(require(IRanges))
		anno$black <- rep(F, nrow(anno))
		for (c in levels(anno$chr)) {
			ar <- IRanges(anno$start[anno$chr == c], anno$end[anno$chr == c])
			br <- IRanges(black$start[black$chr == c], black$end[black$chr == c])
			anno$black[anno$chr == c] <- ar %over% br
		}
	}
	list(cov=cov, anno=anno)
}

normalize <- function(cnv, excludeXYM=T, trimOutliers=T, useRefCols=NULL) {
	sexchromosomes <- if(excludeXYM) c("X","Y","MT", "chrX", "chrY", "chrM")

	cov <- cnv$cov
	trim <- rep(F, nrow(cov))
	# if there is a (set of) reference columns use those to deterime the outliers, 
	# otherwise use the rowMedians
	if(trimOutliers){ 
		ref <- if (is.null(useRefCols)) apply(cov,1,median) 
		       else rowMeans(cov[,useRefCols, drop=F])
		meancov <- mean(ref, trim=.1)
		meansd <- sd.trim(ref, trim=.01)
		trim <- (ref < (meancov-3*meansd)) | (ref > (meancov+3*meansd))
		cat("trimming", sum(trim), "from", length(trim), "bins\n")
	}
	
	cs <- colSums(cov[!(cnv$anno[,1] %in% sexchromosomes | trim),])

	norm <- t(t(cov) / cs)
	list(cov=norm, anno=cnv$anno)
}

normalizeSampleGC <- function(x, gc, maxpoints=15000, plot=F) {
	x <- 2^x
	valid <- is.finite(x) & !is.na(gc) #& x > 0
	use <- which(valid)
	if(length(use) > maxpoints) use <- sample(use, maxpoints)

	df <- data.frame(x=x[use],gc=gc[use])
	fit <- loess(x~gc, data=df)
	normv <- rep(NA, length(x))
	normv[valid] <- predict(fit, data.frame(gc=gc[valid]))
	if(plot) {
		plot(x~gc, data=df)
		lines(normv, col=2)
	}
	log2(x/(normv/median(normv, na.rm=T)))
}

normalizeGC <- function(ratios) {
	m <- ratios$ratios
	for (i in 1:ncol(m))
	m[,i] <- normalizeSampleGC(m[,i],ratios$anno$X5_pct_gc)
	list(ratios=m, anno=ratios$anno)
}

makeRatios <- function(cnv, reference=c("columns", "rowmedian", "rowmeans", "colmedian"), refCols=NULL, fixNonFinite=T) {

	methods <- c("columns", "rowmedian", "rowmeans", "colmedian")
	refmethod <- pmatch(reference, methods)[1]
	ref <- if(refmethod == 1) {
				#todo check col numbers/ids on data
				cat("Using column(s)", refCols,"to calculate ratios\n")
				rowMeans(cnv$cov[,refCols,drop=F])
			} else {
				switch(methods[refmethod], rowmedian=apply(cnv$cov,1,median),
				                  rowmeans=rowMeans(cnv$cov),
								  colmedian=matrix(apply(cnv$cov,2,median), nrow=nrow(cnv$cov), ncol=ncol(cnv$cov), byrow=T))
			}
	ratios <- log2(cnv$cov) - log2(ref)
	#remove refCol from coverage if one ref is used (will be 0)
	if(length(refCols) == 1) ratios <- ratios[,-refCols, drop=F]
	if(fixNonFinite) {
		ratios[ratios == -Inf] <- -4 #zero cov in sample
		ratios[ratios == Inf] <-0  #zero cov in reference (ratio unknown)
		ratios[is.na(ratios)] <-0 	#other stuff
	}

	list(ratios=ratios, anno=cnv$anno)
}


tng <- function(df, use, correctmapa=TRUE,  plot=NULL, verbose=T) {
	
	#tests
	if(!is.logical(use) && length(use) ==nrow(df))
		stop("use should be logicval vector with same size as df")
	#df colums?

	if(!is.null(plot)) {
		if(!is.logical(plot)) {
			if(verbose) cat("Plotting to file", plot,"\n")
			png(plot, width=700, height=1400)
			par(mfrow=c(2,1))
			on.exit(dev.off())
			plot <- TRUE
		} else if(plot) {
			par(mfrow=c(2,1))
		}
	}

	#exclude contains the points to exclude in the 
	#fitting (usually sex chromosomes and blacklisted regions)
	# gc fits  also excludes the low mappability data

	#correct gc using double lowess
	gcuse <- (use & !is.na(df$mapa) & df$mapa > .8 & !is.na(df$gc) & df$gc > 0)
	rough <- loess(count ~ gc, data=df, subset=gcuse, span = 0.03)
	i <- seq(0, 1, by = 0.001)
	final <- loess(predict(rough, i) ~ i, span = 0.3)
	normv <- predict(final, df$gc)
	df$countgcloess <- df$count/(normv/median(normv, na.rm=T))

	if(plot) {
		plot(count ~ gc, data=df, subset=gcuse, ylim=quantile(df$count[gcuse], c(0.0001, .999)), xlim=c(.1,.8), pch=".")
		points(count ~ gc, data=df, subset=!gcuse, col=rgb(1,0,0,.3), pch=".")
		lines(i, predict(rough, i), col=3)
		points(df$gc, normv, col=2, pch=".")
	}

	#correct mapa using linear function that intercepts zero
	if(correctmapa) {
	mapause <- (use & !is.na(df$mapa))
	lm(countgcloess~0+mapa, data=df, subset=mapause) ->fll
	if(verbose) print(summary(fll))

	if (plot) {
		plot(countgcloess ~ mapa, data=df, subset=mapause, ylim=quantile(df$countgcloess, c(0.0001, .999), na.rm=T), pch=".")
		points(countgcloess ~ mapa, data=df, subset=!mapause, col=rgb(1,0,0,.3), pch=".")
		abline(0, fll$coef, col=2)
	}

	return(log2(df$countgcloess / (df$mapa * fll$coef)))
	} else {
		#corerct agains median value (exluding sex chr)
		log2(df$countgcloess / median(df$countgcloess[use], na.rm=T))
	}	
}





clusterCor <- function(cnv, excludeXYM=T) {
	sexchromosomes <- if(excludeXYM) c("X","Y","MT", "chrX", "chrY", "chrM")
	plot(hclust(as.dist(1-cor(cnv$ratios[!(cnv$anno[,1] %in% sexchromosomes),], use="pair"))))
}

plotCNV <- function(ratios, sample=1, chr=NULL, ...) {
	require(gtools)
	m <- tapply(ratios$anno[,3], ratios$anno[,1], max)
	#order the chr levels
	ord <- mixedsort(names(m))
	ofs <- c(0,cumsum(m[ord]/1e6))
	names(ofs) <- c(ord, "end")

	pos <- rowMeans(ratios$anno[,2:3])/1e6
	if(is.null(chr)) pos <- pos + ofs[as.character(ratios$anno[,1])]

	y <- ratios$ratios[,sample]

	df <- data.frame(x=pos, y=y)
	subset <- 1:nrow(df)
	if(!is.null(chr)) subset <- which(ratios$anno[,1] == chr)
	plot(y ~ x, data=df, subset=subset,  xlab="Genomic position", ylab="Log2 ratio", ...)
	#highlight spots outside plot region!
	out <- which(y > par("usr")[4] | y < par("usr")[3])
	op <- par("xpd")
	par(xpd=T)
	if (length(out) > 0) points(x=pos[out],ifelse(y[out] < par("usr")[3], par("usr")[3], par("usr")[4]) , pch=4, cex=.7)
	par(xpd=op)

	abline(h=0, col="gray")
	abline(v=ofs, col="lightblue")
	ym <- par("usr")[4] * .9
	text(x=(ofs[1:(length(ofs)-1)] + ofs[2:length(ofs)]) / 2, y=ym, labels=ord, cex=.5)
}

plotAll <- function(ratios, ...) {
	if(ncol(ratios$ratios) > 1) 
		par(mfrow=c(2,2))
	else
		par(mfrow=c(1,1))
	for(i in seq.int(1,ncol(ratios$ratios),4)) {
		for(p in i:min((i+3),ncol(ratios$ratios))) {
			plotCNV(ratios, sample=p, pch=".", main=colnames(ratios$ratios)[p], ...)
		}
		if(dev.capabilities()$locator) readline()
	}
}

writeNexusNormalized <- function(r,anno, filename, path=".") {
	writeto <- paste(path, filename, sep="/")
	
	write.table(
		data.frame(CHROMOSOME=paste0("chr",anno[,1]),
		           CHR_POSITION=floor(rowMeans(anno[,2:3])),
		           RAW_DATAPOINTS=rep(1,nrow(anno)), 
				   POSITION_COUNT=rep(1,nrow(anno)),
				   RATIO_CORRECTED=round(r, digits=5),
				   WINDOW_SIZE=rep(1,nrow(anno))), 
			file=writeto, sep="\t", row.names=F, quote=F)
	cat("Written to", writeto,"\n")
}

writeAllNexusNormalized <- function(data, filenames=NULL, path=".") {
	if(is.null(filenames)) filenames <- paste(colnames(data$ratios), "txt", sep=".")
	print(filenames)
	for(i in 1:ncol(data$ratios))
		writeNexusNormalized(data$ratios[,i], data$anno, filename=filenames[i], path=path)
}

writeAllCN <- function(data, filename, path=".") {
	writeto <- paste(path, filename, sep="/")
	df <- data.frame(SNP=paste(data$anno[,1], data$anno[,2],sep=":"),
	   Chromosome=data$anno[,1],
	   PhysicalPosition=floor(rowMeans(data$anno[,2:3])),
	   data$ratios)
	
	write.table(df,file=writeto, sep="\t", row.names=F, quote=F)
}

toHMMCopy <- function(covdata) {
	require("IRanges")
	ran <- IRanges(start=covdata$anno$start, end=covdata$anno$end)
	sapply(1:ncol(covdata$cov), function(s) {
		RangedData(ran, space=covdata$anno$chr, reads=covdata$cov[,s], gc=covdata$anno$gc, map=covdata$anno$mapa)
	})
}

sd.trim <- function(x, trim=0, na.rm=FALSE, ...)
{
  if(!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
    warning("argument is not numeric or logical: returning NA")
    return(NA_real_)
  }
  if(na.rm) x <- x[!is.na(x)]
  if(!is.numeric(trim) || length(trim) != 1)
    stop("'trim' must be numeric of length one")
  n <- length(x)
  if(trim > 0 && n > 0) {
     if(is.complex(x)) stop("trimmed sd are not defined for complex data")
     if(trim >= 0.5) return(0)
     lo <- floor(n * trim) + 1
     hi <- n + 1 - lo
     x <- sort.int(x, partial = unique(c(lo, hi)))[lo:hi]
  }
  sd(x)
}


removeCommonFix <- function(names, distance=1) {
	l <- strsplit(names,"")

	#clip prefix
	pclip <- 1
	while( length(unique(sapply(l, "[", pclip))) <= distance) {
		pclip <- pclip + 1
	}

	#reverse strings for end clip pos
	l <- lapply(l, rev)
	eclip <- 1
	while( length(unique(sapply(l, "[", eclip))) <= distance) {
		eclip <- eclip + 1
	}

	sapply(names, function(x) substr(x, pclip, nchar(x) - eclip), USE.NAMES=F)

}

fastseg <- function(data) {
	require(fastseg)

	seg <- fastseg(scdata$ratios)
	toDNAcopyObj(seg, data$anno$chr, data$anno$end, data$ratios, colnames(data$ratios))
}
