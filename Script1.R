
library(ggplot2)
library(gtools)
library(reshape2)
library(rjags)
library(runjags)
library(parallel)
library(coda)
library(mvtnorm)
library(rstan)



runjags.options( inits.warning=FALSE, rng.warning=FALSE )
nChainsDefault = 4
runjagsMethodDefault = "parallel"


###############################################################################################


openMyWindow = function( width, height, ... ) {
    tryInfo = try( windows( width=width, height=height, ... ) ) 
  return ( tryInfo );
}    

openGraph = function( width=7, height=7, mag=1.0, ... ) {
  tryInfo = openMyWindow( width=width*mag, height=height*mag, ...) ;
  if ( class(tryInfo)=="try-error" ) {
    lineInput = readline("WARNING: Previous graphics windows will be closed -- too many open windows?\nTO CONTINUE, PRESS <ENTER> IN R CONSOLE.\n")
    graphics.off()
    openMyWindow( width=width*mag, height=height*mag, ...) ;
  }
}

saveGraph = function( file="saveGraphOutput", type="png", ... ) {
  file=paste0(file,".",type) 
  savePlot( file=file, type=type, ... )
}

#------------------------------------------------------------------------------
# Functions for computing limits of HDI's:

HDIofMCMC <- function( sampleVec, credMass=0.95 ) {
  sortedPts <- sort( sampleVec )
  ciIdxInc <- ceiling( credMass * length( sortedPts ) )
  nCIs <- length( sortedPts ) - ciIdxInc
  ciWidth <- rep( 0, nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] <- sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin <- sortedPts[ which.min( ciWidth ) ]
  HDImax <- sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim <- c( HDImin, HDImax )
  return( HDIlim )
}
#------------------------------------------------------------------------------
# Function(s) for plotting properties of mcmc coda objects.

DbdaAcfPlot <- function( codaObject, parName=varnames(codaObject)[1], plColors=NULL ) {
  if ( all( parName != varnames(codaObject) ) ) { 
    stop("parName must be a column name of coda object")
  }
  nChain <- length(codaObject)
  if ( is.null(plColors) ) plColors=1:nChain
  xMat <- NULL
  yMat <- NULL
  for ( cIdx in 1:nChain ) {
    acfInfo <- acf(codaObject[,c(parName)][[cIdx]],plot=FALSE) 
    xMat <- cbind(xMat,acfInfo$lag)
    yMat <- cbind(yMat,acfInfo$acf)
  }
  matplot( xMat, yMat, type="o", pch=20, col=plColors, ylim=c(0,1),
           main="", xlab="Lag", ylab="Autocorrelation" )
  abline(h=0,lty="dashed")
  EffChnLngth <- effectiveSize(codaObject[,c(parName)])
  text( x=max(xMat), y=max(yMat), adj=c(1.0,1.0), cex=1.25,
        labels=paste("ESS =",round(EffChnLngth,1)) )
}

DbdaDensPlot <- function( codaObject, parName=varnames(codaObject)[1], plColors=NULL ) {
  if ( all( parName != varnames(codaObject) ) ) { 
    stop("parName must be a column name of coda object")
  }
  nChain <- length(codaObject) # or nchain(codaObject)
  if ( is.null(plColors) ) plColors=1:nChain
  xMat <- NULL
  yMat <- NULL
  hdiLims <- NULL
  for ( cIdx in 1:nChain ) {
    sampleVec <- codaObject[, c(parName)][[cIdx]]
    densInfo <- density(sampleVec, from = min(sampleVec), to = max(sampleVec)) 
    #densInfo = density(codaObject[,c(parName)][[cIdx]]) 
    xMat <- cbind(xMat,densInfo$x)
    yMat <- cbind(yMat,densInfo$y)
    hdiLims <- cbind(hdiLims,HDIofMCMC(codaObject[,c(parName)][[cIdx]]))
  }
  matplot( xMat, yMat, type="l", col=plColors, 
           main="", xlab="Param. Value", ylab="Density" )
  abline(h=0)
  points( hdiLims[1,], rep(0,nChain), col=plColors, pch="|" )
  points( hdiLims[2,], rep(0,nChain), col=plColors, pch="|" )
  text( mean(hdiLims), 0, "95% HDI", adj=c(0.5,-0.2) )
  EffChnLngth <- effectiveSize(codaObject[,c(parName)])
  MCSE <- sd(as.matrix(codaObject[,c(parName)]))/sqrt(EffChnLngth) 
  text( max(xMat), max(yMat), adj=c(1.0,1.0), cex=1.25,
        paste("MCSE =\n",signif(MCSE,3)) )
}

diagMCMC <- function( codaObject, parName=varnames(codaObject)[1],
                       saveName=NULL, saveType="png" ) {
  DBDAplColors <- c("skyblue","black","royalblue","steelblue")
  openGraph(height=10,width=14)
  par( mar=0.5+c(3,4,1,0), oma=0.1+c(0,0,2,0), mgp=c(2.25,0.7,0), 
       cex.lab=1.5 )
  layout(matrix(1:4,nrow=2))
  
  # traceplot and gelman.plot are from CODA package:
  
  coda::traceplot( codaObject[,c(parName)], main="", ylab="Param. Value",
                   col=DBDAplColors ) 
  tryVal <- try(
    coda::gelman.plot( codaObject[,c(parName)], main="", auto.layout=FALSE, 
                       col=DBDAplColors )
  )  
  # if it runs, gelman.plot returns a list with finite shrink values:
  if ( class(tryVal)=="try-error" ) {
    plot.new() 
    print(paste0("Warning: coda::gelman.plot fails for ",parName))
  } else { 
    if ( class(tryVal)=="list" & !is.finite(tryVal$shrink[1]) ) {
      plot.new() 
      print(paste0("Warning: coda::gelman.plot fails for ",parName))
    }
  }
  DbdaAcfPlot(codaObject,parName,plColors=DBDAplColors)
  DbdaDensPlot(codaObject,parName,plColors=DBDAplColors)
  titleText <- switch(parName,
                      "mu" = expression(mu),
                      "sigma" = expression(sigma),
                      parName)  # default to parName if not matched
  
  mtext(text = titleText, outer = TRUE, adj = c(0.5, 0.5), cex = 2.0)
  if ( !is.null(saveName) ) {
    saveGraph( file=paste0(saveName,"診斷_",parName), type=saveType)
  }
  dev.off()
}

#------------------------------------------------------------------------------
# Functions for summarizing and plotting distribution of a large sample; 
# typically applied to MCMC posterior.

normalize <- function( v ){ return( v / sum(v) ) }

 # loaded by rjags, but redundancy doesn't hurt

summarizePost <- function( paramSampleVec, ROPE=NULL, credMass=0.95 ) {
  
  meanParam <- mean( paramSampleVec )
  medianParam <- median( paramSampleVec )
  dres <- density( paramSampleVec )
  modeParam <- dres$x[which.max(dres$y)]
  mcmcEffSz <- round( effectiveSize( paramSampleVec ), 1 )
  names(mcmcEffSz) <- NULL
  hdiLim <- HDIofMCMC( paramSampleVec, credMass=credMass )

  result <- NA
  if ( !is.null(ROPE) ) {
    pcltRope <- ( sum( paramSampleVec < ROPE[1] ) 
                 / length( paramSampleVec ) )
    pcgtRope <- ( sum( paramSampleVec > ROPE[2] ) 
                 / length( paramSampleVec ) )
    pcinRope <- 1-(pcltRope+pcgtRope) # 改成0-1而非0-100%
    
    if (pcinRope>=0.95) {
      result <- "="
    } else if (pcgtRope>=0.95) {
      result <- ">"
    } else if (pcltRope>=0.95) {
      result <- "<"
    } else {
      result <- "No Decision"
    }
    
  } else { 
    ROPE <- c(NA,NA)
    pcltRope <- NA 
    pcgtRope <- NA 
    pcinRope <- NA 
  }  
  
  
  return( c( Mean=meanParam, Median=medianParam, Mode=modeParam, 
             ESS=mcmcEffSz,
             HDImass=credMass, HDIlow=hdiLim[1], HDIhigh=hdiLim[2], 
             ROPElow=ROPE[1], ROPEhigh=ROPE[2],
             PcntLtROPE=pcltRope, PcntInROPE=pcinRope, PcntGtROPE=pcgtRope, Result=result ) )
}

#------------------------------------------------------------------------------

plotPostMulti <- function( paramSampleVec, cenTend=c("mode","median","mean")[1], 
                          ROPE=NULL, credMass=0.95, HDItextPlace=0.7, 
                          xlab=NULL, xlim=NULL, yaxt=NULL, ylab=NULL, 
                          main=NULL, cex=NULL, cex.lab=NULL,
                          col=NULL, border=NULL, 
                          breaks=NULL, 
                          ... ) {
  # Override defaults of hist function, if not specified by user:
  # (additional arguments "..." are passed to the hist function)
  if ( is.null(xlim) ) xlim <- range( c( ROPE, paramSampleVec ) )
  if ( is.null(cex.lab) ) cex.lab <- 1.5
  if ( is.null(cex) ) cex <- 1.4
  if ( is.null(main) ) main <- ""
  if ( is.null(yaxt) ) yaxt <- "n"
  if ( is.null(ylab) ) ylab <- ""
  if ( is.null(col) ) col <- "#a2c3da"
  if ( is.null(border) ) border <- "white"
  
  # convert coda object to matrix:
  if ( inherits( paramSampleVec, "mcmc.list") )  {
    paramSampleVec <- as.matrix(paramSampleVec)
  }
  
  
  summaryColNames <- c("ESS","mean","median","mode",
                      "hdiMass","hdiLow","hdiHigh",
                      "ROPElow","ROPEhigh","pLtROPE","pInROPE","pGtROPE")
  postSummary <- matrix( NA, nrow=1, ncol=length(summaryColNames), 
                        dimnames=list( c( xlab ), summaryColNames ) )
  
  postSummary[,"ESS"] <- effectiveSize(paramSampleVec)
  
  postSummary[,"mean"] <- mean(paramSampleVec)
  postSummary[,"median"] <- median(paramSampleVec)
  mcmcDensity <- density(paramSampleVec)
  postSummary[,"mode"] <- mcmcDensity$x[which.max(mcmcDensity$y)]
  
  HDI <- HDIofMCMC( paramSampleVec, credMass )
  postSummary[,"hdiMass"]<-credMass
  postSummary[,"hdiLow"]<-HDI[1]
  postSummary[,"hdiHigh"]<-HDI[2]
  
  # Plot histogram.
  ropeCol <- "steelblue"
  if ( is.null(breaks) ) {
    if ( max(paramSampleVec) > min(paramSampleVec) ) {
      breaks <- c( seq( from=min(paramSampleVec), to=max(paramSampleVec),
                       by=(HDI[2]-HDI[1])/18 ), max(paramSampleVec) )
    } else {
      breaks<-c(min(paramSampleVec)-1.0E-6,max(paramSampleVec)+1.0E-6)
      border<-"skyblue"
    }
  }
  par(xpd=NA)
  histinfo <- hist( paramSampleVec, xlab=xlab, yaxt=yaxt, ylab=ylab,
                    freq=F, border=border, col=col,
                    xlim=xlim, main=main, cex=cex, cex.lab=cex.lab,
                    breaks=breaks, ... )
  
  cenTendHt <- 0.9*max(histinfo$density)
  ROPEtextHt <- 1.1*max(histinfo$density)
  
  # Display central tendency:
  mn <- mean(paramSampleVec)
  med <- median(paramSampleVec)
  mcmcDensity <- density(paramSampleVec)
  mo <- mcmcDensity$x[which.max(mcmcDensity$y)]
  if ( cenTend=="mode" ){ 
    text( mo, cenTendHt,
          bquote(mode==.(signif(mo,3))), adj=c(.5,0), cex=cex )
  }
  if ( cenTend=="median" ){ 
    text( med, cenTendHt,
          bquote(median==.(signif(med,3))), adj=c(.5,0), cex=cex 
          #, col=cvCol 
          )
  }
  if ( cenTend=="mean" ){ 
    text( mn, cenTendHt,
          bquote(mean==.(signif(mn,3))), adj=c(.5,0), cex=cex )
  }
  
  # Display the ROPE.
  if ( !is.null( ROPE ) ) {
    pInROPE <- ( sum( paramSampleVec > ROPE[1] & paramSampleVec < ROPE[2] )
                / length( paramSampleVec ) )
    pGtROPE <- ( sum( paramSampleVec >= ROPE[2] ) / length( paramSampleVec ) )
    pLtROPE <- ( sum( paramSampleVec <= ROPE[1] ) / length( paramSampleVec ) )
    lines( c(ROPE[1],ROPE[1]), c(0.96*ROPEtextHt,0), lty="dotted", lwd=3,
           col=ropeCol )
    lines( c(ROPE[2],ROPE[2]), c(0.96*ROPEtextHt,0), lty="dotted", lwd=3,
           col=ropeCol)
    xlimrange <- xlim[2]-xlim[1]
    
    postSummary[,"ROPElow"]<-ROPE[1] 
    postSummary[,"ROPEhigh"]<-ROPE[2] 
    postSummary[,"pLtROPE"]<-pLtROPE
    postSummary[,"pInROPE"]<-pInROPE
    postSummary[,"pGtROPE"]<-pGtROPE
  }
  
  # Display the HDI.
  lines( HDI, c(0,0), lwd=4, lend=1 )
  text( mean(HDI), 0, bquote(.(100*credMass) * "% HDI" ),
        adj=c(.5,-1.7), cex=cex )
  text( HDI[1], 0, bquote(.(signif(HDI[1],3))),
        adj=c(HDItextPlace,-0.5), cex=cex )
  text( HDI[2], 0, bquote(.(signif(HDI[2],3))),
        adj=c(1.0-HDItextPlace,-0.5), cex=cex )
  par(xpd=F)
  #
  return( postSummary )
}
#------------------------------------------------------------------------------




#########################################################################################
gammaShRaFromModeSD <- function( mode , sd ) {
  if ( mode <=0 ) stop("mode must be > 0")
  if ( sd <=0 ) stop("sd must be > 0")
  rate <- ( mode + sqrt( mode^2 + 4 * sd^2 ) ) / ( 2 * sd^2 )
  shape <- 1 + mode * rate
  return( list( shape=shape , rate=rate ) )
}
