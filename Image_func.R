

source("Script1.R")

Sys.setlocale("LC_ALL", "en_US.UTF-8")
options(warn=-1)

yName='accuracy'
xName='method'
graphFileType='png'

#===============================================================================

genMCMC2 <- function(datFrm,
                     yName = "y", xName = "x",
                     saveName = NULL,
                     numSavedSteps = 50000, adaptSteps = 500, burnInSteps = 1000, thinSteps = 1,
                     nChains = 4) {
  
  # Data preparation
  y <- as.numeric(datFrm[,yName])
  x <- as.numeric(as.factor(datFrm[,xName]))
  
  y1 <- y[x==1]
  y2 <- y[x==2]
  
  f <- 10
  diff_y <- y2 - y1
  nte <- 1
  ntr <- f - 1
  
  dataList <- list(
    y = diff_y,
    n = length(diff_y),
    rho = nte / (nte + ntr),
    mu0 = 0,
    k0 = 1000,
    a = 10,
    b = 10
  )
  
  # Stan model
  stanModel <- "
    data {
      int<lower=0> n;
      vector[n] y;
      real mu0;
      real k0;
      real a;
      real b;
      real rho;
    }
    
    parameters {
      real mu;
      real<lower=0> sigma2;
    }
    
    transformed parameters {
      real<lower=0> sigma;
      sigma = sqrt(sigma2);
    }
    
    model {
      matrix[n, n] Cov_mat;
      vector[n] mu_vec;
    
      // Construct the covariance matrix
      for (i in 1:n) {
        for (j in 1:n) {
          Cov_mat[i, j] = 0;
        }
      }
    
      for (blocki in 1:(n / 10)) {
        for (blockj in 1:(n / 10)) {
          for (i in 1:10) {
            for (j in 1:10) {
              if (blocki != blockj) {
                Cov_mat[(blocki-1)*10 + i, (blockj-1)*10 + j] = sigma2 * rho;
              } else {
                if (i == j) {
                  Cov_mat[(blocki-1)*10 + i, (blockj-1)*10 + j] = sigma2;
                } else {
                  Cov_mat[(blocki-1)*10 + i, (blockj-1)*10 + j] = sigma2 * rho;
                }
              }
            }
          }
        }
      }
    
      // Define the mean vector
      for (i in 1:n) {
        mu_vec[i] = mu;
      }
    
      // Likelihood
      y ~ multi_normal(mu_vec, Cov_mat);
    
      // Priors
      mu ~ normal(mu0, sigma*sqrt(k0));
      sigma2 ~ inv_gamma(a, b);
    }
  "
  
  # Compile the Stan model
  stan_model <- stan_model(model_code = stanModel)
  
  # Initial values for the chains
  initsList <- list(
    list(mu = 0, sigma2 = var(diff_y)),
    list(mu = 0, sigma2 = var(diff_y)),
    list(mu = 0, sigma2 = var(diff_y)),
    list(mu = 0, sigma2 = var(diff_y))
  )
  
  # Run Stan
  fit <- sampling(stan_model, 
                  data = dataList, 
                  iter = numSavedSteps + burnInSteps, 
                  warmup = burnInSteps, 
                  chains = nChains, 
                  thin = thinSteps, 
                  init = initsList,
                  cores = nChains,
                  seed = 26,
                  control = list(adapt_delta = 0.95))
  
  # Extract the samples
  codaSamples <- extract(fit)
  
  if (!is.null(saveName)) {
    save(codaSamples, file = paste0(saveName, "mcmc.Rdata"))
  }
  
  return(codaSamples)
}


genMCMC3 <- function(datFrm,
                     yName = "y", xName = "x",
                     saveName = NULL,
                     numSavedSteps = 50000, adaptSteps = 500, burnInSteps = 1000, thinSteps = 1,
                     nChains = 4) {
  
  # Data preparation
  y <- as.numeric(datFrm[,yName])
  x <- as.numeric(as.factor(datFrm[,xName]))
  
  y1 <- y[x==1]
  y2 <- y[x==2]
  
  f <- 10
  diff_y <- y2 - y1
  nte <- 1
  ntr <- f - 1
  
  dataList <- list(
    y = diff_y,
    n = length(diff_y),
    rho = nte / (nte + ntr),
    mu0 = 0,
    k0 = 1000,
    a = 10,
    b = 10
  )
  
  # Stan model
  stanModel <- "
    data {
      int<lower=0> n;
      vector[n] y;
      real mu0;
      real k0;
      real a;
      real b;
      real rho;
    }
    
    parameters {
      real mu;
      real<lower=0> sigma2;
    }
    
    transformed parameters {
      real<lower=0> sigma;
      sigma = sqrt(sigma2);
    }
    
    model {
      matrix[n, n] Cov_mat;
      vector[n] mu_vec;
    
      // Construct the covariance matrix
      for (i in 1:n) {
        for (j in 1:n) {
          Cov_mat[i, j] = 0;
        }
      }
    
      for (blocki in 1:(n / 10)) {
        for (blockj in 1:(n / 10)) {
          for (i in 1:10) {
            for (j in 1:10) {
              if (blocki != blockj) {
                Cov_mat[(blocki-1)*10 + i, (blockj-1)*10 + j] = 0;
              } else {
                if (i == j) {
                  Cov_mat[(blocki-1)*10 + i, (blockj-1)*10 + j] = sigma2;
                } else {
                  Cov_mat[(blocki-1)*10 + i, (blockj-1)*10 + j] = sigma2 * rho;
                }
              }
            }
          }
        }
      }
    
      // Define the mean vector
      for (i in 1:n) {
        mu_vec[i] = mu;
      }
    
      // Likelihood
      y ~ multi_normal(mu_vec, Cov_mat);
    
      // Priors
      mu ~ normal(mu0, sigma*sqrt(k0));
      sigma2 ~ inv_gamma(a, b);
    }
  "
  
  # Compile the Stan model
  stan_model <- stan_model(model_code = stanModel)
  
  # Initial values for the chains
  initsList <- list(
    list(mu = 0, sigma2 = var(diff_y)),
    list(mu = 0, sigma2 = var(diff_y)),
    list(mu = 0, sigma2 = var(diff_y)),
    list(mu = 0, sigma2 = var(diff_y))
  )
  
  # Run Stan
  fit <- sampling(stan_model, 
                  data = dataList, 
                  iter = numSavedSteps + burnInSteps, 
                  warmup = burnInSteps, 
                  chains = nChains, 
                  thin = thinSteps, 
                  init = initsList,
                  cores = nChains,
                  seed = 26,
                  control = list(adapt_delta = 0.95))
  
  # Extract the samples
  codaSamples <- extract(fit)
  
  if (!is.null(saveName)) {
    save(codaSamples, file = paste0(saveName, "mcmc.Rdata"))
  }
  
  return(codaSamples)
}






#===============================================================================

smryMCMC1 <- function(  codaSamples, RopeMuDiff=NULL, saveName=NULL ) {
  summaryInfo <- NULL
  mcmcMat <- as.matrix(codaSamples,chains=TRUE)
  summaryInfo <- rbind( summaryInfo, 
                        "mu" = summarizePost( mcmcMat[,"mu"],ROPE=RopeMuDiff ))
  summaryInfo <- rbind( summaryInfo, 
                        "sigma" = summarizePost( mcmcMat[,"sigma"] ) )
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo, file=paste(saveName,"SummaryInfo.csv",sep=""), fileEncoding = "UTF-8" )
  }
  return( summaryInfo )
}

#===============================================================================

smryMCMC3 <- function(  codaSamples, RopeMuDiff=NULL, saveName=NULL ) {
  summaryInfo <- NULL
  mcmcMat <- as.data.frame(codaSamples,chains=TRUE)
  summaryInfo <- rbind( summaryInfo, 
                        "mu" = summarizePost( mcmcMat[,"mu"],ROPE=RopeMuDiff ))
  summaryInfo <- rbind( summaryInfo, 
                        "sigma" = summarizePost( mcmcMat[,"sigma"] ) )
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo, file=paste(saveName,"SummaryInfo.csv",sep=""), fileEncoding = "UTF-8" )
  }
  return( summaryInfo )
}

#===============================================================================

plotMCMC1 <- function( codaSamples,
                       datFrm,
                       yName="y", xName="x",
                       RopeMuDiff=NULL,
                       showCurve=FALSE, 
                       pairsPlot=FALSE,
                       saveName=NULL, 
                       saveType="png") {
  #-----------------------------------------------------------------------------
  mcmcMat <- as.matrix(codaSamples,chains=TRUE)
  chainLength <- NROW( mcmcMat )
  mu <- mcmcMat[,"mu"]
  sigma <- mcmcMat[,"sigma"]
  #-----------------------------------------------------------------------------
  # mu
  openGraph(width=18.0,height=6.0)
  layout( matrix( c(1,2), nrow=1, byrow=TRUE ) )
  par( mar=c(5,3.5,3.5,3.5), mgp=c(2.25,0.7,0) )
  
  xlim <- range(mu)
  histInfo <- plotPost( mu,  
                        xlim=xlim, 
                        cex.lab = 1.5,
                        showCurve=showCurve,
                        xlab=bquote(mu), 
                        col="#a2c3da" )
  xlim <- range(sigma)
  histInfo <- plotPost( mu2,  
                        xlim=xlim, 
                        cex.lab = 1.5,
                        showCurve=showCurve,
                        xlab=bquote(sigma), 
                        col="#a2c3da" )
  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"參數後驗",sep=""), type=saveType)
  }
}


#===============================================================================

plotDiffMu <- function( codaSamples, RopeMuDiff=NULL, data=data ) {
  #-----------------------------------------------------------------------------
  mcmcMat <- as.matrix(codaSamples,chains=TRUE)
  chainLength <- NROW( mcmcMat )
  mu <- mcmcMat[,"mu"]
  #-----------------------------------------------------------------------------
  # Set up window and layout:
  histInfo  <-  plotPostMulti( mu, 
                               cex.lab = 1.5,
                               cex = 1.4,
                               ROPE=RopeMuDiff,
                               main=paste0("Dataset",data),
                               col="#a2c3da" )
  return(histInfo)
}


#===============================================================================

plotDiffMu3 <- function( codaSamples, RopeMuDiff=NULL, data=data ) {
  #-----------------------------------------------------------------------------
  mcmcMat <- as.data.frame(codaSamples,chains=TRUE)
  chainLength <- NROW( mcmcMat )
  mu <- mcmcMat[,"mu"]
  #-----------------------------------------------------------------------------
  # Set up window and layout:
  histInfo  <-  plotPostMulti( mu, 
                               cex.lab = 1.5,
                               cex = 1.4,
                               ROPE=RopeMuDiff,
                               main=paste0("Dataset",data),
                               col="#a2c3da" )
  return(histInfo)
}

###################################################################################
# Real Data

data_list <- list()
file_names <- list.files("D:/200/220NCCU/222Thesis/coding_new/data", 
                         pattern = ".csv", full.names = TRUE)

for (i in 1:length(file_names)) {
  data_list[[i]] <- read.csv(file_names[i], header = TRUE, row.names = 1,fileEncoding = "UTF-8") %>%
    setNames(gsub("Index\\.", "", names(.))) %>%
    tibble::rownames_to_column("method") %>%
    pivot_longer(
      cols = -method,
      names_to = "Index",
      values_to = "accuracy"
    ) %>%
    mutate(method = factor(method, levels = unique(method))) %>%
    mutate(method = as.numeric(as.integer(method)),dataset = i)
}

ImageData  <-  data.frame(bind_rows(data_list))
ImageData$accuracy  <- (ImageData$accuracy)*100

#######################################################

call_data <- function(dataframe, data_idx, mtd_idx){
  return(dataframe[(dataframe$dataset %in% data_idx) & (dataframe$method %in% mtd_idx),])
}

