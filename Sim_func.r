source("Script0.R")

Sys.setlocale("LC_ALL", "en_US.UTF-8")
options(warn=-1)

yName='value'
xName='group'
graphFileType='png'

#################################################################################

genMCMC <- function(datFrm, yName = "value", xName = "group",
                     prior_params, rho_0, saveName = NULL,
                     numSavedSteps = 50000, adaptSteps = 500, burnInSteps = 1000, thinSteps = 1,
                     nChains = 4) {
  
  y1 <- datFrm[datFrm[[xName]] == 1, yName]
  y2 <- datFrm[datFrm[[xName]] == 2, yName]
  f <- 10
  diff_y <- y2 - y1
  nte <- 1
  ntr <- f - 1
  rho <- nte / (nte + ntr)
  
  dataList <- list(
    y = diff_y,
    n = length(diff_y),
    rho = rho,
    rho2 = ifelse(rho_0,0,rho),
    mu0 = 0,
    k0 = prior_params$k0,
    a = prior_params$a,
    b = prior_params$b
  )
  
  stanModel <- "
    data {
      int<lower=0> n;
      vector[n] y;
      real mu0;
      real k0;
      real a;
      real b;
      real rho;
      real rho2;
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
                Cov_mat[(blocki-1)*10 + i, (blockj-1)*10 + j] = sigma2 * rho2;
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
    
      for (i in 1:n) {
        mu_vec[i] = mu;
      }
    
      y ~ multi_normal(mu_vec, Cov_mat);
      mu ~ normal(mu0, sigma*sqrt(k0));
      sigma2 ~ inv_gamma(a, b);
    }
  "
  
  stan_model <- stan_model(model_code = stanModel)
  
  initsList <- list(
    list(mu = 0, sigma2 = var(diff_y)),
    list(mu = 0, sigma2 = var(diff_y)),
    list(mu = 0, sigma2 = var(diff_y)),
    list(mu = 0, sigma2 = var(diff_y))
  )
  
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
  
  codaSamples <- extract(fit)
  
  if (!is.null(saveName)) {
    save(codaSamples, file = paste0(saveName, "mcmc.Rdata"))
  }
  
  return(codaSamples)
}

#################################################################################


smryMCMC <- function(  codaSamples, RopeMuDiff=NULL, saveName=NULL ) {
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

#################################################################################


plotDiffMu <- function( codaSamples, RopeMuDiff=NULL, data=data ) {
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

#################################################################################

#################################################################################

test1 <- function(dataframe){
  n_data <- length(unique(dataframe$dataset))
  test1_decision <- data.frame(
    dataset = integer(),
    p_value_less = numeric(),
    p_value_equal = numeric(),
    p_value_large = numeric(),
    decision = character(),
    stringsAsFactors = FALSE
  )
  for ( i in 1:n_data ){
    data <- dataframe[dataframe$dataset==i,]
    
    y1 <- data[data$group==1,]$value
    y2 <- data[data$group==2,]$value
    (data[data$group==1,]$index==data[data$group==2,]$index)
    
    f <- 10
    r <- length(unique(data$index))
    n <- f*r
    nte <- 1
    ntr <- f-1
    x <- y2-y1
    
    x_bar <- mean(x)
    s <- sd(x)
    rho <- nte / (nte + ntr)
    se <- sqrt( s^2 * (1/n + rho/(1 - rho)) )
    
    t_stat <- x_bar / se
    df <- n - 1
    p_value_equal <- round(2 * (1 - pt(abs(t_stat), df=df, lower.tail=T)),4)
    p_value_large <- round(pt(t_stat, df=df, lower.tail=T),4)
    p_value_less <- round(1 - pt(t_stat, df=df, lower.tail=T),4)
    
    if (p_value_less<0.05){
      decision <- ">"
    }
    else if (p_value_large<0.05){
      decision <- "<"
    }
    else {decision <- "not reject ="}
    
    test1_decision <- rbind(test1_decision,data.frame(
      dataset = i,
      p_value_less = p_value_less,
      p_value_equal = p_value_equal,
      p_value_large = p_value_large,
      decision = decision,
      stringsAsFactors = FALSE
    ))
  }
  return(test1_decision)
}

# improper prior
correlatedBayesianTtest <- function(diff_a_b,rho,rope_min,rope_max){
  delta <- mean(diff_a_b)
  n <- length(diff_a_b)
  df <- n-1
  stdX <- sd(diff_a_b)
  sp <- sd(diff_a_b)*sqrt(1/n + rho/(1-rho))
  p.left <- pt((rope_min - delta)/sp, df)
  p.rope <- pt((rope_max - delta)/sp, df)-p.left
  results <- list('left'=p.left,'rope'=p.rope,'right'=1-p.left-p.rope)
  return (results)
}
test2 <- function(dataframe,RopeMuDiff){
  n_data <- length(unique(dataframe$dataset))
  
  test2_decision <- data.frame(
    dataset = integer(),
    p_left = numeric(),
    p_rope = numeric(),
    p_right = numeric(),
    decision = character(),
    stringsAsFactors = FALSE
  )
  for ( i in 1:n_data ){
    data <- dataframe[dataframe$dataset==i,]
    a <- data[data$group==1,]$value
    b <- data[data$group==2,]$value
    diff_a_b <- b-a
    rho <- 1/10
    rope_min <- RopeMuDiff[1]
    rope_max <- RopeMuDiff[2]
    result <- correlatedBayesianTtest(diff_a_b,rho,rope_min,rope_max)
    left <-  round(result$left,4)
    rope  <-  round(result$rope,4)
    right  <-  round(result$right,4)
    if(left>0.95){
      decision <- "<" 
    } else if (rope>0.95){
      decision <- "="
    } else if (right>0.95){
      decision <- ">"
    } else {decision <- "No decision"}
    
    test2_decision <- rbind(test2_decision,data.frame(
      dataset = i,
      p_left = left,
      p_rope = rope,
      p_right = right,
      decision = decision,
      stringsAsFactors = FALSE
    ))
  }
  return(test2_decision)
}

# proper prior
test2_2 <- function(dataframe,RopeMuDiff){
  n_data <- length(unique(dataframe$dataset))
  prior_params <- list(k0 = 1000, a = 10, b = 10)
  rho_0 <- FALSE
  
  png(filename = paste0(plotRoot,fileName,"plot.png"), width = 12, height = 6, units = "in", res = 300)
  layout(matrix(1:n_data, nrow = n_data/4, ncol = 4, byrow = TRUE))
  par(mar = c(2, 2, 2, 2), oma = c(1, 1, 1, 1), mgp = c(2, 0.5, 0))
  
  test2_decision <- data.frame(
    dataset = integer(),
    p_left = numeric(),
    p_rope = numeric(),
    p_right = numeric(),
    decision = character(),
    stringsAsFactors = FALSE
  )
  for ( i in 1:n_data ){
    data <- dataframe[dataframe$dataset==i,]
    mcmcCoda  <-  genMCMC( datFrm=data, yName=yName, xName=xName, 
                           prior_params=prior_params,rho_0=rho_0,
                            numSavedSteps=numSavedSteps,
                            saveName=paste0(fileRoot,"mcmcCoda/test_2_2_dataset",i,"_") )
    histInfo <- plotDiffMu(mcmcCoda,RopeMuDiff = RopeMuDiff, data = i)
    summaryInfo  <-  smryMCMC( mcmcCoda, 
                                RopeMuDiff=RopeMuDiff,
                                saveName=paste0(fileRoot,"summaryInfo/test_2_2_dataset",i,"_"))
    result <- summaryInfo["mu","Result"]
    pcltRope <- round(as.numeric(summaryInfo["mu","PcntLtROPE"]),4)
    pcinRope <- round(as.numeric(summaryInfo["mu","PcntInROPE"]),4)
    pcgtRope <- round(as.numeric(summaryInfo["mu","PcntGtROPE"]),4)
    HDI <- paste0("(",round(as.numeric(summaryInfo["mu","HDIlow"]),4),",",round(as.numeric(summaryInfo["mu","HDIhigh"]),4),")")
    
    test2_decision <- rbind(test2_decision,data.frame(
      dataset = i,
      p_left = pcltRope,
      p_rope = pcinRope,
      p_right = pcgtRope,
      decision = result,
      stringsAsFactors = FALSE
    ))
  }
  dev.off()
  return(test2_decision)
}

test3 <- function(dataframe,RopeMuDiff){
  n_data <- length(unique(dataframe$dataset))
  prior_params <- list(k0 = 1000, a = 10, b = 10)
  rho_0 <-  TRUE
  
  png(filename = paste0(plotRoot,fileName,"plot.png"), width = 12, height = 6, units = "in", res = 300)
  layout(matrix(1:n_data, nrow = 4, ncol = 4, byrow = TRUE))
  par(mar = c(2, 2, 2, 2), oma = c(1, 1, 1, 1), mgp = c(2, 0.5, 0))
  
  test3_decision <- data.frame(
    dataset = integer(),
    p_left = numeric(),
    p_rope = numeric(),
    p_right = numeric(),
    decision = character(),
    stringsAsFactors = FALSE
  )
  for ( i in 1:n_data ){
    data <- dataframe[dataframe$dataset==i,]
    
    mcmcCoda  <-  genMCMC( datFrm=data, yName=yName, xName=xName,
                            prior_params=prior_params,rho_0 =rho_0,
                            numSavedSteps=numSavedSteps,
                            saveName=paste0(fileRoot,"mcmcCoda/test_3_dataset",i,"_") )
    histInfo <- plotDiffMu(mcmcCoda,RopeMuDiff = RopeMuDiff, data = i)
    summaryInfo  <-  smryMCMC( mcmcCoda, 
                                RopeMuDiff=RopeMuDiff,
                                saveName=paste0(fileRoot,"summaryInfo/test_3_dataset",i,"_"))
    result <- summaryInfo["mu","Result"]
    pcltRope <- round(as.numeric(summaryInfo["mu","PcntLtROPE"]),4)
    pcinRope <- round(as.numeric(summaryInfo["mu","PcntInROPE"]),4)
    pcgtRope <- round(as.numeric(summaryInfo["mu","PcntGtROPE"]),4)
    HDI <- paste0("(",round(as.numeric(summaryInfo["mu","HDIlow"]),4),",",round(as.numeric(summaryInfo["mu","HDIhigh"]),4),")")
    
    test3_decision <- rbind(test3_decision,data.frame(
      dataset = i,
      p_left = pcltRope,
      p_rope = pcinRope,
      p_right = pcgtRope,
      decision = result,
      stringsAsFactors = FALSE
    ))
  }
  dev.off()
  return(test3_decision)
}






###################################################################################
# Simulation Data

##############################
cube_matrix <- function(f,r,sigma,rho,rho2){
  n <- f*r
  Cov_mat <- matrix(rep(0,n), nrow = n, ncol = n)
  for (blocki in 1:r) {
    for (blockj in 1:r) {
      for (i in 1:f) {
        for (j in 1:f) {
          if (blocki != blockj) {
            Cov_mat[(blocki-1)*f + i, (blockj-1)*f + j]  <-  sigma^2 * rho2;
          } else if (i == j) {
            Cov_mat[(blocki-1)*f + i, (blockj-1)*f + j]  <-  sigma^2;
          } else {
            Cov_mat[(blocki-1)*f + i, (blockj-1)*f + j]  <-  sigma^2 * rho;
          }
        }
      }
    }
  }
  return(Cov_mat)
}

generate_data <- function(params,rho_0){
  distr <- as.character(params[,'distr'])
  mu1 <- params[,'mu1']
  mu2 <- params[,'mu2']
  sigma1 <- params[,'sigma1']
  sigma2 <- params[,'sigma2']
  r <- params[,'r']
  f <- 10
  n <- r*f
  rho <- 1/f
  rho2 <- ifelse(rho_0,0,rho)
  
  if (sigma2 < sigma1){
    mud <- mu1-mu2
    sigmad <- sqrt(sigma1^2-sigma2^2)
    
    mu2_vec <- rep(mu2,n)
    mud_vec <- rep(mud,n)
    sigma2_mat <- cube_matrix(f,r,sigma2,rho,rho2)
    sigmad_mat <- cube_matrix(f,r,sigmad,rho,rho2)
    
    y2 <- do.call(distr,list(n = 1, mean = mu2_vec, sigma = sigma2_mat))
    d <- do.call(distr,list(n = 1, mean = mud_vec, sigma = sigmad_mat))
    y1 <- y2+d
  } else if (sigma2 == sigma1) {
    mu0 <- (mu1+mu2)/2
    mud1 <- mu1-mu0
    mud2 <- mu2-mu0
    sigma0 <- 0.0005
    sigmad <- sqrt(sigma1^2-sigma0^2)
    
    mu0_vec <- rep(mu0,n)
    mud1_vec <- rep(mud1,n)
    mud2_vec <- rep(mud2,n)
    sigma0_mat <- cube_matrix(f,r,sigma0,rho,rho2)
    sigmad_mat <- cube_matrix(f,r,sigmad,rho,rho2)
    
    y0 <- do.call(distr,list(n = 1, mean = mu0_vec, sigma = sigma0_mat))
    d1 <- do.call(distr,list(n = 1, mean = mud1_vec, sigma = sigmad_mat))
    d2 <- do.call(distr,list(n = 1, mean = mud2_vec, sigma = sigmad_mat))
    
    y1 <- y0+d1
    y2 <- y0+d2
  } else {
    mud <- mu2-mu1
    sigmad <- sqrt(sigma2^2-sigma1^2)
    
    mu1_vec <- rep(mu1,n)
    mud_vec <- rep(mud,n)
    sigma1_mat <- cube_matrix(f,r,sigma1,rho,rho2)
    sigmad_mat <- cube_matrix(f,r,sigmad,rho,rho2)
    
    y1 <- do.call(distr,list(n = 1, mean = mu1_vec, sigma = sigma1_mat))
    d <- do.call(distr,list(n = 1, mean = mud_vec, sigma = sigmad_mat))
    y2 <- y1+d
  }
  output <- list(y1=y1, y2=y2)
  return(output)
}

checking <- function(dataframe){
  n_dataset <- length(unique(dataframe$dataset))
  sry <- data.frame(
    dataset = integer(),
    mean1 = numeric(),
    mean2 = numeric(),
    sd1 = numeric(),
    sd2 = numeric(),
    stringsAsFactors = FALSE
  )
  for (i in 1:n_dataset){
    data <- subset(dataframe,dataset==i)
    y1 <- subset(data,group==1)
    y2 <- subset(data,group==2)
    mean1 <- round(mean((y1$value)/100),4)
    mean2 <- round(mean((y2$value)/100),4)
    sd1 <- round(sd((y1$value)/100),4)
    sd2 <- round(sd((y2$value)/100),4)
    sry <- rbind(sry,data.frame(
      dataset = i,
      mean1 = mean1,
      mean2 = mean2,
      sd1 = sd1,
      sd2 = sd2,
      stringsAsFactors = FALSE
    ))
  }
  return(sry)
}

##############################

r_options <- c(3,10)
mu1_options <- c(0.5)
mu2_options <- c(0.5, 0.53)
sigma1_options <- c(0.01, 0.05)
sigma2_options <- c(0.01, 0.05)
distr_options <- c('rmvnorm')

params_list <- expand.grid(
  mu1 = mu1_options,
  mu2 = mu2_options,
  sigma1 = sigma1_options,
  sigma2 = sigma2_options,
  distr = distr_options,
  r = r_options
)

##############################

SimulationData_0 <- data.frame(
  dataset = integer(),
  group = integer(),
  index = integer(),
  value = numeric(),
  stringsAsFactors = FALSE
)
n_data <- nrow(params_list)

set.seed(26)
for (i in 1:n_data) {
  params <- params_list[i,]
  data <- generate_data(params,rho_0=T)
  y1 <- data$y1
  y2 <- data$y2
  
  r <- params$r
  
  SimulationData_0 <- rbind(SimulationData_0, data.frame(
    dataset = i,
    group = rep(1, r*10),
    index = 1:(r*10),
    value = c(y1*100)
  ))
  
  SimulationData_0 <- rbind(SimulationData_0, data.frame(
    dataset = i,
    group = rep(2, r*10),
    index = 1:(r*10),
    value = c(y2*100)
  ))
}

checking(SimulationData_0)



##############################

SimulationData_rho <- data.frame(
  dataset = integer(),
  group = integer(),
  index = integer(),
  value = numeric(),
  stringsAsFactors = FALSE
)
n_data <- nrow(params_list)

set.seed(26)
for (i in 1:n_data) {
  params <- params_list[i,]
  data <- generate_data(params,rho_0=F)
  y1 <- data$y1
  y2 <- data$y2
  
  r <- params$r
  
  SimulationData_rho <- rbind(SimulationData_rho, data.frame(
    dataset = i,
    group = rep(1, r*10),
    index = 1:(r*10),
    value = c(y1*100)
  ))
  
  SimulationData_rho <- rbind(SimulationData_rho, data.frame(
    dataset = i,
    group = rep(2, r*10),
    index = 1:(r*10),
    value = c(y2*100)
  ))
}

checking(SimulationData_rho)


