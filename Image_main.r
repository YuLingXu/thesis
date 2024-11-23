graphics.off()
rm(list=ls())

getwd()
setwd('D:/200/220NCCU/222Thesis/coding_new')

source("Script_ImageData.R")

myDataFrame <- ImageData

sum(myDataFrame$value<0)
sum(myDataFrame$value>100)

numSavedSteps <- 20000

fileRoot <- "5_2影像/"


#########################################################################################################

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
    
    y <- as.numeric(data[,yName])
    x <- as.numeric(as.factor(data[,xName]))
    y1 <- y[x==1]
    y2 <- y[x==2]
    
    f <- 10
    r <- length(unique(data$Index))
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
    
    if (p_value_less<0.05 & p_value_equal<0.05){
      decision <- ">"
    }
    else if (p_value_large<0.05 & p_value_equal<0.05){
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

test2 <- function(dataframe){
  n_data <- length(unique(dataframe$dataset))
  RopeMuDiff = c(-1,1)
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
    
    y <- as.numeric(data[,yName])
    x <- as.numeric(as.factor(data[,xName]))
    y1 <- y[x==1]
    y2 <- y[x==2]

    diff_a_b <- y2-y1
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

test2_2 <- function(dataframe){
  n_data <- length(unique(dataframe$dataset))
  RopeMuDiff = c(-1,1)
  
  png(filename = paste0(fileNameRootRoot,"後驗.png"), width = 12, height = 6, units = "in", res = 300)
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
    mcmcCoda  <-  genMCMC2( datFrm=data, yName=yName, xName=xName,
                            numSavedSteps=numSavedSteps,
                            saveName=paste0(fileRoot,"mcmcCoda/",fileNameRootRoot,"dataset",i,"_") )
    histInfo <- plotDiffMu3(mcmcCoda,RopeMuDiff = RopeMuDiff, data = i)
    summaryInfo  <-  smryMCMC3( mcmcCoda, 
                                RopeMuDiff=RopeMuDiff,
                                saveName=paste0(fileRoot,"summaryInfo/",fileNameRootRoot,"dataset",i,"_"))
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

test3 <- function(dataframe){
  n_data <- length(unique(dataframe$dataset))
  RopeMuDiff = c(-1,1)
  
  png(filename = paste0(fileNameRootRoot,"後驗.png"), width = 12, height = 6, units = "in", res = 300)
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
    mcmcCoda  <-  genMCMC3( datFrm=data, yName=yName, xName=xName,
                            numSavedSteps=numSavedSteps,
                            saveName=paste0(fileRoot,"mcmcCoda/",fileNameRootRoot,"dataset",i,"_") )
    histInfo <- plotDiffMu3(mcmcCoda,RopeMuDiff = RopeMuDiff, data = i)
    summaryInfo  <-  smryMCMC3( mcmcCoda, 
                                RopeMuDiff=RopeMuDiff,
                                saveName=paste0(fileRoot,"summaryInfo/",fileNameRootRoot,"dataset",i,"_"))
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

##########################################################################################

test_for_contrasts <- function(test, myDataFrame, methods){
  contrasts <- t(combn(methods, 2))
  n_contrast <- dim(contrasts)[1]
  
  for (i in 1:n_contrast){
    
    method1 <- contrasts[i,1]
    method2 <- contrasts[i,2]
    contrast <- paste0("mu",method2,"_mu",method1,"_")
    fileNameRootRoot <- paste0(fileName,contrast)
    
    dataframe <- subset(myDataFrame,method==method1|method==method2)
    params <- list(dataframe = dataframe)
    
    test_decision <- do.call(test, params)
    
    write.csv(test_decision, 
              paste0('D:/200/220NCCU/222Thesis/coding_new',fileRoot,fileNameRootRoot,"結果",".csv"), 
              row.names = FALSE, fileEncoding = "UTF-8")
  }
}

##########################################################################################

normality_test <-  function(dataframe){
  normality_results <- data.frame(
    dataset = integer(),
    method = integer(),
    p_value = numeric(),
    is_normal = character(),
    stringsAsFactors = FALSE
  )
  for (i in unique(dataframe$dataset)) {
    data_subset <- subset(dataframe, dataset == i)
    
    for (method in unique(data_subset$method)) {
      method_data <- data_subset[data_subset$method == method, "accuracy"]
      shapiro_test <- shapiro.test(method_data)
      normality_results <- rbind(normality_results, data.frame(
        dataset = i,
        method = method,
        p_value = shapiro_test$p.value,
        is_normal = ifelse(shapiro_test$p.value > 0.05, "Yes", "No"),
        stringsAsFactors = FALSE
      ))
    }
  }
  normality_summary <- dcast(normality_results, dataset ~ method, value.var = "is_normal")
  print(normality_results)
  print(normality_summary)
  return(normality_summary)
}

homogeneity_test <- function(dataframe){
  homogeneity_results <- data.frame(
    dataset = integer(),
    p_value = numeric(),
    is_homogeneous = character(),
    stringsAsFactors = FALSE
  )
  
  for (i in unique(dataframe$dataset)) {
    data_subset <- subset(dataframe, dataset == i)
    bartlett_test <- bartlett.test(accuracy ~ as.factor(method), data = data_subset)
    
    p_value <- bartlett_test$p.value
    is_homogeneous <- ifelse(p_value > 0.05, "Yes", "No")
    homogeneity_results <- rbind(homogeneity_results, data.frame(
      dataset = i,
      p_value = p_value,
      is_homogeneous = is_homogeneous,
      stringsAsFactors = FALSE
    ))
  }
  return(homogeneity_results)
}

##########################################################################################

fileName <- "5_2影像_0資料_"
fileNameRoot <- paste0(fileRoot,fileName)

plotdata <- myDataFrame
plotdata$dataset <- factor(plotdata$dataset)

method_colors <- c("#fdb462","#80b1d3","#b3de69","#ffed6f")

ggplot(plotdata, aes(x = factor(method), y = accuracy, fill = factor(method))) +
  geom_boxplot() +
  coord_flip() +
  labs(x = "Method",
       y = "Accuracy",
       fill = "Method") +
  theme(legend.position = "right",
        text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 16)) +
  scale_x_discrete(labels = NULL) +
  guides(fill = guide_legend(title = "Method",
                             title.position = "top",
                             title.hjust = 0.5,
                             ncol = 1)) +
  scale_fill_discrete(labels = c("Original", "Contrast Stretching", "Histogram Equalization", "Sharpening")) +
  scale_fill_manual(values = method_colors, 
                    labels = c("Original", "Contrast Stretching", "Histogram Equalization", "Sharpening")) +
  guides(fill = guide_legend(reverse = TRUE)) +
  facet_wrap(~ dataset, ncol = 2, dir = "v", switch = "both")

ggsave(filename = paste0(fileNameRoot, "box.", graphFileType), dpi = 300, width = 12, height = 10)


####################################################################################
fileName <- "5_2影像_0檢定_"
fileNameRoot <- paste0(fileRoot,fileName)

##
normality_test_result <- normality_test(myDataFrame)
write.csv(normality_test_result, file = paste0(fileNameRoot,"常態檢定.csv"), row.names = FALSE)

####################################################################################
fileName <- "5_2影像_1檢定1_"
fileNameRoot <- paste0(fileRoot,fileName)

methods <- c(1, 2, 3, 4)
test_for_contrasts(test1, myDataFrame, methods)

####################################################################################
fileName <- "5_2影像_2檢定2_"
fileNameRoot <- paste0(fileRoot,fileName)

methods <- c(1, 2, 3, 4)
test_for_contrasts(test2, myDataFrame, methods)

######
fileName <- "5_2影像_2檢定2_2_"
fileNameRoot <- paste0(fileRoot,fileName)

methods <- c(1, 2, 3, 4)
test_for_contrasts(test2_2, myDataFrame, methods)

######################################################################################
fileName <- "5_2影像_3檢定3_"
fileNameRoot <- paste0(fileRoot,fileName)

methods <- c(1, 2, 3, 4)
test_for_contrasts(test3, myDataFrame, methods)

######################################################################################
fileName <- "5_2影像_4比較_"
fileNameRoot <- paste0(fileRoot,fileName)




##################################################################################

##   end end end end end end end end end end end end end end end end end end    ##
##   end end end end end end end end end end end end end end end end end end    ##
##   end end end end end end end end end end end end end end end end end end    ##
##   end end end end end end end end end end end end end end end end end end    ##
##   end end end end end end end end end end end end end end end end end end    ##
##   end end end end end end end end end end end end end end end end end end    ##

##################################################################################
