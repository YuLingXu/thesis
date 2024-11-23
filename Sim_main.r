graphics.off()
rm(list=ls())

getwd()
setwd("C:/Users/User/Desktop/code")

source("Sim_func.R")

myDataFrame <- SimulationData_0
fileRoot <- "4_Sim/"

sum(myDataFrame$value<0)
sum(myDataFrame$value>100)

numSavedSteps <- 50000

plotRoot <- paste0(fileRoot,"result(plot)/")
csvRoot <- paste0(fileRoot,"result(csv)/")
mcmcCodaRoot <- paste0(fileRoot,"mcmcCoda/")
summaryInfoRoot <- paste0(fileRoot,"summaryInfo/")

##########################################################################################

fileName <- "4_Sim_0資料_"

group_colors <- c("#fdb462","#80b1d3")

plotdata <- myDataFrame
plotdata$dataset <- factor(plotdata$dataset)

plotdata <- plotdata %>%
  group_by(dataset) %>%
  mutate(
    iqr = IQR(value),
    max_value = max(value),
    min_value = min(value),
    adjusted_max = max_value + iqr * 1, # 擴展範圍的倍數可以根據需要調整
    adjusted_min = min_value
  )

ggplot(plotdata, aes(x = factor(group), y = value, fill = factor(group))) +
  geom_boxplot() +
  coord_flip() +
  labs(x = "group", y = "value", fill = "group") +
  theme(
    legend.position = "right", text = element_text(size = 16), axis.title = element_text(size = 16),
    axis.text = element_text(size = 14), strip.text = element_text(size = 16)
  ) +
  scale_x_discrete(labels = NULL) +
  guides(fill = guide_legend(title = "group", title.position = "top", title.hjust = 0.5, ncol = 1)) +
  scale_fill_manual(values = group_colors, labels = c("1", "2")) +
  guides(fill = guide_legend(reverse = TRUE)) +
  facet_wrap(~ dataset, ncol = 2, scales = "free_x", dir = "h", switch = "both") +
  #geom_text(data = dist_labels, aes(x = group, y = Inf, label = label), vjust = 0.5, hjust = 1.2, size = 5, position = position_dodge(width = 0.75)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  geom_blank(aes(y = adjusted_max)) + # 增加透明點以擴展y軸範圍
  geom_blank(aes(y = adjusted_min))   # 增加透明點以擴展y軸範圍

ggsave(filename = paste0(plotRoot, "box.", graphFileType), dpi = 300, width = 12, height = 13)

########################################################################

fileName <- "4_Sim_1檢定1_"

test1_decision <- test1(myDataFrame)
test1_decision

write.csv(test1_decision, 
          paste0(csvRoot,fileName,"結果",".csv"), 
          row.names = FALSE, fileEncoding = "UTF-8")

####################################################################################

fileName <- "4_Sim_2檢定2_"

RopeMuDiff <- c(-1,1)

test2_decision <- test2(myDataFrame,RopeMuDiff)
test2_decision

write.csv(test2_decision, 
          paste0(csvRoot,fileName,"結果",".csv"), 
          row.names = FALSE, fileEncoding = "UTF-8")

######

fileName <- "4_Sim_2檢定2_2_"

RopeMuDiff <- c(-1,1)

test2_2_decision <- test2_2(myDataFrame,RopeMuDiff)
test2_2_decision

write.csv(test2_2_decision, 
          paste0(csvRoot,fileName,"結果",".csv"), 
          row.names = FALSE, fileEncoding = "UTF-8")


######################################################################################

fileName <- "4_Sim_3檢定3_"

RopeMuDiff <- c(-1,1)

test3_decision <- test3(myDataFrame,RopeMuDiff)
test3_decision

write.csv(test3_decision,
          paste0(csvRoot,fileName,"結果",".csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

######################################################################################

fileName <- "4_Sim_4交叉熵_k10"

k0_options <- c(10)
a_options <- c(0.1, 1, 10)
b_options <- c(0.1, 1, 10)

numSavedSteps <- 10000

RopeMuDiff <- c(-1,1)

prior_params_list <- expand.grid(
  k0 = k0_options,
  a = a_options,
  b = b_options
)
prior_params_list <- apply(prior_params_list, 1, as.list)

compare_models <- function(myDataFrame, prior_params_list) {
  results <- data.frame(
    dataset = integer(),
    prior_params = character(),
    model = character(),
    p_left = numeric(),
    p_rope = numeric(),
    p_right = numeric(),
    decision = character(),
    stringsAsFactors = FALSE
  )
  
  mcmc_storage <- list()  # 用來存儲 mcmcCoda
  summary_storage <- list()  # 用來存儲 summaryInfo
  
  for (prior_params in prior_params_list) {
    for (i in 1:8) { #
      
      datFrm <- subset(myDataFrame, dataset == i)
      
      # 模型 5
      mcmc2_samples <- genMCMC( datFrm=datFrm, yName=yName, xName=xName,
                                numSavedSteps=numSavedSteps, 
                                prior_params = prior_params, rho_0=FALSE)
      summaryInfo2  <-  smryMCMC( mcmc2_samples, 
                                   RopeMuDiff=RopeMuDiff)
      result <- summaryInfo2["mu","Result"]
      pcltRope <- round(as.numeric(summaryInfo2["mu","PcntLtROPE"]),4)
      pcinRope <- round(as.numeric(summaryInfo2["mu","PcntInROPE"]),4)
      pcgtRope <- round(as.numeric(summaryInfo2["mu","PcntGtROPE"]),4)
      
      # 儲存 mcmcCoda 和 summaryInfo
      mcmc_storage[[paste0("test2_", i, "_", paste(prior_params, collapse = "_"))]] <- mcmc2_samples
      summary_storage[[paste0("test2_", i, "_", paste(prior_params, collapse = "_"))]] <- summaryInfo2
      
      results <- rbind(results, data.frame(
        dataset = i,
        prior_params = paste(prior_params, collapse = ","),
        model = "test2",
        p_left = pcltRope,
        p_rope = pcinRope,
        p_right = pcgtRope,
        decision = result,
        stringsAsFactors = FALSE
      ))
      
      # 模型 4
      mcmc3_samples <- genMCMC( datFrm=datFrm, yName=yName, xName=xName,
                                numSavedSteps=numSavedSteps, 
                                prior_params = prior_params, rho_0=TRUE)
      summaryInfo3  <-  smryMCMC( mcmc3_samples, 
                                   RopeMuDiff=RopeMuDiff)
      result <- summaryInfo3["mu","Result"]
      pcltRope <- round(as.numeric(summaryInfo3["mu","PcntLtROPE"]),4)
      pcinRope <- round(as.numeric(summaryInfo3["mu","PcntInROPE"]),4)
      pcgtRope <- round(as.numeric(summaryInfo3["mu","PcntGtROPE"]),4)
      
      # 儲存 mcmcCoda 和 summaryInfo
      mcmc_storage[[paste0("test3_", i, "_", paste(prior_params, collapse = "_"))]] <- mcmc3_samples
      summary_storage[[paste0("test3_", i, "_", paste(prior_params, collapse = "_"))]] <- summaryInfo3
      
      results <- rbind(results, data.frame(
        dataset = i,
        prior_params = paste(prior_params, collapse = ","),
        model = "test3",
        p_left = pcltRope,
        p_rope = pcinRope,
        p_right = pcgtRope,
        decision = result,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # 保存 mcmc_storage 和 summary_storage 到文件
  save(mcmc_storage, file = paste0(mcmcCodaRoot,fileName, "mcmc_storage.RData"))
  save(summary_storage, file = paste0(summaryInfoRoot,fileName, "summary_storage.RData"))
  
  return(results)
}

results <- compare_models(myDataFrame, prior_params_list)

print(results)

write.csv(results,
          paste0(csvRoot,fileName,"結果.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

######################################################################################

fileName <- "4_Sim_4交叉熵_k100"

k0_options <- c(100)
a_options <- c(0.1, 1, 10)
b_options <- c(0.1, 1, 10)

numSavedSteps <- 10000

RopeMuDiff <- c(-1,1)

prior_params_list <- expand.grid(
  k0 = k0_options,
  a = a_options,
  b = b_options
)
prior_params_list <- apply(prior_params_list, 1, as.list)

compare_models <- function(myDataFrame, prior_params_list) {
  results <- data.frame(
    dataset = integer(),
    prior_params = character(),
    model = character(),
    p_left = numeric(),
    p_rope = numeric(),
    p_right = numeric(),
    decision = character(),
    stringsAsFactors = FALSE
  )
  
  mcmc_storage <- list()  # 用來存儲 mcmcCoda
  summary_storage <- list()  # 用來存儲 summaryInfo
  
  for (prior_params in prior_params_list) {
    for (i in 1:8) { #
      
      datFrm <- subset(myDataFrame, dataset == i)
      
      # 模型 5
      mcmc2_samples <- genMCMC( datFrm=datFrm, yName=yName, xName=xName,
                                numSavedSteps=numSavedSteps, 
                                prior_params = prior_params, rho_0=FALSE)
      summaryInfo2  <-  smryMCMC( mcmc2_samples, 
                                  RopeMuDiff=RopeMuDiff)
      result <- summaryInfo2["mu","Result"]
      pcltRope <- round(as.numeric(summaryInfo2["mu","PcntLtROPE"]),4)
      pcinRope <- round(as.numeric(summaryInfo2["mu","PcntInROPE"]),4)
      pcgtRope <- round(as.numeric(summaryInfo2["mu","PcntGtROPE"]),4)
      
      # 儲存 mcmcCoda 和 summaryInfo
      mcmc_storage[[paste0("test2_", i, "_", paste(prior_params, collapse = "_"))]] <- mcmc2_samples
      summary_storage[[paste0("test2_", i, "_", paste(prior_params, collapse = "_"))]] <- summaryInfo2
      
      results <- rbind(results, data.frame(
        dataset = i,
        prior_params = paste(prior_params, collapse = ","),
        model = "test2",
        p_left = pcltRope,
        p_rope = pcinRope,
        p_right = pcgtRope,
        decision = result,
        stringsAsFactors = FALSE
      ))
      
      # 模型 4
      mcmc3_samples <- genMCMC( datFrm=datFrm, yName=yName, xName=xName,
                                numSavedSteps=numSavedSteps, 
                                prior_params = prior_params, rho_0=TRUE)
      summaryInfo3  <-  smryMCMC( mcmc3_samples, 
                                  RopeMuDiff=RopeMuDiff)
      result <- summaryInfo3["mu","Result"]
      pcltRope <- round(as.numeric(summaryInfo3["mu","PcntLtROPE"]),4)
      pcinRope <- round(as.numeric(summaryInfo3["mu","PcntInROPE"]),4)
      pcgtRope <- round(as.numeric(summaryInfo3["mu","PcntGtROPE"]),4)
      
      # 儲存 mcmcCoda 和 summaryInfo
      mcmc_storage[[paste0("test3_", i, "_", paste(prior_params, collapse = "_"))]] <- mcmc3_samples
      summary_storage[[paste0("test3_", i, "_", paste(prior_params, collapse = "_"))]] <- summaryInfo3
      
      results <- rbind(results, data.frame(
        dataset = i,
        prior_params = paste(prior_params, collapse = ","),
        model = "test3",
        p_left = pcltRope,
        p_rope = pcinRope,
        p_right = pcgtRope,
        decision = result,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # 保存 mcmc_storage 和 summary_storage 到文件
  save(mcmc_storage, file = paste0(mcmcCodaRoot,fileName, "mcmc_storage.RData"))
  save(summary_storage, file = paste0(summaryInfoRoot,fileName, "summary_storage.RData"))
  
  return(results)
}

results <- compare_models(myDataFrame, prior_params_list)

print(results)

write.csv(results,
          paste0(csvRoot,fileName,"結果.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

######################################################################################

fileName <- "4_Sim_4交叉熵_k1000"

k0_options <- c(1000)
a_options <- c(0.1, 1, 10)
b_options <- c(0.1, 1, 10)

numSavedSteps <- 10000

RopeMuDiff <- c(-1,1)

prior_params_list <- expand.grid(
  k0 = k0_options,
  a = a_options,
  b = b_options
)
prior_params_list <- apply(prior_params_list, 1, as.list)

compare_models <- function(myDataFrame, prior_params_list) {
  results <- data.frame(
    dataset = integer(),
    prior_params = character(),
    model = character(),
    p_left = numeric(),
    p_rope = numeric(),
    p_right = numeric(),
    decision = character(),
    stringsAsFactors = FALSE
  )
  
  mcmc_storage <- list()  # 用來存儲 mcmcCoda
  summary_storage <- list()  # 用來存儲 summaryInfo
  
  for (prior_params in prior_params_list) {
    
    for (i in 1:8) { #
      print(paste0("prior: ",prior_params,", dataset: ",i))
      datFrm <- subset(myDataFrame, dataset == i)
      
      # 模型 5
      mcmc2_samples <- genMCMC( datFrm=datFrm, yName=yName, xName=xName,
                                numSavedSteps=numSavedSteps, 
                                prior_params = prior_params, rho_0=FALSE)
      summaryInfo2  <-  smryMCMC( mcmc2_samples, 
                                  RopeMuDiff=RopeMuDiff)
      result <- summaryInfo2["mu","Result"]
      pcltRope <- round(as.numeric(summaryInfo2["mu","PcntLtROPE"]),4)
      pcinRope <- round(as.numeric(summaryInfo2["mu","PcntInROPE"]),4)
      pcgtRope <- round(as.numeric(summaryInfo2["mu","PcntGtROPE"]),4)
      
      # 儲存 mcmcCoda 和 summaryInfo
      mcmc_storage[[paste0("test2_", i, "_", paste(prior_params, collapse = "_"))]] <- mcmc2_samples
      summary_storage[[paste0("test2_", i, "_", paste(prior_params, collapse = "_"))]] <- summaryInfo2
      
      results <- rbind(results, data.frame(
        dataset = i,
        prior_params = paste(prior_params, collapse = ","),
        model = "test2",
        p_left = pcltRope,
        p_rope = pcinRope,
        p_right = pcgtRope,
        decision = result,
        stringsAsFactors = FALSE
      ))
      
      # 模型 4
      mcmc3_samples <- genMCMC( datFrm=datFrm, yName=yName, xName=xName,
                                numSavedSteps=numSavedSteps, 
                                prior_params = prior_params, rho_0=TRUE)
      summaryInfo3  <-  smryMCMC( mcmc3_samples, 
                                  RopeMuDiff=RopeMuDiff)
      result <- summaryInfo3["mu","Result"]
      pcltRope <- round(as.numeric(summaryInfo3["mu","PcntLtROPE"]),4)
      pcinRope <- round(as.numeric(summaryInfo3["mu","PcntInROPE"]),4)
      pcgtRope <- round(as.numeric(summaryInfo3["mu","PcntGtROPE"]),4)
      
      # 儲存 mcmcCoda 和 summaryInfo
      mcmc_storage[[paste0("test3_", i, "_", paste(prior_params, collapse = "_"))]] <- mcmc3_samples
      summary_storage[[paste0("test3_", i, "_", paste(prior_params, collapse = "_"))]] <- summaryInfo3
      
      results <- rbind(results, data.frame(
        dataset = i,
        prior_params = paste(prior_params, collapse = ","),
        model = "test3",
        p_left = pcltRope,
        p_rope = pcinRope,
        p_right = pcgtRope,
        decision = result,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # 保存 mcmc_storage 和 summary_storage 到文件
  save(mcmc_storage, file = paste0(mcmcCodaRoot,fileName, "mcmc_storage.RData"))
  save(summary_storage, file = paste0(summaryInfoRoot,fileName, "summary_storage.RData"))
  
  return(results)
}

results <- compare_models(myDataFrame, prior_params_list)

print(results)

write.csv(results,
          paste0(csvRoot,fileName,"結果.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")





##################################################################################

##   end end end end end end end end end end end end end end end end end end    ##
##   end end end end end end end end end end end end end end end end end end    ##
##   end end end end end end end end end end end end end end end end end end    ##
##   end end end end end end end end end end end end end end end end end end    ##
##   end end end end end end end end end end end end end end end end end end    ##
##   end end end end end end end end end end end end end end end end end end    ##

##################################################################################


