# imports
library("tidyverse")
library("ggplot2")
library("plotly")
library("corrplot")
library("Hmisc")
library("ggcorrplot")

data <- read.csv("data\\full_data.csv")[,-2]
print.data.frame(data)
data <- data %>% 
  group_by(Iteration) %>% 
  summarise(across(everything(), mean))
data[] <- lapply(data, function(x) ifelse(x > 1 , 1 , 0 ))
data <- as.data.frame(colSums(data))
colnames(data) <- c("Value")
ggplot(data, aes(x=row.names(data), y=Value)) +
  geom_bar(stat="identity")


data <- data %>% 
  group_by(Step, Iteration) %>% 
  summarise(across(everything(), sum))
data0 <- data[data$Iteration == 0,]
data1 <- data[data$Iteration == 1,]
data0 <- data0[,-1]
data0 <- data0[,-1]
data1 <- data1[,-1]
data1 <- data1[,-1]
data0 = data0[, colMeans(data0) > 1]
data1 = data1[, colMeans(data1) > 1]

res <- round(cor(data0),2)
print(res)
res[,order(colnames(res))]

corrplot(res, method = "color", type = "upper",
     tl.cex=0.43, order = 'alphabet', tl.col = "black")

res2 <- rcorr(as.matrix(data0))

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

res3 <- flattenCorrMatrix(res2$r, res2$P)
res3[order(-res3$cor),]
print(tibble(data0[,"H1N3"]), n=50)



corr_simple <- function(data=df,sig=0.4){
  #convert data to numeric in order to run correlations
  #convert to factor first to keep the integrity of the data - each value will become a number rather than turn into NA
  df_cor <- data %>% mutate_if(is.character, as.factor)
  df_cor <- df_cor %>% mutate_if(is.factor, as.numeric)
  #run a correlation and drop the insignificant ones
  corr <- cor(df_cor)
  #prepare to drop duplicates and correlations of 1     
  corr[lower.tri(corr,diag=TRUE)] <- NA 
  #drop perfect correlations
  corr[corr == 1] <- NA
  #turn into a 3-column table
  corr <- as.data.frame(as.table(corr))
  #remove the NA values from above 
  corr <- na.omit(corr) 
  #select significant values  
  corr <- subset(corr, abs(Freq) > sig) 
  #sort by highest correlation
  corr <- corr[order(-abs(corr$Freq)),] 
  #print table
  print(corr)
  #turn corr back into matrix in order to plot with corrplot
  mtx_corr <- reshape2::acast(corr, Var1~Var2, value.var="Freq")
  print(mtx_corr)
  #plot correlations visually
  corrplot(mtx_corr, is.corr=FALSE, tl.col="black", na.label=" ", tl.cex = 0.6, method = "square")
}
corr_simple(data0)
corr_simple(data1)
