# imports
library("tidyverse")
library("gganimate")
library("ggplot2")
library("gifski")
library("extrafont")
loadfonts(device = "win")

#data preprocessing
data <- read.csv("data\\data4.csv")
y <- data[,-2]
z <- y[,-2]
data1 <- y %>% 
  group_by(Step, Species) %>% 
  summarise(across(everything(), sum))
data2 <- z %>% 
  group_by(Step) %>% 
  summarise(across(everything(), sum))
data1 <- as.data.frame(data1)
data2 <- as.data.frame(data2)
data1 <- data1[, colSums(data1 != 0) > 0]
data2 <- data2[, colSums(data2 != 0) > 0]
data1[,1] = rep(0:119, each=4)
data2[,1] = c(0:119)
long_df <- data1 %>% gather(Key, Value, -Step, -Species)
long_df_2 <- data2 %>% gather(Key, Value, -Step)
long_df_2 = long_df_2 %>%
  group_by(Step)%>%      
  mutate(rank = rank(-Value, ties.method = 'first'), #row_number(),
         label = paste0(" ", Value)) %>%
  group_by(Value)
long_df_2 = long_df_2[,-3]
total <- merge(long_df,long_df_2,by=c("Step","Key"))

plot = ggplot(data=total, aes(x=Step, y=Value, group = Key, colour = Key)) +
  geom_line() +
  theme(legend.position = "none")
print(plot)
