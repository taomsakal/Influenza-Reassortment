library("tidyverse")
library("gganimate")
library("ggplot2")
data <- read.csv("data\\data.csv")
numeric_data <- data[, -c(2:3)]
y <- filter(numeric_data, Step == 0)
data1 <- as.matrix(colSums(y))
data1 <- t(data1)
data1 <- as.data.frame(data1)
for (i in 1:9) {
    data1 <- rbind(data1, colSums(filter(numeric_data, Step == i)))
}
data1 <- as.data.frame(data1[, -1])
data1 <- as.data.frame(t(data1))
data1 <- data1[rowSums(data1[])>0,]
data1 <- as.data.frame(t(data1))
rownames(data1)<-c(0:9)
data1 <- tibble::rownames_to_column(data1, "step")
long_df <- data1 %>% gather(Key, Value, -step)
long_df = long_df %>% 
  group_by(step)%>%      
  mutate(rank = rank(-Value, ties.method = 'first')) %>%
  group_by(Value)
print.data.frame(long_df)
anim = ggplot(data=long_df, mapping=aes(x=rank, group=Key, fill=Key), fill=colnames(data1))+
        geom_tile(aes(y = (Value/2) + 25,
                      height = Value,
                      width = 0.7)) +
        geom_text(aes(y = 0, label = Key, vjust = 0.2, hjust = 0.5, size = 1, angle = 90)) +
        theme(legend.position = "none",
              axis.line=element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              panel.background=element_blank(),
              panel.border=element_blank(),
              plot.background=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              plot.margin = margin(0, 0, 00, 0, "mm")) +
              transition_states(step, transition_length = 4, state_length = 1)+
              ease_aes('sine-in-out')

library(gifski)
library(png)
animate(anim, nframes = 350,fps = 25,  width = 1200, height = 1000, 
        renderer = gifski_renderer("strains.gif"))



