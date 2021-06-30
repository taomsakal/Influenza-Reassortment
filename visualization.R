library("tidyverse")
library("gganimate")
library("ggplot2")
library("gifski")
data <- read.csv("data\\data.csv")
y <- filter(numeric_data, Step == 0)
data1 <- as.matrix(colSums(y))
data1 <- t(data1)
data1 <- as.data.frame(data1)
for (i in 1:9) {
    data1 <- rbind(data1, colSums(filter(numeric_data, Step == i)))
}
data1 <- data1[,colSums(data1)>0]
data1[,1] = c(0:9)
long_df <- data1 %>% gather(Key, Value, -Step)
long_df = long_df %>%
  group_by(Step)%>%      
  mutate(rank = rank(-Value, ties.method = 'first'),
         label = paste0("       ", Value)) %>%
  group_by(Value)
print.data.frame(long_df)
anim = ggplot(data=long_df, mapping=aes(x=rank, group=Key, fill=Key))+
        geom_tile(aes(y = (Value/2) + 12,
                      height = Value,
                      width = 0.7)) +
        geom_text(aes(y = 10, label = Key, vjust = 0.2, hjust = 1.1, size = 1, angle = 90)) +
        geom_text(aes(y = Value, label = label, hjust = 0, size = 1, angle = 90)) +
        labs(y = "Population Infected", x="Steps") +
        ggtitle("Virus Strains Over Time") +
        theme(legend.position = "none",
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks=element_blank(),
              axis.title.y=element_text(vjust=-7, hjust=0.5, size=25, family = "Franklin Gothic Demi"),
              axis.title.x=element_text(vjust=2.5, hjust=0.5, size=25, family = "Franklin Gothic Demi"),
              panel.background=element_blank(),
              panel.border=element_blank(),
              plot.title = element_text(size=50, vjust=-2, hjust=0.5, family = "Rockwell"),
              plot.background=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              plot.margin = margin(0, 0, 0, 0, "mm")) +
        transition_states(Step, transition_length = 1, state_length = 1) +
        ease_aes('sine-in-out')

animate(anim, nframes = 600,fps = 40,  width = 1200, height = 850, 
        renderer = gifski_renderer("strains.gif"))

