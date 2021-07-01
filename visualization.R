library("tidyverse")
library("gganimate")
library("ggplot2")
library("gifski")
library("extrafont")
loadfonts(device = "win")
data <- read.csv("data\\data.csv")
y <- data[,-2]
z <- y[,-2]
data1 <- y %>% 
  group_by(Step, Species) %>% 
  summarise(across(everything(), sum))
data2 <- z %>% 
  group_by(Step) %>% 
  summarise(across(everything(), sum))
data1 = as.data.frame(data1)
data2 = as.data.frame(data2)
data1 <- data1[,colSums(data1)>0]
data2 <- data2[,colSums(data2)>0]
data1[,1] = rep(0:9, each=4)
data2[,1] = c(0:9)
long_df <- data1 %>% gather(Key, Value, -Step,Species)
long_df_2 <- data2 %>% gather(Key, Value, -Step,Species)
print.data.frame(data2)
long_df = long_df %>%
  group_by(Step)%>%      
  mutate(rank = rank(-Value, ties.method = 'first'),
         label = paste0(" ", Value)) %>%
  group_by(Value)
print.data.frame(long_df)
anim = ggplot(data=long_df, mapping=aes(x=rank, y=Value, fill=Key, frame = Step))+
        geom_bar(stat = "identity") +
        geom_text(aes(y = 0, label = Key, vjust = 0, hjust = 1.1, size = 1, angle = 90)) +
        geom_text(aes(y = Value, label = label, hjust = 0, size = 1, angle = 90)) +
        labs(y = "Population Infected", x="Steps") +
        ggtitle("Virus Strains Over Time") +
        theme(legend.position = "none",
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks=element_blank(),
              axis.title.y=element_text(vjust=-7, hjust=0.5, size=25, family = "Franklin Gothic Demi"),
              axis.title.x=element_text(vjust=1, hjust=0.5, size=25, family = "Franklin Gothic Demi"),
              panel.background=element_blank(),
              panel.border=element_blank(),
              plot.title = element_text(size=50, vjust=-2, hjust=0.5, family = "Rockwell"),
              plot.background=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              plot.margin = margin(0, 0, 0, 0, "mm")) +
        transition_states(Step, transition_length = 1, state_length = 1) +
        ease_aes('sine-in-out')

animate(anim, nframes = 10,fps = 1,  width = 1200, height = 850, 
        renderer = gifski_renderer("strains.gif"))

