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
data1 <- as.data.frame(data1)
data2 <- as.data.frame(data2)
data1 <- data1[, colSums(data1 != 0) > 0]
data2 <- data2[, colSums(data2 != 0) > 0]
data1[,1] = rep(0:9, each=4)
data2[,1] = c(0:9)
long_df <- data1 %>% gather(Key, Value, -Step, -Species)
long_df_2 <- data2 %>% gather(Key, Value, -Step)
print.data.frame(long_df_2)
long_df_2 = long_df_2 %>%
  group_by(Step)%>%      
  mutate(rank = rank(-Value, ties.method = 'first'),
         label = paste0(" ", Value)) %>%
  group_by(Value)
long_df_2 = long_df_2[,-3]
total <- merge(long_df,long_df_2,by=c("Step","Key"))
print.data.frame(total)
anim = ggplot(data=total, mapping=aes(x=rank, y=Value, fill=Species, frame = Step))+
        geom_bar(stat = "identity", position = "stack") +
        geom_text(aes(y = -1, label = Key, vjust = 0.25, hjust = 1, angle = 90, family = "Arial Narrow"), size = 6) +
        geom_text(aes(y = as.integer(sub('.', '', label)), label = label, hjust = 0, vjust = 0.25, angle = 90, family = "Arial Narrow"), size = 6) +
        labs(y = "Population Infected", x="Steps") +
        ggtitle("Virus Strains Over Time") +
        guides(size = FALSE) +
        theme(legend.position = c(0.9, 0.8),
              legend.key.size = unit(1, 'cm'), #change legend key size
              legend.key.height = unit(1, 'cm'), #change legend key height
              legend.key.width = unit(1, 'cm'), #change legend key width
              legend.title = element_text(size=16, family = "Franklin Gothic Demi"), #change legend title font size
              legend.text = element_text(size=14, family = "Franklin Gothic Demi Cond"),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks=element_blank(),
              axis.title.y=element_text(vjust=-7, hjust=0.5, size=25, family = "Franklin Gothic Demi Cond"),
              axis.title.x=element_text(vjust=5, hjust=0.5, size=25, family = "Franklin Gothic Demi Cond"),
              panel.background=element_blank(),
              panel.border=element_blank(),
              plot.title = element_text(size=50, vjust=-2, hjust=0.5, family = "Franklin Gothic Heavy"),
              plot.background=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              plot.margin = margin(0, 0, 0, 0, "mm")) +
              scale_y_continuous(expand = c(.1, .1)) +
        transition_states(Step, transition_length = 1, state_length = 1) +
        ease_aes('sine-in-out')

animate(anim, nframes = 525,fps = 35,  width = 1200, height = 900, 
        renderer = gifski_renderer("animated_bar.gif"))

