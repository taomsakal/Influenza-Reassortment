# imports
library("tidyverse")
library("gganimate")
library("ggplot2")
library("gifski")
library("extrafont")
library("grid")
library("ggalt")
font_import()
library(extrafont)
loadfonts(device="win")  
windowsFonts(x=windowsFont("Times New Roman"))
#data preprocessing
data <- read.csv("data\\data9.csv")
print.data.frame(data)

data_bar <- data[data$Iteration == 0,]
data_bar <- data[data$Step > 20,]
data_bar <- data_bar[,-2]
data_bar <- data_bar[,-2]
data_bar <- data_bar %>% 
  group_by(Step) %>% 
  summarise(across(everything(), sum))
data_bar
data_bar <- data_bar[, colMeans(data_bar) > 35]
long_df <- data_bar %>% gather(Key, Value, -Step)
long_df
total

plot = ggplot(data=long_df, aes(x=Step, y=Value, group = Key, color=Key)) +
  geom_xspline(spline_shape=2, size=0.75) +
  scale_colour_grey(start = 0.3, end = .75) + 
  theme_classic() + 
  ylab(label = "Individuals Infected") + 
  guides(color=guide_legend(title="Strains")) +
  theme(plot.margin = unit(c(1,3,1,1), "lines"),  text=element_text(family="Times New Roman", size=12), legend.position = "none")
#geom_text(data = subset(long_df, Step == 400), aes(label = Key, colour = Key, x = Inf, y = Value), hjust = -.1, size=2.5) +

gt <- ggplotGrob(plot)
gt$layout$clip[gt$layout$name == "panel"] <- "off"

png("lineplot.png", 7, 4, units = 'in', res = 300,)
grid.draw(gt)
dev.off()

