# imports
library("tidyverse")
library("ggplot2")
library("plotly")
library("corrplot")
library("Hmisc")
library("igraph")
library("ggraph")
library("Rcpp")
library("cowplot")

data <- read.csv("data\\data9.csv")
print.data.frame(data)

#making barplot
data_bar <- data[data$Iteration == 0,]
data_bar
data_bar <- data_bar[,-2]

data_bar <- data_bar %>% 
  group_by(Step) %>% 
  summarise(across(everything(), sum))
data_bar <- data_bar[, colSums(data_bar != 0) > 0]
data_bar
long_df <- data_bar %>% gather(Key, Value, -Step, -Species)
long_df
total <- long_df[long_df$Step == 375,]
ggplot(data=total, mapping=aes(x=Key, y=Value, fill=Species))+
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(y = -1, label = Key, vjust = 0.05, hjust = 0.87, angle = 90), size = 2.5) +
  labs(y = "Population Infected", x="Steps") +
  ggtitle("Step 350: Virus Strains") +
  theme(legend.position = c(0.9, 0.8),
        legend.key.size = unit(0.5, "cm"), #change legend key size
        legend.key.height = unit(0.5, "cm"), #change legend key height
        legend.key.width = unit(0.5, "cm"), #change legend key width
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=6),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.y=element_text(vjust=-7, hjust=0.5, size=10),
        axis.title.x=element_text(vjust=2, hjust=0.5, size=10),
        panel.background=element_blank(),
        panel.border=element_blank(),
        plot.title = element_text(size=20, vjust=-2, hjust=0.5),
        plot.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.margin = margin(0, 0, 0, 0, "mm"))+
        scale_y_continuous(expand = c(.1, .1))

#makign correlation matrices
data <- read.csv("data\\data9.csv")
#data <- data[data$Species == "Human",]
data <- data[,-3]
data <- data %>% 
  group_by(Step, Iteration) %>% 
  summarise(across(everything(), sum))

data <- data[data$Step > 220,]
data <- data[,colSums(data!= 0) > 0]
data <- data[,-1]
data <- data[,-1]

data0 <- data[data$Iteration == 0,]
data <- data[data$Step > 250,]
data0 <- data0[,-1]
data0 <- data0[,-1]
data0 = data0[, colMeans(data0) > 1]

res <- cor(as.matrix(data))
corrplot(res,order = 'alphabet', method = 'color', type='upper',tl.cex=0.6)

flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  = (cormat[ut])
  )
}

res3 <- flattenCorrMatrix(res)
res3[order(-res3$cor),]
print(tibble(data0[,"H1N3"]), n=50)

#Makign Adjacency Networks

res[ res < 0.8 ] <- 0
diag(res) <- 0
graph <- graph.adjacency(mode = "undirected", res, weighted=TRUE)
V(graph)$name <- colnames(res)
Isolated = which(degree(graph)==0)
graph = delete.vertices(graph, Isolated)
l <- layout.fruchterman.reingold(graph, niter=5000, area=vcount(graph)^4*10)
l2 <- 100*layout_with_fr(graph)


gg = ggraph(graph, layout=l) + 
  geom_edge_link(color="gray50", width=0.6) +
  geom_node_point(color="black") +
  geom_node_text(aes(label = name), size=2, repel = TRUE)

pdf("plot.pdf", 5, 2.5)
#corrplot(res,order = 'alphabet', method = 'color', type='upper',tl.cex=0.6)
plot(gg)
dev.off()

#boxplots
data <- read.csv("data\\full_data.csv")
print.data.frame(data)
data_bar <- data[data$Iteration == 0,]
data_bar <- data_bar[,-2] #or-3
data_bar <- data_bar %>% 
  group_by(Step, Species) %>% 
  summarise(across(everything(), sum))
data_bar <- data_bar[, colSums(data_bar != 0) > 0]
data_bar
long_df <- data_bar %>% gather(Key, Value, -Step, -Species)
plot1 <- ggplot(data = long_df[which(long_df$Value>0),], aes(x=Species, y=Value)) + 
  geom_boxplot(aes(fill=Species)) +
  scale_fill_grey() + theme_classic() + theme(legend.position = "none", text=element_text(family="Times New Roman", size=11)) +
  ylab("Number of Individuals a Strain Infects")

long_df_2 <- long_df %>% 
  group_by(Step, Species) %>% 
  summarise(counts = sum(Value > 0, na.rm = TRUE))
long_df_2
plot2 <- ggplot(data = long_df_2, aes(x=Species, y=counts)) + 
  geom_boxplot(aes(fill=Species)) +
  scale_fill_grey() + theme_classic() + theme(legend.position = "none", text=element_text(family="Times New Roman", size=11)) +
  ylab("Number of Strains in circulation")

png("boxplot.png", 7, 4, units = 'in', res = 300,)
plot_grid(plot1, plot2, labels = "AUTO")
dev.off()
