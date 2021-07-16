# imports
library("tidyverse")
library("ggplot2")
library("shiny")
library("plotly")
library(Cairo)
options(shiny.usecairo=TRUE)

#data preprocessing
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
data1[,1] = rep(0:99, each=4)
data2[,1] = c(0:99)
long_df <- data1 %>% gather(Key, Value, -Step, -Species)
long_df_2 <- data2 %>% gather(Key, Value, -Step)
long_df_2 = long_df_2 %>%
  group_by(Step)%>%      
  mutate(rank = rank(-Value, ties.method = 'first'), #row_number(),
         label = Value)
long_df_2 = long_df_2[,-3]
total <- merge(long_df,long_df_2,by=c("Step","Key"))

#Shiny App

#UI
ui <- fluidPage(
  mainPanel(plotlyOutput("plot2")))

#Server
server <- function(input,output, session){
  output$plot2<-renderPlotly({
    data <- total
    b <- length(unique(data[["Key"]]))-1
    max = max(data$label)
    print(max)
    plot_ly(data[data$Species == "Bird",], x = ~Key, y = ~Value, type = 'bar', name='Bird', frame= ~Step) %>% 
    add_trace(data = data[data$Species == "Poultry",], y = ~Value, name = 'Poultry') %>% 
    add_trace(data = data[data$Species == "Pig",],y = ~Value, name = 'Pig') %>% 
    add_trace(data = data[data$Species == "Human",],y = ~Value, name = 'Human') %>%
    layout(yaxis = list(title = 'Count',range = c(0,max+1)), barmode = 'stack', xaxis = list(title = 'Strains', tickvals = 0:b + 1/b, tickfont = list(size = 6), autorange = TRUE)) %>%
    animation_opts(
        1000, redraw = FALSE
    )
    
    })
}

#Launch
shinyApp(ui=ui, server=server)
