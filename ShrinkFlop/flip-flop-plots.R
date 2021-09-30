library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(viridis)

plot_flipflop <- function(df, plot_title, show_legend = T, point_size = 1, text_size = 12, show_plot = T) {
  
  mytheme <- theme(legend.position="top",
                  legend.text = element_text(family = "Helvetica", size = rel(1)), 
                   axis.title = element_text(family = "Helvetica", size = rel(1)), 
                   axis.text = element_text(family = "Helvetica", size = rel(1)), 
                   axis.line = element_line(size = 1,colour = "black"), 
                   axis.ticks = element_line(colour="black",size = rel(1)),
                   panel.grid.major = element_line(colour="grey90",size = rel(0.5)), 
                   panel.grid.minor = element_blank(), 
                   panel.background = element_rect(fill = "grey98"), 
                   legend.key = element_rect(fill = "grey98"), 
                   legend.title = element_text(family = "Helvetica", size = rel(1)), 
                   plot.title = element_text(face = "bold", size = rel(1.25),family = "Helvetica"))
  
  sc <- scale_color_viridis(discrete = F, option = "plasma", begin = 0, end = 0.95)
  
 #df2 = data.frame(data.frame(x=df$Regularizer,y=df$Zhou))
  
  p1 <- ggplot(data = data.frame(x=df$Regularizer, y = df$"Regularized Sinkhorn"),color="red") +
    aes(x,y)+
    geom_line(data = data.frame(data.frame(x=df$Regularizer,y=df$Zhou)),color="green")+
    geom_point(data = data.frame(data.frame(x=df$Regularizer,y=df$Zhou)))+

    geom_point(size = point_size) +
    geom_line(color="red")+
    labs(x = "Regularization", y = "Relative error", color = "Class") +
    #theme(legend.position="top")+
    mytheme+
    scale_x_continuous(trans='log2') + scale_y_continuous(trans='log2') 
    
    #aes(x = X1, y = X2, color = y) +
    
  if (show_plot) {
    print(p1)
  }
}
