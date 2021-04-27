library(reshape2)
library(scales)

#ideally this would be able to output a list of plots, one for each error metric. 

SimplePlot <- function(results){
  mytheme <- theme(legend.text = element_text(family = "Helvetica", size = rel(1)), 
                 #legend.position = "top",
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

  df = melt(results, id.vars = "Regularizer",variable.name = "Algorithm", value.name = "Error")
  #df = melt(results, id.vars = "Regularizer",variable.name = "Algorithm")
  lp <- ggplot(data=df, aes(x=Regularizer, y=Error, group=Algorithm, color=Algorithm)) + geom_line() + geom_point() + scale_x_log10(breaks = 10^(-10:10),labels = trans_format("log10", math_format(10^.x)))+ scale_y_log10(breaks = 10^(-10:10),labels = trans_format("log10", math_format(10^.x)))+ mytheme
  
  return(lp)}

