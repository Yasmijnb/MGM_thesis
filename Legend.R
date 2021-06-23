# Create legend for network

################################################################################

library(ggplot2)

################################################################################

# Change a dataframe to comply with the network

# Node colour
mtcars$gear<- rep(c('Bacterial gene','Clinical','Cytokine','Human gene'), nrow(mtcars)/4)
colnames(mtcars)[which(colnames(mtcars)=='gear')] <- 'Colour' 
mtcars$Colour <- as.factor(mtcars$Colour)

# Node shape
mtcars$carb<- rep(c('Continous','Discrete'), nrow(mtcars)/2)
colnames(mtcars)[which(colnames(mtcars)=='carb')] <- 'Shape' 
mtcars$Shape <- as.factor(mtcars$Shape)

# Edge colour
mtcars$qsec<- rep(c('Positive','Negative'), nrow(mtcars)/2)
colnames(mtcars)[which(colnames(mtcars)=='qsec')] <- 'Correlation' 
mtcars$Correlation<-as.factor(mtcars$Correlation)

################################################################################

# Make the network for the nodes
p <- ggplot(data = mtcars, 
            aes(x=mpg, y=wt, color=Colour, shape=Shape))+
  geom_point()
p + scale_color_manual(values=c('blue','magenta','orange','cyan')) +
  scale_shape_manual(values=c(19,15))

################################################################################

# Make the network for the edges
q <- ggplot(data = mtcars, 
            aes(x=mpg, y=wt, color=Correlation))+
  geom_line()
q + scale_color_manual(values=c('red','darkgreen'))
