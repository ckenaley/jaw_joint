
library (ggplot2)
YP501 <- read.csv("DataforLabs/JAW_YP5_01.csv",stringsAsFactors=T)
YP502 <- read.csv("DataforLabs/JAW_YP5_02.csv",stringsAsFactors=T)
YP503 <- read.csv("DataforLabs/JAW_YP5_03.csv",stringsAsFactors=T)
YP504 <- read.csv("DataforLabs/JAW_YP5_04.csv",stringsAsFactors=T)
YP505 <- read.csv("DataforLabs/JAW_YP5_05.csv",stringsAsFactors=T)



#create scatterplot of x1 vs. y1
plot(YP501$ms, YP501$g, col='red', pch=19, 
     xlab='time in ms', ylab='grams', main='Grams as a fct of time YP5')

YP502$g
#overlay scatterplot of x2 vs. y2
points(YP502$ms, YP502$g, col='blue', pch=19)

points(YP503$ms, YP503$g, col='burlywood', pch=19)

points(YP504$ms, YP504$g, col='chartreuse', pch=19)
points(YP505$ms, YP504$g, col='lightsalmon4', pch=19)


#add legend
legend(1, 1, legend=c('Trial 1', 'Trial 2', 'Trial 3','Trial 4','Trial 5'), pch=c(19, 19,19, 19,19 ), col=c('red', 'blue','burlywood', 'chartreuse','lightsalmon4'))
