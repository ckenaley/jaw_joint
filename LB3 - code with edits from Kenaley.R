library (ggplot2)
library(tidyr)
YP501 <- read.csv("DataforLabs/JAW_LB3_01.csv",stringsAsFactors=T)
YP502 <- read.csv("DataforLabs/JAW_LB3_02.csv",stringsAsFactors=T)
YP503 <- read.csv("DataforLabs/JAW_LB3_03.csv",stringsAsFactors=T)
YP504 <- read.csv("DataforLabs/JAW_LB3_04.csv",stringsAsFactors=T)
YP505 <- read.csv("DataforLabs/JAW_LB3_05.csv",stringsAsFactors=T)

#create scatterplot of x1 vs. y1
plot(YP501$ms, YP501$g, col='red', pch=19, 
     xlab='time in ms', ylab='grams', main='Grams as a fct of time LB3')

YP502$g
#overlay scatterplot of x2 vs. y2
points(YP502$ms, YP502$g, col='blue', pch=19)

points(YP503$ms, YP503$g, col='green', pch=19)

points(YP504$ms, YP504$g, col='yellow', pch=19)


points(YP505$ms, YP505$g, col='purple', pch=19)

#add legend
legend(1,3, legend=c('Trial 1','Trial 2', 'Trial 3', ' Trial 4','Trial 5'), pch=c(19, 19,19,19,19), col=c('red', 'blue','green','yellow','purple')


lmb3_05 <- read_csv("DataforLabs/Jaw_LB3_05.csv")

lmb3_05 %>% 
  ggplot(aes(ms,g))+geom_point()



f <- list.files("DataforLabs",pattern = "JAW",full.names = T)

dat <- read_csv(f,id = "file")
