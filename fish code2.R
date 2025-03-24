#install.packages(c("Rcpp","readr"))
library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)
lmb3_05 <- read_csv("DataForLabs/JAW_LB3_05.csv")
lmb3_05
lmb3_05 %>% 
  ggplot(aes(ms,g))+geom_point()
f <- list.files("DataForLabs",pattern = "JAW",full.names = T)
f
dat <- read_csv(f,id = "file")
CP5_01 <- read_csv("DataForLabs/JAW_CP1_01 (1).csv")
CP5_01
LB1_02 <- read_csv("DataForLabs/JAW_LB1_02.csv")
#LB1_02[nrow(LB1_02)+1,]=c(0,0,200,18,0.279)
#LB1_02<-LB1_02%>%add_row(0="0",0.000="0",183.00="00",0.368="0.368")
