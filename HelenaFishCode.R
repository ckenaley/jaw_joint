library(tidyverse) 

f <- list.files(pattern="JAW") 

dat <- read_csv(f,id="file") %>% 
  mutate(fish=gsub("JAW_(.*)_.*.csv","\\1",file),
         species=gsub("(\\w*)\\d","\\1",fish),
         trial=gsub("JAW_.*_(.*).csv","\\1",file)
  ) %>% 
  group_by(file) %>% 
  filter(ms<ms[which.max(g)]) %>% 
  mutate(g=g-first(g),
         n=1:n(),
         ang=(n/max(n)*30),
         ang=plyr::round_any(ang,0.5)
  ) %>% 
  group_by(fish,species,trial,ang) %>% 
  summarize(g=mean(g)) %>% 
  arrange(ang,fish,species,trial)




dat %>% 
  ggplot(aes(ang,g,col=trial))+
  geom_point()+facet_wrap(.~fish)

#dat_sum <- 

se <- function(x) {sd(x) / sqrt(sum(!is.na(x)))}
dat%>% 
  group_by(ang,fish,species) %>% 
  filter(g!=0) %>% 
  summarize(g=mean(g),se=se(g),n=n()) %>% view()


