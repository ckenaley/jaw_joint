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
  mutate(f=g*0.0098066500286389) %>% 
arrange(ang,fish,species,trial)
  



dat %>% 
  ggplot(aes(ang,f,col=trial))+
  geom_point()+facet_wrap(.~fish)

#dat_sum <- 
  
se <- function(x) {sd(x) / sqrt(sum(!is.na(x)))}
dat_sum <-   dat%>% 
  group_by(ang,species,fish) %>% 
  filter(g!=0) %>% 
  summarize(se=se(f),f=mean(f),n=n()) 

w <- dat_sum %>% 
  ggplot(aes(ang,f,col=species))+
  geom_errorbar(aes(x=ang,ymin=f-se,ymax=f+se))+
  geom_point()+
  xlab("Angular displacemnt")+
  ylab("Force (N)")+
  scale_x_continuous(labels = ~ paste0(.x, "\u00B0"))+
    geom_text(x=20,y=0.045,label="1.36e-02 N cm",col="black")+
  theme_classic(12)
  
ggsave("jaw_curves.png",p,device = "png")
  
  expand_grid(dat_sum,jaw=c(15,30)) %>% 
    mutate(Nm=jaw/10*f) %>%
   group_by(ang,species,fish) %>% 
    summarise(T_max=max(Nm)) %>% 
    filter(T_max>0) %>% 
    pull(T_max) %>% 
    range()
