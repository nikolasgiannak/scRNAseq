library(ggplot2) 
library(tidyr)
library(dplyr)

head(ToothGrowth)
str(ToothGrowth)

# assign colors
wisteria <- c("grey65", "burlywood3", "khaki2", "plum1", "lightcyan2", "cornflowerblue", "slateblue3")

ToothGrowth %>% 
  ggplot(aes(x = dose, y = len)) +
  geom_point(aes(color = supp)) +
  scale_color_manual(values = wisteria[c(1, 7)]) +
  labs(x = "dose",              # the labs option allow you to change axis and legend labels. 
       y = "teeth length",
       color = "supplement")

data.frame(
  x = c("group1", "group2"),
  y = c(10, 20)
) %>% 
  ggplot(aes(x = x, y = y)) +
  geom_bar(stat = "identity", aes(fill = x)) +
  ggtitle("fill by x")

data.frame(
  x = c("group1", "group2"),
  y = c(10, 20)
) %>% 
  ggplot(aes(x = x, y = y)) +
  geom_bar(stat = "identity", aes(color = x)) +
  ggtitle("color by x")


ToothGrowth %>% 
  ggplot(aes(x = dose, y = len)) +
  stat_summary(geom = "line",  # draw a line
               fun = mean, # summarise by averaging, thus mean
               aes(group = supp, color = supp), # one line per supp, color by supp
               size = 1.2) + 
  geom_point(aes(fill = supp), 
             position = position_jitter(0.1, seed = 666),    
             alpha = 0.8,             
             size = 3,
             shape = 21, 
             color = "black") +    
  scale_fill_manual(values = wisteria[c(1, 7)]) +  
  labs(x = "dose",              
       y = "teeth length",
       fill = "supplement")


ToothGrowth %>% 
  ggplot(aes(x = dose, y = len)) +
  stat_summary(geom = "line",  
               fun = mean,  
               aes(group = supp, color = supp),  
               size = 1.2) + 
  geom_point(aes(fill = supp), 
             position = position_jitter(0.1, seed = 666),    
             alpha = 0.8,             
             size = 3,
             shape = 21, 
             color = "black") +    
  scale_fill_manual(values = wisteria[c(1, 7)]) +  
  scale_color_manual(values = wisteria[c(1, 7)]) + #add scale_color_manual here  
  labs(x = "dose",              
       y = "teeth length",
       fill = "supplement",
       color = "supplement") #change the color legend title here 




ToothGrowth %>% 
  ggplot(aes(x = dose, y = len)) +
  stat_summary(geom = "line",  
               fun = mean,  
               aes(group = supp, color = supp),  
               size = 1.2) + 
  stat_summary(geom = "errorbar",  # draw errorbars (another option is geom = "linerange")
               fun.data = mean_se, # mean_se specifies the errorbar to be based on standard error (se)
               aes(group = supp), # one errorbar per dose per supp
               width = 0.1) +   # make the bar narrower 
  geom_point(aes(fill = supp), 
             position = position_jitter(0.1, seed = 666),    
             alpha = 0.8,             
             size = 3,
             shape = 21, 
             color = "black") +    
  scale_fill_manual(values = wisteria[c(1, 7)]) +  
  scale_color_manual(values = wisteria[c(1, 7)]) + 
  labs(x = "dose",              
       y = "teeth length",
       fill = "supplement",
       color = "supplement")




ToothGrowth %>% 
  ggplot(aes(x = dose, y = len)) +
  stat_summary(geom = "line",  
               fun = mean,  
               aes(group = supp, color = supp),  
               size = 1.2) + 
  stat_summary(geom = "ribbon",  #draw ribbons around line
               fun.data = mean_se, 
               aes(group = supp,
                   fill = supp), 
               alpha = 0.5) +   #make ribbons more transparent 
  geom_point(aes(fill = supp), 
             position = position_jitter(0.1, seed = 666),    
             alpha = 0.8,             
             size = 3,
             shape = 21, 
             color = "black") +    
  scale_fill_manual(values = wisteria[c(1, 7)]) +  
  scale_color_manual(values = wisteria[c(1, 7)]) + 
  labs(x = "dose",              
       y = "teeth length",
       fill = "supplement",
       color = "supplement")




ToothGrowth %>% 
  ggplot(aes(x = dose, y = len)) +
  stat_summary(geom = "line",  
               fun = mean,  
               aes(group = supp, color = supp),  
               size = 1.2) + 
  stat_summary(geom = "ribbon",   
               fun.data = mean_se, 
               aes(group = supp,
                   fill = supp), 
               alpha = 0.5) +   
  geom_point(aes(fill = supp), 
             position = position_jitter(0.1, seed = 666),    
             alpha = 0.8,             
             size = 3,
             shape = 21, 
             color = "black") +    
  scale_fill_manual(values = wisteria[c(1, 7)]) +  
  scale_color_manual(values = wisteria[c(1, 7)]) + 
  scale_x_continuous(breaks = c(0.5, 1, 2)) +  
  labs(x = "dose",              
       y = "teeth length",
       fill = "supplement",
       color = "supplement") +
  guides(color = "none") +
  theme_minimal() +
  theme(
    axis.line = element_line(size = 1.2), # add axis line back
    text = element_text(size = 12, color = "black", face = "bold"), # make text clearer and darker
    axis.text = element_text(size = 12, color = "black", face = "bold"),
    legend.position = c(0.8, 0.2) # move legend to inside graph area
  ) 




ToothGrowth %>% 
  ggplot(aes(x = dose, y = len)) +
  facet_grid(. ~ supp) +  
  stat_summary(geom = "line",  
               fun = mean,  
               aes(group = supp, color = supp),  
               size = 1.2) + 
  stat_summary(geom = "ribbon",   
               fun.data = mean_se, 
               aes(group = supp,
                   fill = supp), 
               alpha = 0.5) +   
  geom_point(aes(fill = supp), 
             position = position_jitter(0.1, seed = 666),    
             alpha = 0.8,             
             size = 3,
             shape = 21, 
             color = "black") +    
  scale_fill_manual(values = wisteria[c(6, 7)]) +  
  scale_color_manual(values = wisteria[c(6, 7)]) + 
  scale_x_continuous(breaks = c(0.5, 1, 2)) +  
  labs(x = "dose",              
       y = "teeth length",
       fill = "supplement",
       color = "supplement") +
  guides(color = "none") +
  guides(fill = guide_legend(nrow = 1, ncol = 2)) +  
  theme_minimal() +
  theme(
    axis.line = element_line(size = 1.2), # add axis line back
    text = element_text(size = 12, color = "black", face = "bold"), # make text clearer and darker
    axis.text = element_text(size = 12, color = "black", face = "bold"),
    legend.position = "none", # move legend to inside graph area
    panel.spacing = unit(1.2, "lines") #make 1.2 lines of spacing between the subplots (panels)
  ) 



ToothGrowth %>% 
  ggplot(aes(x = supp, y = len)) +
  facet_grid(. ~ dose, switch = "x") + #switch = "x" puts x panel labels to the bottom instead of top    
  geom_point(aes(fill = supp), 
             position = position_jitter(0.1, seed = 666),    
             alpha = 0.8,             
             size = 3,
             shape = 21, 
             color = "black") + 
  stat_summary(geom = "errorbar", 
               width = 0.2, size = 1,
               aes(group = supp),  
               fun.data = mean_se) +     
  scale_fill_manual(values = wisteria[c(1, 7)]) +  
  scale_color_manual(values = wisteria[c(1, 7)]) + 
  scale_x_discrete(labels = NULL) +  
  labs(x = "dose",  #now x is dose again               
       y = "teeth length",
       fill = "supplement",
       color = "supplement") +
  guides(color = "none") +
  guides(fill = guide_legend(nrow = 1, ncol = 2)) +  
  theme_minimal() +
  theme(
    axis.line = element_line(size = 1.2), # add axis line back
    text = element_text(size = 12, color = "black", face = "plain"), # make text clearer and darker
    axis.text = element_text(size = 12, color = "black", face = "bold"),
    legend.position = "bottom", # move legend to inside graph area
    panel.spacing = unit(1.2, "lines"), #make 1.2 lines of spacing between the subplots (panels)
    strip.placement = "outside"
  ) 