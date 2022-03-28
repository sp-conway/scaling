# pairwise.R
# limited function for generating unique combos of angle and radius in turtle/protractor stimuli
# Sean Conway
# last modified Feb. 2022

# libraries needed
library(dplyr)
library(tidyr)
library(ggplot2)

pairwise <- function(angle, radius, plot=F){
  # all angle-radius combos
  stim <- crossing(angle, radius) %>%
    mutate(stimulus=1:n()) 
  
  # all combos of unique stimuli
  combos <- t(combn(1:nrow(stim),2))
  
  # "x" stimuli 
  x <- tibble(
    stimulus=combos[,1],
    which="x",
    trial=1:nrow(combos)
  ) %>%
    left_join(stim,by="stimulus")
  
  # "y" stimuli
  y <- tibble(
    stimulus=combos[,2],
    which="y",
    trial=1:nrow(combos)
  ) %>%
    left_join(stim,by="stimulus")
  
  # bind and pivot x & y stim
  all <- bind_rows(x,y) %>%
    select(-stimulus) %>%
    pivot_wider(id_cols=trial,
                names_from=which,
                values_from=c(angle,radius))
  
  # plot stimulus space 
  if(plot){ 
    # axis limits
    alim <- c(floor(min(stim$angle)),
              ceiling(max(stim$angle)))
    rlim <- c(floor(min(stim$radius)),
              ceiling(max(stim$radius)))
    # plotting code
    pl <- stim %>%
      ggplot(aes(angle,radius))+
      geom_point(alpha=.6,col="darkblue",size=3.5)+
      ggthemes::theme_few()+
      scale_x_continuous(limits=alim,breaks=round(seq(alim[1],alim[2],length.out=4)))+
      scale_y_continuous(limits=rlim,breaks=round(seq(rlim[1],rlim[2],length.out=4)))+
      labs(title="stimulus space")+
      theme(plot.title=element_text(hjust=0.5))
    return(list(all, pl))
  }else{
    return(all)
  }
}

# values used in scaling study
angle <- c(10, 17, 25, 33, 40)
radius <- c(30,  57, 113,  85, 140)

# all combos
pairwise(angle,radius)

# same function but now plot them!
pairwise(angle,radius,plot=T)

stim <- crossing(angle,radius)

write_csv(stim, here("data","stim.csv"))
