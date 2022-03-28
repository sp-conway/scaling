# turtle_MDS.R
# Sean Conway
# Last modified March 2022

# setup =====================================================================
rm(list=ls())

# libraries
library(here)
library(dplyr)
library(ggplot2)
library(tidyr)
library(fs)
library(purrr)
library(patchwork)

# controls
save_plot <- T

# The MASS package is used for running MDS, but loading it masks some dplyr functions I need
if("MASS" %in% (.packages())){
  detach("package:MASS",unload=T)
}

# data cleaning =====================================================================

# Read in and clean raw data files
files <- dir_ls(here("data","turtles_Spring2022","raw"))
dat <- map_dfr(files, data.table::fread) %>%
  as_tibble() %>%
  select(-c(internal_node_id,
            success,
            timeout,
            failed_images,
            failed_audio,
            failed_video,
            view_history,
            stimulus,
            question_order,
            trial_type)) 

dat_clean <- dat %>%
  filter(trialType=="rating") %>%
  mutate(across(c(rt,response),as.numeric)) %>%
  mutate(check=case_when(
    angle_a==angle_b & radius_a==radius_b ~ TRUE,
    TRUE~FALSE
  ),
  dissim=(1+max(response)-response) # convert similarity to dissimilarity by subtracting response from 
  ) %>% 
  unite("turtle_a",c(angle_a,radius_a),remove=F,sep="_") %>%
  unite("turtle_b",c(angle_b,radius_b),remove=F,sep="_")


# Read in stimulus values
stim <- read_csv(here("data","turtles_Spring2022","stim.csv"))

# angle-radius values as character string
ang_rad <- as.character(unite(stim,"ang_rad",c(angle,radius),sep="_")$ang_rad)

# number of unique stimuli
n_stim <- length(ang_rad)

# Dissimilarity matrix =====================================================================
# some of this adapted code from Andrew Cohen (2019)
dissim_mat <- matrix(0,nrow=n_stim,ncol=n_stim)
s1_num <- 1
for (s1 in ang_rad) {
  s2_num<-1
  for(s2 in ang_rad) {
    if (s1 == s2) {
      resp<-0
    } else {
      resp <- dat_clean[(dat_clean$turtle_a==s1 & dat_clean$turtle_b==s2) |
                        (dat_clean$turtle_a==s2 & dat_clean$turtle_b==s1), ]$dissim
    }
    dissim_mat[s1_num, s2_num]<-mean(resp)
    dissim_mat[s2_num, s1_num]<-mean(resp)
    s2_num <- s2_num+1
  }
  s1_num<-s1_num+1
} 


# MDS =====================================================================
fit <- MASS::isoMDS(dissim_mat, k=2)
fit

MDS_points <- tibble(
  angle=fit$points[,1],
  radius=fit$points[,2]
)


# PLOTTING =====================================================================
MDS_plot <- ggplot(MDS_points, aes(angle,radius))+
  geom_point(alpha=.8,col="#00BDD0",size=3.5)+
  coord_fixed(xlim=c(-4,4),ylim=c(-4,4))+
  labs(title="MDS Solution")+
  ggthemes::theme_few()+
  theme(plot.title=element_text(hjust=0.5))

orig_stim_plot <- ggplot(stim,aes(angle,radius))+
  geom_point(alpha=.8,col="#00BDD0",size=3.5)+
  coord_cartesian(xlim=c(10,40),ylim=c(30,140))+
  labs(title="Raw Stimulus Values")+
  ggthemes::theme_few()+
  theme(plot.title=element_text(hjust=0.5))

all_plot <- orig_stim_plot/MDS_plot
all_plot

if(save_plot){
  ggsave(here("R","plots","MDS_Plots.pdf"),all_plot)
}


