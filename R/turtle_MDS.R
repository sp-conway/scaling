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
library(readr)

# controls
save_plot <- T

# The MASS package is used for running MDS, but loading it masks some dplyr functions I need
if("MASS" %in% (.packages())){
  detach("package:MASS",unload=T)
}

# data cleaning =====================================================================

# Read in and clean raw data files
files <- dir_ls(here("data","turtles_Spring2022","raw"))
dat <- map_dfr(files, ~data.table::fread(.x) %>% as_tibble()) %>% 
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

n_subs <- length(unique(dat$sub_n))
cat("\n===============\n",n_subs,"subjects\n===============\n")

# clean data for MDS stuff
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
stim <- readr::read_csv(here("data","turtles_Spring2022","stim.csv")) %>% 
  mutate(which=1:n()) 

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
mds_fit <- MASS::isoMDS(dissim_mat)
mds_fit

MDS_points <- tibble(
  angle=mds_fit$points[,1],
  radius=mds_fit$points[,2]
)

# Get mds matrix (original)
mds_mat <- matrix(c(MDS_points$angle,
                    MDS_points$radius),
                  ncol=2)

# Functions for transforming 
# Transform function
mds_transform <- function(mds_mat, s, theta, t_x, t_y) {
  
  # Reflection matrix
  ref_mat <- matrix(c(0,1,1,0), nrow=2, ncol=2, byrow=TRUE)
  
  # Rotation matrix
  rot_mat <- matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)), nrow=2, ncol=2, byrow=TRUE)
  
  # Reflection 1 
  mds_mat_Ref <- mds_mat %*% ref_mat
  
  # Reflection 2
  mds_mat_Ref2 <- matrix(c(mds_mat_Ref[,1],-mds_mat_Ref[,2]),ncol=2)
  
  # Rotation
  mds_mat_r <- mds_mat_Ref2 %*% rot_mat
  
  # Scaling
  mds_mat_s_r <- mds_mat_r*s
  
  # Translate
  mds_mat_t_s_r <- mds_mat_s_r + c(t_x, t_y)
  
  # Return transformed matrix
  return(mds_mat_t_s_r)
  
}

# get ss
cost <- function(params, stim_mat, mds_mat){
  s <- params[1]
  theta <- params[2]
  t_x <- params[3]
  t_y <- params[4]
  mds_transformed <- mds_transform(mds_mat, s, theta, t_x, t_y)
  ss <- sum((stim_mat-mds_transformed)^2)
  return(ss)
}

# starting params
s_params <- c(0, 0, 1, 0)

# fit model
fit <- optim(par=s_params, 
             fn=cost, 
             method='L-BFGS-B',
             lower=c(-100, -100, -10, -2*pi), 
             upper=c( 100,  100,  10,  2*pi),
             stim_mat=matrix(c(stim$angle,stim$radius),ncol=2),
             mds_mat=mds_mat)

# PLOTTING =====================================================================
MDS_mat_transf <- mds_transform(mds_mat, fit$par[1], fit$par[2], fit$par[3], fit$par[4])
MDS_points_transf <- tibble(
  angle=MDS_mat_transf[,1],
  radius=MDS_mat_transf[,2],
  which=1:nrow(MDS_mat_transf)
)

MDS_plot_orig <- MDS_points %>%
  mutate(which=1:n()) %>%
  ggplot(aes(angle,radius))+
  geom_text(aes(label=which),alpha=.8,size=3)+
  coord_fixed(xlim=c(-5,5),ylim=c(-5,5))+
  labs(title="MDS Solution original")+
  ggthemes::theme_few()+
  theme(plot.title=element_text(hjust=0.5))

MDS_plot_transf <- MDS_points_transf %>%
  mutate(which=1:n()) %>%
  ggplot(aes(angle,radius))+
  geom_text(aes(label=which),alpha=.8,size=3)+
  coord_cartesian(xlim=range(MDS_points_transf$angle),
                  ylim=range(MDS_points_transf$radius))+
  labs(title="MDS Solution transformed")+
  ggthemes::theme_few()+
  theme(plot.title=element_text(hjust=0.5))

orig_stim_plot <- stim %>%
  ggplot(aes(angle,radius))+
  geom_text(aes(label=which),alpha=.8,size=3)+
  coord_cartesian(xlim=range(stim$angle),
                  ylim=range(stim$radius))+
  labs(title="Raw Stimulus Values")+
  ggthemes::theme_few()+
  theme(plot.title=element_text(hjust=0.5))

all_plot <- orig_stim_plot|(MDS_plot_orig/MDS_plot_transf)
all_plot

# Save results ================================================================================
# save plot
if(save_plot){
  ggsave(here("R","plots","MDS_Plots.pdf"),all_plot,
         width=8,height=8)
}

# save transformed mds coords
write_csv(MDS_points_transf, file=here("R","mds_results","mds_turtles.csv"))
write_csv(stim,file=here("R","mds_results","orig_turtles.csv"))


