library(tidyverse)
library(nloptr)


# load functions

source(here::here("scripts","fcn","harvest.R"))
source(here::here("scripts","fcn","growth.R"))
source(here::here("scripts","fcn","pi.R"))
source(here::here("scripts","fcn","Payoff_post.R"))
source(here::here("scripts","fcn","Payoff_pre.R"))
source(here::here("scripts","fcn","ut_post.R"))
source(here::here("scripts","fcn","vfi_post_fcn.R"))
source(here::here("scripts","fcn","vfi_pre_fcn.R"))
source(here::here("scripts","fcn","ut.R"))
findNearest<-function(x,y){
  index<-which.min(abs(x-y))
  return(index)
}



# Set up parameters



param_grid<-expand_grid(sizex=50,n_z=13,cshape=1.3,delta=0.96,r=0.8,K=1000,p=5,c=250,mean_theta=0,sigma_theta=c(0.15,0.3),gamma=seq(0,1000,by=250),coverage=0,ut_par=c(paste0("0.5","_","power"),
                                                                                                                                                                       paste0("2","_","power"))) |> 
  separate(ut_par,into=c("a","ut_mod"),sep="_") |> 
  mutate(a=as.numeric(a))


vfi_post_out<-pmap(param_grid,vfi_post_fcn)

save(file=here::here("data","vfi_post_out_insurance_1-20-26.RData"),vfi_post_out)


## Run VFI pre-shock

vfi_pre_out<-pmap(param_grid,vfi_pre_fcn)

save(file=here::here("data","vfi_pre_out_insurance_1-29-26.RData"),vfi_pre_out)

