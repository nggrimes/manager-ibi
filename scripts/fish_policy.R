library(tidyverse)

source(here::here("scripts","fcn","fisher_ut.R"))
source(here::here("scripts","fcn","discretize_normal_shock.R"))
source(here::here("scripts","fcn","payoff_fcn.R"))
source(here::here("scripts","fcn","vfi_pre_fish_fcn.R"))
source(here::here("scripts","fcn","utility.R"))
source(here::here("scripts","fcn","growth2.R"))
source(here::here("scripts","fcn","insurance_payout.R"))
source(here::here("scripts","fcn","payoff_post_fcn.R"))
source(here::here("scripts","fcn","vfi_post_fish_fcn.R"))




param_grid<-expand_grid(n_shocks=13,mu_w=0,sigma_w=c(0.2),m=c(-.8,-.5),trigger=c(-0.01,-0.15),n_biomass=35,K=100,r=0.8,c_cost=15,pay=c(0,1),beta=0.95,ut_par=c(paste0("0.3","_","cara"),
                                                                                                                                                                                paste0("0.03","_","cara"))) |> 
  separate(ut_par,into=c("gamma","ut_mod"),sep="_") |> 
  mutate(gamma=as.numeric(gamma))


vfi_pre_out<-pmap(param_grid,vfi_pre_fish_fcn)


vfi_post_out<-pmap(param_grid,vfi_post_fish_fcn,.progress=TRUE)


save_date<-Sys.Date()
save_name<-paste0('vfi_post_out_',save_date,'.Rdata')

save(vfi_post_out,file=here::here('data',save_name))

save_name<-paste0('vfi_pre_out_',save_date,'.Rdata')

save(vfi_pre_out,file=here::here('data',save_name))
