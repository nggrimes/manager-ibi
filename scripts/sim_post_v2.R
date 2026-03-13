###sim with the post runs
library(tidyverse)
library(fields)

## Simulation parameters
n_years<-30
n_sims<-1000
b0<-50    # Choose between small and K
gamma=0.03
m=-0.5
trigger=-.01

load(here::here("data","vfi_post_out_2026-03-11.Rdata"))


## Load functions
source(here::here("scripts","fcn","fisher_ut.R"))
source(here::here("scripts","fcn","discretize_normal_shock.R"))
source(here::here("scripts","fcn","utility.R"))
source(here::here("scripts","fcn","growth2.R"))
source(here::here("scripts","fcn","insurance_payout.R"))

noi <- bind_rows(vfi_post_out) |> 
  filter(parameters$gamma==gamma) |> 
  filter(parameters$m==m) |> 
  filter(parameters$trigger==trigger)|>
  filter(parameters$pay==0)

ins<-bind_rows(vfi_post_out) |> 
  filter(parameters$gamma==gamma) |> 
  filter(parameters$m==m) |> 
  filter(parameters$trigger==trigger)|>
  filter(parameters$pay==1)

parameters<-unique(ins$parameters)

b<-unique(noi$conv$b)

noi_f<-noi$conv |> 
  select(-c(Vstar,cstar)) |> 
  pivot_wider(names_from=w,values_from=fstar) |> 
  select(-b) |> 
  as.matrix()

noi_c<-noi$conv |> 
  select(-c(Vstar,fstar)) |> 
  pivot_wider(names_from=w,values_from=cstar) |> 
  select(-b) |> 
  as.matrix()

ins_f<-ins$conv |> 
  select(-c(Vstar,cstar)) |> 
  pivot_wider(names_from=w,values_from=fstar) |> 
  select(-b) |> 
  as.matrix()

ins_c<-ins$conv |> 
  select(-c(Vstar,fstar)) |> 
  pivot_wider(names_from=w,values_from=cstar) |> 
  select(-b) |> 
  as.matrix()

# make 30 year random weather draws 500 times
set.seed(42)

shock_process <- discretize_normal_shock(parameters$n_shocks, parameters$mu_w, parameters$sigma_w)
w_grid <- shock_process$w_grid
w_probs <- shock_process$probs


prem_rate<-sum(pmax(trigger - w_grid, 0) * w_probs) * (1 + m)

rnd_weather <- replicate(n_sims, rnorm(n_years, 0, parameters$sigma_w))

rnd_weather[rnd_weather < min(w_grid)] <- min(w_grid)+0.0001
rnd_weather[rnd_weather > max(w_grid)] <- max(w_grid)-0.0001

final_ut<-rep(NA,n_sims)
final_b<-rep(NA,n_sims)

final_ut_pub_ins<-rep(NA,n_sims)


final_ut_pub_no<-rep(NA,n_sims)
final_b_pub_no<-rep(NA,n_sims)


for(j in 1:ncol(rnd_weather)){
  
  # make a data frame to store the results
  sim <- data.frame(year = 1:n_years, 
                    ins_ut = rep(NA, n_years), 
                    pub_ut = rep(NA, n_years),
                    pub_ins_ut = rep(NA, n_years),
                    ins_b= rep(NA, n_years),
                    pub_b = rep(NA, n_years))
  
  # set the initial population size
  sim$ins_b[1] <- b0
  sim$pub_b[1] <- b0
  
  # Interpolate
  f_ins <- interp.surface(list(x = b,
                               y = w_grid,
                               z = ins_f), 
                          loc = cbind(b0, rnd_weather[1,j]))
  
  f_noi<-interp.surface(list(x = b,
                             y = w_grid,
                             z = noi_f), 
                        loc = cbind(b0, rnd_weather[1,j]))
  
  c_ins <- interp.surface(list(x = b,
                               y = w_grid,
                               z = ins_c), 
                          loc = cbind(b0, rnd_weather[1,j]))
  
  c_noi<-interp.surface(list(x = b,
                             y = w_grid,
                             z = noi_c), 
                        loc = cbind(b0, rnd_weather[1,j]))
  
  
  ins_pay<-map_dbl(.x=shock,
                   .f=~insurance_payout(shock=.x,
                                        coverage=c_ins,
                                        trigger=trigger))
  ins_prem<-c_ins*prem_rate
  
  pub_pay<-map_dbl(.x=shock,
                   .f=~insurance_payout(shock=.x,
                                        coverage=c_noi,
                                        trigger=trigger))
  pub_prem<-c_noi*prem_rate
  
  pi_ins<-map2_dbl(.x=shock,.y=ins_pay,~f_ins-c_cost*f_ins/((1+.x)*sim$ins_b[1])+.y-ins_prem)
  
  sim$ins_ut[1]<-sum(map_dbl(pi_ins,~utility(.x,gamma,boost=0,mod=parameters$ut_mod))*w_probs)
  
  pi_pub_ins<-map2_dbl(.x=shock,.y=pub_pay,~f_noi-c_cost*f_noi/((1+.x)*sim$pub_b[1])+.y-ins_prem)
  
  sim$pub_ins_ut[1]<-sum(map_dbl(pi_pub_ins,~utility(.x,gamma,boost=0,mod=parameters$ut_mod))*w_probs)
  
  
  pi_pub<-map2_dbl(.x=shock,.y=0,~f_noi-c_cost*f_noi/((1+.x)*sim$pub_b[1])+.y-0)
  
  sim$pub_ut[1]<-sum(map_dbl(pi_pub,~utility(.x,gamma,boost=0,mod=parameters$ut_mod))*w_probs)
  
  
  # loop over the years
  for(i in 2:n_years){

    #print(paste0("Sim ",j," Year ",i))
    
    sim$ins_b[i] <- growth(max(sim$ins_b[i-1]*(1+rnd_weather[i-1,j])-f_ins,parameters$K/1000),parameters$r,parameters$K)
    sim$pub_b[i] <- growth(max(sim$pub_b[i-1]*(1+rnd_weather[i-1,j])-f_noi,parameters$K/1000),parameters$r,parameters$K)
    
    # just make sure we don't go outside our spline bounds
    if(sim$ins_b[i]<min(b)) {
      sim$ins_b[i]<-min(b)+0.0001
    } 
    
    
    if(sim$pub_b[i]<min(b)) {
      sim$pub_b[i]<-min(b)+0.0001
    } 
    
    f_ins <- interp.surface(list(x = b,
                                 y = w_grid,
                                 z = ins_f), loc = cbind(sim$ins_b[i], rnd_weather[i,j]))
    
    f_noi<-interp.surface(list(x = b,
                               y = w_grid,
                               z = noi_f), 
                          loc = cbind(sim$pub_b[i], rnd_weather[i,j]))
    
    c_ins <- interp.surface(list(x = b,
                                 y = w_grid,
                                 z = ins_c), 
                            loc = cbind(sim$ins_b[i], rnd_weather[i,j]))
    
    c_noi<-interp.surface(list(x = b,
                               y = w_grid,
                               z = noi_c), 
                          loc = cbind(sim$pub_b[i], rnd_weather[i,j]))
    
    # get insurance utility
  
    
    
    ins_pay<-map_dbl(.x=shock,
                     .f=~insurance_payout(shock=.x,
                                          coverage=c_ins,
                                          trigger=trigger))
    ins_prem<-c_ins*prem_rate
    
    pub_pay<-map_dbl(.x=shock,
                     .f=~insurance_payout(shock=.x,
                                          coverage=c_noi,
                                          trigger=trigger))
    pub_prem<-c_noi*prem_rate
    
    pi_ins<-map2_dbl(.x=shock,.y=ins_pay,~f_ins-c_cost*f_ins/((1+.x)*sim$ins_b[i])+.y-ins_prem)
    
    sim$ins_ut[i]<-sum(map_dbl(pi_ins,~utility(.x,gamma,boost=0,mod=parameters$ut_mod))*w_probs)
    
    pi_pub_ins<-map2_dbl(.x=shock,.y=pub_pay,~f_noi-c_cost*f_noi/((1+.x)*sim$pub_b[i])+.y-pub_prem)
    
    sim$pub_ins_ut[i]<-sum(map_dbl(pi_pub_ins,~utility(.x,gamma,boost=0,mod=parameters$ut_mod))*w_probs)
    
    
    pi_pub<-map2_dbl(.x=shock,.y=0,~f_noi-c_cost*f_noi/((1+.x)*sim$pub_b[i])+.y-0)
    
    sim$pub_ut[i]<-sum(map_dbl(pi_pub,~utility(.x,gamma,boost=0,mod=parameters$ut_mod))*w_probs)
    
  }
  
  # discount the utilities
  sim$ins_ut_disc <- sim$ins_ut * parameters$beta^(0:(n_years-1))
  sim$pub_ut_disc <- sim$pub_ut * parameters$beta^(0:(n_years-1))
  sim$pub_ins_ut_disc<-sim$pub_ins_ut * parameters$beta^(0:(n_years-1))
  
  # store the final biomass and summed discounted utilities in storage vec
  
  
  final_ut[j]<-sum(sim$ins_ut_disc)
  final_b[j]<-sim$ins_b[n_years]
  
  final_ut_pub_no[j]<-sum(sim$pub_ut_disc)
  final_b_pub_no[j]<-sim$pub_b[n_years]
  
  final_ut_pub_ins[j]<-sum(sim$pub_ins_ut_disc)
  
  
}


df<-data.frame(final_ut, final_b, final_ut_pub_no, final_b_pub_no,final_ut_pub_ins,sim=1:n_sims,b0=b0) |> 
  drop_na() |> 
  pivot_longer(cols=c(final_ut, final_b, final_ut_pub_no, final_b_pub_no,final_ut_pub_ins), names_to='variable', values_to='value')


df |> filter(variable %in% c('final_ut', 'final_ut_pub_no','final_ut_pub_ins')) |> 
  ggplot(aes(x=value, fill=variable,color=variable)) +
  geom_density(alpha=0.33,linewidth=2.5) +
  theme_minimal() +
  scale_color_manual(name="",values=c("#003660","#79A540",'#900C3F'),guide='none')+
  scale_fill_manual(name="",values=c("#003660","#79A540",'#900C3F'),labels=c("final_ut"="Insurance","final_ut_pub_no"="No Insurance","final_ut_pub_ins"="Fisher Insurance"))+
  theme_classic()+
  theme(legend.text=element_text(size=24))+
  theme(legend.title =element_text(size=28))+
  theme(axis.text =element_text(size=22))+
  theme(axis.title = element_text(size=26))+
  labs(title='', x='Sum of discounted utility', y='Density')+
  scale_y_continuous(expand = c(0,0))+
  guides(fill = guide_legend(byrow=TRUE, override.aes = list(shape = NA, color = NA, fill = c("#003660", "#79A540",'#900C3F'))))+
  theme(legend.spacing.y = unit(10, "cm"),
        axis.text = element_text(size=txt_size,family = 'serif'),
        axis.title = element_text(size=txt_size+3,family='serif'),
        legend.text = element_text(size=txt_size+3,family='serif'),
        legend.position = c(0.15, 0.8))

ggsave(here::here('fig','post_shock_ut_b_ra03_m5.png'),dpi=300,width=7,height=7)


df |> filter(variable %in% c('final_b', 'final_b_pub_no')) |> 
  ggplot(aes(x=value, fill=variable, color=variable)) +
  geom_density(aes(y = after_stat(density * length(df$value))),alpha=0.33,linewidth=2.5) +
  theme_minimal() +
  scale_color_manual(name="",values=c("#003660","#79A540"),guide='none')+
  scale_fill_manual(name="",values=c("#003660","#79A540"),labels=c("final_b"="Insurance","final_b_pub_no"="No Insurance"))+
  theme_classic()+
  theme(legend.text=element_text(size=24))+
  theme(legend.title =element_text(size=28))+
  theme(axis.text =element_text(size=22))+
  theme(axis.title = element_text(size=26))+
  labs(title='', x='Final Biomass', y='Density')+
  scale_y_continuous(expand = c(0,0))+
  guides(fill = guide_legend(byrow=TRUE, override.aes = list(shape = NA, color = NA, fill = c("#003660", "#79A540"))))+
  theme(legend.spacing.y = unit(10, "cm"),
        axis.text = element_text(size=txt_size,family = 'serif'),
        axis.title = element_text(size=txt_size+3,family='serif'),
        legend.text = element_text(size=txt_size+3,family='serif'),
        legend.position = c(0.25, 0.8))


ggsave(here::here('fig','post_shock_sim_b_ra03_m5.png'),dpi=300,width=7,height=7)

