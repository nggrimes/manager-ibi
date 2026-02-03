###sim with the post runs
library(tidyverse)
library(fields)

## Simulation parameters
n_years<-30
n_sims<-500
b0<-500    # Choose between small and K


load(here::here("data","vfi_post_out_insurance_1-14-26.Rdata"))

## Load functions
source(here::here("scripts","fcn","harvest.R"))
source(here::here("scripts","fcn","growth.R"))
source(here::here("scripts","fcn","pi.R"))
source(here::here("scripts","fcn","ut_post.R"))

noi <- bind_rows(vfi_post_out) |> 
  filter(parameters$gamma==0) |> 
  filter(parameters$sigma_theta==0.2) |> 
  filter(parameters$ut_mod=='log')

ins<-bind_rows(vfi_post_out) |> 
  filter(parameters$gamma==250) |> 
  filter(parameters$sigma_theta==0.2)|> 
  filter(parameters$ut_mod=='log')

parameters<-unique(ins$parameters)

b<-unique(noi$conv$b)
theta<-unique(noi$conv$theta)

noi_pol<-noi$conv |> 
  select(-Vstar) |> 
  pivot_wider(names_from=theta,values_from=fstar) |> 
  select(-b) |> 
  as.matrix()

ins_pol<-ins$conv |> 
  select(-Vstar) |> 
  pivot_wider(names_from=theta,values_from=fstar) |> 
  select(-b) |> 
  as.matrix()



# make 30 year random weather draws 500 times
set.seed(42)

rnd_weather <- replicate(n_sims, rnorm(n_years, 0, parameters$sigma_theta))

rnd_weather[rnd_weather < min(theta)] <- min(theta)+0.0001
rnd_weather[rnd_weather > max(theta)] <- max(theta)-0.0001

final_ut<-rep(NA,n_sims)
final_b<-rep(NA,n_sims)

final_ut_no<-rep(NA,n_sims)
final_b_no<-rep(NA,n_sims)

for(j in 1:ncol(rnd_weather)){
  

  # make a data frame to store the results
  sim <- data.frame(year = 1:n_years, 
                    ins_ut = rep(NA, n_years), 
                    no_ut = rep(NA, n_years),
                    ins_b= rep(NA, n_years),
                    no_b = rep(NA, n_years))
  
  # set the initial population size
  sim$ins_b[1] <- b0
  sim$no_b[1] <- b0
  
  interp_obj <- list(x = b,
                     y = theta,
                     z = ins_pol)  # keep as matrix
  
  # Interpolate
  f_ins <- interp.surface(interp_obj, loc = cbind(b0, rnd_weather[1,j]))
  
  f_noi<-interp.surface(list(x = b,
                             y = theta,
                             z = noi_pol), 
                        loc = cbind(b0, rnd_weather[1,j]))
  
  # get insurance utility
  temp_pi<-pi(f_ins,
              sim$ins_b[1],
              rnd_weather[1,j],
              parameters$p,
              parameters$c,
              parameters$cshape,
              parameters$trigger,
              parameters$premium,
              parameters$gamma)
  
  sim$ins_ut[1]<-ut_post(temp_pi,
                         ra=parameters$a,
                         ut_mod=parameters$ut_mod)
  
  # get no insurance utility
  temp_pi<-pi(f_noi,
              sim$ins_b[1],
              rnd_weather[1,j],
              parameters$p,
              parameters$c,
              parameters$cshape,
              parameters$trigger,
              premium=0,
              gamma=0)
  
  
  sim$no_ut[1]<-ut_post(temp_pi,
                         ra=parameters$a,
                         ut_mod=parameters$ut_mod)  


  # loop over the years
  for(i in 2:n_years){

    #print(paste0("Sim ",j," Year ",i))
    
    sim$ins_b[i] <- growth(sim$ins_b[i-1],f_ins,rnd_weather[i-1,j],parameters$r,parameters$K)
    sim$no_b[i] <- growth(sim$no_b[i-1],f_noi,rnd_weather[i-1,j],parameters$r,parameters$K)
    
    
    # just make we don't go outside our spline bounds
    if(sim$ins_b[i]<min(b)) {
      sim$ins_b[i]<-10
    } 

      
    if(sim$no_b[i]<min(b)) {
        sim$no_b[i]<-10
      } 
    f_ins<-interp.surface(list(x = b,
                               y = theta,
                               z = ins_pol), 
                          loc = cbind(sim$ins_b[i], rnd_weather[i,j]))
    
    f_noi<-interp.surface(list(x = b,
                               y = theta,
                               z = noi_pol), 
                          loc = cbind(sim$no_b[i], rnd_weather[i,j]))

    # get insurance utility
    temp_pi<-pi(f_ins,
                sim$ins_b[i],
                rnd_weather[i,j],
                parameters$p,
                parameters$c,
                parameters$cshape,
                parameters$trigger,
                parameters$premium,
                parameters$gamma)
    
    sim$ins_ut[i]<-ut_post(temp_pi,
                           ra=parameters$a,
                           ut_mod=parameters$ut_mod)
    
    # get no insurance utility
    temp_pi<-pi(f_noi,
                sim$ins_b[i],
                rnd_weather[i,j],
                parameters$p,
                parameters$c,
                parameters$cshape,
                parameters$trigger,
                premium=0,
                gamma=0)
    
    
    sim$no_ut[i]<-ut_post(temp_pi,
                           ra=parameters$a,
                           ut_mod=parameters$ut_mod) 
    # get the next year's population size
    
    
  }
  
  # discount the utilities
  sim$ins_ut <- sim$ins_ut * parameters$delta^(0:(n_years-1))
  sim$no_ut <- sim$no_ut * parameters$delta^(0:(n_years-1))
  
  # store the final biomass and summed discounted utilities in storage vec
  
  final_ut[j] <- sum(sim$ins_ut)
  final_b[j] <- sim$ins_b[n_years]
  
  final_ut_no[j] <- sum(sim$no_ut)
  
  final_b_no[j] <- sim$no_b[n_years]
  
  
}


df<-data.frame(final_ut, final_b, final_ut_no, final_b_no,sim=1:n_sims,b0=b0) |> 
  drop_na() |> 
  pivot_longer(cols=c(final_ut, final_b, final_ut_no, final_b_no), names_to='variable', values_to='value')


df |> filter(variable %in% c('final_ut', 'final_ut_no')) |> 
  ggplot(aes(x=value, fill=variable,color=variable)) +
  geom_density(alpha=0.33,linewidth=2.5) +
  theme_minimal() +
  scale_color_manual(name="",values=c("#003660","#79A540"),guide='none')+
  scale_fill_manual(name="",values=c("#003660","#79A540"),labels=c("final_ut"="Insurance","final_ut_no"="No Insurance"))+
  theme_classic()+
  theme(legend.text=element_text(size=24))+
  theme(legend.title =element_text(size=28))+
  theme(axis.text =element_text(size=22))+
  theme(axis.title = element_text(size=26))+
  labs(title='', x='Sum of discounted utility', y='Density')+
  scale_y_continuous(expand = c(0,0))+
  guides(fill = guide_legend(byrow=TRUE, override.aes = list(shape = NA, color = NA, fill = c("#003660", "#79A540"))))+
  theme(legend.spacing.y = unit(10, "cm"))

df |> filter(variable %in% c('final_b', 'final_b_no')) |> 
  ggplot(aes(x=value, fill=variable, color=variable)) +
  geom_density(aes(y = after_stat(density * length(df$value))),alpha=0.33,linewidth=2.5) +
  theme_minimal() +
  scale_color_manual(name="",values=c("#003660","#79A540"),guide='none')+
  scale_fill_manual(name="",values=c("#003660","#79A540"),labels=c("Insurance","No Insurance"))+
  theme_classic()+
  theme(legend.text=element_text(size=24))+
  theme(legend.title =element_text(size=28))+
  theme(axis.text =element_text(size=22))+
  theme(axis.title = element_text(size=26))+
  labs(title='', x='Final Biomass', y='Density')+
  scale_y_continuous(expand = c(0,0))+
  guides(fill = guide_legend(byrow=TRUE, override.aes = list(shape = NA, color = NA, fill = c("#003660", "#79A540"))))+
  theme(legend.spacing.y = unit(10, "cm"))




