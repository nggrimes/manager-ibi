### Run a forward simulation using the harvest control rules

library(tidyverse)

## Simulation parameters
n_years<-30
n_sims<-500
b0<-500    # Choose between small and K


## Load functions
source(here::here("scripts","fcn","harvest.R"))
source(here::here("scripts","fcn","growth.R"))
source(here::here("scripts","fcn","pi.R"))
source(here::here("scripts","fcn","ut_post.R"))

load(here::here('data','vfi_pre_out_insurance_1-14-26.Rdata'))


noi <- bind_rows(vfi_pre_out) |> 
  filter(parameters$gamma==0) |> 
  filter(parameters$sigma_theta==0.5) |> 
  filter(parameters$ut_mod=='cara')

ins<-bind_rows(vfi_pre_out) |> 
  filter(parameters$gamma==250) |> 
  filter(parameters$sigma_theta==0.5)|> 
  filter(parameters$ut_mod=='cara')

parameters<-unique(ins$parameters)

b<-unique(noi$conv$b)


noi_pol<-noi$conv |> 
  select(-Vstar)

ins_pol<-ins$conv |> 
  select(-Vstar)



# make 30 year random weather draws 500 times
set.seed(42)


rnd_weather <- replicate(n_sims, rnorm(n_years, 0, parameters$sigma_theta))

#replace any values of rnd_Weather less than -1 with -0.99

rnd_weather[rnd_weather < -1] <- -0.99

final_ut<-rep(NA,500)
final_b<-rep(NA,500)

final_ut_no<-rep(NA,500)
final_b_no<-rep(NA,500)

final_ut_rn<-rep(NA,500)
final_b_rn<-rep(NA,500)

for(j in 1:ncol(rnd_weather)){
  
  # make a data frame to store the results
  sim <- data.frame(year = 1:30, 
                    ins_ut = rep(NA, 30), 
                    no_ut = rep(NA, 30),
                    ins_b= rep(NA, 30),
                    no_b = rep(NA, 30),
                    rn_ut = rep(NA, 30),
                    rn_b = rep(NA, 30))
  
  # set the initial population size
  sim$ins_b[1] <- b0
  sim$no_b[1] <- b0
  sim$rn_b[1] <- b0
  
  f_ins<-spline(x=b,y=ins_pol$fstar,xout=sim$ins_b[1])$y
  f_no<-spline(x=b,y=noi_pol$fstar,xout=sim$no_b[1])$y

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

    
    sim$ins_b[i] <- growth(sim$ins_b[i-1],f_ins,rnd_weather[i-1,j],parameters$r,parameters$K)
    sim$no_b[i] <- growth(sim$no_b[i-1],f_noi,rnd_weather[i-1,j],parameters$r,parameters$K)
    
    
    f_ins<-spline(x=b,y=ins_pol$fstar,xout=sim$ins_b[i])$y
    f_no<-spline(x=b,y=noi_pol$fstar,xout=sim$no_b[i])$y
    
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

# Run a t.test to compare if the distributions are different
t.test(final_ut ~ final_ut_no, data=df |> pivot_wider(names_from='variable', values_from='value'))




plot_df<-rbind(no_pol_conv %>% 
                 mutate(ins='No Insurance'),
               ins_pol_conv %>% 
                 mutate(ins='Insurance'),
               rn_pol_conv %>% 
                 mutate(ins='Risk Neutral'))




ggplot(plot_df |> filter(ins!='Risk Neutral'),aes(x=b,y=pol_opt,color=ins))+
  geom_line(linewidth=2.5)+
  scale_color_manual(name="",values=c("#003660","#79A540"))+
  theme_classic()+
  labs(y='Optimal Harvest Rate',x='Biomass')+
  theme(legend.text=element_text(size=24))+
  theme(legend.title =element_text(size=28))+
  theme(axis.text =element_text(size=22))+
  theme(axis.title = element_text(size=26))


ggplot(plot_df |> filter(ins!='Insurance'),aes(x=b,y=pol_opt,color=ins))+
  geom_line(linewidth=2.5)+
  scale_color_manual(name="",values=c("#79A540","#900C3F"))+
  theme_classic()+
  labs(y='Optimal Harvest Rate',x='Biomass')+
  theme(legend.text=element_text(size=24))+
  theme(legend.title =element_text(size=28))+
  theme(axis.text =element_text(size=22))+
  theme(axis.title = element_text(size=26))

ggplot(plot_df |> filter(ins=='Risk Neutral'),aes(x=b,y=pol_opt,color=ins))+
  geom_line(linewidth=2.5)+
  scale_color_manual(name="",values=c("#900C3F"))+
  theme_classic()+
  labs(y='Optimal Harvest Rate',x='Biomass')+
  theme(legend.text=element_text(size=24))+
  theme(legend.title =element_text(size=28))+
  theme(axis.text =element_text(size=22))+
  theme(axis.title = element_text(size=26))

