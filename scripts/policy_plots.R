library(tidyverse)

#load runs here
load(here::here("data","vfi_post_out_insurance_1-9-26.Rdata"))

post_df<-bind_rows(vfi_post_out)

sig_choice<-unique(post_df$parameters$sigma_theta)[2]

## Get realized theta_prop for post evaluation

shock_process <- discretize_normal_shock(15, 0, sig_choice)
w_grid <- shock_process$w_grid
w_probs <- shock_process$probs
theta_prop<-w_probs

mean_theta=0
sigma_theta=sig_choice
n_z=13
n_sims=10000

set.seed(42)

#Get distribution and bins of
theta_raw<-rnorm(n_sims,mean=mean_theta,sd=sigma_theta)
theta_dist<-theta_raw[which(theta_raw>-1)]

thetahist<-hist(theta_dist,breaks=seq(min(theta_dist),max(theta_dist),l=n_z+1),plot=FALSE)
## Make bins classification
bins=theta_dist
for(k in 2:length(thetahist$breaks)){
  index=which(theta_dist>thetahist$breaks[k-1]&theta_dist<=thetahist$breaks[k])
  bins[index]=k-1
}

theta_df<-data.frame(theta_dist=theta_dist,bins=bins)

means<- theta_df |> 
  group_by(bins) |> 
  summarize(mean=mean(theta_dist))

#elminate first observation
theta_prop=thetahist$counts/n_sims
theta_val=means$mean

# Summarize best converged gamma for each type of model

summary_post<-post_df |> 
  group_by(parameters$gamma,parameters$ut_mod,conv$b,parameters$sigma_theta) |> 
  summarize(mean=weighted.mean(conv$Vstar,w_probs)) |>
  ungroup() |> 
  rename(gamma="parameters$gamma",
         ut_mod="parameters$ut_mod",
         b="conv$b",
         mean=mean) |> 
  group_by(gamma,ut_mod) |>
  summarize(v=mean(mean)) |> 
  group_by(ut_mod) |> 
  slice_max(v)



pre_df<-bind_rows(vfi_pre_out)

summary_pre<-pre_df |> 
  group_by(parameters$gamma,parameters$ut_mod) |>
  summarize(v=mean(conv$Vstar)) |> 
  rename(gamma="parameters$gamma",
         ut_mod="parameters$ut_mod",
         v=v) |>
  group_by(ut_mod) |> 
  slice_max(v)


policy_post<-post_df |> 
  filter(parameters$gamma==0.25 & parameters$ut_mod=='power' & parameters$sigma_theta==sig_choice  & parameters$a==0.5 & parameters$trigger==-0.1)

policy_pre<-pre_df |>
  filter(parameters$gamma==0.25 & parameters$ut_mod=='power' & parameters$sigma_theta==sig_choice  & parameters$a==0.5 & parameters$trigger==-0.1)

policy_post_noi<-post_df |> 
  filter(parameters$gamma==0 & parameters$ut_mod=='power' & parameters$sigma_theta==sig_choice  & parameters$a==0.5 & parameters$trigger==-0.1)

policy_pre_noi<-pre_df |> 
  filter(parameters$gamma==0 & parameters$ut_mod=='power' & parameters$sigma_theta==sig_choice  & parameters$a==0.5 & parameters$trigger==-0.1)

vals_rescaled <- scales::rescale(w_grid)

# diverging palette
pal_fun <- scales::col_numeric(
  palette = c("red", "grey80", "blue"),
  domain  = range(w_grid)
)

# named vector for scale_color_manual
cols <- setNames(pal_fun(w_grid), levels(as.factor(w_grid)))

# pal_fun <- col_numeric(
#   palette = c("red", "grey80", "blue"),
#   domain  = c(min(theta), 0, max(vals))
# )

trig<-unique(policy_post$parameters$trigger)
trig_index<-which(w_grid==trig)
trig_end<-length(w_grid)

two_line_noi<-policy_post_noi %>% 
  mutate(pay=ifelse(conv$theta>trig,'no','yes')) |> 
  group_by(conv$b,pay) %>%
  mutate(custom_weight=case_when(
    pay=="yes"~theta_prop[1:trig_index][row_number()],
    pay=="no"~theta_prop[(trig_index+1):trig_end][row_number()]
  )) |> 
  summarize(fstar=weighted.mean(conv$fstar,w=custom_weight)) |>
  ungroup() |> 
  mutate(ins='No Insurance') |> 
  rename(b='conv$b')

two_line_i<-policy_post |> 
  mutate(pay=ifelse(conv$theta>trig,'no','yes')) |> 
  group_by(conv$b,pay) %>%
  mutate(custom_weight=case_when(
    pay=="yes"~theta_prop[1:trig_index][row_number()],
    pay=="no"~theta_prop[(trig_index+1):trig_end][row_number()]
  )) |> 
  summarize(fstar=weighted.mean(conv$fstar,w=custom_weight,na.rm=TRUE)) |>
  ungroup() |> 
  mutate(ins='Insurance') |> 
  rename(b='conv$b')


rbind(two_line_noi,two_line_i) |> 
  ggplot(aes(x=b,y=fstar,color=pay))+
  geom_line(aes(linetype = ins),linewidth=1.5)+
  ylim(0,1)+
  scale_color_manual(name='Environmental\nShock',values=c('no'='#047C90','yes'='#78A540'),labels=c('yes'='Negative Shock','no'='Positive Shock'))+
  theme_minimal()+
  labs(title=paste0('Post','- ','Risk Aversion:',policy_post$parameters$ut_mod[1],' ','Sigma:',sig_choice))

p_noi_df<-policy_pre_noi |> 
  mutate(ins='No Insurance') 

p_df<-policy_pre |> 
  mutate(ins="Insurance") 

rbind(p_noi_df,p_df) |> 
  ggplot(aes(x=conv$b,y=conv$fstar,linetype=ins))+
  geom_line(linewidth=1.5)+
  ylim(0,1)+
  theme_minimal()+
  labs(title=paste0('Pre','- ','Risk Aversion:',policy_pre$parameters$ut_mod[1],' ','Sigma:',sig_choice))


