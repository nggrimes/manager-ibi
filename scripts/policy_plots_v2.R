# policy_plots v2

library(tidyverse)

# parameters
m=0.05
gamma=0.03
trigger=-0.01


#load runs here
load(here::here("data","vfi_post_out_2026-03-11.Rdata"))
load(here::here("data","vfi_pre_out_2026-03-11.Rdata"))

source(here::here("scripts","fcn","discretize_normal_shock.R"))

post_df<-bind_rows(vfi_post_out)

sig_choice<-unique(post_df$parameters$sigma_w)[1]

## Get realized theta_prop for post evaluation

shock_process <- discretize_normal_shock(post_df$parameters$n_shocks[1], 0, sig_choice)
w_grid <- shock_process$w_grid
w_probs <- shock_process$probs



pre_df<-bind_rows(vfi_pre_out)



policy_post<-post_df |> 
  filter(parameters$gamma==gamma & parameters$ut_mod=='cara' & parameters$pay==1  & parameters$m==m & parameters$trigger==trigger)

policy_pre<-pre_df |>
  filter(parameters$gamma==gamma & parameters$ut_mod=='cara' & parameters$pay==1  & parameters$m==m & parameters$trigger==trigger)

policy_post_noi<-post_df |> 
  filter(parameters$gamma==gamma& parameters$ut_mod=='cara' & parameters$pay==0  & parameters$m==m & parameters$trigger==trigger)

policy_pre_noi<-pre_df |> 
  filter(parameters$gamma==gamma& parameters$ut_mod=='cara' & parameters$pay==0  & parameters$m==m & parameters$trigger==trigger)

vals_rescaled <- scales::rescale(w_grid)

# diverging palette
pal_fun <- scales::col_numeric(
  palette = c("red", "grey80", "blue"),
  domain  = range(w_grid)
)

# named vector for scale_color_manual
cols <- setNames(pal_fun(w_grid), levels(as.factor(w_grid)))

## text sizes

txt_size<-12

trig<-unique(policy_post$parameters$trigger)
trig_index<-which(w_grid<trig)
trig_end<-length(w_grid)
index_end<-length(trig_index)

two_line_noi<-policy_post_noi %>% 
  mutate(payout=ifelse(conv$w>trig,'no','yes')) |> 
  group_by(conv$b,payout) %>%
  mutate(custom_weight=case_when(
    payout=="yes"~w_probs[trig_index][row_number()],
    payout=="no"~w_probs[(index_end+1):trig_end][row_number()]
  )) |> 
  summarize(fstar=weighted.mean(conv$fstar,w=custom_weight)) |>
  ungroup() |> 
  mutate(ins='No Insurance') |> 
  rename(b='conv$b')

two_line_i<-policy_post |> 
  mutate(payout=ifelse(conv$w>trig,'no','yes')) |> 
  group_by(conv$b,payout) %>%
  mutate(custom_weight=case_when(
    payout=="yes"~w_probs[trig_index][row_number()],
    payout=="no"~w_probs[(index_end+1):trig_end][row_number()]
  )) |> 
  summarize(fstar=weighted.mean(conv$fstar,w=custom_weight)) |>
  ungroup() |> 
  mutate(ins='Insurance') |> 
  rename(b='conv$b')


rbind(two_line_noi,two_line_i) |> 
  ggplot(aes(x=b,y=fstar,color=payout))+
  geom_line(aes(linetype = ins),linewidth=1.5)+
  scale_color_manual(name='',values=c('no'='#047C90','yes'='#78A540'),labels=c('yes'='Negative Shock','no'='Positive Shock'))+
  theme_classic()+
  labs(title='',y="Optimal Quota\n",x='Observed Biomass',linetype='')+
  scale_y_continuous(expand=c(0,0))+
  theme(axis.text = element_text(size=txt_size,family = 'serif'),
        axis.title = element_text(size=txt_size+3,family='serif'),
        legend.text = element_text(size=txt_size+3,family='serif'),
        legend.position = c(0.25, 0.8))

ggsave(here::here('fig','post_shock_ra03_m05.png'),dpi=300,width=7,height=7)

## % diff in harvest
pct_diff<-two_line_i|>mutate(diff=(two_line_i$fstar-two_line_noi$fstar)/two_line_noi$fstar)

pct_diff|>
  filter(fstar>5)|>
  ggplot(aes(x=b,y=diff,color=payout))+
  geom_line(linewidth=1.5)+
  scale_color_manual(name='',values=c('no'='#047C90','yes'='#78A540'),labels=c('yes'='Negative Shock','no'='Positive Shock'))+
  scale_y_continuous(labels=scales::percent,expand=c(0,0))+
  labs(y='Percent difference in Quota',x='Observed Biomass')+
  theme_classic()+
  theme(axis.text = element_text(size=txt_size,family = 'serif'),
        axis.title = element_text(size=txt_size+3,family='serif'),
        legend.text = element_text(size=txt_size+3,family='serif'),
        legend.position = c(0.8, 0.85))

ggsave(here::here('fig','post_shock_diff_ra03_m05.png'),dpi=300,width=7,height=7)

p_noi_df<-policy_pre_noi |> 
  mutate(ins='No Insurance') 

p_df<-policy_pre |> 
  mutate(ins="Insurance") 

rbind(p_noi_df,p_df) |> 
  ggplot(aes(x=conv$b,y=conv$fstar,linetype=ins))+
  geom_line(linewidth=1.5)+
  theme_classic()+
  labs(title='',y="Optimal Quota\n",x='Observed Biomass',linetype='')+
  scale_y_continuous(expand=c(0,0))+
  theme(axis.text = element_text(size=txt_size,family = 'serif'),
        axis.title = element_text(size=txt_size+3,family='serif'),
        legend.text = element_text(size=txt_size+3,family='serif'),
        legend.position = c(0.8, 0.25))


ggsave(here::here('fig','pre_shock_ra03_m05.png'),dpi=300,width=7,height=7)


### View coverage choices

# ex-ante

ex_cov<-pre_df|>
  filter(parameters$ut_mod=='cara' & parameters$trigger==trigger & parameters$m!= -0.8 )|>
mutate(m= case_when(parameters$m==-0.2 & parameters$gamma==0.003~ -0.5,
                    TRUE~ parameters$m))|>
  filter(m!= -0.2)

new_names <- c("0.003" = "Low Risk Aversion", "0.03" = "Mild Risk Aversion", "0.3"="High Risk Aversion")

ggplot(ex_cov)+
  geom_line(aes(x=conv$b,y=conv$cstar,color=factor(m)),linewidth=1.5)+
  facet_wrap(~parameters$gamma,labeller=as_labeller(new_names))+
  theme_bw()+
  scale_color_manual(values=c('#047C90','#78A540','#900C3F'))+
  labs(x='Observed Biomass\n',y='Optimal Insurance Coverage',color='Premium Loading Factor (m)')+
  scale_y_continuous(expand=c(0,0),limits = c(0,105))+
  theme(axis.text = element_text(size=txt_size,family = 'serif'),
        axis.title = element_text(size=txt_size+1,family='serif'),
        legend.text = element_text(size=txt_size+1,family='serif'),
        legend.title = element_text(size=txt_size+2,family='serif'),
        strip.text.x = element_text(size = txt_size+2,family='serif'),
        legend.position = "bottom",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))

ggsave(here::here('fig','pre_shock_coverage.png'),dpi=300,width=7,height=7)


post_cov<-  post_df|>
  filter(parameters$ut_mod=='cara' & parameters$trigger==trigger & parameters$m!= -0.8 & parameters$pay==0 )|>
  mutate(m= case_when(parameters$m==-0.2 & parameters$gamma==0.003~ -0.5,
                      TRUE~ parameters$m))|>
  filter(m!= -0.2)|>
  group_by(parameters$gamma,parameters$ut_mod,conv$b,parameters$sigma_theta,m) |> 
  summarize(mean=weighted.mean(conv$cstar,w_probs))|>
  rename('gamma'='parameters$gamma',
         'b'='conv$b')

ggplot(post_cov)+
  geom_line(aes(x=b,y=mean,color=factor(m)),linewidth=1.5)+
  facet_wrap(~gamma,labeller=as_labeller(new_names))+
  theme_bw()+
  scale_color_manual(values=c('#047C90','#78A540','#900C3F'))+
  labs(x='Observed Biomass\n',y='Optimal Insurance Coverage',color='Premium Loading Factor (m)')+
  scale_y_continuous(expand=c(0,0),limits = c(0,105))+
  theme(axis.text = element_text(size=txt_size,family = 'serif'),
        axis.title = element_text(size=txt_size+1,family='serif'),
        legend.text = element_text(size=txt_size+1,family='serif'),
        legend.title = element_text(size=txt_size+2,family='serif'),
        strip.text.x = element_text(size = txt_size+2,family='serif'),
        legend.position = "bottom",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))

ggsave(here::here('fig','post_shock_coverage.png'),dpi=300,width=7,height=7)
