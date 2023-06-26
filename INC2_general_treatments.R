## SOIL INCUBATION CALCULATIONS, JAN 2023
## Author: Michelle S. Wang, michelle.s.wang.th@dartmouth.edu

# Load packages + functions
library(tidyverse)
# library(SoilR)
library(FME)
library(ggpubr)

# Read in data
data <- read.csv("IRGA_Measurements.csv", stringsAsFactors = FALSE, header = TRUE) # scan in document formatted like example
last_day <- max(data$Day, na.rm = TRUE) # [days] final day of measurement for this datasheet

############################################################################
# GENERAL CALCULATIONS #######################################################
# Constants
R <- 82.05746   # [mL*atm/(K*mol)]

# Room Parameters
Pr <- .98 # [atm]
Tr <- 22 + 273  # [K]

# Jar/Soil Parameters
Vjar_P <- 473.176 - 46  # [mL] pint jar - filled sample cup, from Google Sheet 'Incubation Initializations <- Bulk Density' 

# n [mol] air inside jar
n_P <- (Pr*Vjar_P) / (R*Tr)     # [mol] Palouse

# Moles/Mass of C inside jar
molmass_C <- 12.011*10^3  # [mg/mol] molar mass of C
data_C <- data %>%
  mutate(moles_C_P = C_ppm*n_P/(10^6)) %>%  # [mol] moles of C in air in jar
  mutate(mass_C_P = moles_C_P*molmass_C)    # [mg] mg of C in air in jar

# Removed air inside syringe
Vrem <- 30      # [mL] CO2 rich air removed from jar
Trem <- 25 + 273  # [K] temp of air removed since in incubator
nrem <- (Pr*Vrem) / (R*Trem)  # [mol] moles of air removed from jar

# CLEAN DATA #######################################################
# Flux
num_labs <- c('DASE_C', 'DASE_O', 'AD_S', 'POET_S', 'NREL_S', 'AD_N', "POET_N", 'NREL_N', 'CS_N', 'GWC16', 'GWC20')
names(num_labs) <- c('1', '2', '3', '4','5', '6', '7', '8', '9', 'S1', 'S2')

data_all <- data_C %>%
  filter(!Sample %in% c('CO2 FREE', '2008', '2%')) %>% # filters out controls 
  separate(Sample, c("Num", "Lett"), sep=cumsum(c(1,1)), remove = FALSE) %>%  # separates out Palouse/Vershire soil and treatments
  mutate(Date.Time = as.POSIXct(Date.Time, format = '%m/%d/%y %H:%M')) %>% # converts Date.Time from characters to date-time format
  mutate(Date = as.Date(Date)) %>%
  group_by(Sample, Flush) %>%  # group by flush 
  arrange(Date.Time) %>%  # arrange in ascending order
  mutate(time_diff = as.numeric(Date.Time - lag(Date.Time, default = first(Date.Time)), units = 'hours')) %>% # [hours] find time difference in flush groups 
  mutate(mass_diff = as.numeric(mass_C_P - lag(mass_C_P, default = first(mass_C_P)))) %>% # [mg] find mass_C difference in flush groups for Palouse and Vershire
  mutate(rem_moles_C = C_ppm*nrem/(10^6)) %>%  # [mol] moles of C in removed air
  mutate(rem_mass_C = rem_moles_C*molmass_C)  %>% # [mg] mg of C in removed air
  mutate(adj_mass_diff = as.numeric(ifelse(time_diff != '0', as.numeric(mass_diff + lag(rem_mass_C, default = first(rem_mass_C))), '0'))) %>% # [mg] find adjusted mass by including removal mass_C difference in flush groups
  filter(!time_diff == '0') %>% # delete used values
  mutate(flux = adj_mass_diff/time_diff) #%>% # this depends on prev. line being right
  #filter_if(~is.numeric(.), all_vars(!is.infinite(.))) # keeps the "last" day of a flux measurement ie. gets rid of the "first" day of each session, that's what we graph


# Respired
data_resp <- data_all %>% 
  group_by(Sample) %>%
  ungroup(Flush) %>%  # ungroup Flush but keep groups by Sample
  #select(Flush, Sample, Date.Time, flux) %>%  # clean it up
  arrange(Date.Time) %>% # rearrange in ascending order
  mutate(time_hours = (Date.Time - lag(Date.Time, k = 1))) %>%     # time difference btwn flux measurements in days
  #mutate(time_hours = (Date.Time - lag(Date.Time, k = 1))*24) %>% # time difference btwn flux measurements in hours
  mutate(C_resp = .5*(time_hours)*(flux+lag(flux))) %>% # [mg] trapezoidal area calculation to get C respired
  drop_na(C_resp) %>% # drops rows w/ NAs which arise from the first trapezoid area measurement
  mutate(C_resp_cum = cumsum(as.numeric(C_resp))) %>% # [mg] cumulatively add together trapezoids
  #mutate(invC_resp_cum = total - C_resp_cum)this doesn't work, but it could be used to generate the inv figure lee thinks abt
  mutate(time = as.numeric(Date.Time - first(Date.Time), units = 'days')) # calculate time difference from first in group in [days]

#write.csv(data_resp, file="respdata1.csv", row.names = FALSE)

stats_resp <- data_resp %>% # output averages plotted in RESP graphs
  group_by(Sample, Num) %>% 
  summarise(max_C_resp_cum = max(C_resp_cum)) %>%
  #group_by(Num) %>%   # comment this in/out if you want it broken up to replicates or not
  summarise(mean_C_resp_cum = mean(max_C_resp_cum)) #%>% # [mg] 
  #mutate(Name = num_labs)    # comment this in/out if you want it broken up to replicates or not

stats_resp2 <- data_resp %>% # output averages plotted in RESP graphs
  group_by(Sample, Num) %>% 
  summarise(max_C_resp_cum = max(C_resp_cum)) %>%
  group_by(Num) %>%   # comment this in/out if you want it broken up to replicates or not
  summarise(mean_C_resp_cum = mean(max_C_resp_cum), stdev = sd(max_C_resp_cum)) # [mg C] 
  #mutate(Num = num_labs)    # comment this in/out if you want it broken up to replicates or not

print(paste("The incubation period currently spans", last_day, "days!"))
#write.csv(stats_resp2, file = 'INC2summary_cumCresp.csv', row.names = FALSE) # CHECK THAT THIS IS CORRECT NAME


####################################################################
# STATISTICS #######################################################

# 2 WAY ANOVA for INC1, test if treatment and soil type have an effect on mean C resp/fraction of C retained by soil/fraction of C retained by residue by end of incubation

# recode Num to factors
thirteenC_data <- read.csv("INC2_2wayanova_13C.csv", stringsAsFactors = FALSE, header = TRUE) # scan in document formatted like example
soil_C_data <- read.csv('summary_SOILcumCresp.csv', stringsAsFactors = FALSE, header = TRUE)
stats_resp <- rbind(stats_resp, soil_C_data)

twowayanova_data <- merge(stats_resp, thirteenC_data, by = 'Sample') # if this excludes the soil data, check to make sure the 'Sample' column for both sheets is labelled correctly

twowayanova_data$Num <- factor(twowayanova_data$Num,
                               levels = names(num_labs),
                               labels = num_labs)

## Open vs. Closed, paired t-test
oc_data <- twowayanova_data %>%
  filter(Num == 'DASE_O' | Num == 'DASE_C')

t.test(mean_C_resp_cum ~ Valve, data = oc_data, paired = TRUE)
t.test(fr ~ Valve, data = oc_data, paired = TRUE)
# We see that for both 13C and inc data, p>.05 or that there is an insignificant difference between Open v. Closed.

## Dosage, 2-way ANOVA
dose_data <- twowayanova_data %>%
  filter(Sub == 'AD' | Sub == 'NREL' | Sub == 'POET') %>% # balanced design since same number of observations per treatment
  mutate(Sub = factor(Sub, levels = c('AD', 'NREL', 'POET'), labels = c("AD2","HLFB2", "HLFB3"))) %>%
  mutate(Dose = factor(Dose, levels = c('S', 'N'), labels = c('Standard', 'Reduced'))) 

dose_summary <- dose_data %>%
  group_by(Num, Sub, Dose) %>%
  summarize(act_mean_C_resp_cum = mean(mean_C_resp_cum)) 

dose_plot1 <- ggboxplot(dose_data, x = 'Sub', y = 'mean_C_resp_cum', color = 'Dose', # boxplot shows various treatments and how they compare to each other + soil type
          xlab = 'Treatment',
          ylab = 'Mean C Respired [mg]')  
dose_plot1
ggsave("dose_plot1.png", plot = dose_plot1, width = 15, height = 15, units = "cm")

  
dose_plot2 <- ggline(dose_data, x = "Sub", y = "mean_C_resp_cum", color = "Dose",
       add = c("mean_se", "dotplot"),
       palette = c("#00AFBB", "#E7B800"),
       xlab = 'Treatment',
       ylab = 'Mean C Respired [mg]')
dose_plot2
ggsave("dose_plot2.png", plot = dose_plot2, width = 15, height = 15, units = "cm")


ggline(dose_data, x = "Sub", y = "fr", color = "Dose",
       add = c("mean_se", "dotplot"),
       palette = c("#00AFBB", "#E7B800"))

res.aov_dose <- aov(mean_C_resp_cum ~ Dose * Sub, data = dose_data) # test interaction btwn Num and Typ
summary(res.aov_dose)

res.aov_dose2 <- aov(fr ~ Dose * Sub, data = dose_data) # test interaction btwn Num and Typ
summary(res.aov_dose2)

# Tukey-Kramer maybe
TukeyHSD(res.aov_dose)  # unclear if this is taking into account unbalanced design, I think it's using Tukey Kramer

## INC1 V INC2, t-test
INC2_DASE_C_data <- twowayanova_data %>%
  #filter(Sub == 'DASE' & Valve == 'C') %>%
  filter(Sub == 'DASE' & Valve == 'C') %>%
  select(Num, mean_C_resp_cum, fr) %>%
  mutate(Inc = '2')

INC1_anovadata <- read.csv("INC1_twowayanova_data.csv", stringsAsFactors = FALSE, header = TRUE) # scan in document formatted like example
INC1_DASE_data <- INC1_anovadata %>%
  filter(Typ == 'P' & Num == 'DASE HLFB') %>%
  select(Num, mean_C_resp_cum, fr) %>%
  mutate(Inc = '1')

DASE_data <- rbind(INC2_DASE_C_data, INC1_DASE_data)

t.test(mean_C_resp_cum ~ Inc, data = DASE_data, paired = FALSE)
t.test(fr ~ Inc, data = DASE_data, paired = FALSE)

# Check SOIL controls
INC2_soil_data <- twowayanova_data %>%
  filter(Sub == 'SOIL') %>%
  select(Num, mean_C_resp_cum) %>%
  mutate(Inc = '2')

INC1_soil_data <- INC1_anovadata %>%
  filter(Typ == 'P' & Num == 'Soil Control') %>%
  select(Num, mean_C_resp_cum) %>%
  mutate(Inc = '1')

SOIL_data <- rbind(INC1_soil_data, INC2_soil_data)

t.test(mean_C_resp_cum ~ Inc, data = SOIL_data, paired = FALSE)

############################################################################
# PLOTTING #######################################################

# Theme and Labels
theme_C <- theme_light() + 
  theme(panel.grid.minor = element_blank(), 
  text = element_text(size = 30), #for facetwrapped plots
  strip.background = element_rect(color="black", fill="#93C5FF", size=1.5, linetype="solid"),
  legend.position = "none",
  plot.title = element_text(hjust = 0.5)
  ) 

# CHANGE THESE DATES FOR YOUR GRAPHING PLEASURE! 
end_date = '2023-03-09 14:00' # !!!! CHANGE THIS TO EXTEND GRAPH !!!
  #'2022-12-13 9:15' <- this is for 42 days
# start_date = #"2022-11-07 7:30" # <- this is for ignoring the initial 6 day bump
start_date = "2022-11-01 19:10" # ACTUAL FIRST MEASUREMENT
lims <- as.POSIXct(strptime(c(start_date, end_date), format = "%Y-%m-%d %H:%M"))  

# FLUX: Mean and SE Each
p1P<- ggplot(data_all, aes(x=Date.Time, y=flux)) +
  geom_point(aes(size = .8)) +
  scale_x_datetime(limits = lims) +
  stat_summary(fun.data = "mean_se", colour = "red", size = .8) +
  facet_wrap(~Num, labeller = labeller(Num = num_labs)) + 
  #facet_wrap(~Num, scales = 'free', labeller = labeller(Num = num_labs)) + # free scale bc 1 is so small 
  theme_C +
  #scale_y_continuous(limits=c(0,.35)) +  # sets all plots start at 0 go to .3
  labs(x = '', y = 'Carbon Flux [mg/hr]', title = 'Carbon Flux Evolution in Various Treatments') 
p1P

#ggsave("flux_mean&se.png", plot = p1P, width = 60, height = 20, units = "cm")  # change this accordingly

# RESPIRED: Mean and SE Each
p2P <- ggplot(data_resp, aes(x=Date.Time, y=C_resp_cum)) +
  geom_point(aes(size = .8)) +
  scale_x_datetime(limits = lims) +
  stat_summary(fun.data = "mean_se", colour = "red", size = .8) +
  facet_wrap(~Num, labeller = labeller(Num = num_labs)) + # NON FREE SCALE
  ##facet_wrap(~Num, scales = 'free', labeller = labeller(Num = num_labs)) + # free scale bc 1 is so small 
  ##geom_vline(xintercept = as.POSIXct(as.Date(c('2021-03-22', '2021-04-22'))), linetype = 'dashed', color = 'blue', size = 2) +  # when water was added, comment this out for no lines 
  theme_C +
  #scale_y_continuous(limits=c(0,500)) +  # sets all plots start at 0 go to unique maxes for each 
  labs(x = '', y = 'Cumulative Carbon Respired [mg]', title = 'Cumulative Carbon Respired in Various Treatments') 
p2P

#ggsave("resp_mean&se.png", plot = p2P, width = 60, height = 20, units = "cm")


# GRAPHS OF RESPIRED ONLY OF THE NEW RATIOS 
edit_data_resp <- data_resp %>%
  filter(Num == '6' | Num == '7'| Num == '8'| Num == '9') 

edit_p2P <- ggplot(edit_data_resp, aes(x=Date.Time, y=C_resp_cum)) +
  geom_point(aes(size = .8)) +
  scale_x_datetime(limits = lims) +
  stat_summary(fun.data = "mean_se", colour = "red", size = .8) +
  facet_wrap(~Num, labeller = labeller(Num = num_labs)) +
  theme_C +
  labs(x = '', y = 'Cumulative Carbon Respired [mg]', title = 'Cumulative Carbon Respired in 50% Residue Dosage Treatments') 
edit_p2P
ggsave("newratio_resp_mean&se.png", plot = edit_p2P, width = 60, height = 20, units = "cm")

# GRAPHS OF RESPIRED ONLY OF OLD RATIOS
edit_data_resp <- data_resp %>%
  filter(Num == '3' | Num == '4'| Num == '5')

edit_p2P <- ggplot(edit_data_resp, aes(x=Date.Time, y=C_resp_cum)) +
  geom_point(aes(size = .8)) +
  scale_x_datetime(limits = lims) +
  stat_summary(fun.data = "mean_se", colour = "red", size = .8) +
  facet_wrap(~Num, labeller = labeller(Num = num_labs)) +
  theme_C +
  labs(x = '', y = 'Cumulative Carbon Respired [mg]', title = 'Cumulative Carbon Respired in Normal Residue Dosage Treatments') 
edit_p2P
ggsave("oldratio_resp_mean&se.png", plot = edit_p2P, width = 60, height = 20, units = "cm")

# GRAPHS OF DASE O/C
OCdata_resp <- data_resp %>%
  filter(Num == '1' | Num == '2')

OC_p <- ggplot(OCdata_resp, aes(x=Date.Time, y=C_resp_cum)) +
  geom_point(aes(size = .8)) +
  scale_x_datetime(limits = lims) +
  stat_summary(fun.data = "mean_se", colour = "red", size = .8) +
  facet_wrap(~Num, labeller = labeller(Num = num_labs)) +
  theme_C +
  labs(x = '', y = 'Cumulative Carbon Respired [mg]', title = 'Cumulative Carbon Respired in DASE O/C Treatments') 
OC_p
ggsave("DASEOC_resp_mean&se.png", plot = OC_p, width = 60, height = 20, units = "cm")

# OC_p2 <- ggplot(OCdata_resp, aes(x=Date.Time, y=C_resp_cum)) +
#   geom_point(aes(size = .8, col = Num)) +
#   scale_x_datetime(limits = lims) +
#   stat_summary(fun.data = "mean_se", colour = "red", size = .8) +
#   #facet_wrap(~Num, labeller = labeller(Num = num_labs)) +
#   theme_C +
#   labs(x = '', y = 'Cumulative Carbon Respired [mg]', title = 'Cumulative Carbon Respired in DASE O/C Treatments') +
#   legend
# OC_p2
# ggsave("sameDASEOC_resp_mean&se.png", plot = OC_p2, width = 60, height = 20, units = "cm")

#lumped figure w/ geom smooth of C respired
theme_lump <- theme_light() + 
  theme(panel.grid.minor = element_blank(), 
        text = element_text(size = 30), #for facetwrapped plots
        strip.background = element_rect(color="black", fill="#93C5FF", size=1.5, linetype="solid"),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)
  ) 

data_resp_old <- data_resp %>%
  filter(Num == '2' |Num == '3' | Num == '4'| Num == '5')

lumped1 <- ggplot(data_resp_old, aes(x=Date.Time, y=C_resp_cum)) +
  geom_smooth(aes(color = Num), se = TRUE) +
  scale_x_datetime(limits = lims) +
  theme_lump +
  scale_color_manual("Treatments", labels = c("DASE1", "AD2", "DASE2", "DASE3"), values = c("2", "3", "4", "5")) +
  scale_y_continuous(limits=c(0,250)) +  # sets all plots start at 0 go to unique maxes for each 
  labs(x = '', y = 'Cumulative Carbon Respired [mg]', title = 'Cumulative Carbon Respired in INCUBATION 2 of HLFB Amended Palouse Soil') 
lumped1

ggsave("Plumped_scale_mean&se.png", plot = lumped_P, width = 60, height = 20, units = "cm")

#C retained throughout Incubation 2
datainitC <- read.csv("justinitC.csv", stringsAsFactors = FALSE, header = TRUE) # scan in document formatted like example
data3 <- left_join(data_resp, datainitC, by = 'Sample') 
data3 <- data3 %>%
  mutate(invC_resp_cum = init_C*1000 - C_resp_cum)  %>%
  filter(Num == '2' |Num == '3' | Num == '4'| Num == '5') %>%
  mutate(invC_resp_cum_ADJ = case_when(Num == '3' ~ .5*invC_resp_cum,
                                       Num == '2' | Num == '4' ~ .35*invC_resp_cum,
                                       Num == '5' ~ .35*invC_resp_cum)) %>%
  mutate(ID = case_when(Num == '2' ~ 'DASE1_2', 
                        Num == '3' ~ 'AD2',
                        Num == '4' ~ 'DASE2',
                        Num == '5' ~ 'DASE3'))

lumped_2 <- ggplot(data3, aes(x=Date.Time, y=invC_resp_cum)) +
  geom_smooth(aes(color = Num), se = TRUE) +
  scale_x_datetime(limits = lims) +
  theme_lump +
  scale_color_manual("Treatments", labels = c("DASE1", "AD2", "DASE2", "DASE3"), values = c("2", "3", "4", "5")) +
  #scale_y_continuous(limits=c(0,250)) +  # sets all plots start at 0 go to unique maxes for each 
  labs(x = '', y = 'Cumulative Carbon Retained [mg]', title = 'Carbon Retained in Treatments Throughout INCUBATION 2 of HLFB Amended Palouse Soil') 
lumped_2

ggsave("Cret4res.png", plot = lumped_2, width = 60, height = 20, units = "cm")
 

lumped_3 <- ggplot(data3, aes(x=Date.Time, y=invC_resp_cum_ADJ)) +
  geom_smooth(aes(color = Num), se = TRUE) +
  scale_x_datetime(limits = lims) +
  theme_lump +
  scale_color_manual("Treatments", labels = c("DASE1", "AD2", "DASE2", "DASE3"), values = c("2", "3", "4", "5")) +
  #scale_y_continuous(limits=c(0,250)) +  # sets all plots start at 0 go to unique maxes for each 
  labs(x = '', y = 'Cumulative Carbon Retained [mg]', title = 'Carbon Retained in Treatments Throughout INCUBATION 2 of HLFB Amended Palouse Soil') 
lumped_3

ggsave("Cret4resADJ.png", plot = lumped_3, width = 60, height = 20, units = "cm")

# combine
inc2data <- read.csv("INC2_invC_resp.csv", stringsAsFactors = FALSE, header = TRUE) # scan in document formatted like example

retdata1 <- data3 %>%
  select(Sample, ID, time, invC_resp_cum, invC_resp_cum_ADJ) 

retdata2 <- inc2data %>%
  filter(time < 135) %>%
  select(Sample, ID, time, invC_resp_cum, invC_resp_cum_ADJ) 

retdata_comb <- rbind(data_frame(retdata1), data_frame(retdata2))

lumped_4 <- ggplot(retdata_comb, aes(x=time, y=invC_resp_cum)) +
  geom_smooth(aes(color = ID), se = TRUE) +
  #scale_x_datetime(limits = lims) +
  theme_lump +
  #scale_color_manual("Treatments", labels = c("DASE1", "AD2", "DASE2", "DASE3"), values = c("2", "3", "4", "5")) +
  #scale_y_continuous(limits=c(0,250)) +  # sets all plots start at 0 go to unique maxes for each 
  labs(x = 'Time [days]', y = 'Cumulative Carbon Retained [mg]', title = 'Carbon Retained in 135 Day Incubations of HLFB Amended Palouse Soil') 
lumped_4

ggsave("totalCret4res.png", plot = lumped_4, width = 60, height = 20, units = "cm")

lumped_5 <- ggplot(retdata_comb, aes(x=time, y=invC_resp_cum_ADJ)) +
  geom_smooth(aes(color = ID), se = TRUE) +
  #scale_x_datetime(limits = lims) +
  theme_lump +
  #scale_color_manual("Treatments", labels = c("DASE1", "AD2", "DASE2", "DASE3"), values = c("2", "3", "4", "5")) +
  #scale_y_continuous(limits=c(0,250)) +  # sets all plots start at 0 go to unique maxes for each 
  labs(x = 'Time [days]', y = 'Cumulative Carbon Retained [mg]', title = 'Carbon Retained in 135 Day Incubations of HLFB Amended Palouse Soil') 
lumped_5

ggsave("ADJtotalCret4res.png", plot = lumped_5, width = 60, height = 20, units = "cm")

#Show Initial C
pV_init <- ggplot(data_resp, aes(x=Date.Time, y=C_resp_cum)) +
  geom_smooth(aes(size = .8)) +
  scale_x_datetime(limits = lims) +
  #stat_summary(fun.data = "mean_se", colour = "red", size = .8) +
  facet_wrap(~Num, labeller = labeller(Num = num_labs)) + # NON FREE SCALE
  #geom_vline(xintercept = as.POSIXct(as.Date(c('2021-03-22', '2021-04-22'))), linetype = 'dashed', color = 'blue', size = 2) +
  #facet_wrap(~Num, scales = 'free', labeller = labeller(Num = num_labs)) + # free scale bc 1 is so small 
  ##geom_vline(xintercept = as.POSIXct(as.Date(c('2021-03-22', '2021-04-22'))), linetype = 'dashed', color = 'blue', size = 2) +  # when water was added, comment this out for no lines 
  theme_C +
  scale_y_continuous(limits=c(0,NA)) +  # sets all plots start at 0 go to unique maxes for each 
  labs(x = '', y = 'Cumulative Carbon Respired [mg]', title = 'Cumulative Carbon in 267 Day Incubation of HLFB Amended Vershire Soil') 
pV_init

# Calculate C retained as percentage of residue C and total treatment C 
datainitC <- read.csv("justinitC.csv", stringsAsFactors = FALSE, header = TRUE) # scan in document formatted like example
data3 <- left_join(data_resp, datainitC, by = 'Sample') 
data3 <- merge(data3, soil_resp, by = c('Typ', 'Flush'))  # matches soil resp. to each measurement at a time point 
data3 <- data3 %>%
  mutate(init_totC = init_C*1000) %>%         # [mg C] initial C (soil+res) in each treatment on average
  mutate(init_resC = init_resC*1000) %>%      # [mg C] initial C (res) in each treatment on average
  mutate(Ctot_ret = 100*(init_totC - C_resp_cum)/init_totC) %>% # [% total C] C retained from total treatment
  mutate(Cres_ret = 100*(init_resC-(C_resp_cum-mean_soil_cum))/init_resC) %>% # [% residue C] C retained from residue in each treatment
  group_by(Num, Typ, Date.Time) %>%
  summarize(meanCtot_ret = mean(Ctot_ret), meanCres_ret = mean(Cres_ret)) %>% # [% total C] average of prev. calculations per treatment on specific days
  ungroup() %>% # necessary to add row after
  add_row(Typ = 'P', Num = c('1', '2', '3', '4', '5'), Date.Time = as.POSIXct('2021-11-01 15:00:00'), meanCres_ret = 100, meanCtot_ret = 100) %>% # add initial anchor point of 100% for all treatments (when incubation began)
  add_row(Typ = 'V', Num = c('1', '2', '3', '4', '5'), Date.Time = as.POSIXct('2021-11-01 15:00:00'), meanCres_ret = 100, meanCtot_ret = 100) 

data3P <- data3 %>%
  filter(Typ == 'P') %>%
  mutate(Ctot_Label = round(ifelse(Date.Time == max(Date.Time), meanCtot_ret, NA), 0)) %>%  # add labels to last point of each line
  mutate(Cres_Label = round(ifelse(Date.Time == max(Date.Time), meanCres_ret, NA), 0))

data3V <- data3 %>%
  filter(Typ == 'V') %>%
  mutate(Ctot_Label = round(ifelse(Date.Time == max(Date.Time), meanCtot_ret, NA), 0)) %>%  # add labels to last point of each line
  mutate(Cres_Label = round(ifelse(Date.Time == max(Date.Time), meanCres_ret, NA), 0))

# Graph C retained graphs 
# C retained of only residue graphs
Cret_P <- ggplot(data3P, aes(x=Date.Time, y=meanCres_ret)) +
  geom_line(aes(color = Num), size = .5) +
  geom_point(size = .25, color = 'black') + 
  scale_x_datetime(date_breaks = '1 month', labels = date_format("%b")) +
  ylim(35, 100) +
  theme_lump +
  scale_color_manual("Treatments", labels = c("Soil Control", "CS", "AD", "C-CBP", "DASE"), values = c("1", "2", "3", "4", "5")) +
  #scale_y_continuous(limits=c(0,500)) +  # sets all plots start at 0 go to 500
  labs(x = '', y = 'Carbon Retained in Residue \n [% of Initial Residue C]', title = 'Palouse Soil Incubations') +   # Carbon Retained in Residue in \n 267 Day Incubation of Biofuel Residues in Palouse Soil
  geom_label_repel(aes(label = Cres_Label), min.segment.length = 0, size = 2, force = 2.1, direction = 'y', hjust = 'left', label.padding = unit(0.1, "lines"), na.rm = TRUE)  # labels last point with final percentage of each line 
Cret_P

# C retained of total treatment graphs
Cret_P2 <- ggplot(data3P, aes(x=Date.Time, y=meanCtot_ret)) +
  geom_line(aes(color = Num), size = .5) +
  geom_point(size = .25, color = 'black') + 
  scale_x_datetime(date_breaks = '1 month', labels = date_format("%b")) +
  ylim(60, 100) +
  #scale_x_datetime(limits = lims) +
  theme_lump +
  scale_color_manual("Treatments", labels = c("Soil Control", "CS", "AD", "C-CBP", "DASE"), values = c("1", "2", "3", "4", "5")) +
  #scale_y_continuous(limits=c(0,500)) +  # sets all plots start at 0 go to 500
  labs(x = '', y = 'Carbon Retained in Treatment \n [% of Initial Treatment C]', title = 'Palouse Soil Incubations') + # 'Carbon Retained in Treatment in \n 267 Day Incubation of Biofuel Residues in Palouse Soil' 
  geom_label_repel(aes(label = Ctot_Label), min.segment.length = 0 , size = 2, force = .6, direction = 'y', hjust = 'left', label.padding = unit(0.1, "lines"), na.rm = TRUE)  # labels last point with final percentage of each line 
Cret_P2

# C retained of only residue graphs
Cret_V <- ggplot(data3V, aes(x=Date.Time, y=meanCres_ret)) +
  geom_line(aes(color = Num), size = .5) +
  geom_point(size = .25, color = 'black') + 
  scale_x_datetime(date_breaks = '1 month', labels = date_format("%b")) +
  ylim(35, 100) +
  #scale_x_datetime(limits = lims) +
  theme_lump +
  #geom_label_repel(aes(), nudge_x = 1, na.rm = TRUE) +  # CHANGE THIS SO THE LABEL WORKS, PAGE IS SAVED IN GOOGLE
  scale_color_manual("Treatments", labels = c("Soil Control", "CS", "AD", "C-CBP", "DASE"), values = c("1", "2", "3", "4", "5")) +
  #scale_y_continuous(limits=c(0,500)) +  # sets all plots start at 0 go to 500
  labs(x = '', y = 'Carbon Retained in Residue \n [% of Initial Residue C]', title = 'Vershire Soil Incubations') +   # Carbon Retained in Residue in \n 267 Day Incubation of Biofuel Residues in Vershire Soil
  geom_label_repel(aes(label = Cres_Label), min.segment.length = 0, size = 2, force = .5, direction = 'y', hjust = 'left', label.padding = unit(0.1, "lines"), na.rm = TRUE)  # labels last point with final percentage of each line 
Cret_V

# C retained of total treatment graphs
Cret_V2 <- ggplot(data3V, aes(x=Date.Time, y=meanCtot_ret)) +
  geom_line(aes(color = Num), size = .5) +
  geom_point(size = .25, color = 'black') +  
  scale_x_datetime(date_breaks = '1 month', labels = date_format("%b")) +
  ylim(60, 100) +
  #scale_x_datetime(limits = lims) +
  theme_lump +
  scale_color_manual("Treatments", labels = c("Soil Control", "CS", "AD", "C-CBP", "DASE"), values = c("1", "2", "3", "4", "5")) +
  #scale_y_continuous(limits=c(0,500)) +  # sets all plots start at 0 go to 500
  labs(x = '', y = 'Carbon Retained in Treatment \n [% of Initial Treatment C]', title = 'Vershire Soil Incubations') + # 'Carbon Retained in Treatment in \n 267 Day Incubation of Biofuel Residues in Vershire Soil'
  geom_label_repel(aes(label = Ctot_Label), min.segment.length = 0, size = 2, force = .6, direction = 'y', hjust = 'left', label.padding = unit(0.1, "lines"), na.rm = TRUE)  # labels last point with final percentage of each line 
Cret_V2

ggsave("CretP_scale_mean&se.png", plot = Cret_P, width = 60, height = 20, units = "cm")
ggsave("CretV_scale_mean&se.png", plot = Cret_V, width = 60, height = 20, units = "cm")



#NOT READY
#################################
# FLUX: Corn Stover Lumped
data_corn_flux <- data_all %>%
  filter(Num != 1)
p_corn_flux <- ggplot(data_corn_flux, aes(x = Date.Time, y = flux, group = Num)) +
  stat_summary(fun.data = "mean_se", aes(group = Num,color = Num)) +
  theme_C +
  labs(x = '', y = 'Carbon Flux [mg/hr]', title = 'Carbon Flux Evolution in Corn Stover \n Treatments Over an 89 Day Incubation') 
p_corn_flux
ggsave("corn_flux.png", plot = p_corn_flux, width = 15, height = 15, units = "cm")

# RESPIRED: Corn Stover Lumped
data_corn_resp <- data_resp %>%
  filter(Num != 1)
p_corn_resp <- ggplot(data_corn_resp, aes(x = Date.Time, y = C_resp_cum, group = Num)) +
  stat_summary(fun.data = "mean_se", aes(group = Num,color = Num)) +
  theme_C +
  labs(x = '', y = 'Cumulative Carbon Respired [mg]', title = 'Cumulative Carbon Respired in Corn Stover \n Treatments Over an 89 Day Incubation') 
p_corn_resp
ggsave("corn_resp.png", plot = p_corn_resp, width = 15, height = 15, units = "cm")

############################################################################
# SOIL MODELLING #######################################################
# Based off of https://www.bgc-jena.mpg.de/TEE/optimization/2015/12/09/Fractions-Incubations/
# Context from https://escholarship.org/uc/item/9h72f7hk

# Clean data for modelling
data_mod <- data_resp %>%
  ungroup(Sample) %>%   # now, not grouped as anything
  select(c('time','Num', 'C_resp_cum'))  %>% # select these columns for ease
  group_by(Num, time) %>%  #  
  summarize(cummCO2 = mean(C_resp_cum)) # sd gives an error for some reason: Stderr = sd(C_resp_cum))   # [mg] amount of carbon respired cumulatively, not in terms of mg C/g soil
  #summarize(cummCO2 = mean(C_resp_cum)/50, Stderr = sd(C_resp_cum/50)) %>%  # /50 so it's in [g C/g soil] since we start w/ ~50g soil, summarizing by all incubations def. loses precision since it's not a rate, it's an absolute amount?, but also it's based off of rate anyways
write.csv(data_mod, file = 'INC2data_mod.csv')

############################################################################
# EXTRA PLOTTING #######################################################

# FLUX: All Plots
plot_all_flux <- ggplot(data = data_all) +
  geom_point(aes(x = Date.Time, y = flux)) +
  labs(x = 'Time', y = 'Flux [mg/hr]', title = 'Soil Incubation Experiments: C Flux v. Time') +
  facet_wrap(~Sample, nrow = 3) +
  theme_light()
plot_all_flux

# Plot 1 FLUX separately
data_1 <- data_all %>%
  filter(Num == '1')
plot_1_flux <- ggplot(data = data_1) +
  geom_point(aes(x = Date.Time, y = flux)) +
  labs(x = 'Time', y = 'Carbon Flux [mg/hr]', title = 'Incubation of (1) Soil Controls: Carbon Flux v. Time') +
  facet_wrap(~Sample, nrow = 1) +
  theme_light()
plot_1_flux

# Plot 2 and 3 FLUX separately
data_23 <- data_all %>%
  filter(Num != '1')
plot_23_flux <- ggplot(data = data_23) +
  geom_point(aes(x = Date.Time, y = flux)) +
  labs(x = 'Time', y = 'Carbon Flux [mg/hr]', title = 'Incubation of (2) 500 um and (3) 8500 um Corn Stover: Carbon Flux v. Time') +
  facet_wrap(~Sample, nrow = 2) +
  theme_light()
plot_23_flux

# Plot LUMPED FLUX v Time [1,2,3 w/ geom_smooth in one image, same scale]
new_labels <- c('1' = 'No Residue Control','2' = '500 um Corn Stover', '3' = '8500 um Corn Stover')
plot_lumped_flux <- ggplot(data_all, aes(x = Date.Time, y = flux)) +
  geom_point() +
  geom_smooth(span = 0.8) +
  labs(x = 'Time', y = 'Carbon Flux [mg/hr]', title = 'Carbon Flux v. Time') +
  facet_wrap(~Num, nrow = 3, labeller = labeller(Num = new_labels)) +
  theme_light()
plot_lumped_flux

# Plot 1 LUMPED FLUX separately
data_lumped1 <- data_all %>%
  filter(Num == '1')
plot_lumped1_flux <- ggplot(data_lumped1, aes(x = Date.Time, y = flux)) +
  geom_point() +
  geom_smooth(span = 0.8) +
  labs(x = 'Time', y = 'Carbon Flux [mg/hr]', title = 'Incubation of (1) Soil Controls: Cumulative C Respired v. Time') +
  facet_wrap(~Num, nrow = 1, labeller = labeller(Num = new_labels)) +
  theme_light()
plot_lumped1_flux

# Plot 2 and 3 LUMPED FLUX separately
data_lumped23 <- data_all %>%
  filter(Num != '1')
plot_lumped23_flux <- ggplot(data_lumped23, aes(x = Date.Time, y = flux)) +
  geom_point() +
  geom_smooth(span = 0.8) +
  labs(x = 'Time', y = 'Carbon Flux [mg/hr]', title = 'Incubation of (2) 500 um and (3) 8500 um Corn Stover: Carbon Flux v. Time') +
  facet_wrap(~Num, nrow = 3, labeller = labeller(Num = new_labels)) +
  theme_light()
plot_lumped23_flux

# Plot 1 RESP separately
data_1_resp <- data_resp %>%
  filter(Num == '1')
plot_1_resp <- ggplot(data = data_1_resp) +
  geom_point(aes(x = Date.Time, y = C_resp_cum)) +
  labs(x = 'Time', y = 'Cumulative C Respired [mg]', title = 'Incubation of (1) Soil Controls: Cumulative C Respired v. Time') +
  facet_wrap(~Sample, nrow = 1) +
  theme_light()
plot_1_resp

# Plot 2 and 3 RESP separately
data_23_resp <- data_resp %>%
  filter(Num != '1')
plot_23_resp <- ggplot(data = data_23_resp) +
  geom_point(aes(x = Date.Time, y = C_resp_cum)) +
  labs(x = 'Time', y = 'Cumulative C Respired [mg]', title = 'Incubation of (2) 500 um and (3) 8500 um Corn Stover: Cumulative C Respired v. Time') +
  facet_wrap(~Sample, nrow = 2) +
  theme_light()
plot_23_resp

#Adjusted for controls
ggplot(data_resp, aes(x = Date.Time, y = C_resp_cum, group = Num, col = Num)) +
  geom_point()

# Plot LUMPED RESP. v Time [1,2,3 w/ geom_smooth in one image, same scale]
plot_lumped_resp <- ggplot(data_resp, aes(x = Date.Time, y = C_resp_cum)) +
  geom_point() +
  geom_smooth(span = 0.8) +
  labs(x = 'Time', y = 'Cumulative C Respired [mg]', title = 'Cumulative C Respired v. Time') +
  facet_wrap(~Num, nrow = 3, labeller = labeller(Num = new_labels)) +
  theme_light()
plot_lumped_resp

# Plot LUMPED 1 separately
data_lumped1_resp <- data_resp %>%
  filter(Num == '1')
plot_lumped1_resp <- ggplot(data_lumped1_resp, aes(x = Date.Time, y = C_resp_cum)) +
  geom_point() +
  geom_smooth(span = 0.8) +
  labs(x = 'Time', y = 'Carbon Flux [mg/hr]', title = 'Incubation of (1) Soil Controls: Cumulative Carbon Respired v. Time') +
  facet_wrap(~Num, nrow = 1, labeller = labeller(Num = new_labels)) +
  theme_light()
plot_lumped1_resp

# Plot LUMPED 2 and 3 separately
data_lumped23_resp <- data_resp %>%
  filter(Num != '1')
plot_lumped23_resp <- ggplot(data_lumped23_resp, aes(x = Date.Time, y = C_resp_cum)) +
  geom_point() +
  geom_smooth(span = 0.8) +
  labs(x = 'Time', y = 'Carbon Flux [mg/hr]', title = 'Incubation of (2) 500 um and (3) 8500 um Corn Stover: Cumulative Carbon Respired v. Time') +
  facet_wrap(~Num, nrow = 1, labeller = labeller(Num = new_labels)) +
  theme_light()
plot_lumped23_resp
