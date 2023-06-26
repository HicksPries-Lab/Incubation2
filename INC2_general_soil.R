## SOIL INCUBATION CALCULATIONS, SEPT 2022
## Author: Michelle S. Wang, michelle.s.wang.th@dartmouth.edu

# Load packages + functions
library(tidyverse)
# library(SoilR)
library(FME)

# Read in data
data <- read.csv("IRGASoil_Measurements.csv", stringsAsFactors = FALSE, header = TRUE) # scan in document formatted like example
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
data_all <- data_C %>%
  #filter(!Sample %in% c('CO2 FREE', '2008', '2%')) %>% # filters out controls 
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
  mutate(flux = adj_mass_diff/time_diff) # %>% # [mg/hr] %>% # this depends on prev. line being right
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

write.csv(data_resp, file="respsoildata1.csv", row.names = FALSE)

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
  summarise(mean_C_resp_cum = mean(max_C_resp_cum), stdev = sd(max_C_resp_cum)) #%>% # [mg] 
#mutate(Name = num_labs)    # comment this in/out if you want it broken up to replicates or not


print(paste("The incubation period currently spans", last_day, "days!"))
#write.csv(stats_resp2, file = 'INC2summary_SOILcumCresp.csv', row.names = FALSE) # CHECK THAT THIS IS CORRECT NAME


####################################################################
# STATISTICS #######################################################

# 2 WAY ANOVA for INC1, test if treatment and soil type have an effect on mean C resp/fraction of C retained by soil/fraction of C retained by residue by end of incubation


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
num_labs <- c('.16GWC', '.2GWC')
names(num_labs) <- c('1', '2')

lims <- as.POSIXct(strptime(c("2022-11-14 9:00","2023-03-27 11:00"), format = "%Y-%m-%d %H:%M"))  # !!!! CHANGE THIS TO EXTEND GRAPH !!!

# FLUX: Mean and SE Each
p1P<- ggplot(data_all, aes(x=Date.Time, y=flux)) +
  geom_point(aes(size = .8)) +
  scale_x_datetime(limits = lims) +
  stat_summary(fun.data = "mean_se", colour = "red", size = .8) +
  facet_wrap(~Num, labeller = labeller(Num = num_labs)) + 
  #facet_wrap(~Num, scales = 'free', labeller = labeller(Num = num_labs)) + # free scale bc 1 is so small 
  theme_C +
  #scale_y_continuous(limits=c(0,.35)) +  # sets all plots start at 0 go to .3
  labs(x = '', y = 'Carbon Flux [mg/hr]', title = 'Carbon Flux Evolution in Soils with Differing GWC Treatments') 
p1P

ggsave("JustsoilFlux_mean&se.png", plot = p1P, width = 60, height = 20, units = "cm")  # change this accordingly

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
  labs(x = '', y = 'Cumulative Carbon Respired [mg]', title = 'Cumulative Carbon Respired in Soils with Differing GWC Treatments') 
p2P

ggsave("JustsoilResp_mean&se.png", plot = p2P, width = 60, height = 20, units = "cm")


# OLD OLD OLD doesn't work

#lumped figure w/ geom smooth of C respired
theme_lump <- theme_light() + 
  theme(panel.grid.minor = element_blank(), 
        text = element_text(size = 30), #for facetwrapped plots
        strip.background = element_rect(color="black", fill="#93C5FF", size=1.5, linetype="solid"),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)
  ) 

lumped_P <- ggplot(data_resp_P, aes(x=Date.Time, y=C_resp_cum)) +
  geom_smooth(aes(color = Num), se = TRUE) +
  scale_x_datetime(limits = lims) +
  theme_lump +
  scale_color_manual("Treatments", labels = c("Soil Control", "CS", "AD", "C-CBP", "DASE"), values = c("1", "2", "3", "4", "5")) +
  scale_y_continuous(limits=c(0,500)) +  # sets all plots start at 0 go to unique maxes for each 
  labs(x = '', y = 'Cumulative Carbon Respired [mg]', title = 'Cumulative Carbon Respired in 267 Day Incubation of HLFB Amended Palouse Soil') 
lumped_P

lumped_V <- ggplot(data_resp_V, aes(x=Date.Time, y=C_resp_cum)) +
  geom_smooth(aes(color = Num), se = TRUE) +
  scale_x_datetime(limits = lims) +
  theme_lump +
  scale_color_manual("Treatments", labels = c("Soil Control", "CS", "AD", "C-CBP", "DASE"), values = c("1", "2", "3", "4", "5")) +
  scale_y_continuous(limits=c(0,500)) +  # sets all plots start at 0 go to unique maxes for each 
  labs(x = '', y = 'Cumulative Carbon Respired [mg]', title = 'Cumulative Carbon Respired in 267 Day Incubation of HLFB Amended Vershire Soil') 
lumped_V

ggsave("Plumped_scale_mean&se.png", plot = lumped_P, width = 60, height = 20, units = "cm")
ggsave("Vlumped_scale_mean&se.png", plot = lumped_V, width = 60, height = 20, units = "cm")

#Lee's graph
datainitC <- read.csv("justinitC.csv", stringsAsFactors = FALSE, header = TRUE) # scan in document formatted like example
data3 <- left_join(data_resp, datainitC, by = 'Sample') 
data3 <- data3 %>%
  mutate(invC_resp_cum = init_C*1000 - C_resp_cum) 

data3P <- data3 %>%
  filter(Typ == 'P')
data3V <- data3 %>%
  filter(Typ == 'V')

lumped_P2 <- ggplot(data3P, aes(x=Date.Time, y=invC_resp_cum)) +
  geom_smooth(aes(color = Num), se = TRUE) +
  scale_x_datetime(limits = lims) +
  theme_lump +
  scale_color_manual("Treatments", labels = c("Soil Control", "CS", "AD", "C-CBP", "DASE"), values = c("1", "2", "3", "4", "5")) +
  #scale_y_continuous(limits=c(0,500)) +  # sets all plots start at 0 go to 500
  labs(x = '', y = 'Carbon Retained in Treatment [mg]', title = 'Carbon Retained in Treatment in 267 Day Incubation of HLFB Amended Palouse Soil') 
lumped_P2

lumped_V2 <- ggplot(data3V, aes(x=Date.Time, y=invC_resp_cum)) +
  geom_smooth(aes(color = Num), se = TRUE) +
  scale_x_datetime(limits = lims) +
  theme_lump +
  scale_color_manual("Treatments", labels = c("Soil Control", "CS", "AD", "C-CBP", "DASE"), values = c("1", "2", "3", "4", "5")) +
  #scale_y_continuous(limits=c(0,500)) +  # sets all plots start at 0 go to 500
  labs(x = '', y = 'Carbon Retained in Treatment [mg]', title = 'Carbon Retained in Treatment in 267 Day Incubation of HLFB Amended Vershire Soil') 
lumped_V2

ggsave("Pinitlumped_scale_mean&se.png", plot = lumped_P2, width = 60, height = 20, units = "cm")
ggsave("Vinitlumped_scale_mean&se.png", plot = lumped_V2, width = 60, height = 20, units = "cm")


#Show Initial C
pV_init <- ggplot(data_resp_V, aes(x=Date.Time, y=C_resp_cum)) +
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

############################################################################
# SOIL MODELLING #######################################################
# Based off of https://www.bgc-jena.mpg.de/TEE/optimization/2015/12/09/Fractions-Incubations/
# Context from https://escholarship.org/uc/item/9h72f7hk

# Clean data for modelling
data_mod <- data_resp %>%
  ungroup(Sample) %>%   # now, not grouped as anything
  select(c('time','Num', 'C_resp_cum'))  %>% # select these columns for ease
  group_by(Num, time) %>%  
  mutate(Num = ifelse(Num == '1','10','11')) %>%
  summarize(cummCO2 = mean(C_resp_cum)) # sd gives an error for some reason: Stderr = sd(C_resp_cum))   # [mg] amount of carbon respired cumulatively, not in terms of mg C/g soil
#summarize(cummCO2 = mean(C_resp_cum)/50, Stderr = sd(C_resp_cum/50)) %>%  # /50 so it's in [g C/g soil] since we start w/ ~50g soil, summarizing by all incubations def. loses precision since it's not a rate, it's an absolute amount?, but also it's based off of rate anyways
write.csv(data_mod, file = 'INC2_dataSOIL_mod.csv')


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
