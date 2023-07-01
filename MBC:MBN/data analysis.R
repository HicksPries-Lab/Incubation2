#Research Project Figures

#start: load every time
library(tidyverse)
library(dplyr)
library(colourpicker)
setwd("/Users/audreyadamchak/Downloads/independent research project/Project data/analysis/R/csv files")


Biomass = read.csv("MB5.9.csv", header=TRUE)
Available = read.csv("av5.9.csv", header=TRUE)
Lignin = read.csv("lignin.csv", header=TRUE)
Cresp = read.csv("Cresp.csv", header=TRUE)

Data <- Available %>% 
  full_join(.,Biomass, by=c("Sample", "Treatment")) %>%
  filter(Treatment != '20GWC') 


Data <- Available %>% 
  full_join(.,Biomass, by=c("Sample", "Treatment")) %>%
  filter(Treatment != '20GWC')%>% 
  mutate(Treatment = as.factor(Treatment),
        Treatment = factor(Treatment, 
                            levels = c("Control", "AD", "CS", "DASE", "NREL", "POET"), 
                            labels = c("Control", "AD2", "CS2", "HLFB1", "HLFB2", "HLFB3")))

Data$Treatment = relevel(Data$Treatment, "CS2")
Data$Treatment = relevel(Data$Treatment, "Control") 

CB <- Biomass %>% 
  full_join(.,Cresp, by=c("Sample", "Treatment"))%>%
  filter(Treatment != '20GWC')%>% 
  mutate(Treatment = as.factor(Treatment),
         Treatment = factor(Treatment, 
                            levels = c("Control", "AD", "CS", "DASE", "NREL", "POET"), 
                            labels = c("Control", "AD2", "CS2", "HLFB1", "HLFB2", "HLFB3")))


#BOXPLOTS

ggplot(Data) + 
  geom_boxplot(aes(x = Treatment, y = MBC, fill = Treatment), show.legend = FALSE) +
  ggtitle('Microbial Biomass Carbon') +
  scale_fill_manual(values=c("#269909F5", "#F0E442", "#E09B1B", "#87CEFA","#6495ED","#1671B3")) +
  theme_classic() +
  xlab("Treatment")+
  ylab("MBC content (mg C/g dry soil)") +
  theme_classic() +
  theme(axis.text.x=element_text(size = 12, colour = "black"),
        axis.text.y=element_text(size = 12, colour = "black"),
        axis.title.x=element_text(size = 14, colour = "black"),
        axis.title.y=element_text(size = 14, colour = "black"), 
        plot.title = element_text(color="black", size=18, hjust = 0.5)) +
  scale_y_continuous(expand = c(0, 0)) 

ggplot(Data) + 
  geom_boxplot(aes(x = Treatment, y = MBN, fill = Treatment), show.legend = FALSE)+
  ggtitle('Microbial Biomass Nitrogen') +
  xlab("Treatment")+
  ylab("MBN content (mg N/g dry soil)") +
  scale_fill_manual(values=c("#269909F5", "#F0E442", "#E09B1B", "#87CEFA","#6495ED","#1671B3")) +
  theme_classic() +
  theme(axis.text.x=element_text(size = 12, colour = "black"),
        axis.text.y=element_text(size = 12, colour = "black"),
        axis.title.x=element_text(size = 14, colour = "black"),
        axis.title.y=element_text(size = 14, colour = "black"), 
        plot.title = element_text(color="black", size=18, hjust = 0.5))+
scale_y_continuous(expand = c(0, 0), limits = c(0, 0.37)) 

ggplot(Data) + 
  geom_boxplot(aes(x = Treatment, y = TOC, fill = Treatment), show.legend = FALSE)+
  scale_fill_manual(values=c("#269909F5", "#F0E442", "#E09B1B", "#87CEFA","#6495ED","#1671B3")) +
  theme_classic() +
  ggtitle('Total Organic Carbon') +
  xlab("Treatment")+
  ylab("TOC (mg C/g dry soil)") +
  theme_classic() +
  theme(axis.text.x=element_text(size = 12, colour = "black"),
        axis.text.y=element_text(size = 12, colour = "black"),
        axis.title.x=element_text(size = 14, colour = "black"),
        axis.title.y=element_text(size = 14, colour = "black"), 
        plot.title = element_text(color="black", size=18, hjust = 0.5))+
  scale_y_continuous(expand = c(0, 0), limits = c(0.1, 0.31))

ggplot(Data, (aes(x = Treatment, y = TN, fill = Treatment))) + 
  geom_boxplot(show.legend = FALSE)+
  ggtitle('Total Nitrogen') +
  xlab("Treatment")+
  ylab("TN (mg N/g dry soil)") +
  scale_fill_manual(values=c("#269909F5", "#F0E442", "#E09B1B", "#87CEFA","#6495ED","#1671B3")) +
  theme_classic() +
  theme(axis.text.x=element_text(size = 12, colour = "black"),
        axis.text.y=element_text(size = 12, colour = "black"),
        axis.title.x=element_text(size = 14, colour = "black"),
        axis.title.y=element_text(size = 14, colour = "black"), 
        plot.title = element_text(color="black", size=18, hjust = 0.5))+
  scale_y_continuous(expand = c(0, 0)) 

#Residue Characteristics

Biomass <- Biomass %>% filter(Treatment != '.20GWC') 

Biomass <- Biomass %>% filter(Treatment != '.16GWC') 

Test <- Lignin %>% 
  full_join(.,Biomass, by=c("Treatment")) %>% 
  mutate(Treatment = as.factor(Treatment),
         Treatment = factor(Treatment, 
                            levels = c("Control", "AD", "CS", "DASE", "NREL", "POET"), 
                            labels = c("Control", "AD2", "CS2", "HLFB1", "HLFB2", "HLFB3"))) %>% 
  filter(Treatment != 'Control') 
"#59ACFF", "#95B2F0","#1671B3"
"#87CEFA","#6495ED","#1671B3"

ggplot(Test) + 
  geom_jitter(aes(x = X.Lignin, y = MBC, color = Treatment), size=3) +
  #geom_smooth(aes(x = X.Lignin, y = MBC)) +
  scale_color_manual(values=c("#E09B1B","#FFD700", "#87CEFA","#6495ED","#1671B3")) +
  theme_classic() +
  ggtitle("MBC vs Lignin") +
  xlab("Lignin")+
  ylab("MBC content (mg C/g dry soil)") +
  theme_classic() +
  theme(axis.text.x=element_text(size = 12, colour = "black"),
        axis.text.y=element_text(size = 12, colour = "black"),
        axis.title.x=element_text(size = 14, colour = "black"),
        axis.title.y=element_text(size = 14, colour = "black"), 
        plot.title = element_text(color="black", size=16, hjust = 0.5))+
  scale_y_continuous(expand = c(0, 0), limits = c(0.1, 0.35)) 


ggplot(Test) + 
  geom_jitter(aes(x = Lignin.N, y = MBC, color = Treatment), , size=3) +
  geom_smooth(aes(x = Lignin.N, y = MBC), method="lm", color="black") +
  scale_color_manual(values=c("#E09B1B","#FFD700", "#87CEFA","#6495ED","#1671B3")) +
  geom_text(data = annotation5, aes (x=x, y=y, label=label), 
            color= "black", size = 4) +
  theme_classic() +
  ggtitle("MBC vs Lignin:N") +
  xlab("Lignin:N")+
  ylab("MBC content (mg C/g dry soil)") +
  theme_classic() +
  theme(axis.text.x=element_text(size = 12, colour = "black"),
        axis.text.y=element_text(size = 12, colour = "black"),
        axis.title.x=element_text(size = 14, colour = "black"),
        axis.title.y=element_text(size = 14, colour = "black"), 
        plot.title = element_text(color="black", size=16, hjust = 0.5))+
  scale_y_continuous(expand = c(0, 0), limits = c(0.1, 0.35)) 

ggplot(Test) + 
  geom_point(aes(x = X.C, y = MBC, color = Treatment), size=3) +
  geom_smooth(aes(x = X.C, y = MBC), color="black", method="lm", formula=y~poly(x,2)) +
  scale_color_manual(values=c("#E09B1B","#FFD700", "#87CEFA","#6495ED","#1671B3")) +
  geom_text(data = annotation3, aes (x=x, y=y, label=label), 
            color= "black", size = 4) +
  theme_classic() +
  ggtitle("MBC vs % Carbon in Substrate") +
  xlab("% Substrate Carbon")+
  ylab("MBC content (mg C/g dry soil)") +
  theme_classic() +
  theme(axis.text.x=element_text(size = 12, colour = "black"),
        axis.text.y=element_text(size = 12, colour = "black"),
        axis.title.x=element_text(size = 14, colour = "black"),
        axis.title.y=element_text(size = 14, colour = "black"), 
        plot.title = element_text(color="black", size=16, hjust = 0.5))



ggplot(Test) + 
  geom_jitter(aes(x = X.N, y = MBC, color = Treatment), size=3) +
  geom_smooth(aes(x = X.N, y = MBC), color="black", method="lm", formula=y~poly(x,2)) +
  scale_color_manual(values=c("#E09B1B","#FFD700", "#87CEFA","#6495ED","#1671B3")) +
  geom_text(data = annotation4, aes (x=x, y=y, label=label), 
            color= "black", size = 4) +
  theme_classic() +
  ggtitle("MBC vs % Nitrogen in Substrate") +
  xlab("% Substrate Nitrogen")+
  ylab("MBC content (mg C/g dry soil)") +
  theme_classic() +
  theme(axis.text.x=element_text(size = 12, colour = "black"),
        axis.text.y=element_text(size = 12, colour = "black"),
        axis.title.x=element_text(size = 14, colour = "black"),
        axis.title.y=element_text(size = 14, colour = "black"), 
        plot.title = element_text(color="black", size=16, hjust = 0.5))



ggplot(Test) + 
  geom_jitter(aes(x = C.N, y = MBC, color = Treatment), size=3) +
  geom_smooth(aes(x = C.N, y = MBC), color="black", method="lm") +
  scale_color_manual(values=c("#E09B1B","#FFD700", "#87CEFA","#6495ED","#1671B3")) +
  geom_text(data = annotation, aes (x=x, y=y, label=label), 
            color= "black", size = 4) +
  theme_classic() +
  ggtitle("MBC vs Substrate C:N") +
  xlab("Substrate C:N")+
  ylab("MBC content (mg C/g dry soil)") +
  theme_classic() +
  theme(axis.text.x=element_text(size = 12, colour = "black"),
        axis.text.y=element_text(size = 12, colour = "black"),
        axis.title.x=element_text(size = 14, colour = "black"),
        axis.title.y=element_text(size = 14, colour = "black"), 
        plot.title = element_text(color="black", size=16, hjust = 0.5))

ggplot(CB) + 
  geom_jitter(aes(y = X.Cresp, x = MBC, color = Treatment), size=3) +
  geom_smooth(aes(y = X.Cresp, x = MBC), color="black", method="lm") +
  scale_color_manual(values=c("#269909F5", "#E09B1B","#FFD700", "#87CEFA","#6495ED","#1671B3")) +
  geom_text(data = annotation2, aes (x=x, y=y, label=label), 
            color= "black", size = 4) +
  theme_classic() +
  ggtitle("MBC vs % Carbon Respired") +
  ylab("% C Respired")+
  xlab("MBC content (mg C/g dry soil)") +
  theme_classic() +
  theme(axis.text.x=element_text(size = 12, colour = "black"),
        axis.text.y=element_text(size = 12, colour = "black"),
        axis.title.x=element_text(size = 14, colour = "black"),
        axis.title.y=element_text(size = 14, colour = "black"),
        plot.title = element_text(color="black", size=16, hjust = 0.5)) 


##STATISTICAL SIGNIFICANCE

#anova test
a.MBC <- aov(MBC ~ Treatment, data = Data)
sa.MBC <- summary(a.MBC)

a.MBN <- aov(MBN ~ Treatment, data = Data)
sa.MBN <- summary(a.MBN)

a.TOC <- aov(TOC ~ Treatment, data = Data)
summary(a.TOC)

a.TN <- aov(TN ~ Treatment, data = Data)
summary(a.TN)


#Tukey Test
t.MBC <- TukeyHSD(a.MBC, conf.level=.95)

t.MBN <- TukeyHSD(a.MBN, conf.level=.95)

t.TOC <- TukeyHSD(a.TOC, conf.level=.95)

t.TN <- TukeyHSD(a.TN, conf.level=.95)

#graph
plot(TukeyHSD(a.MBC, conf.level=.95), las = 2)

plot(TukeyHSD(a.MBN, conf.level=.95), las = 2)

plot(TukeyHSD(a.TOC, conf.level=.95), las = 2)

plot(TukeyHSD(a.TN, conf.level=.95), las = 2)
#if it does not cross 0 than the groups are significantly different

#regressions
modCN = lm(MBC~C.N, data = Test)
summary(modCN)
#P value: 4.976e-05 Adjusted R-squared:  0.4973 highly significant! add to graph

annotation <- data.frame(
  x= c(50, 50), 
  y= c(.15, .13),
  label = c("R-squared = 0.50", "P < 0.0001")
)

modCresp = lm(MBC~X.Cresp, data = CB)
summary(modCresp)
#p-value: 3.267e-06 R-squared:  0.5552 

annotation2 <- data.frame(
  x= c(.29, .29), 
  y= c(.1, .08),
  label = c("R-squared = 0.55", "P < 0.00001")
)

modC = lm(MBC~X.C, data = Test)
summary(modC)
#p-value: 0.0003237, Adjusted R-squared:  0.4121

annotation3 <- data.frame(
  x= c(48, 48), 
  y= c(.15, .13),
  label = c("R-squared = 0.41", "P < 0.001")
)

modN = lm(MBC~X.N, data = Test)
summary(modN)
#p-value: 0.0009139, Adjusted R-squared:  0.3595 

annotation4 <- data.frame(
  x= c(2.5, 2.5), 
  y= c(.35, .33),
  label = c("R-squared = 0.36", "P < 0.001")
)


modlignin = lm(MBC~X.Lignin, data = Test)
summary(modlignin)
#p-value: 0.3484 Adjusted R-squared:  -0.0035

modligninN = lm(MBC~Lignin.N, data = Test)
summary(modligninN)
#p-value: 0.007596 Adjusted R-squared:  0.2396
Lignin:N

annotation5 <- data.frame(
  x= c(35, 35), 
  y= c(.15, .13),
  label = c("R-squared = 0.24", "P < 0.01")
)
