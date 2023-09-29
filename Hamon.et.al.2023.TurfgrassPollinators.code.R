#comparing abundance of insect families/guilds 
# Read in the pan trap abundance by bout data set
alldat = read.csv('sampling.by.date.csv') 
#load sampling dates matched with unique site/date IDs
samplecodes = read.csv('sampling.bout.codes.csv')
#load guild data
guildlist = read.csv('guild.designations.csv')

#load packages
library(tidyverse)
library(ggplot2)
library(lme4) #for constructing mixed models
library(lmerTest) #for displaying p-values
library(car) #for Anova function
library(emmeans) #for lmer posthoc test
library(splusTimeDate) #for fixing date designations
#for biodiversity 
library(BiodiversityR) # also loads vegan

#merge standardized dates on the basis of unique site/date IDs
alldat2 <- merge(alldat, samplecodes, by = c("Date", "Site"))

#filter out Date/Site combinations that are invalid
alldat3 <- alldat2[!(alldat2$toss=="Y"),]

#convert abundance column to numeric
alldat3$Abundance <- as.numeric(alldat3$Abundance) 

#Filter out May and June collections to make sampling time period the same across years.
alldat3 <- subset(alldat3, !(Month %in% c(5, 6)))

#make date.clean read as a date
alldat3$date.clean <- as.Date(alldat3$date.clean, "%m/%d/%Y")

#subset families with 30 or more total individuals across samples
#used for initial family selection
goodfamilies <- alldat2 %>%
  group_by(Family) %>%
  summarise(number = sum(Abundance)) %>%
  filter(number > 30) 

alldat4 <- alldat3 %>%
  right_join(goodfamilies, by = c("Family"))

#create dataframe of number of individuals per family per sampling bout
familyboutinitial <- alldat4 %>%
  group_by(Family, sampling.bout, .drop=FALSE) %>%
  summarise(number=sum(Abundance))
familybout <- familyboutinitial %>% spread(sampling.bout, number, fill = 0) %>% 
  gather(sampling.bout, number, -Family)

#pull just sampling.bout, year, and month designations
myvars <- c("sampling.bout", "Year", "Month", "site.clean", "date.clean", "jd", "jd.year", "jd.adjusted")
samplecodes2 <- samplecodes[myvars]
samplecodes2 <- na.omit(samplecodes2)
samplecodes3 <- samplecodes2 %>% distinct()

#again, clean and optimize data
#merge year designations back into familybout table
familybout2 <- merge(familybout, samplecodes3, by = c("sampling.bout"), all.x = TRUE, all.y = FALSE)
#make date.clean column read as a date
familybout2$date.clean <- as.Date(familybout2$date.clean, "%m/%d/%Y")
#make year a factor
familybout2[,'Year']<-as.character(familybout2[,'Year'])

##transforming data to approach normality
#create a column with sqrt transformed number
familybout2$sqrt.number<-sqrt((familybout2$number)+1)

#testing the effect of year*jd on FAMILY abundance
df = familybout2[(familybout2$Family=="Aphididae"),] #pull out a single family
mod1 <- lmer(sqrt.number~Year*jd + (1|site.clean), data = df)
summary(mod1)
#Use type 3 if interaction is present
Anova(mod1, type=3, test.statistic=c("F"))
emmeans(mod1, specs = pairwise ~ Year, adjust = "tukey")
#visualize model estimates
#load palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

df$fit<-predict(mod1)
Scoliidaeplot<- ggplot(df,aes(jd, fit, linetype=Year, shape=Year)) +
  theme_classic()+
  geom_smooth(method="lm", size=0.8, aes(color=Year)) +
  geom_point(size = 2.4, aes(color=Year))+
  scale_color_manual(values=cbbPalette)+
  labs(x ="Julian day", y =expression(sqrt("Abundance + 1")))+
  xlim(190,303)+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=18,family="serif", colour = "black"),
        axis.text.y = element_text(size=18,family="serif", colour = "black"),
        legend.position = "none")
Scoliidaeplot

#merge alldat2 (i.e. all families) with guild designations
guilddat <- merge(alldat3, guildlist, by = c("Family"), all.x = TRUE, all.y = FALSE)

#select guilds with 30 or more total individuals across samples
#used for initial data exploration
goodguilds <- guilddat %>%
  group_by(Guild) %>%
  summarise(number = sum(Abundance)) %>%
  filter(number > 30) 
guilddat2 <- guilddat %>%
  right_join(goodguilds, by = c("Guild"))

#create dataframe of number of individuals per guild per sampling bout
guildboutinitial <- guilddat2 %>%
  group_by(Guild, sampling.bout, .drop=FALSE) %>%
  summarise(number=sum(Abundance))
guildbout <- guildboutinitial %>% spread(sampling.bout, number, fill = 0) %>% 
  gather(sampling.bout, number, -Guild)

#again, clean and optimize data
#merge dates back into guildbout
guildbout2 <- merge(guildbout, samplecodes3, by = c("sampling.bout"), all.x = TRUE, all.y = FALSE)
#make date.clean column read as a date
guildbout2$date.clean <- as.Date(guildbout2$date.clean, "%m/%d/%Y")
#make year a factor
guildbout2[,'Year']<-as.character(guildbout2[,'Year'])
#create a column with the sqrt transformed number 
guildbout2$sqrt.number <- sqrt((guildbout2$number)+1)

#testing the effect of year*jd on GUILD abundance
df = guildbout2[(guildbout2$Guild=="bees"),] #pull out a single guild
mod1 <- lmer(sqrt.number~Year*jd + (1|site.clean), data = df)
summary(mod1)
#Use type 3 if interaction is present
Anova(mod1, type=3, test.statistic=c("F"))
emmeans(mod1, specs = pairwise ~ Year, adjust = "tukey")
#visualize model estimates
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

df$fit<-predict(mod1)
butterflyplot<- ggplot(df,aes(jd, fit, linetype=Year, shape=Year)) +
  theme_classic()+
  geom_smooth(method="lm", size=0.8, aes(color=Year)) +
  geom_point(size = 2.4, aes(color=Year))+
  scale_color_manual(values=cbbPalette)+
  labs(x ="Julian day", y =expression(sqrt("Abundance + 1")))+
  xlim(190,303)+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x=element_text(size=18,family="serif", colour = "black"),
        axis.text.y=element_text(size=18,family="serif", colour = "black"),
        legend.position = "none")
butterflyplot

#species accumulation curve by year
#convert dataset to community dataset
#subset the rows we need
myvars <- c("sampling.bout", "Family", "Abundance")
commdat <- alldat4[myvars]

library(labdsv) #package for running matrify, to convert to matrix
matcomm<-matrify(commdat)

#create environmental data
myvars2 <- c("sampling.bout", "Year")
envdat <- alldat4[myvars2]
envdat2<-envdat %>% 
  distinct(sampling.bout, .keep_all=TRUE)

#convert year to factor
envdat2[,'Year']<-as.character(envdat2[,'Year'])

#acccum comp result
accum.1 <- accumcomp(matcomm, y = envdat2, factor = 'Year', method = 'exact',
                     conditioned = FALSE, plotit = FALSE)
accum.long1 <- accumcomp.long(accum.1, ci=NA, label.freq=5)

#theme for species accumulation plot
BioR.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  text = element_text(size = 16),
  axis.text = element_text(size = 16, colour = "gray25"),
  axis.title = element_text(size = 18, colour = "gray25"),
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 16),
  legend.key = element_blank())

#species accumulation plot
plotgg2 <- ggplot(data=accum.long1, aes(x = Sites, y = Richness, ymax = UPR, ymin = LWR)) + 
  scale_x_continuous(expand=c(0, 1), sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_color_manual(values = c("2018" = "#000000",
                                "2019"="#E69F00",
                                "2020"="#56B4E9"))+
  scale_fill_manual(values = c("2018" = "#000000",
                               "2019"="#E69F00",
                               "2020"="#56B4E9"))+
  geom_line(aes(colour=Grouping), size=2) +
  geom_point(data=subset(accum.long1, labelit==TRUE), 
             aes(colour=Grouping, shape=Grouping), size=5) +
  geom_ribbon(aes(colour=Grouping, fill = Grouping), alpha=0.2) + 
  BioR.theme +
  labs(x = "No. sampling events", y = "Family diversity", colour = "Grouping", shape = "Grouping")+
  theme(legend.position="none", axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)))
plotgg2







