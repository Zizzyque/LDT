#load packages
library(tidyverse) #wrangling & tidying
library(lme4) #mixed models
library(lmerTest) #model output viewer
library(emmeans) #estimated marginal means
library(performance) #checking model assumptions
library(buildmer) #step-wise regression
library(visdat) #check for missing data

#sort out scientific notation
options(scipen=999,digits=3)

#load in dataset
ldt <- read_csv("mixed_models.csv")

#annoymise
ldt_tidied <- ldt %>%
  select(-UNIQUE_ID)

#rename variables
ldt_tidied <- ldt_tidied %>%
  rename(Participant = par_id,
         Age = p_age,
         ProbeWord = p_word,
         ProbeType = type_of_probe,
         ProbeColour = colour,
         PrimeFreq = prime_f,
         ProbeFreq = probe_f,
         LDTAcc = acc,
         RT = reaction_time,
         StimID = stim,
         CosineSimilarity = cos)

#renaming levels
ldt_tidied <- ldt_tidied %>%
  mutate(ProbeColour = case_when(ProbeColour == "r" ~ "RED",
                                 ProbeColour == "b" ~ "BLACK"))
#removing irrelevant trials
ldt_tidied <- ldt_tidied %>%
  #remove NONCE trials
  filter(ProbeType !="NONCE") %>%
  #remove inaccurate trials
  filter(!LDTAcc == 0)

#factorising
ldt_tidied <- ldt_tidied %>%
  mutate(ProbeColour = factor(ProbeColour), #IV1
         ProbeType = factor(ProbeType), #IV2
         Participant = factor(Participant), #Participant ranfx
         StimID = factor(StimID)) #Item ranfx

#missing data check
vis_miss(ldt_tidied)

#summary stats
ldt_tidied %>%
  group_by(ProbeType,ProbeColour) %>% #group by condition
  summarise(MeanRT = mean(RT), #calculate means
            SD = sd(RT)) %>% #calculate SDs
  arrange(MeanRT) #arrange fast to slow

#datavis
ldt_tidied %>%
  ggplot(aes(x = ProbeType:ProbeColour, #x AXIS
             y = log(RT), #y axis
             colour = ProbeColour))+ #colour codes levels of ProbeColour
  geom_boxplot()+ #build boxplot
  guides(colour = "none") + stat_summary(fun.data = "mean_cl_boot",
               colour = "black")+
  theme_minimal()+
  theme(text = element_text(size = 12))+
  labs(x = "ProbeType * ProbeColour",
       y = "RT(ms.)")+
  ggtitle("Effect of ProbeType * ProbeColour Interaction on RT")+
  theme(plot.title = element_text(hjust=0.5))+
  scale_colour_manual(values = c("black","red"))

#contrast coding
#contrast codes for ProbeColour
contrasts(ldt_tidied$ProbeColour) <- matrix(c(.5,-.5))
#contrast codes for PT
contrasts(ldt_tidied$ProbeType) <- matrix(c(.5,-.5))

#build linear mixed effects model (LMM)
#maximal model(have everything in it)
ldt.mod1 <-lmer(RT ~ ProbeType * ProbeColour + #main fx & ixn
                  (1 + ProbeType * ProbeColour | Participant) + #ranfx
                  (1 + ProbeType * ProbeColour | StimID), #ranfx
                data = ldt_tidied) #dataset

summary(ldt.mod1) #singular fit

#reduce complexity of ranfx
ldt.mod2 <-lmer(RT ~ ProbeType * ProbeColour + #main fx & ixn
                  (1 | Participant) + #dropped ixn
                  (1 + ProbeType * ProbeColour | StimID), #ranfx
                data = ldt_tidied) #dataset

summary(ldt.mod2)
#要把variance最小的删去；
ldt.mod3 <-lmer(RT ~ ProbeType * ProbeColour + #main fx & ixn
                  (1 | Participant) + #dropped ixn
                  (1 + ProbeType | StimID), #ranfx
                data = ldt_tidied) #dataset
summary(ldt.mod3)

#只剩下probetype了，删去
ldt.mod4 <-lmer(RT ~ ProbeType * ProbeColour + #main fx & ixn
                  (1 | Participant) + #dropped ixn
                  (1 | StimID), #ranfx
                data = ldt_tidied) #dataset
summary(ldt.mod4)

ldt.mod5 <- lmer(RT~ProbeType * ProbeColour +
                   TrialN+PrimeFreq+ProbeFreq+
                   (1+ProbeType*ProbeColour|Participant)+
                   (1+ProbeType*ProbeColour|StimID),
                 data=ldt_tidied)
summary(ldt.mod5)

ldt.mod6 <- lmer(RT ~ ProbeType * ProbeColour +
                   TrialN+PrimeFreq+ProbeFreq+
                   (1|Participant)+
                   (1+ProbeType*ProbeColour|StimID),
                 data=ldt_tidied)
summary(ldt.mod6)


ldt.mod7 <- lmer(RT ~ ProbeType * ProbeColour +
                   TrialN+PrimeFreq+ProbeFreq+
                   (1|Participant)+
                   (1+ProbeType|StimID),
                 data=ldt_tidied)
summary(ldt.mod7)

ldt.mod8 <- lmer(RT ~ ProbeType * ProbeColour +
                   TrialN+PrimeFreq+ProbeFreq+
                   (1|Participant)+
                   (1|StimID),
                 data=ldt_tidied)
summary(ldt.mod8)


#model comparison - Likelihood Ratio Test (LRT)
anova(ldt.mod4, ldt.mod7)

ldt.mod0 <- lmer(RT ~ 1+ #null no IVs
       (1 | Participant)+ #drop ixn
       (1 | StimID),
     data=ldt_tidied)

ldt.mod9 <- lmer(RT ~ CosineSimilarity * ProbeColour + #main fx & ixn
                    TrialN + PrimeFreq + ProbeFreq+ #covariance
                    (1 + CosineSimilarity * ProbeColour| Participant)+ 
                    (1 + CosineSimilarity * ProbeColour| StimID),
                  data=ldt_tidied)
summary(ldt.mod9)

ldt.mod10 <- lmer(RT ~ CosineSimilarity * ProbeColour + #main fx & ixn
                    TrialN + PrimeFreq + ProbeFreq+ #covariance
                   (1 + CosineSimilarity | Participant)+ 
                   (1 + CosineSimilarity | StimID),
                 data=ldt_tidied)
summary(ldt.mod10)

ldt.mod11 <- lmer(RT ~ CosineSimilarity * ProbeColour + #main fx & ixn
                    TrialN + PrimeFreq + ProbeFreq+ #covariance
                    (1 | Participant)+ 
                    (1 | StimID),
                  data=ldt_tidied)
summary(ldt.mod11)
#estimated marginal means for "ALL" possible contrasts

library(emmeans)
ldt_emms_all <- emmeans(ldt.mod4,
                        pairwise ~ ProbeType * ProbeColour,#pairwise comparisons
                        adjust = "bonferroni") #adjust p-values
  
#assume:
#REL-RED < UNREL-RED
#REL-RED < UNREL-BLACK
#UNREL-RED < UNREL-BLACK

#specific contrasts
ldt_emms_spec <- emmeans(ldt.mod4,#model
                        pairwise ~ ProbeType * ProbeColour,#pairwise comparisons
                        adjust = "none") #do not adjust p-values
#correct for contrasts of interest
ldt_emms_all_bonf <- ldt_emms_all$contrasts %>%
  data.frame() %>%
  slice(4:6) %>%
  mutate(p.value = p.value*3)


  
  
  



