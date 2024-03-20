#load packages
library(tidyverse) #for wrangling and tidying
library(lme4) #for LMMs
library(lmerTest) #for better model output
library(emmeans) #to calculate Estimated Marginal Means (EMMs)
library(performance) #check model assumptions
library(buildmer) #for step-wise regression
library(visdat) #check for missing data

#sort out scientific notation
options(scipen = 999, digits = 3)

#load data
ldt <- read_csv("mixed_models.csv")

#wrangling and tidying
ldt_tidied <- ldt %>%
  #renaming
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
         CosineSimilarity = cos) %>%
  #renaming levels
  mutate(ProbeColour = case_when(ProbeColour == "r" ~ "RED",
                                 ProbeColour == "b" ~ "BLACK"))

#removing irrelevant trials
#remove nonce trials
ldt_tidied <- ldt_tidied %>%
  filter(!ProbeType == "NONCE") %>%
  #remove incorrect trials (standard for LDT analyses)
  filter(!LDTAcc == 0)

#factorising
ldt_tidied <- ldt_tidied %>%
  mutate(ProbeType = factor(ProbeType)) %>% #IV1
  mutate(ProbeColour = factor(ProbeColour)) %>% #IV2
  mutate(Participant = factor(Participant)) %>% #Participant ranfx.
  mutate(Prime = factor(StimID)) #Item ranfx.

#check for missing data
vis_miss(ldt_tided)

#summary stats.
ldt_tidied %>%
  group_by(ProbeType, ProbeColour) %>% #group by condition
  summarise(Mean = mean(RT), SD = sd(RT)) %>% #calculate mean and SD
  arrange(Mean) #arrange fast to slow

#anonymise
ldt_tidied <- ldt_tidied %>%
  select(-UNIQUE_ID) #remove unique identifier column leaving just Participant no.

#dataviz.
ldt_tidied %>%
  #Reorder conditions from fastest to slowest RT, colour by condition
  ggplot(aes(x = ProbeType:ProbeColour, y = log(RT), colour = ProbeColour)) +
  #Produce violin plot
  geom_boxplot() +
  #Add some jitter to avoid overlap
  #No legend as not needed
  guides(colour = "none") + 
  #Provide CIs and mean point
  stat_summary(fun.data = "mean_cl_boot", colour = "black") + 
  #No grey background
  theme_minimal() + 
  #Adjust text size to 12pt
  theme(text = element_text(size = 12)) + 
  #Add axis labels
  labs(x = "ProbeType x ProbeColour", y = "RT(ms.)") + 
  #Add chart title
  ggtitle("Effect of ProbeType x ProbeColour Interaction on RT") +
  #Center chart title
  theme(plot.title = element_text(hjust = 0.5)) +
  #Assign different colors to specific levels of ProbeColour
  scale_colour_manual(values = c("black", "red"))  # Adjust the colors as needed

#contrast coding
contrasts(ldt_tidied$ProbeColour) <- matrix(c(.5, -.5))
contrasts(ldt_tidied$ProbeType) <- matrix(c(.5, -.5))

#building the model
ldt.mod1 <- lmer(RT ~ ProbeType * ProbeColour + #main fx & ixn
                    (1 + ProbeType * ProbeColour | Participant) + #ranfx
                    (1 + ProbeType * ProbeColour | StimID), #ranfx
                  data = ldt_tidied) #dataset

summary(ldt.mod1) #singular fit

#singular fit/convergence issues
#reduce ranefx structure as per literature (Matuschek et al. 2017:
#"Balancing Type I error and power in linear mixed models")
ldt.mod2 <- lmer(RT ~ ProbeType * ProbeColour + 
                    (1 | Participant) +
                    (1 + ProbeType * ProbeColour | StimID),
                  data = ldt_tidied)

summary(ldt.mod2) #singular fit

#repeat process
ldt.mod3 <- lmer(RT ~ ProbeType * ProbeColour + 
                   (1 | Participant) +
                   (1 + ProbeType | StimID),
                 data = ldt_tidied)

summary(ldt.mod3) #singular fit

#repeat process
ldt.mod4 <- lmer(RT ~ ProbeType * ProbeColour + 
                   (1 | Participant) +
                   (1 | StimID),
                 data = ldt_tidied)

summary(ldt.mod4) #OK

#check model assumptions - residuals roughly normal?
check_model(ldt.mod4)

#Estimated Marginal Means (EMMs)

#assume we predicted:
#REL-RED < all other conditions
#UNREL-BLACK > all other conditions
#UNREL-RED similar to REL-BLACK

#EMMs for all possible contrasts
ldt_emms_all <- emmeans(ldt.mod4, #the model
        pairwise #pairwise comparisons
        ~ ProbeType * ProbeColour, #comparisons among conditions
        adjust = "bonferroni") #assuming we predicted differences across the board

#assume we predicted:
#REL-RED < UNREL-RED
#REL-RED < UNREL-BLACK
#UNREL-RED < UNREL-BLACK

#EMMs for specific contrasts
#Bonferroni adjustment = p x no. contrasts of interest
ldt_emms_all_bonf <- ldt_emms_all$contrasts %>% #pull out contrasts table
  data.frame() %>% #make into data frame
  slice(4:6) %>% #extract rows of interest
  mutate(p.value = p.value * 3) #multiply p-vals by number of contrasts

#### --- USING BUILDMER --- ###

#using buildmer to do stepwise regression 'out of the box'
ldt.blmr1 <- buildmer(RT ~ ProbeType * ProbeColour + 
                        (1 + ProbeType * ProbeColour | Participant) +
                        (1 + ProbeType * ProbeColour | Prime),
                      data = ldt_tidied)
summary(ldt.blmr1)

#using buildmer just to fit appropriate ranfx structure
ldt.blmr2 <- buildmer(RT ~ 
                        (1 + ProbeType * ProbeColour | Participant) +
                        (1 + ProbeType * ProbeColour | Prime),
                      data = ldt_tidied,
                      buildmerControl = list(include = ~ ProbeType * ProbeColour))
summary(ldt.blmr2)
