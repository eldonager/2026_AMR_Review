
library(tidyverse)

AMR_clean<-read.csv("AMR_clean.csv")



AMR_clean<- AMR_clean %>%
  dplyr::mutate(no_isolates_susceptible=no_isolate-no_isolates_resistant)

AMR_clean<- AMR_clean %>% 
  dplyr::select(-no_isolates_intermediate)


AMRdf<- AMR_clean%>%
  dplyr::select(-sampling_start_year)%>%
  dplyr::mutate(Year=sampling_end_year)%>%
  dplyr::select(-sampling_end_year)%>%
  dplyr::select(doi,country, region, Year, species,species_2,strain,no_isolates_resistant,
                no_isolates_susceptible,no_isolate,
                antimicrobial, antibiotic_class, ast_method)




data1<- AMRdf%>%
  dplyr::group_by(doi,country, region, Year, species,species_2,strain,no_isolates_resistant,
                  no_isolates_susceptible,no_isolate,
                  antimicrobial, antibiotic_class, ast_method)%>%
  dplyr::mutate(MeanRes=sum(no_isolates_resistant)/sum(no_isolate))


#data1%>%
# count(MeanRes)%>%
#View()


data1 <- data1[!is.na(data1$Year), ]

data1%>%
  na.omit()
na.omit(data1)

data1 <- data1[data1$MeanRes != 0, ]
# Assuming data1 is your data frame


data1$Year <- as.factor(data1$Year)

data1$Year <- factor(data1$Year)

data1$country <- factor(data1$country)
data1$region <- factor(data1$region)
data1$species <- factor(data1$species)
data1$species_2 <- factor(data1$species_2)
data1$strain <- factor(data1$strain)
data1$antimicrobial <- factor(data1$antimicrobial)
#data1$who_classification <- factor(data1$who_classification)
data1$ast_method <- factor(data1$ast_method)
data1$antibiotic_class<- factor(data1$antibiotic_class)
data1$no_isolates_resistant<- as.numeric(data1$no_isolates_resistant)
data1$no_isolates_susceptible<-as.numeric(data1$no_isolates_susceptible)
#data1$doi <- factor(data1$doi)

#library(lme4)
# Convert doi to numeric index
data1$doi <- as.numeric(as.factor(data1$doi))

data2<- aggregate(MeanRes~region+Year+species+antimicrobial+MeanRes+no_isolate, data = data1, FUN = mean)

##seeing whats happening-nothing much
data2%>%
ggplot(aes(Year, MeanRes, group = region)) +
  geom_line(alpha = 1/3)


##Nesting the dataset

by_region <- data2 %>%
  group_by( region) %>% 
  
  nest()
by_region


##EA
by_region$data[[1]]


##Fitting the model

region_model <- function(data2) {
  glm(MeanRes ~ Year+species+antimicrobial, data = data2, weights = no_isolate)
}

##The df is in list-use purrr::map() to apply country_model to each element
  
models <- map(by_region$data, region_model)


by_region1 <- by_region %>%
  mutate(model = map(data, region_model))
by_region1


#Add residuals
library(modelr)
by_region2 <- by_region1 %>%
  mutate(
    resids = map2(data, model, add_residuals)
  )
by_region2

##unnest

resids <- unnest(by_region2, resids)
resids


p10<-resids %>%
  ggplot(aes(Year, resid)) +
  geom_line(aes(group = region), alpha = 1 / 3) +
  geom_smooth(se = FALSE)+
  facet_wrap(~region)

p10



by_region2 %>%
  mutate(glance = map(model, broom::glance)) %>%
  unnest(glance)


##inspecting models
glance <- by_region2 %>%
  mutate(glance = map(model, broom::glance)) %>%
  unnest(glance, .drop = TRUE)##add .drop=TRUE to drop some variables
glance

##r squared values
glance %>%
  arrange(r.squared)





