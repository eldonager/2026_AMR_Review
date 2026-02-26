#Mapping AST Methods
AMR_clean<-read.csv("AMR_clean.csv")

AMR_clean%>%
  count(doi,species)%>%
  unique()


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



data2<- aggregate(MeanRes~Year+doi+antibiotic_class+country+region+species+no_isolate, data=data1, FUN=mean)
# Convert doi to numeric index
data2$doi <- as.numeric(as.factor(data2$doi))


model1 <- glm(cbind(no_isolates_resistant/no_isolate) ~ Year + country + region + species +antimicrobial+strain+
                ast_method+ (1 | doi),
              family = quasibinomial, weights=no_isolate,data = data1)

#summary(model1)

model2<-glm(cbind(no_isolates_resistant/no_isolate) ~ Year + country + region + species +antimicrobial+strain+
              (1 | doi),
            family=quasibinomial, weights=no_isolate, data = data1)
#summary(model2)


model3<-glm(cbind(no_isolates_resistant/no_isolate) ~ Year + country + region + species +antimicrobial+ (1 | doi),
            family=quasibinomial, weights=no_isolate,
            data = data1)
#summary(model3)

model4<-glm(cbind(no_isolates_resistant/no_isolate) ~ Year + country + region + species + (1 | doi),
            family=quasibinomial, weights=no_isolate, data = data1)
#summary(model4)

model5<-glm(cbind(no_isolates_resistant/no_isolate) ~ Year + country + region + (1 | doi),
            family=quasibinomial, weights=no_isolate, data = data1)
#summary(model5)


model6<-glm(cbind(no_isolates_resistant/no_isolate) ~ Year + country + (1 | doi),
            
            family=quasibinomial, weights=no_isolate, data = data1)
#summary(model6)

model7<- glm(cbind(no_isolates_resistant/no_isolate) ~ Year +  (1 | doi),
             
             family=quasibinomial, weights=no_isolate, data = data1)

#model7




##Diagnostics using boot package
diag_plots<-glm.diag.plots(model7, glmdiag = glm.diag(model7), subset = NULL,
                           iden = FALSE, labels = NULL, ret = FALSE)
#glm.diag.plots(model3, glmdiag = glm.diag(model3), subset = NULL,
#                            iden = FALSE, labels = NULL, ret = FALSE)
# diag_plots
# png("model7_diagnostics_highres.png", 
#     width = 10, height = 8, 
#     units = "in", 
#     res = 300, 
#     type = "cairo")

# Generate the plot
# glm.diag.plots(model7, glmdiag = glm.diag(model7), subset = NULL,
#                iden = FALSE, labels = NULL, ret = FALSE)
# 
# # Close the device
# dev.off()


data2$Year <- as.numeric(as.character(data2$Year))

model <- glm(MeanRes ~ Year, data = data2, family = quasibinomial, weights = no_isolate)

new_preds <- predict(model, newdata = data.frame(Year = seq(min(data2$Year), max(data2$Year), by = 1)), 
                     type = "response", se.fit = TRUE)


pred_df <- data.frame(
  Year = seq(min(data2$Year), max(data2$Year), by = 1),
  fit = new_preds$fit,
  se.fit = new_preds$se.fit
)



fit1<-predict(model,type = "response", se.fit = T)


df1<-as.data.frame(fit1)

df2<- cbind(data2, df1)

df3<- df2%>%
  mutate(lower = fit - 1.96 * se.fit)%>%
  mutate(upper = fit+ 1.96 * se.fit)

df4<- aggregate(MeanRes~doi+Year+region, data = df3, FUN = mean)


##overal trends
library(ggplot2)

trendplot <- ggplot() +
  geom_point(data = df4, aes(x = Year, y = MeanRes),
             shape = 21, color = "black", alpha = 0.4, size = 3) +
  geom_line(data = pred_df, aes(x = Year, y = fit), color = "blue", linewidth = 1) +
  geom_ribbon(data = pred_df, aes(x = Year, ymin = fit - 1.96*se.fit, ymax = fit + 1.96*se.fit),
              fill = "gray70", alpha = 0.3) +
  geom_boxplot(data = data2, aes(x = Year, y = MeanRes, group = Year),
               outlier.shape = NA, fill = "skyblue", color = "black", alpha = 0.2) +
  theme_Publication() +
  labs(x = "Year", y = "Mean Resistance (%)") + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_x_continuous(breaks = seq(min(df4$Year), max(df4$Year), by = 2)) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 14),  # Increased size for x-axis text
    axis.text.y = element_text(size = 14),  # Increased size for y-axis text
    axis.title.x = element_text(size = 16),  # Increased size for x-axis label
    axis.title.y = element_text(size = 16)   # Increased size for y-axis label
  )

trendplot

# ggsave(filename = "trendplotppt.tiff",
#        plot = trendplot,
#        device = "tiff",
#        width = 10,
#        height = 8,
#        units = "in",
#        dpi = 700,
#        compression = "lzw")

#regional trends

#east african trends

ea<- df4%>%
  filter(region=="Eastern Africa")

eadata2<- data2%>%
  filter(region=="Eastern Africa")

eamodel <- glm(MeanRes ~ Year, data = eadata2, family = quasibinomial, weights = no_isolate)

eanew_preds <- predict(eamodel, newdata = data.frame(Year = seq(min(eadata2$Year), max(eadata2$Year), by = 1)), 
                       type = "response", se.fit = TRUE)


eapred_df <- data.frame(
  Year = seq(min(eadata2$Year), max(eadata2$Year), by = 1),
  fit = eanew_preds$fit,
  se.fit = eanew_preds$se.fit
)



eafit1<-predict(eamodel,type = "response", se.fit = T)
print(eafit1)

eadf1<-as.data.frame(eafit1)

eadf2<- cbind(eadata2, eadf1)

eadf3<- eadf2%>%
  mutate(lower = fit - 1.96 * se.fit)%>%
  mutate(upper = fit+ 1.96 * se.fit)

eadf4<- aggregate(MeanRes~doi+Year+region, data = eadf3, FUN = mean)


eaplot <- ggplot() +
  geom_point(data = eadf4, aes(x = Year, y = MeanRes),
             shape = 21, color = "black", alpha = 0.4, size = 3) +
  geom_line(data = eapred_df, aes(x = Year, y = fit), color = "blue", linewidth = 1) +
  geom_ribbon(data = eapred_df, aes(x = Year, ymin = fit - 1.96*se.fit, ymax = fit + 1.96*se.fit),
              fill = "gray70", alpha = 0.3) +
  geom_boxplot(data = eadata2, aes(x = Year, y = MeanRes, group = Year),
               outlier.shape = NA, fill = "skyblue", color = "black", alpha = 0.2) +
  theme_Publication() +
  labs(x = "Year", y = "Mean Resistance (%)", subtitle = "East Africa") + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_x_continuous(breaks = seq(min(eadf4$Year), max(eadf4$Year), by = 2)) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 14),  # Increased size for x-axis text
    axis.text.y = element_text(size = 14),  # Increased size for y-axis text
    axis.title.x = element_text(size = 16),  # Increased size for x-axis label
    axis.title.y = element_text(size = 16)   # Increased size for y-axis label
  )



#North Africa plot

na<- df4%>%
  filter(region=="Northern Africa")

nadata2<- data2%>%
  filter(region=="Northern Africa")

namodel <- glm(MeanRes ~ Year, data = nadata2, family = quasibinomial, weights = no_isolate)

nanew_preds <- predict(namodel, newdata = data.frame(Year = seq(min(nadata2$Year), max(nadata2$Year), by = 1)), 
                       type = "response", se.fit = TRUE)


napred_df <- data.frame(
  Year = seq(min(nadata2$Year), max(nadata2$Year), by = 1),
  fit = nanew_preds$fit,
  se.fit = nanew_preds$se.fit
)



nafit1<-predict(namodel,type = "response", se.fit = T)


nadf1<-as.data.frame(nafit1)

nadf2<- cbind(nadata2, nadf1)

nadf3<- nadf2%>%
  mutate(lower = fit - 1.96 * se.fit)%>%
  mutate(upper = fit+ 1.96 * se.fit)

nadf4<- aggregate(MeanRes~doi+Year+region, data = nadf3, FUN = mean)


naplot <- ggplot() +
  geom_point(data = nadf4, aes(x = Year, y = MeanRes),
             shape = 21, color = "black", alpha = 0.4, size = 3) +
  geom_line(data = napred_df, aes(x = Year, y = fit), color = "blue", linewidth = 1) +
  geom_ribbon(data = napred_df, aes(x = Year, ymin = fit - 1.96*se.fit, ymax = fit + 1.96*se.fit),
              fill = "gray70", alpha = 0.3) +
  geom_boxplot(data = nadata2, aes(x = Year, y = MeanRes, group = Year),
               outlier.shape = NA, fill = "skyblue", color = "black", alpha = 0.2) +
  theme_Publication() +
  labs(x = "Year", y = "Mean Resistance (%)", subtitle = "North Africa") + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_x_continuous(breaks = seq(min(nadf4$Year), max(nadf4$Year), by = 2)) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 14),  # Increased size for x-axis text
    axis.text.y = element_text(size = 14),  # Increased size for y-axis text
    axis.title.x = element_text(size = 16),  # Increased size for x-axis label
    axis.title.y = element_text(size = 16)   # Increased size for y-axis label
  )


#West africa

wa<- df4%>%
  filter(region=="West Africa")

wadata2<- data2%>%
  filter(region=="West Africa")

wamodel <- glm(MeanRes ~ Year, data = wadata2, family = quasibinomial, weights = no_isolate)

wanew_preds <- predict(wamodel, newdata = data.frame(Year = seq(min(wadata2$Year), max(wadata2$Year), by = 1)), 
                       type = "response", se.fit = TRUE)


wapred_df <- data.frame(
  Year = seq(min(wadata2$Year), max(wadata2$Year), by = 1),
  fit = wanew_preds$fit,
  se.fit = wanew_preds$se.fit
)



wafit1<-predict(wamodel,type = "response", se.fit = T)


wadf1<-as.data.frame(wafit1)

wadf2<- cbind(wadata2, wadf1)

wadf3<- wadf2%>%
  mutate(lower = fit - 1.96 * se.fit)%>%
  mutate(upper = fit+ 1.96 * se.fit)

wadf4<- aggregate(MeanRes~doi+Year+region, data = wadf3, FUN = mean)


waplot <- ggplot() +
  geom_point(data = wadf4, aes(x = Year, y = MeanRes),
             shape = 21, color = "black", alpha = 0.4, size = 3) +
  geom_line(data = wapred_df, aes(x = Year, y = fit), color = "blue", linewidth = 1) +
  geom_ribbon(data = wapred_df, aes(x = Year, ymin = fit - 1.96*se.fit, ymax = fit + 1.96*se.fit),
              fill = "gray70", alpha = 0.3) +
  geom_boxplot(data = wadata2, aes(x = Year, y = MeanRes, group = Year),
               outlier.shape = NA, fill = "skyblue", color = "black", alpha = 0.2) +
  theme_Publication() +
  labs(x = "Year", y = "Mean Resistance (%)", subtitle = "West Africa") + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_x_continuous(breaks = seq(min(wadf4$Year), max(wadf4$Year), by = 2)) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 14),  # Increased size for x-axis text
    axis.text.y = element_text(size = 14),  # Increased size for y-axis text
    axis.title.x = element_text(size = 16),  # Increased size for x-axis label
    axis.title.y = element_text(size = 16)   # Increased size for y-axis label
  )






#south africa

sa<- df4%>%
  filter(region=="Southern Africa")

sadata2<- data2%>%
  filter(region=="Southern Africa")

samodel <- glm(MeanRes ~ Year, data = sadata2, family = quasibinomial, weights = no_isolate)

sanew_preds <- predict(samodel, newdata = data.frame(Year = seq(min(sadata2$Year), max(sadata2$Year), by = 1)), 
                       type = "response", se.fit = TRUE)


sapred_df <- data.frame(
  Year = seq(min(sadata2$Year), max(sadata2$Year), by = 1),
  fit = sanew_preds$fit,
  se.fit = sanew_preds$se.fit
)



safit1<-predict(samodel,type = "response", se.fit = T)


sadf1<-as.data.frame(safit1)

sadf2<- cbind(sadata2, sadf1)

sadf3<- sadf2%>%
  mutate(lower = fit - 1.96 * se.fit)%>%
  mutate(upper = fit+ 1.96 * se.fit)

sadf4<- aggregate(MeanRes~doi+Year+region, data = sadf3, FUN = mean)


saplot <- ggplot() +
  geom_point(data = sadf4, aes(x = Year, y = MeanRes),
             shape = 21, color = "black", alpha = 0.4, size = 3) +
  geom_line(data = sapred_df, aes(x = Year, y = fit), color = "blue", linewidth = 1) +
  geom_ribbon(data = sapred_df, aes(x = Year, ymin = fit - 1.96*se.fit, ymax = fit + 1.96*se.fit),
              fill = "gray70", alpha = 0.3) +
  geom_boxplot(data = sadata2, aes(x = Year, y = MeanRes, group = Year),
               outlier.shape = NA, fill = "skyblue", color = "black", alpha = 0.2) +
  theme_Publication() +
  labs(x = "Year", y = "Mean Resistance (%)",
       subtitle = "Southern Africa") + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_x_continuous(breaks = seq(min(sadf4$Year), max(sadf4$Year), by = 2)) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 14),  # Increased size for x-axis text
    axis.text.y = element_text(size = 14),  # Increased size for y-axis text
    axis.title.x = element_text(size = 16),  # Increased size for x-axis label
    axis.title.y = element_text(size = 16)   # Increased size for y-axis label
  )
saplot


cbndp<-ggarrange(eaplot , waplot, naplot , saplot,
          labels = c("(a)", "(b)", "(c)", "(d)"),
          common.legend = TRUE)
# 
# ggsave(filename = "allregionscbndtrends.tiff",
#        plot = cbndp,
#        device = "tiff",
#        width = 10,
#      height = 8,
#        units = "in",
#        dpi = 700,
#        compression = "lzw")