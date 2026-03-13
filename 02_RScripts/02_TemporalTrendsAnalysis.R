###########
#This script is used to generate descriptive statistics and temporal trends of AMR data. 

###########

library(ggpubr)

library(boot)
library(lme4)
library(boot)
library(ggthemes)
library(patchwork)
library(tidyverse)


#code to generate publication quality figures
theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}



library(tidyverse)
library(boot)
library(ggthemes)
library(patchwork)



# Load data
AMR_clean <- read.csv("AMR_clean.csv", stringsAsFactors = FALSE)

# Data wrangling
AMR_clean <- AMR_clean %>%
  mutate(no_isolates_susceptible = no_isolate - no_isolates_resistant) %>%
  select(-no_isolates_intermediate, -sampling_start_year) %>%
  rename(Year = sampling_end_year) %>%
  filter(
    !is.na(Year),
    !is.na(no_isolates_resistant),
    !is.na(no_isolates_susceptible),
    no_isolate > 0
  )

# Aggregate counts per study-antibiotic unit then compute MeanRes
data1 <- AMR_clean %>%
  mutate(doi_id = as.numeric(as.factor(doi))) %>%
  group_by(doi_id, region, Year, species, species_2,
           strain, antimicrobial, antibiotic_class, ast_method) %>%
  summarise(
    no_isolates_resistant   = sum(no_isolates_resistant,   na.rm = TRUE),
    no_isolates_susceptible = sum(no_isolates_susceptible, na.rm = TRUE),
    no_isolate              = sum(no_isolate,              na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    MeanRes          = no_isolates_resistant / no_isolate,
    Year_f           = factor(Year),
    region           = factor(region),
    species          = factor(species),
    species_2        = factor(species_2),
    strain           = factor(strain),
    antimicrobial    = factor(antimicrobial),
    antibiotic_class = factor(antibiotic_class),
    ast_method       = factor(ast_method)
  ) %>%
  filter(!is.na(MeanRes))

# Aggregate to study-year-region level for trend models and plots
data2 <- data1 %>%
  group_by(doi_id, Year, region) %>%
  summarise(
    MeanRes                 = mean(MeanRes,                na.rm = TRUE),
    no_isolates_resistant   = sum(no_isolates_resistant,   na.rm = TRUE),
    no_isolates_susceptible = sum(no_isolates_susceptible, na.rm = TRUE),
    no_isolate              = sum(no_isolate,              na.rm = TRUE),
    .groups = "drop"
  )

###########################################################################
# Models fitted on data1 with Year as factor
# Purpose: identify which covariates drive resistance 
# quasibinomial for overdispersion; weights = log10(no_isolate+1). Quasibinomial is chosen because 
#AMR data is overdispersed by nature
###########################################################################

model1 <- glm(MeanRes ~ Year_f + region + species + antimicrobial + strain + ast_method,
              family = quasibinomial, weights = log10(no_isolate + 1), data = data1)
summary(model1)

model2 <- glm(MeanRes ~ Year_f + region + species + antimicrobial + strain,
              family = quasibinomial, weights = log10(no_isolate + 1), data = data1)
summary(model2)

model3 <- glm(MeanRes ~ Year_f + region + species + antimicrobial,
              family = quasibinomial, weights = log10(no_isolate + 1), data = data1)
summary(model3)

model4 <- glm(MeanRes ~ Year_f + region + species,
              family = quasibinomial, weights = log10(no_isolate + 1), data = data1)
summary(model4)

model5 <- glm(MeanRes ~ Year_f + region,
              family = quasibinomial, weights = log10(no_isolate + 1), data = data1)
summary(model5)

model6 <- glm(MeanRes ~ Year_f + species,
              family = quasibinomial, weights = log10(no_isolate + 1), data = data1)
summary(model6)

model7 <- glm(MeanRes ~ Year_f,
              family = quasibinomial, weights = log10(no_isolate + 1), data = data1)
summary(model7)

# QAIC for model selection - lower = better fit
# QAIC = -2*(logLik/dispersion) + 2K
dispersion1 <- summary(model1)$dispersion
QAIC1 <- -2 * (-0.5 * model1$deviance / dispersion1) + 2 * length(coef(model1))
print(paste("Model 1 QAIC:", round(QAIC1, 2)))

dispersion2 <- summary(model2)$dispersion
QAIC2 <- -2 * (-0.5 * model2$deviance / dispersion2) + 2 * length(coef(model2))
print(paste("Model 2 QAIC:", round(QAIC2, 2)))

dispersion3 <- summary(model3)$dispersion
QAIC3 <- -2 * (-0.5 * model3$deviance / dispersion3) + 2 * length(coef(model3))
print(paste("Model 3 QAIC:", round(QAIC3, 2)))

dispersion4 <- summary(model4)$dispersion
QAIC4 <- -2 * (-0.5 * model4$deviance / dispersion4) + 2 * length(coef(model4))
print(paste("Model 4 QAIC:", round(QAIC4, 2)))

dispersion5 <- summary(model5)$dispersion
QAIC5 <- -2 * (-0.5 * model5$deviance / dispersion5) + 2 * length(coef(model5))
print(paste("Model 5 QAIC:", round(QAIC5, 2)))

dispersion6 <- summary(model6)$dispersion
QAIC6 <- -2 * (-0.5 * model6$deviance / dispersion6) + 2 * length(coef(model6))
print(paste("Model 6 QAIC:", round(QAIC6, 2)))

dispersion7 <- summary(model7)$dispersion
QAIC7 <- -2 * (-0.5 * model7$deviance / dispersion7) + 2 * length(coef(model7))
print(paste("Model 7 QAIC:", round(QAIC7, 2)))

###############################################################
########################################################
#Model 1 has the lowest QAIC (2540.5) — it is the best fitting model.
#This means that all covariates matter: region, species, antimicrobial, strain, and AST method 
#all contribute meaningfully to explaining resistance variation. 
####################################################################
#######################################################################
##################################


# Diagnostic plots for best model 
par(mfrow = c(2, 2))
plot(fitted(model1), residuals(model1, type = "pearson"),
     xlab = "Fitted values", ylab = "Pearson residuals",
     main = "Residuals vs Fitted", pch = 16, col = "steelblue")
abline(h = 0, lty = 2, col = "red")

qqnorm(residuals(model1, type = "pearson"), main = "Normal Q-Q",
       pch = 16, col = "steelblue")
qqline(residuals(model1, type = "pearson"), col = "red")

plot(fitted(model1), sqrt(abs(residuals(model1, type = "pearson"))),
     xlab = "Fitted values", ylab = "sqrt(|Pearson residuals|)",
     main = "Scale-Location", pch = 16, col = "steelblue")

plot(hatvalues(model1), residuals(model1, type = "pearson"),
     xlab = "Leverage", ylab = "Pearson residuals",
     main = "Leverage vs Residuals", pch = 16, col = "steelblue")
abline(h = 0, lty = 2, col = "red")
par(mfrow = c(1, 1))




########################
##########################
# Diagnostic plots for the best-fitting variable selection model (model 1, QAIC = 2540.5).
# The Residuals vs Fitted plot shows a funnel-shaped pattern with greater spread at 
# low fitted values narrowing at higher values, which is expected for AMR
# data where most studies report low resistance.
# The Normal Q-Q plot shows moderate deviation from the reference line 
# at both tails, indicating right-skewed residuals; this reflects the skewed distribution
# of resistance proportions in the dataset where many studies report low resistance 
# and few report very high resistance. Since quasibinomial regression does not 
# assume normally distributed residuals this does not invalidate the model. 
# The Scale-Location plot shows decreasing variance as fitted values increase, 
# a pattern inherent to outcomes bounded between 0 and 1 that the quasibinomial 
# variance function partially accommodates. The Leverage vs Residuals plot shows 
# that the majority of observations have low leverage and cluster near zero residuals, 
# indicating they do not exert undue influence on model estimates; 
# a small number of high-residual observations are present but their low leverage
# means they do not meaningfully distort the fitted coefficients. 
# Taken together, these diagnostics indicate no major violations of
# model assumptions and support the validity of model 1 for inference.
################################################################################
#################################################################################
################################





###########################################################################
# PART 2: TEMPORAL TREND
# Fitted on data2 (study-year-region aggregated) with numeric Year
# Purpose: test whether AMR prevalence increases or decreases over time
# Coefficient of Year: positive = increasing, negative = decreasing
# p-value of Year: significance of trend
###########################################################################

model <- glm(MeanRes ~ Year,
             family  = quasibinomial,
             weights = log10(no_isolate + 1),
             data    = data2)
summary(model)
print(summary(model)$coefficients["Year", ])

###############################
# coefficient = 0.075, p < 0.001, meaning AMR is significantly increasing over time. 
# The positive coefficient on Year means that for every one year increase, 
# mean resistance increases on the log-odds scale by 0.075. 
#########################################
##############################################################




new_preds <- predict(model,
                     newdata = data.frame(Year = seq(min(data2$Year), max(data2$Year), by = 1)),
                     type = "response", se.fit = TRUE)

pred_df <- data.frame(
  Year   = seq(min(data2$Year), max(data2$Year), by = 1),
  fit    = new_preds$fit,
  se.fit = new_preds$se.fit
)

fit1 <- predict(model, type = "response", se.fit = TRUE)
df1  <- as.data.frame(fit1)
df2  <- cbind(data2, df1)

df3 <- df2 %>%
  mutate(lower = fit - 1.96 * se.fit,
         upper = fit + 1.96 * se.fit)

df4 <- aggregate(MeanRes ~ doi_id + Year + region, data = df3, FUN = mean)

# overall temporal trend
trendplot <- ggplot() +
  geom_point(data = df4, aes(x = Year, y = MeanRes),
             shape = 21, color = "black", alpha = 0.4, size = 3) +
  geom_line(data = pred_df, aes(x = Year, y = fit),
            color = "blue", linewidth = 1) +
  geom_ribbon(data = pred_df, aes(x = Year,
                                  ymin = fit - 1.96 * se.fit,
                                  ymax = fit + 1.96 * se.fit),
              fill = "gray70", alpha = 0.3) +
  geom_boxplot(data = data2, aes(x = Year, y = MeanRes, group = Year),
               outlier.shape = NA, fill = "skyblue", color = "black", alpha = 0.2) +
  theme_Publication() +
  labs(x = "Year", y = "Mean Resistance") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_x_continuous(breaks = seq(min(df4$Year), max(df4$Year), by = 2)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

trendplot

###########################################################################
# PART 3: REGIONAL TRENDS
# MeanRes ~ Year model fitted separately per region on data2
###########################################################################

# Eastern Africa
eadata2 <- data2 %>% filter(region == "Eastern Africa")

eamodel <- glm(MeanRes ~ Year,
               family  = quasibinomial,
               weights = log10(no_isolate + 1),
               data    = eadata2)
summary(eamodel)
print(summary(eamodel)$coefficients["Year", ])

eanew_preds <- predict(eamodel,
                       newdata = data.frame(Year = seq(min(eadata2$Year), max(eadata2$Year), by = 1)),
                       type = "response", se.fit = TRUE)

eapred_df <- data.frame(
  Year   = seq(min(eadata2$Year), max(eadata2$Year), by = 1),
  fit    = eanew_preds$fit,
  se.fit = eanew_preds$se.fit
)

eafit1 <- predict(eamodel, type = "response", se.fit = TRUE)
eadf1  <- as.data.frame(eafit1)
eadf2  <- cbind(eadata2, eadf1)
eadf3  <- eadf2 %>% mutate(lower = fit - 1.96 * se.fit, upper = fit + 1.96 * se.fit)
eadf4  <- aggregate(MeanRes ~ doi_id + Year + region, data = eadf3, FUN = mean)

eaplot <- ggplot() +
  geom_point(data = eadf4, aes(x = Year, y = MeanRes),
             shape = 21, color = "black", alpha = 0.4, size = 3) +
  geom_line(data = eapred_df, aes(x = Year, y = fit),
            color = "blue", linewidth = 1) +
  geom_ribbon(data = eapred_df, aes(x = Year,
                                    ymin = fit - 1.96 * se.fit,
                                    ymax = fit + 1.96 * se.fit),
              fill = "gray70", alpha = 0.3) +
  geom_boxplot(data = eadata2, aes(x = Year, y = MeanRes, group = Year),
               outlier.shape = NA, fill = "skyblue", color = "black", alpha = 0.2) +
  theme_Publication() +
  labs(x = "Year", y = "Mean Resistance", subtitle = "East Africa") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_x_continuous(breaks = seq(min(eadf4$Year), max(eadf4$Year), by = 2)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

eaplot


# Northern Africa
nadata2 <- data2 %>% filter(region == "Northern Africa")

namodel <- glm(MeanRes ~ Year,
               family  = quasibinomial,
               weights = log10(no_isolate + 1),
               data    = nadata2)
summary(namodel)
print(summary(namodel)$coefficients["Year", ])

nanew_preds <- predict(namodel,
                       newdata = data.frame(Year = seq(min(nadata2$Year), max(nadata2$Year), by = 1)),
                       type = "response", se.fit = TRUE)

napred_df <- data.frame(
  Year   = seq(min(nadata2$Year), max(nadata2$Year), by = 1),
  fit    = nanew_preds$fit,
  se.fit = nanew_preds$se.fit
)

nafit1 <- predict(namodel, type = "response", se.fit = TRUE)
nadf1  <- as.data.frame(nafit1)
nadf2  <- cbind(nadata2, nadf1)
nadf3  <- nadf2 %>% mutate(lower = fit - 1.96 * se.fit, upper = fit + 1.96 * se.fit)
nadf4  <- aggregate(MeanRes ~ doi_id + Year + region, data = nadf3, FUN = mean)

naplot <- ggplot() +
  geom_point(data = nadf4, aes(x = Year, y = MeanRes),
             shape = 21, color = "black", alpha = 0.4, size = 3) +
  geom_line(data = napred_df, aes(x = Year, y = fit),
            color = "blue", linewidth = 1) +
  geom_ribbon(data = napred_df, aes(x = Year,
                                    ymin = fit - 1.96 * se.fit,
                                    ymax = fit + 1.96 * se.fit),
              fill = "gray70", alpha = 0.3) +
  geom_boxplot(data = nadata2, aes(x = Year, y = MeanRes, group = Year),
               outlier.shape = NA, fill = "skyblue", color = "black", alpha = 0.2) +
  theme_Publication() +
  labs(x = "Year", y = "Mean Resistance", subtitle = "North Africa") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_x_continuous(breaks = seq(min(nadf4$Year), max(nadf4$Year), by = 2)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

naplot


# West Africa
wadata2 <- data2 %>% filter(region == "West Africa")

wamodel <- glm(MeanRes ~ Year,
               family  = quasibinomial,
               weights = log10(no_isolate + 1),
               data    = wadata2)
summary(wamodel)
print(summary(wamodel)$coefficients["Year", ])

wanew_preds <- predict(wamodel,
                       newdata = data.frame(Year = seq(min(wadata2$Year), max(wadata2$Year), by = 1)),
                       type = "response", se.fit = TRUE)

wapred_df <- data.frame(
  Year   = seq(min(wadata2$Year), max(wadata2$Year), by = 1),
  fit    = wanew_preds$fit,
  se.fit = wanew_preds$se.fit
)

wafit1 <- predict(wamodel, type = "response", se.fit = TRUE)
wadf1  <- as.data.frame(wafit1)
wadf2  <- cbind(wadata2, wadf1)
wadf3  <- wadf2 %>% mutate(lower = fit - 1.96 * se.fit, upper = fit + 1.96 * se.fit)
wadf4  <- aggregate(MeanRes ~ doi_id + Year + region, data = wadf3, FUN = mean)

waplot <- ggplot() +
  geom_point(data = wadf4, aes(x = Year, y = MeanRes),
             shape = 21, color = "black", alpha = 0.4, size = 3) +
  geom_line(data = wapred_df, aes(x = Year, y = fit),
            color = "blue", linewidth = 1) +
  geom_ribbon(data = wapred_df, aes(x = Year,
                                    ymin = fit - 1.96 * se.fit,
                                    ymax = fit + 1.96 * se.fit),
              fill = "gray70", alpha = 0.3) +
  geom_boxplot(data = wadata2, aes(x = Year, y = MeanRes, group = Year),
               outlier.shape = NA, fill = "skyblue", color = "black", alpha = 0.2) +
  theme_Publication() +
  labs(x = "Year", y = "Mean Resistance", subtitle = "West Africa") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_x_continuous(breaks = seq(min(wadf4$Year), max(wadf4$Year), by = 2)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

waplot


# Southern Africa
sadata2 <- data2 %>% filter(region == "Southern Africa")

samodel <- glm(MeanRes ~ Year,
               family  = quasibinomial,
               weights = log10(no_isolate + 1),
               data    = sadata2)
summary(samodel)
print(summary(samodel)$coefficients["Year", ])

sanew_preds <- predict(samodel,
                       newdata = data.frame(Year = seq(min(sadata2$Year), max(sadata2$Year), by = 1)),
                       type = "response", se.fit = TRUE)

sapred_df <- data.frame(
  Year   = seq(min(sadata2$Year), max(sadata2$Year), by = 1),
  fit    = sanew_preds$fit,
  se.fit = sanew_preds$se.fit
)

safit1 <- predict(samodel, type = "response", se.fit = TRUE)
sadf1  <- as.data.frame(safit1)
sadf2  <- cbind(sadata2, sadf1)
sadf3  <- sadf2 %>% mutate(lower = fit - 1.96 * se.fit, upper = fit + 1.96 * se.fit)
sadf4  <- aggregate(MeanRes ~ doi_id + Year + region, data = sadf3, FUN = mean)


sapred_df <- data.frame(
  Year   = seq(min(sadata2$Year), max(sadata2$Year), by = 1),
  fit    = sanew_preds$fit,
  se.fit = sanew_preds$se.fit
) %>%
  mutate(
    lower = pmax(0, fit - 1.96 * se.fit),
    upper = pmin(1, fit + 1.96 * se.fit)
  )


head(sapred_df)


saplot <- ggplot() +
  geom_point(data = sadf4, aes(x = Year, y = MeanRes),
             shape = 21, color = "black", alpha = 0.4, size = 3) +
  geom_line(data = sapred_df, aes(x = Year, y = fit),
            color = "blue", linewidth = 1) +
  geom_ribbon(data = sapred_df, aes(x = Year, ymin = lower, ymax = upper),
              fill = "gray70", alpha = 0.3) +
  geom_boxplot(data = sadata2, aes(x = Year, y = MeanRes, group = Year),
               outlier.shape = NA, fill = "skyblue", color = "black", alpha = 0.2) +
  theme_Publication() +
  labs(x = "Year", y = "Mean Resistance", subtitle = "Southern Africa") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_x_continuous(breaks = seq(min(sadf4$Year), max(sadf4$Year), by = 2)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

saplot

figure5 <- (eaplot | naplot) / (waplot | saplot) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  ) &
  theme(plot.tag = element_text(face = "bold", size = 14))

figure5



# Save figures
ggsave("trendplot.tiff", plot = trendplot, device = "tiff",
       width = 10, height = 8, units = "in", dpi = 700, compression = "lzw")

ggsave("figure5_regional_trends.tiff", plot = figure5, device = "tiff",
       width = 14, height = 12, units = "in", dpi = 700, compression = "lzw")



# #Loading data
# AMR_clean<-read.csv("AMR_clean.csv")
# 
# 
# #data wrangling
# AMR_clean<- AMR_clean %>%
#   dplyr::mutate(no_isolates_susceptible=no_isolate-no_isolates_resistant)
# 
# AMR_clean<- AMR_clean %>% 
#   dplyr::select(-no_isolates_intermediate)
# 
# 
# AMRdf<- AMR_clean%>%
#   dplyr::select(-sampling_start_year)%>%
#   dplyr::mutate(Year=sampling_end_year)%>%
# dplyr::select(-sampling_end_year)%>%
#   dplyr::select(doi,country, region, Year, species,species_2,strain,no_isolates_resistant,
#                 no_isolates_susceptible,no_isolate,
#                 antimicrobial, antibiotic_class, ast_method)
#   
# 
# 
# 
# data1<- AMRdf%>%
#  dplyr::group_by(doi,country, region, Year, species,species_2,strain,no_isolates_resistant,
#                  no_isolates_susceptible,no_isolate,
#                 antimicrobial, antibiotic_class, ast_method)%>%
#   dplyr::mutate(MeanRes=sum(no_isolates_resistant)/sum(no_isolate))
#   
# 
# #data1%>%
#  # count(MeanRes)%>%
#   #View()
# 
# 
# data1 <- data1[!is.na(data1$Year), ]
# 
# data1%>%
#   na.omit()
# na.omit(data1)
# 
# data1 <- data1[data1$MeanRes != 0, ]
# # Assuming data1 is your data frame
# 
# 
# data1$Year <- as.factor(data1$Year)
# 
# data1$Year <- factor(data1$Year)
# 
# data1$country <- factor(data1$country)
# data1$region <- factor(data1$region)
# data1$species <- factor(data1$species)
# data1$species_2 <- factor(data1$species_2)
# data1$strain <- factor(data1$strain)
# data1$antimicrobial <- factor(data1$antimicrobial)
# #data1$who_classification <- factor(data1$who_classification)
# data1$ast_method <- factor(data1$ast_method)
# data1$antibiotic_class<- factor(data1$antibiotic_class)
# data1$no_isolates_resistant<- as.numeric(data1$no_isolates_resistant)
# data1$no_isolates_susceptible<-as.numeric(data1$no_isolates_susceptible)
# #data1$doi <- factor(data1$doi)
# 
# #library(lme4)
# # Convert doi to numeric index
# data1$doi <- as.numeric(as.factor(data1$doi))
# 
# 
# #aggregate data
# data2<- aggregate(MeanRes~Year+doi+antibiotic_class+country+region+species+no_isolate, data=data1, FUN=mean)
# # Convert doi to numeric index
# data2$doi <- as.numeric(as.factor(data2$doi))
# 
# 
# 
# #fitting the glm models starting with the most complex model
# #study id(doi) is added as a random effect
# #The model is fitted using a quasi-binomial distribution to account for overdispersion
# 
# model1 <- glm(cbind(no_isolates_resistant/no_isolate) ~ Year + country + region + species +antimicrobial+strain+
#                 ast_method+ (1 | doi),
#               family = quasibinomial, weights=no_isolate,data = data1)
# 
# summary(model1)
# 
# model2<-glm(cbind(no_isolates_resistant/no_isolate) ~ Year + country + region + species +antimicrobial+strain+
#               (1 | doi),
#             family=quasibinomial, weights=no_isolate, data = data1)
# summary(model2)
# 
# 
# model3<-glm(cbind(no_isolates_resistant/no_isolate) ~ Year + country + region + species +antimicrobial+ (1 | doi),
#             family=quasibinomial, weights=no_isolate,
#             data = data1)
# summary(model3)
# 
# model4<-glm(cbind(no_isolates_resistant/no_isolate) ~ Year + country + region + species + (1 | doi),
#             family=quasibinomial, weights=no_isolate, data = data1)
# summary(model4)
# 
# model5<-glm(cbind(no_isolates_resistant/no_isolate) ~ Year + country + region + (1 | doi),
#             family=quasibinomial, weights=no_isolate, data = data1)
# summary(model5)
# 
# 
# model6<-glm(cbind(no_isolates_resistant/no_isolate) ~ Year + country + (1 | doi),
#             
#             family=quasibinomial, weights=no_isolate, data = data1)
# summary(model6)
# 
# model7<- glm(cbind(no_isolates_resistant/no_isolate) ~ Year +  (1 | doi),
#              
#              family=quasibinomial, weights=no_isolate, data = data1)
# 
# model7
# 
# 
# #Model selection using QAIC
# #Note that we use QAIC instead of AIC because we are using quasi-binomial distribution which accounts for overdispersion
# #The QAIC is calculated as -2*logLik + 2*K where K is the number of parameters in the model
# #The model with the lowest QAIC is the best model
# #we calculated all the QAIC values for the models using the following code
# #example provided below is for model 7 which has the lowest QAIC value
# #you can replace model7 with the model of interest
# dispersion <- summary(model7)$dispersion
# logLik_quasi <- -0.5 * model7$deviance / dispersion
# K <- length(coef(model7))  # Number of parameters
# QAIC <- -2 * logLik_quasi + 2 * K
# print(QAIC)
# QAIC(model7)
# 
# 
# 
# 
# #We use the glm.diag.plots function from the boot package to generate diagnostic plots for the model
# #The diagnostic plots are used to check the assumptions of the model
# #note that this code was used to produce the diagnostic plots (Figure 2) in the manuscript
# diag_plots<-glm.diag.plots(model7, glmdiag = glm.diag(model7), subset = NULL,
#                iden = FALSE, labels = NULL, ret = FALSE)
# #glm.diag.plots(model3, glmdiag = glm.diag(model3), subset = NULL,
# #                            iden = FALSE, labels = NULL, ret = FALSE)
# # diag_plots
# 
# # Save the diagnostic plots to a file
# # png("model7_diagnostics_highres.png", 
# #     width = 10, height = 8, 
# #     units = "in", 
# #     res = 300, 
# #     type = "cairo")
# 
# # Generate the plot
# # glm.diag.plots(model7, glmdiag = glm.diag(model7), subset = NULL,
# #                iden = FALSE, labels = NULL, ret = FALSE)
# # 
# # # Close the device
# # dev.off()
# 
# 
# data2$Year <- as.numeric(as.character(data2$Year))
# 
# model <- glm(MeanRes ~ Year, data = data2, family = quasibinomial, weights = no_isolate)
# 
# #predicting the model
# new_preds <- predict(model, newdata = data.frame(Year = seq(min(data2$Year), max(data2$Year), by = 1)), 
#                      type = "response", se.fit = TRUE)
# 
# #create a data frame with the predicted values
# pred_df <- data.frame(
#   Year = seq(min(data2$Year), max(data2$Year), by = 1),
#   fit = new_preds$fit,
#   se.fit = new_preds$se.fit
# )
# 
# 
# fit1<-predict(model,type = "response", se.fit = T)
# 
# 
# df1<-as.data.frame(fit1)
# 
# df2<- cbind(data2, df1)
# 
# #calculate the confidence intervals
# df3<- df2%>%
#   mutate(lower = fit - 1.96 * se.fit)%>%
#   mutate(upper = fit+ 1.96 * se.fit)
# 
# df4<- aggregate(MeanRes~doi+Year+region, data = df3, FUN = mean)
# 
# 
# ##overal trends
# library(ggplot2)
# 
# 
# ##This is the code used to generate figure 4 in the manuscript (overall temporal trends)
# trendplot <- ggplot() +
#   geom_point(data = df4, aes(x = Year, y = MeanRes),
#              shape = 21, color = "black", alpha = 0.4, size = 3) +
#   geom_line(data = pred_df, aes(x = Year, y = fit), color = "blue", linewidth = 1) +
#   geom_ribbon(data = pred_df, aes(x = Year, ymin = fit - 1.96*se.fit, ymax = fit + 1.96*se.fit),
#               fill = "gray70", alpha = 0.3) +
#   geom_boxplot(data = data2, aes(x = Year, y = MeanRes, group = Year),
#                outlier.shape = NA, fill = "skyblue", color = "black", alpha = 0.2) +
#   theme_Publication() +
#   labs(x = "Year", y = "Mean Resistance") + 
#   scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
#   scale_x_continuous(breaks = seq(min(df4$Year), max(df4$Year), by = 2)) +
#   theme(axis.text.x = element_text(angle = 60, hjust = 1))
# 
# trendplot
# 
# # ggsave(filename = "trendplot.tiff",
# #        plot = trendplot,
# #        device = "tiff",
# #        width = 10,
# #        height = 8,
# #        units = "in",
# #        dpi = 700,
# #        compression = "lzw")
# 
# #regional trends
# 
# 
# 
# ##we now generate plots for each regionand generate figure5
# 
# #east african trends
# 
# ea<- df4%>%
#   filter(region=="Eastern Africa")
# 
# eadata2<- data2%>%
#   filter(region=="Eastern Africa")
# 
# eamodel <- glm(MeanRes ~ Year, data = eadata2, family = quasibinomial, weights = no_isolate)
# 
# eanew_preds <- predict(eamodel, newdata = data.frame(Year = seq(min(eadata2$Year), max(eadata2$Year), by = 1)), 
#                      type = "response", se.fit = TRUE)
# 
# 
# eapred_df <- data.frame(
#   Year = seq(min(eadata2$Year), max(eadata2$Year), by = 1),
#   fit = eanew_preds$fit,
#   se.fit = eanew_preds$se.fit
# )
# 
# 
# eafit1<-predict(eamodel,type = "response", se.fit = T)
# print(eafit1)
# 
# eadf1<-as.data.frame(eafit1)
# 
# eadf2<- cbind(eadata2, eadf1)
# 
# eadf3<- eadf2%>%
#   mutate(lower = fit - 1.96 * se.fit)%>%
#   mutate(upper = fit+ 1.96 * se.fit)
# 
# eadf4<- aggregate(MeanRes~doi+Year+region, data = eadf3, FUN = mean)
# 
# 
# eaplot <- ggplot() +
#   geom_point(data = eadf4, aes(x = Year, y = MeanRes),
#              shape = 21, color = "black", alpha = 0.4, size = 3) +
#   geom_line(data = eapred_df, aes(x = Year, y = fit), color = "blue", linewidth = 1) +
#   geom_ribbon(data = eapred_df, aes(x = Year, ymin = fit - 1.96*se.fit, ymax = fit + 1.96*se.fit),
#               fill = "gray70", alpha = 0.3) +
#   geom_boxplot(data = eadata2, aes(x = Year, y = MeanRes, group = Year),
#                outlier.shape = NA, fill = "skyblue", color = "black", alpha = 0.2) +
#   theme_Publication() +
#   labs(x = "Year", y = "Mean Resistance", subtitle = "East Africa") + 
#   scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
#   scale_x_continuous(breaks = seq(min(eadf4$Year), max(eadf4$Year), by = 2)) +
#   theme(axis.text.x = element_text(angle = 60, hjust = 1))
# eaplot
# 
# 
# #North Africa plot
# 
# na<- df4%>%
#   filter(region=="Northern Africa")
# 
# nadata2<- data2%>%
#   filter(region=="Northern Africa")
# 
# namodel <- glm(MeanRes ~ Year, data = nadata2, family = quasibinomial, weights = no_isolate)
# 
# nanew_preds <- predict(namodel, newdata = data.frame(Year = seq(min(nadata2$Year), max(nadata2$Year), by = 1)), 
#                        type = "response", se.fit = TRUE)
# 
# 
# napred_df <- data.frame(
#   Year = seq(min(nadata2$Year), max(nadata2$Year), by = 1),
#   fit = nanew_preds$fit,
#   se.fit = nanew_preds$se.fit
# )
# 
# 
# 
# nafit1<-predict(namodel,type = "response", se.fit = T)
# 
# 
# nadf1<-as.data.frame(nafit1)
# 
# nadf2<- cbind(nadata2, nadf1)
# 
# nadf3<- nadf2%>%
#   mutate(lower = fit - 1.96 * se.fit)%>%
#   mutate(upper = fit+ 1.96 * se.fit)
# 
# nadf4<- aggregate(MeanRes~doi+Year+region, data = nadf3, FUN = mean)
# 
# 
# naplot <- ggplot() +
#   geom_point(data = nadf4, aes(x = Year, y = MeanRes),
#              shape = 21, color = "black", alpha = 0.4, size = 3) +
#   geom_line(data = napred_df, aes(x = Year, y = fit), color = "blue", linewidth = 1) +
#   geom_ribbon(data = napred_df, aes(x = Year, ymin = fit - 1.96*se.fit, ymax = fit + 1.96*se.fit),
#               fill = "gray70", alpha = 0.3) +
#   geom_boxplot(data = nadata2, aes(x = Year, y = MeanRes, group = Year),
#                outlier.shape = NA, fill = "skyblue", color = "black", alpha = 0.2) +
#   theme_Publication() +
#   labs(x = "Year", y = "Mean Resistance", subtitle = "North Africa") + 
#   scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
#   scale_x_continuous(breaks = seq(min(nadf4$Year), max(nadf4$Year), by = 2)) +
#   theme(axis.text.x = element_text(angle = 60, hjust = 1))
# naplot
# 
# #West africa
# 
# wa<- df4%>%
#   filter(region=="West Africa")
# 
# wadata2<- data2%>%
#   filter(region=="West Africa")
# 
# wamodel <- glm(MeanRes ~ Year, data = wadata2, family = quasibinomial, weights = no_isolate)
# 
# wanew_preds <- predict(wamodel, newdata = data.frame(Year = seq(min(wadata2$Year), max(wadata2$Year), by = 1)), 
#                        type = "response", se.fit = TRUE)
# 
# 
# wapred_df <- data.frame(
#   Year = seq(min(wadata2$Year), max(wadata2$Year), by = 1),
#   fit = wanew_preds$fit,
#   se.fit = wanew_preds$se.fit
# )
# 
# 
# 
# wafit1<-predict(wamodel,type = "response", se.fit = T)
# 
# 
# wadf1<-as.data.frame(wafit1)
# 
# wadf2<- cbind(wadata2, wadf1)
# 
# wadf3<- wadf2%>%
#   mutate(lower = fit - 1.96 * se.fit)%>%
#   mutate(upper = fit+ 1.96 * se.fit)
# 
# wadf4<- aggregate(MeanRes~doi+Year+region, data = wadf3, FUN = mean)
# 
# 
# waplot <- ggplot() +
#   geom_point(data = wadf4, aes(x = Year, y = MeanRes),
#              shape = 21, color = "black", alpha = 0.4, size = 3) +
#   geom_line(data = wapred_df, aes(x = Year, y = fit), color = "blue", linewidth = 1) +
#   geom_ribbon(data = wapred_df, aes(x = Year, ymin = fit - 1.96*se.fit, ymax = fit + 1.96*se.fit),
#               fill = "gray70", alpha = 0.3) +
#   geom_boxplot(data = wadata2, aes(x = Year, y = MeanRes, group = Year),
#                outlier.shape = NA, fill = "skyblue", color = "black", alpha = 0.2) +
#   theme_Publication() +
#   labs(x = "Year", y = "Mean Resistance", subtitle = "West Africa") + 
#   scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
#   scale_x_continuous(breaks = seq(min(wadf4$Year), max(wadf4$Year), by = 2)) +
#   theme(axis.text.x = element_text(angle = 60, hjust = 1))
# waplot
# 
# 
# 
# 
# 
# #south africa
# 
# sa<- df4%>%
#   filter(region=="Southern Africa")
# 
# sadata2<- data2%>%
#   filter(region=="Southern Africa")
# 
# samodel <- glm(MeanRes ~ Year, data = sadata2, family = quasibinomial, weights = no_isolate)
# 
# sanew_preds <- predict(samodel, newdata = data.frame(Year = seq(min(sadata2$Year), max(sadata2$Year), by = 1)), 
#                        type = "response", se.fit = TRUE)
# 
# 
# sapred_df <- data.frame(
#   Year = seq(min(sadata2$Year), max(sadata2$Year), by = 1),
#   fit = sanew_preds$fit,
#   se.fit = sanew_preds$se.fit
# )
# 
# 
# 
# safit1<-predict(samodel,type = "response", se.fit = T)
# 
# 
# sadf1<-as.data.frame(safit1)
# 
# sadf2<- cbind(sadata2, sadf1)
# 
# sadf3<- sadf2%>%
#   mutate(lower = fit - 1.96 * se.fit)%>%
#   mutate(upper = fit+ 1.96 * se.fit)
# 
# sadf4<- aggregate(MeanRes~doi+Year+region, data = sadf3, FUN = mean)
# 
# 
# saplot <- ggplot() +
#   geom_point(data = sadf4, aes(x = Year, y = MeanRes),
#              shape = 21, color = "black", alpha = 0.4, size = 3) +
#   geom_line(data = sapred_df, aes(x = Year, y = fit), color = "blue", linewidth = 1) +
#   geom_ribbon(data = sapred_df, aes(x = Year, ymin = fit - 1.96*se.fit, ymax = fit + 1.96*se.fit),
#               fill = "gray70", alpha = 0.3) +
#   geom_boxplot(data = sadata2, aes(x = Year, y = MeanRes, group = Year),
#                outlier.shape = NA, fill = "skyblue", color = "black", alpha = 0.2) +
#   theme_Publication() +
#   labs(x = "Year", y = "Mean Resistance",
#        subtitle = "Southern Africa") + 
#   scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
#   scale_x_continuous(breaks = seq(min(sadf4$Year), max(sadf4$Year), by = 2)) +
#   theme(axis.text.x = element_text(angle = 60, hjust = 1))
# saplot
# 
# 
# 
# cp<-  ggarrange(eaplot , waplot, naplot , saplot,
#                 labels = c("(a)", "(b)", "(c)", "(d)"),
#                 common.legend = TRUE)
# 
# 
# yleft <- textGrob("Percent resistance", rot = 90, gp = gpar(fontsize = 20))
# bottom <- textGrob("Antimicrobial", gp = gpar(fontsize = 20))
# 
# #tiff("cp.tiff", width=2300, height=2000, res=300)
# 
# dev.off()
# library(patchwork)
# cp.region_plot <- (ea.p | wa.p)/ (na.p | sa.p)+
#   plot_annotation(tag_levels = "a",
#                   tag_prefix = "(",
#                   tag_suffix = ")")&
#   theme(plot.tag = element_text(face = "bold"),
#         legend.position = "top")
# 
# cp.region_plot 
# 
# # ggsave(filename = "cp.tiff",
# #        plot = cp,
# #        device = "tiff",
# #        width = 10,
# #        height = 8,
# #        units = "in",
# #        dpi = 700,
# #        compression = "lzw")
# 
# # p1 <- ggboxplot(df4, x = "Year", y = "MeanRes", fill = "MeanRes",
# #                    xlab = "Year", ylab = "Percent resistance") +
# # geom_boxplot(fill = "skyblue2") +
# # geom_jitter(alpha = 0.4) +
# # rotate_x_text(60) +
# #  ggpubr::font("xlab", size = 17)+
# #  ggpubr::font("ylab", size = 17)+
# # ggpubr::font("xy.text", size = 17)+
# # ggpubr::font("legend.title",size = 17)+
# # ggpubr::font("legend.text", size = 17)
# # 
# # 
# # # 
# # p2<-p1 +
# #   geom_smooth(aes(as.numeric(Year),
# #           (MeanRes)),
# #          method = "glm",se = T, colour = "maroon",lty=2,
# #       method.args = list(family = "quasibinomial")) +
# #   theme_Publication() 
#     #geom_line(aes(y = lower), lty = 2) +
#     #geom_line(aes(y = upper), lty = 2) +
#   #facet_wrap(~region)
# 
# # p2<-p2 + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))+
# #   rotate_x_text(60)
# # 
# # 
# # p2
# 
# 
# 
# df4.1<-aggregate(MeanRes~Year+region+fit, data = df3, FUN = mean)
# 
# 
# ##Fisher exact test used to test for differences in resistance between years
# 
# ftest<-fisher.test(table(df4.1$Year, df4.1$MeanRes))
# ftest <- fisher.test(table(df4$Year, df4$MeanRes), workspace=2e8)
# 
# # Using simulate.p.value to handle large tables
# ftest <- fisher.test(table(df4.1$Year, df4.1$MeanRes), simulate.p.value=TRUE, B=1e5)
# 
# ftest$p.value
# 
# #East Africa
# ea<- df4.1%>%
#   filter(region=="Eastern Africa")
# ea_ftest<-fisher.test(table(ea$Year, ea$MeanRes))
# ea_ftest <- fisher.test(table(ea$Year, ea$MeanRes), simulate.p.value=TRUE, B=1e5)
# 
# ea_ftest$p.value
# 
# 
# #central africa
# ca<- df4.1%>%
#   filter(region=="Central Africa")
# ca_ftest<-fisher.test(table(ca$Year, ca$MeanRes))
# ca_ftest <- fisher.test(table(ca$Year, ca$MeanRes), simulate.p.value=TRUE, B=1e5)
# ca_ftest$p.value
# 
# #Northern Africa
# na<- df4.1%>%
#   filter(region=="Northern Africa")
# na_ftest<-fisher.test(table(na$Year, na$MeanRes), simulate.p.value=TRUE, B=1e5)
# na_ftest$p.value
# 
# #Southern Africa
# 
# sa<- df4.1%>%
#   
#   filter(region=="Southern Africa")
# 
# sa_ftest<-fisher.test(table(sa$Year, sa$MeanRes), simulate.p.value=TRUE, B=1e5)
# sa_ftest$p.value
# #West Africa
# 
# wa<- df4.1%>%
#   filter(region=="West Africa")
# 
# wa_ftest <- fisher.test(table(wa$Year, wa$MeanRes), simulate.p.value=TRUE, B=1e5)
# 
# wa_ftest$p.value
# 
# 
# # 
# # df_overall<-aggregate(MeanRes~Year+region, data = df4, FUN = mean)
# # 
# # ea2<- aggregate(fit~Year, data = df3, FUN = mean)
# # 
# # #Number of studies included in temporal trend analysis
# # df4.2<-aggregate((no_isolates_resistant/no_isolate)~doi+Year+region, data = AMRdf, FUN = mean) %>%
# #   na.omit()
# # 
# # ea2<- df4.2%>%
# #   filter(region=="Eastern Africa")
# # 
# # ea2%>%
# #   count(doi)
# # 
# # wa2<- df4.2%>%
# #   filter(region=="West Africa")
# # 
# # wa2%>%
# #   count(doi)
# # 
# # na2<- df4.2%>%
# #   filter(region=="Northern Africa")
# # 
# # na2%>%
# #   count(doi)
# # 
# # sa2<- df4.2%>%
# #   filter(region=="Southern Africa")
# # 
# # sa2%>%
# #   count(doi)
# # 
# # ca2<- df4.2%>%
# #   
# #   filter(region=="Central Africa")
# # 
# # ca2%>%
# #   count(doi)
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
