###########
#Descriptive
###########
 library(dplyr)
library(ggplot2)
 library(cowplot)
library(ggpubr)
library(kimisc)
 library(AICcmodavg)
library(Metrics)
library(detectseparation)
library(brglm2)
library(car)
library(boot)#diagnostic plots
library(lme4)##glm
library(Matrix)
library(glmmTMB)
library(sjPlot) #for plotting lmer and glmer mods
library(sjmisc) 
library(effects)
 library(sjstats)
library(tidyverse)



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

#Mapping AST Methods
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
  dplyr::mutate(MeanRes=(no_isolates_resistant)/(no_isolate))
  

#data1%>%
 # count(MeanRes)%>%
  #View()


data1 <- data1[!is.na(data1$Year), ]

data1%>%
  na.omit()
na.omit(data1)

data1 <- data1[data1$MeanRes != 0, ]
# data1 is my data frame


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

#data2<- aggregate(MeanRes~Year+doi+antibiotic_class+country+region+species+no_isolate, data=data1, FUN=mean)




model1 <- glm(cbind(no_isolates_resistant/no_isolate) ~ Year + country + region + species +antimicrobial+strain+
                ast_method+ (1 | doi),
                family = quasibinomial, weights=no_isolate,data = data1)

summary(model1)

model2<-glm(cbind(no_isolates_resistant/no_isolate) ~ Year + country + region + species +antimicrobial+strain+
              (1 | doi),
            family=quasibinomial, weights=no_isolate, data = data1)
summary(model2)


model3<-glm(cbind(no_isolates_resistant/no_isolate) ~ Year + country + region + species +antimicrobial+ (1 | doi),
            family=quasibinomial, weights=no_isolate,
            data = data1)
summary(model3)

model4<-glm(cbind(no_isolates_resistant/no_isolate) ~ Year + country + region + species + (1 | doi),
            family=binomial, weights=no_isolate, data = data1)
summary(model4)

model5<-glm(cbind(no_isolates_resistant/no_isolate) ~ Year + country + region + (1 | doi),
            family=binomial, weights=no_isolate, data = data1)
summary(model5)


model6<-glm(cbind(no_isolates_resistant/no_isolate) ~ Year + country + (1 | doi),
            
            family=binomial, weights=no_isolate, data = data1)
summary(model6)

model7<- glm(cbind(no_isolates_resistant/no_isolate) ~ Year +  (1 | doi),
             
             family=quasibinomial, weights=no_isolate, data = data1)

model7

#model7<- glm(cbind(no_isolates_resistant/no_isolate) ~ Year,
             
             #family=binomial, weights=no_isolate, data = data1)

##Diagnostics using DHArma
glm.diag.plots(model7, glmdiag = glm.diag(model7), subset = NULL,
               iden = FALSE, labels = NULL, ret = FALSE)
diag_plots<-glm.diag.plots(model3, glmdiag = glm.diag(model3), subset = NULL,
                           iden = FALSE, labels = NULL, ret = FALSE)
diag_plots



#After trying inspecting models 1-7, model 7 is the best model

data1$Year <- as.numeric(as.character(data1$Year))

new_preds <- predict(model7, newdata = data.frame(Year = seq(min(data1$Year), max(data1$Year), by = 1)), 
                     type = "response", se.fit = TRUE)


pred_df <- data.frame(
  Year = seq(min(data1$Year), max(data1$Year), by = 1),
  fit = new_preds$fit,
  se.fit = new_preds$se.fit
)


trendplot <- ggplot() +
  geom_point(data = data1, aes(x = Year, y = MeanRes),
             shape = 21, color = "black", alpha = 0.5, size = 3) +
  geom_line(data = pred_df, aes(x = Year, y = fit), color = "red", linewidth = 1) +
  geom_ribbon(data = pred_df, aes(x = Year, ymin = fit - 1.96*se.fit, ymax = fit + 1.96*se.fit),
              fill = "blue", alpha = 0.3) +
  geom_boxplot(data = data2, aes(x = Year, y = MeanRes, group = Year),
               outlier.shape = NA, fill = "skyblue2", color = "black", alpha = 0.2) +
  theme_Publication() 






data2<- aggregate(MeanRes~Year+doi+antibiotic_class+country+region+species+no_isolate, data=data1, FUN=mean)





vars=names(data1)[14]#mean resistance

fits = lapply(vars, function(x) {
  glm(substitute(i ~ Year, list(i = as.name(x))), 
      data = data2, 
      family = quasibinomial, 
      no_isolate)
})

lapply(fits, summary)

#Extracting Coefficients and P-values:

#extract coefficient and p value
sapply(fits, function(f) summary(f)$coefficients[,c(1,4)])[c(2,4),]



#Predicting the model
mypreds = lapply(fits, function(f) predict(f,data = data2, type = "response", se.fit = T))



##trend plot

trendplot =  ggplot(data = data2) +
  ggplot2::geom_point(mapping = aes(x=Year, y = MeanRes),
              shape=21, fill="skyblue", col="black", alpha=0.5, size=3, stroke=0.4,
              position = position_jitter(width = 0.2, height = 0))+
  geom_boxplot(mapping = aes(x=Year, y = MeanRes, group = Year),
                outlier.shape = NA, fill="grey 70", col="black", alpha=0.2)+
  geom_ribbon(mapping = aes(x=Year, y=MeanRes, 
                            ymin=mypreds[[2]]$fits- 1.96*mypreds[[2]]$se.fit,
                            ymax=mypreds[[2]]$fits+ 1.96*mypreds[[2]]$se.fit),
              fill="blue", alpha=0.3) +
  geom_line(mapping = aes(x=Year, y=mypreds[[1]]$fits), col="red") 







fit1<-predict(model4,type = "response", se.fit = T)
print(fit1)

df1<-as.data.frame(fit1)

df2<- cbind(data1,df1)

df3<- df2%>%
  mutate(lower = fit - 1.96 * se.fit)%>%
  mutate(upper = fit + 1.96 * se.fit)

df4<- aggregate(MeanRes~doi+Year+region, data = df3, FUN = mean)


#df4<- aggregate(MeanRes~Year+country+region+species+antibiotic_class+fit+lower+upper, data = df3, FUN = mean)

# Convert 'Year' to a factor (discrete variable)

# Ensure Year is treated as a factor if it represents a categorical variable
#
# 

# df4$upper <- ifelse(df4$upper < 0, 0, ifelse(df4$upper > 1, 1, df4$upper))
# df4$lower <- ifelse(df4$lower < 0, 0, ifelse(df4$lower > 1, 1, df4$lower))


#tiff("p2.tiff", width=2300, height=2000, res=300)

p1 <- ggboxplot(df4, x = "Year", y = "MeanRes", fill = "MeanRes",
                   xlab = "Year", ylab = "Percent resistance") +
geom_boxplot(fill = "skyblue2") +
geom_jitter(alpha = 0.4) +
rotate_x_text(60) +
 ggpubr::font("xlab", size = 17)+
 ggpubr::font("ylab", size = 17)+
ggpubr::font("xy.text", size = 17)+
ggpubr::font("legend.title",size = 17)+
ggpubr::font("legend.text", size = 17)


# 
p2<-p1 +
  geom_smooth(aes(as.numeric(Year),
          (MeanRes)),
         method = "glm",se = T, colour = "maroon",lty=2,
      method.args = list(family = "quasibinomial")) +
  theme_Publication() 
    #geom_line(aes(y = lower), lty = 2) +
    #geom_line(aes(y = upper), lty = 2) +
  #facet_wrap(~region)

p2<-p2 + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))+
  rotate_x_text(60)


p2


#tiff("p2.tiff", width=2300, height=2000, res=300)
#dev.off()



#By region
# 
#subsettin to regions

ea<- df4%>%
  filter(region=="Eastern Africa")
na<- df4%>%
  filter(region=="Northern Africa")
sa<- df4%>%
  filter(region=="Southern Africa")
wa<- df4%>%
  filter(region=="West Africa")


#EA
# Load necessary libraries
library(ggplot2)
library(ggpubr)

# Convert Year to numeric in the data frame
ea$Year <- as.numeric(as.character(ea$Year))

# Create the plot using ggplot directly
ea.p <- ggplot(ea, aes(x = Year, y = MeanRes)) +
  geom_boxplot(aes(group = Year), fill = "skyblue2") +
  geom_jitter(alpha = 0.4) +
  geom_smooth(method = "glm", se = TRUE, colour = "maroon", lty = 2,
              method.args = list(family = "quasibinomial")) +
  labs(x = "Year", y = "Percent resistance", subtitle = "East Africa") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  ggpubr::font("xlab", size = 17) +
  ggpubr::font("ylab", size = 17) +
  ggpubr::font("xy.text", size = 17) +
  ggpubr::font("legend.title", size = 17) +
  ggpubr::font("legend.text", size = 17) +
  scale_x_continuous(breaks = seq(2004, 2022, by = 2), limits = c(2004, 2022))+
  theme_Publication()

#adjust y axis scale to 100

ea.p<-ea.p + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))+
  rotate_x_text(60)






#NA
# Load necessary libraries


# Convert Year to numeric in the data frame
na$Year <- as.numeric(as.character(na$Year))

# Create the plot using ggplot directly
na.p <- ggplot(na, aes(x = Year, y = MeanRes)) +
  geom_boxplot(aes(group = Year), fill = "skyblue2") +
  geom_jitter(alpha = 0.4) +
  geom_smooth(method = "glm", se = TRUE, colour = "maroon", lty = 2,
              method.args = list(family = "quasibinomial")) +
  labs(x = "Year", y = "Percent resistance", subtitle = "North Africa") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  ggpubr::font("xlab", size = 17) +
  ggpubr::font("ylab", size = 17) +
  ggpubr::font("xy.text", size = 17) +
  ggpubr::font("legend.title", size = 17) +
  ggpubr::font("legend.text", size = 17) +
  scale_x_continuous(breaks = seq(2005, 2020, by = 2), limits = c(2005, 2020))+
  theme_Publication()


#adjust y axis scale to 100

na.p<-na.p + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))+
  rotate_x_text(60)






#SA
# Convert Year to numeric in the data frame
sa$Year <- as.numeric(as.character(sa$Year))

# Create the plot using ggplot directly
sa.p <- ggplot(sa, aes(x = Year, y = MeanRes)) +
  geom_boxplot(aes(group = Year), fill = "skyblue2") +
  geom_jitter(alpha = 0.4) +
  geom_smooth(method = "glm", se = TRUE, colour = "maroon", lty = 2,
              method.args = list(family = "quasibinomial")) +
  labs(x = "Year", y = "Percent resistance", subtitle = "Southern Africa") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  ggpubr::font("xlab", size = 17) +
  ggpubr::font("ylab", size = 17) +
  ggpubr::font("xy.text", size = 17) +
  ggpubr::font("legend.title", size = 17) +
  ggpubr::font("legend.text", size = 17) +
  scale_x_continuous(breaks = seq(2005, 2017, by = 2), limits = c(2005, 2017))+
  theme_Publication()

#adjust y axis scale to 100

sa.p<-sa.p + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))+
  rotate_x_text(60)


#WA

# Convert Year to numeric in the data frame
wa$Year <- as.numeric(as.character(wa$Year))

# Create the plot using ggplot directly
wa.p <- ggplot(wa, aes(x = Year, y = MeanRes)) +
  geom_boxplot(aes(group = Year), fill = "skyblue2") +
  geom_jitter(alpha = 0.4) +
  geom_smooth(method = "glm", se = TRUE, colour = "maroon", lty = 2,
              method.args = list(family = "quasibinomial")) +
  labs(x = "Year", y = "Percent resistance", subtitle = "West Africa") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  ggpubr::font("xlab", size = 17) +
  ggpubr::font("ylab", size = 17) +
  ggpubr::font("xy.text", size = 17) +
  ggpubr::font("legend.title", size = 17) +
  ggpubr::font("legend.text", size = 17) +
  scale_x_continuous(breaks = seq(2002, 2020, by = 2), limits = c(2002, 2020))+
  theme_Publication()


#adjust y axis scale to 100

wa.p<-wa.p + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))+
  rotate_x_text(60)




cp<-  ggarrange(ea.p , wa.p, na.p , sa.p,
                labels = c("(a)", "(b)", "(c)", "(d)"),
                common.legend = TRUE)


yleft <- textGrob("Percent resistance", rot = 90, gp = gpar(fontsize = 20))
bottom <- textGrob("Antimicrobial", gp = gpar(fontsize = 20))

#tiff("cp.region_plot.tiff", width=2300, height=2000, res=300)
library(patchwork)
cp.region_plot <- (ea.p | wa.p)/ (na.p | sa.p)+
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")")&
  theme(plot.tag = element_text(face = "bold"),
        legend.position = "top")

cp.region_plot 
# plot.comb <- patchwork::patchworkGrob(cp.region_plot )
# 
# f.p<- grid.arrange(plot.comb, ncol = 1, nrow = 1,left = yleft, bottom = bottom)
#   
# f.p


#dev.off()
#tiff("cp.region_plot.tiff", width=2300, height=2000, res=300)
#dev.off()




# Assuming 'model4' is your GLM object
# Create effect plot

# library(jtools)
# effect_plot(model4, pred = "Year", interval = TRUE)
# 
# 
# 
#library(broom)
#fit1<-tidy(model4, conf.int = TRUE) %>%
#   mutate(model = "Model general")
# 

#   install.packages("ggeffects")
#   library (ggeffect)
# y_hat <- ggeffects::ggpredict(model4, terms = c("Year"))
# 
# 
# Trendplot =  ggplot(data = df4) +
#   ggplot2::geom_point(mapping = aes(x=Year, y = MeanRes),
#              shape=21, fill="skyblue", col="black", alpha=0.5, size=3, stroke=0.4,
#              position = position_jitter(width = 0.2, height = 0))+
#  geom_boxplot(mapping = aes(x=Year, y = MeanRes, group = Year),
#                outlier.shape = NA, fill="grey 70", col="black", alpha=0.2)+ 
#   geom_smooth(method = "glm", se = F, colour = "brown") +
#   geom_line(aes(y = lower), lty = 2) +
#   geom_line(aes(y = upper), lty = 2) 
#   
# 
# data10<-aggregate(MeanRes~Year, data=data1,FUN=mean)
# 
# ggplot(data = df4) +
#   ggplot2::geom_point(mapping = aes(x = Year, y = MeanRes),
#                       shape = 21, fill = "skyblue", col = "black", alpha = 0.5, size = 3, stroke = 0.4,
#                       position = position_jitter(width = 0.2, height = 0)) +
#   geom_boxplot(mapping = aes(x = Year, y = MeanRes, group = Year),
#                outlier.shape = NA, fill = "grey70", col = "black", alpha = 0.2) + 
#   
#   
#   
# 
#   geom_ribbon(data = df4, aes(x=Year, y = lower, ymax = upper),
#              fill = "blue", alpha = 0.3) +
#   geom_line(data = df4, aes(y = fit))
# 
# 
# 
# 
# ggplot(df4, aes(x = Year, y = MeanRes)) +
#   geom_boxplot(position = position_jitter(width = 0.05, height = 0.05)) +
#   geom_ribbon(data = df4, aes(x=Year, y = lower, ymax = upper),
#               fill = "blue", alpha = 0.3) +
#   geom_line(data = df4, aes(y = fit)) 
#   
#   
# geom_ribbon(data= x_urch, aes(x=c.urchinden, ymin=lower, ymax=upper), alpha= 0.3, fill="blue")
#   
#   
# geom_abline(intercept = y_hat[[1, 2]], slope = y_hat[[19, 2]], 
#             color = "blue")+
#   
#   geom_ribbon(data= y_hat, aes(ymin = conf.low, ymax = conf.high), alpha = 0.3)
#   
#   
#   geom_abline(intercept = y_hat1[[1, 4]], slope = y_hat1[[19, 4]], 
#               color = "green")+#lower
#   geom_abline(intercept = y_hat1[[1, 5]], slope = y_hat1[[19, 5]], 
#               color = "yellow")#upper
#   
#   
# y_hat1<- as.data.frame(y_hat)%>%
#   dplyr::mutate(lower=predicted-1.96*std.error)%>%
#   dplyr::mutate(upper=predicted+1.96*std.error)
#   
# 
# ggplot(data = y_hat, aes(x = x, y = predicted)) +
#   # plot the fitted line
#   geom_line(color = "red")
# 
# 
# 
# 
# 
# plot(model3)
# 
# fit<-predict(model3, interval = "cofidence", level=.95)
# 
# new_data<- cbind(data1, fit)
# 
# 
# 
# pd = data.frame(data)
# 
# 
# 
# 
# # use type = "response" for probability-scale predictions    
# preds = predict(model8, newdata = pd, type = "response", se.fit = TRUE)
# pd$fit = preds$fit
# pd$se = preds$se.fit
# 
# glm_model<-glm(am~.,data = mtcars)  
# glm_model    
# 
# library(sjPlot)    
# plot_model(model3, vline.color = "red")    
# plot_model(glm_model, show.values = TRUE, value.offset = .3)   
# 
# sjPlot::plot_model(model3, type = "eff", terms = "Year")
# 
# 
# jtools::effect_plot(model3, pred = Year, interval = TRUE, y.label = "% Resistance")
# 
# 
# 
# require(ggiraph)
# require(ggiraphExtra)
# require(plyr)
# ggPredict(model3,se=TRUE)
# 
# jtools::effect_plot(model3, pred = Year, interval = TRUE, plot.trend= TRUE, 
#             jitter = 0.05)
# 

#overal trends


df4.1<-aggregate(MeanRes~Year+region+fit, data = df3, FUN = mean)
  
  


#Fisher's exact test

ftest<-fisher.test(table(df4.1$Year, df4.1$MeanRes))
ftest <- fisher.test(table(df4$Year, df4$MeanRes), workspace=2e8)

# Using simulate.p.value to handle large tables
ftest <- fisher.test(table(df4.1$Year, df4.1$MeanRes), simulate.p.value=TRUE, B=1e5)

ftest$p.value

#East Africa
ea<- df4.1%>%
  filter(region=="Eastern Africa")
ea_ftest<-fisher.test(table(ea$Year, ea$MeanRes))
ea_ftest <- fisher.test(table(ea$Year, ea$MeanRes), simulate.p.value=TRUE, B=1e5)

ea_ftest$p.value


#central africa
ca<- df4.1%>%
  filter(region=="Central Africa")
ca_ftest<-fisher.test(table(ca$Year, ca$MeanRes))
ca_ftest <- fisher.test(table(ca$Year, ca$MeanRes), simulate.p.value=TRUE, B=1e5)
ca_ftest$p.value

#Northern Africa
na<- df4.1%>%
  filter(region=="Northern Africa")
na_ftest<-fisher.test(table(na$Year, na$MeanRes), simulate.p.value=TRUE, B=1e5)
na_ftest$p.value

#Southern Africa

sa<- df4.1%>%
  
  filter(region=="Southern Africa")

sa_ftest<-fisher.test(table(sa$Year, sa$MeanRes), simulate.p.value=TRUE, B=1e5)
sa_ftest$p.value
#West Africa

wa<- df4.1%>%
  filter(region=="West Africa")

wa_ftest <- fisher.test(table(wa$Year, wa$MeanRes), simulate.p.value=TRUE, B=1e5)

wa_ftest$p.value



df_overall<-aggregate(MeanRes~Year+region, data = df4, FUN = mean)

ea2<- aggregate(fit~Year, data = df3, FUN = mean)

#Number of studies included in temporal trend analysis
df4.2<-aggregate((no_isolates_resistant/no_isolate)~doi+Year+region, data = AMRdf, FUN = mean) %>%
  na.omit()

ea2<- df4.2%>%
  filter(region=="Eastern Africa")

ea2%>%
  count(doi)

wa2<- df4.2%>%
  filter(region=="West Africa")

wa2%>%
  count(doi)

na2<- df4.2%>%
  filter(region=="Northern Africa")

na2%>%
  count(doi)

sa2<- df4.2%>%
  filter(region=="Southern Africa")

sa2%>%
  count(doi)

ca2<- df4.2%>%
  
  filter(region=="Central Africa")

ca2%>%
  count(doi)










