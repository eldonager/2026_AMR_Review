library (tidyverse)
library(janitor)
library(ggpubr)
library(patchwork)
library(grid)
library(gridExtra)





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





AMR_clean <- read.csv("AMR_clean.csv")

AMR_clean %>%
  count(no_isolates_susceptible)%>%
  View()

colnames(AMR_clean)


AMR_clean<- AMR_clean %>%
  dplyr::select(doi,country,region, iso_3, y_coordinate,x_coordinate,sampling_start_year,
                sampling_end_year, species, pathogen,sal_prevalence, no_sample,
                antimicrobial,antimicrobial_compound,antibiotic_class,
                who_classification, no_isolate,no_isolates_resistant, 
                no_isolates_intermediate,no_isolates_susceptible,mdr_percentage)                                      
                                  
                            
  AMR_clean<- AMR_clean %>%
    dplyr::mutate(no_isolates_susceptible=no_isolate-no_isolates_resistant)
  
  AMR_clean<- AMR_clean %>% 
    dplyr::select(-no_isolates_intermediate)
  
 
  ##Pooled prevalence for salmonella by antibiotic compounds- bargraphs
  #aggregate all resistant isolates and all NIsolates by speciesag, compound and region
  SalRes = aggregate(no_isolates_resistant ~ antimicrobial_compound + species + region, data = AMR_clean, FUN = sum)
  SalAll = aggregate(no_isolate ~ antimicrobial_compound + species + region, data = AMR_clean, FUN = sum)
  
  #and divide ResIso/NIsolates to return true mean
  SalMean = round((SalRes$no_isolates_resistant /SalAll$no_isolate)*100, digits = 0)
  SalMeandf = as.data.frame(cbind(SalRes$antimicrobial_compound, SalRes$species, SalRes$region, SalMean, SalAll$no_isolate))
  colnames(SalMeandf) = c("Compound","Species","region","Mean","NIsolates")
  SalMeandf$Mean = as.numeric(as.character(SalMeandf$Mean))
  SalMeandf$NIsolates = as.numeric(as.character(SalMeandf$NIsolates))
  
  #restrict analysis to bacteria-drug pairings where NIsolates > =10
  SalMeandf = SalMeandf[which(SalMeandf$NIsolates >= 10),]
  
  hist(SalMeandf$Mean, main = "Distribution of Mean Prevalence", xlab = "Mean Prevalence")
  
  
  # Compute the 95% CI of proportion where x = p_hat and y = n two tailed z = 1.96
  CI.function <- function(x,y) {
    x + c(-1.96,1.96)*sqrt(x*(1-x)/y)}
  
  #95% CI
  CIsal = as.data.frame(t(mapply(CI.function, SalMeandf$Mean/100,  SalMeandf$NIsolates)))
  SalMeandf  = cbind(SalMeandf, round(CIsal*100, digits = 0))
  colnames(SalMeandf) = c("Compound","Species","Region","Mean","NIsolates","CILow","CIHigh")
  SalMeandf$CILow[SalMeandf$CILow < 0] = 0
  SalMeandf$CIHigh[SalMeandf$CIHigh >100] = 100
  
  col_index <- 4  # Assuming you want to check the second column
  
  # Select rows where the specified column has a value other than zero
  SalMeandf <- SalMeandf[SalMeandf[, col_index] != 0, ]
  
  #East African Bar plots
  sal.ea<-SalMeandf %>%
    dplyr::filter(Compound %in% c('ERY', 'STR', 'TET', 'AMP', 'NIT', 'KAN', 
                                  'NAL', 'SOX', 'AMC' 
                                  #'CHL', 'COT', 'CIP',
                                  #'CFL', 'GEN', 'NEO', 'TMP'
                                  ), Region == "Eastern Africa")
  
  sal.ea$Compound <- reorder(sal.ea$Compound, -sal.ea$Mean)
  
  sal.ea.plot <- ggbarplot(sal.ea, "Compound", "Mean", 
                                 fill = "Species", position = position_dodge(0.7),
              subtitle = "Eastern Africa",# Cattle, n = 3,523, 
              # Chicken, n = 5,906, Pigs = n = 3,252", 
                                 xlab = FALSE, ylab = FALSE,
                                 legend.title = "Species",
                                 font.subtitle = c(12))+
    rotate_x_text(90)+
    geom_linerange(aes(group = Species, ymax = CIHigh, ymin = CILow),
                   position = position_dodge(width = 0.7))+
    font("xy.text", size = 19)+
    font("legend.title",size = 19)+
    font("legend.text", size = 19)+
    theme_Publication()+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  
  sum_by_speciesEA <- sal.ea %>%
    group_by(Species) %>%
    summarize(total_isolates = sum(NIsolates))
  
  
  sal.ea.plot <- ggpar(sal.ea.plot, ylim = c(0, 100))
  
  
  #West African Bar plots
  sal.wa<-SalMeandf %>%
    dplyr::filter(Compound %in% c('ERY', 'AMX', 'AMP', 'AMC', 'COT', 'CTX', 
                                  'SMX' #'TET', 'MIN', 'STR', 'NAL',
                                  #'GEN', 'CAZ', 'CIP', 'CHL'
                                  ), Region == "West Africa")
  
  
  sal.wa$Compound <- reorder(sal.wa$Compound, -sal.wa$Mean)
 
  
  
  
  
  species <- c("Cattle", "Chicken", "Environment", "Pigs")
  compounds <- c("AMP", "AMX", "ERY", "AMC", "COT", "CTX", "SMX")
  
  # Create an empty data frame
  data1 <- data.frame(Compound = character(),
                     Species = character(),
                     Region = character(),
                     Mean = numeric(),
                     NIsolates = numeric(),
                     CILow = numeric(),
                     CIHigh = numeric(),
                     stringsAsFactors = FALSE)
  
  # Populate the data frame with all combinations of species and compounds
  for (sp in species) {
    for (cp in compounds) {
      data1 <- rbind(data1, data.frame(Compound = cp, Species = sp, Region = "West Africa", Mean = 0, NIsolates = 0, CILow = 0, CIHigh = 0))
    }
  }
  
 
  
  
  sal.wa<-rbind(sal.wa, data1)
  
  
  
  
  
  
  sum_by_speciesWA <- sal.wa %>%
    group_by(Species) %>%
    summarize(total_isolates = sum(NIsolates)) 
  
  sal.wa.plot <- ggbarplot(sal.wa, "Compound", "Mean", 
                           fill = "Species", position = position_dodge(0.7),
  subtitle = "West Africa", #Cattle n = 1,643, Chicken n = 19,595, 
  # Environment n = 1,395, Pig n = 2,676, Sheep n = 32",
                           xlab = FALSE, ylab = FALSE,
                           legend.title = "Species",
                           font.subtitle = c(12))+
    rotate_x_text(90)+
    geom_linerange(aes(group = Species, ymax = CIHigh, ymin = CILow),
                   position = position_dodge(width = 0.7))+
    font("xy.text", size = 19)+
    font("legend.title",size = 19)+
    font("legend.text", size = 19)+
    theme_Publication()+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  
  sal.wa.plot <- ggpar(sal.wa.plot, ylim = c(0, 100))
  
  
  #North African Bar plots
  sal.na<-SalMeandf %>%
    dplyr::filter(Compound %in% c( 'AMX', 'AMP', 'AMC', 'COT', 'CTX', 
                                  'SMX', 'TET', 'TMP' 
                                  #'MIN', 'STR', 'NAL',
                                  #'GEN', 'CIP', 'CHL'
                                  ), Region == "Northern Africa")
  
  

  
  additional_data <- data.frame(
    Compound = c("SMX", "AMX", "CTX"),
    Species = c("Turkey", "Turkey",  "Cattle"),
    Mean = c(0.01, 0.01,  0.01))
  
  
  additional_data <- additional_data %>%
    mutate(Region = "Northern Africa",
           NIsolates = 78,
           CILow = 0,
           CIHigh = 0)
  
  # Combine the original data with the additional data
  sal.na <- rbind(sal.na, additional_data)
  
  sal.na%>%
    mutate(Mean1=Mean)
  
  sal.na$Compound <- reorder(sal.na$Compound, -sal.na$Mean)
  
  sum_by_speciesNA <- sal.na %>%
    group_by(Species) %>%
    summarize(total_isolates = sum(NIsolates)) 
  
  sal.na.plot <- ggbarplot(sal.na, "Compound", "Mean", 
                           fill = "Species", position = position_dodge(0.7),
      subtitle = "Northern Africa",# Cattle n = 1,438,
      # Chicken n = 7,154, Turkey n = 744",
                           xlab = FALSE, ylab = FALSE,
                           legend.title = "Species",
                           font.subtitle = c(12))+
    rotate_x_text(90)+
    geom_linerange(aes(group = Species, ymax = CIHigh, ymin = CILow),
                   position = position_dodge(width = 0.7))+
    font("xy.text", size = 19)+
    font("legend.title",size = 19)+
    font("legend.text", size = 19)+
    theme_Publication()+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  
  sal.na.plot 
  
  
  sal.na.plot <- ggpar(sal.na.plot, ylim = c(0, 100))
  
  #South African Bar plots
  sal.sa<-SalMeandf %>%
    dplyr::filter(Compound %in% c('SOX', 'ERY', 'AMX', 'TET', 'STR', 'CTX', 
                                  'CFX', 'AMP', 'AMC', 'COT' #'CHL', 'CFL',
                                  #'NOR', 'KAN', 'CIP', 'NAL'
                                  ), Region == "Southern Africa")
  
  sal.sa$Compound <- reorder(sal.sa$Compound, -sal.sa$Mean)
  
  
  
  # Define the species and compounds
  species <- c("Cattle", "Chicken", "Goats", "Pigs", "Sheep")
  compounds <- c("SOX", "ERY", "AMX", "TET", "STR", "CTX", "CFX", "AMP", "AMC", "COT")
  
  # Create an empty data frame
  data <- data.frame(Compound = character(),
                     Species = character(),
                     Region = character(),
                     Mean = numeric(),
                     NIsolates = numeric(),
                     CILow = numeric(),
                     CIHigh = numeric(),
                     stringsAsFactors = FALSE)
  
  # Populate the data frame with all combinations of species and compounds
  for (sp in species) {
    for (cp in compounds) {
      data <- rbind(data, data.frame(Compound = cp, Species = sp, Region = "Southern Africa", Mean = 0, NIsolates = 0, CILow = 0, CIHigh = 0))
    }
  }
  
  # Print the data frame
  print(data)
  
  sal.sa<-rbind(sal.sa, data) 
  
  
  
  
  
  sum_by_speciesSA <- sal.sa %>%
    group_by(Species) %>%
    summarize(total_isolates = sum(NIsolates))
  
  
  sal.sa.plot <- ggbarplot(sal.sa, "Compound", "Mean", 
                           fill = "Species", position = position_dodge(0.7),
  subtitle = "Southern Africa",
  #Cattle n = 2,239, Chicken n = 7,639,
  # Goats n = 612, Pigs n = 1,068, Sheep n = 264",
  #                         xlab = FALSE, ylab = FALSE,
                           legend.title = "Species",
                           legend.position = "top",
                           font.subtitle = c(12))+
    rotate_x_text(90)+
    geom_linerange(aes(group = Species, ymax = CIHigh, ymin = CILow),
                   position = position_dodge(width = 0.7))+
    font("xy.text", size = 19)+
    font("legend.title",size = 19)+
    font("legend.text", size = 19)+
    theme_Publication()+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  
  sal.sa.plot <- ggpar(sal.sa.plot, ylim = c(0, 100))
  
  
p1<-  ggarrange(sal.ea.plot , sal.wa.plot, sal.na.plot , sal.sa.plot,
            labels = c("(a)", "(b)", "(c)", "(d)"))
  
yleft <- textGrob("Percent resistance", rot = 90, gp = gpar(fontsize = 20))
bottom <- textGrob("Antimicrobial", gp = gpar(fontsize = 20))


region_plot <- (sal.ea.plot | sal.wa.plot)/ (sal.na.plot | sal.sa.plot)+
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")")&
  theme(plot.tag = element_text(face = "bold"),
        legend.position = "top")


plot1 <- patchwork::patchworkGrob(region_plot)
#tiff("p1.tiff", width=3500, height=2000, res=300)
p1 <- grid.arrange(plot1, ncol = 1, nrow = 1,left = yleft, bottom = bottom)
p1


ggsave(filename = "p1.tiff",
       plot = p1,
       device = "tiff",
       width = 10,
       height = 8,
       units = "in",
       dpi = 700,
       compression = "lzw")



#dev.off()

  