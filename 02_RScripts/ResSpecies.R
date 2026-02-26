
library (tidyverse)
library(janitor)
library(ggpubr)
library(wesanderson)
library(RColorBrewer)
library("ggsci")
library(patchwork)


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
            axis.title = element_text(face = "bold", size = rel(1)),
            axis.title.y = element_text(angle=90, vjust = 2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "",
            legend.direction = "horizontal",
            legend.key.size = unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_blank(),  # Remove legend title
            legend.text = element_blank(),   # Optionally, remove legend text
            plot.margin = unit(c(10,5,5,5), "mm"),
            strip.background = element_rect(colour="#f0f0f0", fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
}

scale_fill_Publication <- function(...) {
  library(scales)
  discrete_scale("fill", "Publication", manual_pal(values = c("#386cb0", "#fdb462", "#7fc97f", "#ef3b2c", "#662506", "#a6cee3", "#fb9a99", "#984ea3", "#ffff33")), ...)
}

scale_colour_Publication <- function(...) {
  library(scales)
  discrete_scale("colour", "Publication", manual_pal(values = c("#386cb0", "#fdb462", "#7fc97f", "#ef3b2c", "#662506", "#a6cee3", "#fb9a99", "#984ea3", "#ffff33")), ...)
}




AMR_clean <- read.csv("AMR_clean.csv")



colnames(AMR_clean)


AMR_clean<- AMR_clean %>%
  dplyr::select(doi,country,region, iso_3, y_coordinate,x_coordinate,sampling_start_year,
                sampling_end_year, species, pathogen,sal_prevalence,
                antimicrobial_compound,antibiotic_class,
                who_classification, no_isolate,no_isolates_resistant, 
                no_isolates_intermediate,no_isolates_susceptible,mdr_percentage)                                      


AMR_clean<- AMR_clean %>%
  dplyr::mutate(no_isolates_susceptible=no_isolate-no_isolates_resistant)


#marge all cephalosporin groups to one
AMR_clean$antibiotic_class[AMR_clean$antibiotic_class=="3rd generation Cephalosporin"] <- "Cephalosporins"
AMR_clean$antibiotic_class[AMR_clean$antibiotic_class=="1st generation Cephalosporin"] <- "Cephalosporins"
AMR_clean$antibiotic_class[AMR_clean$antibiotic_class=="2nd generation Cephalosporin"] <- "Cephalosporins"
AMR_clean$antibiotic_class[AMR_clean$antibiotic_class=="4th generation Cephalosporin"] <- "Cephalosporins"
#Merge amphenicol groups to flouroquinolones


AMR_clean$antibiotic_class[AMR_clean$antibiotic_class=="Amphenicols"] <- "Fluoroquinolones"

AMR_clean<- AMR_clean %>% 
  dplyr::select(-no_isolates_intermediate)

unique(AMR_clean$antibiotic_class)
AMR_tet<- AMR_clean %>%
  filter(antibiotic_class=="Tetracycline")

unique(AMR_tet$antibiotic_class)
AMR_pen<- AMR_clean %>%
  filter(antibiotic_class=="Penicillin " )
unique(AMR_pen$antibiotic_class)
AMR_ceph<- AMR_clean %>%
  filter(antibiotic_class=="Cephalosporins")

unique(AMR_ceph$antibiotic_class)
AMR_amino<- AMR_clean %>%
  filter(antibiotic_class=="Aminoglycosides")

unique(AMR_amino$antibiotic_class)
AMR_floro<- AMR_clean %>%
  filter(antibiotic_class=="Fluoroquinolones")

unique(AMR_floro$antibiotic_class)
AMR_sulf<- AMR_clean %>%
  filter(antibiotic_class=="Sulfonamides")

unique(AMR_sulf$antibiotic_class)



#Significant difference between median resistance levels for tetracycline by antibiotic class

amr_sig <- aggregate((no_isolates_resistant/no_isolate)~antibiotic_class+species, data = AMR_clean, FUN = median)


amr_sig<- amr_sig %>%
  filter(antibiotic_class %in% c("Tetracycline", "Penicillin " , "Cephalosporins",
                                 "Aminoglycosides", "Fluoroquinolones", "Sulfonamides"))


# t test comparing penicillin and tetracycline resistance levels

# Rename the column using base R
colnames(amr_sig)[colnames(amr_sig) == "(no_isolates_resistant/no_isolate)"] <- "median resistance prevalence"

# t test comparing penicillin and tetracycline resistance levels using median resistance prevalence

t.test(amr_sig$`median resistance prevalence`[amr_sig$antibiotic_class == "Penicillin "], 
       amr_sig$`median resistance prevalence`[amr_sig$antibiotic_class == "Tetracycline"], 
       alternative = "two.sided", paired = FALSE, var.equal = FALSE, conf.level = 0.95)
#t test comparing tetracyclin and cephalosporin resistance levels using median resistance prevalence

t.test(amr_sig$`median resistance prevalence`[amr_sig$antibiotic_class == "Tetracycline"], 
       amr_sig$`median resistance prevalence`[amr_sig$antibiotic_class == "Cephalosporins"], 
       alternative = "two.sided", paired = FALSE, var.equal = FALSE, conf.level = 0.95)

#t test comparing tetracyclin and aminoglycosides resistance levels using median resistance prevalence

t.test(amr_sig$`median resistance prevalence`[amr_sig$antibiotic_class == "Tetracycline"], 
       amr_sig$`median resistance prevalence`[amr_sig$antibiotic_class == "Aminoglycosides"], 
       alternative = "two.sided", paired = FALSE, var.equal = FALSE, conf.level = 0.95)

#t test comparing tetracyclin and fluoroquinolones resistance levels using median resistance prevalence

t.test(amr_sig$`median resistance prevalence`[amr_sig$antibiotic_class == "Tetracycline"], 
       amr_sig$`median resistance prevalence`[amr_sig$antibiotic_class == "Fluoroquinolones"], 
       alternative = "two.sided", paired = FALSE, var.equal = FALSE, conf.level = 0.95)

#t test comparing tetracyclin and sulfonamides resistance levels using median resistance prevalence

t.test(amr_sig$`median resistance prevalence`[amr_sig$antibiotic_class == "Tetracycline"], 
       amr_sig$`median resistance prevalence`[amr_sig$antibiotic_class == "Sulfonamides"], 
       alternative = "two.sided", paired = FALSE, var.equal = FALSE, conf.level = 0.95)



#AMR resistance levels for tetracycline by species

##Pooled prevalence for salmonella by antibiotic compounds- bargraphs
#aggregate all resistant isolates and all NIsolates by speciesag, compound and region
tetSalRes = aggregate(no_isolates_resistant ~  antimicrobial_compound+species+country , data = AMR_tet, FUN = sum)
tetSalAll = aggregate(no_isolate ~  antimicrobial_compound+species+country , data = AMR_tet, FUN = sum)

#and divide ResIso/NIsolates to return true mean
tetSalMean = round((tetSalRes$no_isolates_resistant /tetSalAll$no_isolate)*100, digits = 0)
tetSalMeandf = as.data.frame(cbind(tetSalRes$antimicrobial_compound,  tetSalRes$species, tetSalMean, tetSalAll$no_isolate, tetSalRes$country))
colnames(tetSalMeandf) = c("antimicrobial_compound", "Species","Mean","NIsolates", "country")
tetSalMeandf$Mean = as.numeric(as.character(tetSalMeandf$Mean))
tetSalMeandf$NIsolates = as.numeric(as.character(tetSalMeandf$NIsolates))


tetSalMeandf$Mean <- as.numeric(as.character(tetSalMeandf$Mean))


tetdf1 <- tetSalMeandf %>%
  group_by(Species) %>%
  select(
         "country",
         "antimicrobial_compound",
         "Mean",
         "NIsolates")


tetdf2 <-aggregate(Mean~Species+antimicrobial_compound, data = tetdf1, FUN=mean)


#changing percentage resistance into numeric class and rounding off  
tetdf2$Mean <- as.numeric(as.character(tetdf2$Mean))

tetdf2$Mean <- round(tetdf2$Mean, digits = 0)

tetdf2<- tetdf2[which(tetdf2$Mean>= 1),]

#tiff("cp.species_plot.tiff", width=2300, height=2000, res=300)
##species barplot
tetspecies_plot<-ggplot(tetdf2, aes(x = Species, y = Mean, fill = Species))+ 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA
  ) +
  
  scale_fill_brewer(palette = "Paired") +
  theme(legend.position = "top")+
  coord_cartesian(xlim = c(1.2, NA), clip = "off")

#theme_pubr
tetpp1<- tetspecies_plot+
  theme_pubr(
    base_size = 17,
    base_family = "",
    border = FALSE,
    margin = TRUE,
    legend = c("top"),
    x.text.angle = 0
  )

#labelling the plot
tet<- tetpp1 + 
  labs(
    x = "Host", 
    y = "Percent resistance", 
    fill = "Host", 
    color = "Host", 
    subtitle = "Tetracycline"
  ) +
  scale_color_discrete(name = "Host") +
  scale_x_discrete(
    name = "Host", 
    limits = c("Cattle", "Chicken", "Goats", "Pigs", "Sheep", "Turkey", "Environment")
  ) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "none"
  ) +
  theme_Publication()

#call the plot
tet<-tet+ scale_y_continuous(limits = c(0, 100))+
  rotate_x_text(60)

tet

#AMR resistance levels for penicillin by species

##Pooled prevalence for salmonella by antibiotic compounds- bargraphs
#aggregate all resistant isolates and all NIsolates by speciesag, compound and region
penSalRes = aggregate(no_isolates_resistant ~  antimicrobial_compound+species+country , data = AMR_pen, FUN = sum)
penSalAll = aggregate(no_isolate ~  antimicrobial_compound+species+country , data = AMR_pen, FUN = sum)
  
#and divide ResIso/NIsolates to return true mean
penSalMean = round((penSalRes$no_isolates_resistant /penSalAll$no_isolate)*100, digits = 0)
penSalMeandf = as.data.frame(cbind(penSalRes$antimicrobial_compound,  penSalRes$species, penSalMean, penSalAll$no_isolate, penSalRes$country))
colnames(penSalMeandf) = c("antimicrobial_compound", "Species","Mean","NIsolates", "country")
penSalMeandf$Mean = as.numeric(as.character(penSalMeandf$Mean))
penSalMeandf$NIsolates = as.numeric(as.character(penSalMeandf$NIsolates))

penSalMeandf$Mean <- as.numeric(as.character(penSalMeandf$Mean))


#group by species
pendf1 <- penSalMeandf %>%
  group_by(Species) %>%
  select(
    "country",
    "antimicrobial_compound",
    "Mean",
    "NIsolates")


##Aggregate mean resistance levels for penicillin by species
pendf2 <-aggregate(Mean~Species+antimicrobial_compound, data = pendf1, FUN=mean)

#changing percentage resistance into numeric class and rounding off 

pendf2$Mean <- as.numeric(as.character(pendf2$Mean))

pendf2$Mean <- round(pendf2$Mean, digits = 0)

pendf2<- pendf2[which(pendf2$Mean>= 1),]

##species barplot
penspecies_plot<-ggplot(pendf2, aes(x = Species, y = Mean, fill = Species))+ 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA
  ) +
  
  scale_fill_brewer(palette = "Paired") +
  theme(legend.position = "top")+
  coord_cartesian(xlim = c(1.2, NA), clip = "off")

#theme_pubr
penpp1<- penspecies_plot+
  theme_pubr(
    base_size = 17,
    base_family = "",
    border = FALSE,
    margin = TRUE,
    legend = c("top"),
    x.text.angle = 0
  )

#labelling the plot
pen<- penpp1 + 
  labs(
    x = "Host", 
    y = "Percent resistance", 
    fill = "Host", 
    color = "Host", 
    subtitle = "Penicillin"
  ) +
  scale_color_discrete(name = "Host") +
  scale_x_discrete(
    name = "Host", 
    limits = c("Cattle", "Chicken", "Goats", "Pigs", "Sheep", "Turkey", "Environment")
  ) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "none"
  ) +
  theme_Publication()

#call the plot

pen<-pen+ scale_y_continuous(limits = c(0, 100))+
  rotate_x_text(60)

pen

#AMR resistance levels for cephalosporins by species

##Pooled prevalence for salmonella by antibiotic compounds- bargraphs
#aggregate all resistant isolates and all NIsolates by speciesag, compound and region
cephSalRes = aggregate(no_isolates_resistant ~  antimicrobial_compound+species+country , data = AMR_ceph, FUN = sum)
cephSalAll = aggregate(no_isolate ~  antimicrobial_compound+species+country , data = AMR_ceph, FUN = sum)

#and divide ResIso/NIsolates to return true mean
cephSalMean = round((cephSalRes$no_isolates_resistant /cephSalAll$no_isolate)*100, digits = 0)
cephSalMeandf = as.data.frame(cbind(cephSalRes$antimicrobial_compound,  cephSalRes$species, cephSalMean, cephSalAll$no_isolate, cephSalRes$country))
colnames(cephSalMeandf) = c("antimicrobial_compound", "Species","Mean","NIsolates", "country")
cephSalMeandf$Mean = as.numeric(as.character(cephSalMeandf$Mean))
cephSalMeandf$NIsolates = as.numeric(as.character(cephSalMeandf$NIsolates))
  
cephSalMeandf$Mean <- as.numeric(as.character(cephSalMeandf$Mean))


#group by species
cephdf1 <- cephSalMeandf %>%
  group_by(Species) %>%
  select(
    "country",
    "antimicrobial_compound",
    "Mean",
    "NIsolates")


##Aggregate mean resistance levels for cephalosporins by species
cephdf2 <-aggregate(Mean~Species+antimicrobial_compound, data = cephdf1, FUN=mean)

#changing percentage resistance into numeric class and rounding off

cephdf2$Mean <- as.numeric(as.character(cephdf2$Mean))

cephdf2$Mean <- round(cephdf2$Mean, digits = 0)

cephdf2<- cephdf2[which(cephdf2$Mean>= 1),]

##species barplot
cephspecies_plot<-ggplot(cephdf2, aes(x = Species, y = Mean, fill = Species))+ 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA
  ) +
  
  scale_fill_brewer(palette = "Paired") +
  theme(legend.position = "top")+
  coord_cartesian(xlim = c(1.2, NA), clip = "off")

#theme_pubr
cephpp1<- cephspecies_plot+
  theme_pubr(
    base_size = 17,
    base_family = "",
    border = FALSE,
    margin = TRUE,
    legend = c("top"),
    x.text.angle = 0
  )
  
#labelling the plot
ceph<- cephpp1 + 
  labs(
    x = "Host", 
    y = "Percent resistance", 
    fill = "Host", 
    color = "Host", 
    subtitle = "Cephalosporins"
  ) +
  scale_color_discrete(name = "Host") +
  scale_x_discrete(
    name = "Host", 
    limits = c("Cattle", "Chicken", "Goats", "Pigs", "Sheep", "Turkey", "Environment")
  ) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "none"
  ) +
  theme_Publication()

#call the plot

ceph<-ceph+ scale_y_continuous(limits = c(0, 100))+
  rotate_x_text(60)

ceph


#AMR resistance levels for aminoglycosides by species

##Pooled prevalence for salmonella by antibiotic compounds- bargraphs
#aggregate all resistant isolates and all NIsolates by speciesag, compound and region
aminoSalRes = aggregate(no_isolates_resistant ~  antimicrobial_compound+species+country , data = AMR_amino, FUN = sum)
aminoSalAll = aggregate(no_isolate ~  antimicrobial_compound+species+country , data = AMR_amino, FUN = sum)

#and divide ResIso/NIsolates to return true mean
aminoSalMean = round((aminoSalRes$no_isolates_resistant /aminoSalAll$no_isolate)*100, digits = 0)
aminoSalMeandf = as.data.frame(cbind(aminoSalRes$antimicrobial_compound,  aminoSalRes$species, aminoSalMean, aminoSalAll$no_isolate, aminoSalRes$country))
colnames(aminoSalMeandf) = c("antimicrobial_compound", "Species","Mean","NIsolates", "country")
aminoSalMeandf$Mean = as.numeric(as.character(aminoSalMeandf$Mean))
aminoSalMeandf$NIsolates = as.numeric(as.character(aminoSalMeandf$NIsolates))

aminoSalMeandf$Mean <- as.numeric(as.character(aminoSalMeandf$Mean))


#group by species
aminodf1 <- aminoSalMeandf %>%
  group_by(Species) %>%
  select(
    "country",
    "antimicrobial_compound",
    "Mean",
    "NIsolates")


##Aggregate mean resistance levels for aminoglycosides by species
aminodf2 <-aggregate(Mean~Species+antimicrobial_compound, data = aminodf1, FUN=mean)

#changing percentage resistance into numeric class and rounding off

aminodf2$Mean <- as.numeric(as.character(aminodf2$Mean))

aminodf2$Mean <- round(aminodf2$Mean, digits = 0)

aminodf2<- aminodf2[which(aminodf2$Mean>= 1),]

##species barplot
aminospecies_plot<-ggplot(aminodf2, aes(x = Species, y = Mean, fill = Species))+ 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA
  ) +
  
  scale_fill_brewer(palette = "Paired") +
  theme(legend.position = "top")+
  coord_cartesian(xlim = c(1.2, NA), clip = "off")

#theme_pubr
aminopp1<- aminospecies_plot+
  theme_pubr(
    base_size = 17,
    base_family = "",
    border = FALSE,
    margin = TRUE,
    legend = c("top"),
    x.text.angle = 0
  )

#labelling the plot

amino<- aminopp1 + 
  labs(
    x = "Host", 
    y = "Percent resistance", 
    fill = "Host", 
    color = "Host", 
    subtitle = "Aminoglycosides"
  ) +
  scale_color_discrete(name = "Host") +
  scale_x_discrete(
    name = "Host", 
    limits = c("Cattle", "Chicken", "Goats", "Pigs", "Sheep", "Turkey", "Environment")
  ) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "none"
  ) +
  theme_Publication()

#call the plot

amino<-amino+ scale_y_continuous(limits = c(0, 100))+
  rotate_x_text(60)

amino


#AMR resistance levels for fluoroquinolones by species

##Pooled prevalence for salmonella by antibiotic compounds- bargraphs
#aggregate all resistant isolates and all NIsolates by speciesag, compound and region
floroSalRes = aggregate(no_isolates_resistant ~  antimicrobial_compound+species+country , data = AMR_floro, FUN = sum)
floroSalAll = aggregate(no_isolate ~  antimicrobial_compound+species+country , data = AMR_floro, FUN = sum)

#and divide ResIso/NIsolates to return true mean
floroSalMean = round((floroSalRes$no_isolates_resistant /floroSalAll$no_isolate)*100, digits = 0)
floroSalMeandf = as.data.frame(cbind(floroSalRes$antimicrobial_compound,  floroSalRes$species, floroSalMean, floroSalAll$no_isolate, floroSalRes$country))
colnames(floroSalMeandf) = c("antimicrobial_compound", "Species","Mean","NIsolates", "country")
floroSalMeandf$Mean = as.numeric(as.character(floroSalMeandf$Mean))
floroSalMeandf$NIsolates = as.numeric(as.character(floroSalMeandf$NIsolates))

floroSalMeandf$Mean <- as.numeric(as.character(floroSalMeandf$Mean))


#group by species
florodf1 <- floroSalMeandf %>%
  group_by(Species) %>%
  select(
    "country",
    "antimicrobial_compound",
    "Mean",
    "NIsolates")


##Aggregate mean resistance levels for fluoroquinolones by species
florodf2 <-aggregate(Mean~Species+antimicrobial_compound, data = florodf1, FUN=mean)

#changing percentage resistance into numeric class and rounding off

florodf2$Mean <- as.numeric(as.character(florodf2$Mean))

florodf2$Mean <- round(florodf2$Mean, digits = 0)

florodf2<- florodf2[which(florodf2$Mean>= 1),]

##species barplot
florospecies_plot<-ggplot(florodf2, aes(x = Species, y = Mean, fill = Species))+ 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA
  ) +
  
  scale_fill_brewer(palette = "Paired") +
  theme(legend.position = "top")+
  coord_cartesian(xlim = c(1.2, NA), clip = "off")

#theme_pubr
floropp1<- florospecies_plot+
  theme_pubr(
    base_size = 17,
    base_family = "",
    border = FALSE,
    margin = TRUE,
    legend = c("top"),
    x.text.angle = 0
  )

#labelling the plot

floro<- floropp1 + 
  labs(
    x = "Host", 
    y = "Percent resistance", 
    fill = "Host", 
    color = "Host", 
    subtitle = "Fluoroquinolones"
  ) +
  scale_color_discrete(name = "Host") +
  scale_x_discrete(
    name = "Host", 
    limits = c("Cattle", "Chicken", "Goats", "Pigs", "Sheep", "Turkey", "Environment")
  ) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "none"
  ) +
  theme_Publication()

#call the plot

floro<-floro+ scale_y_continuous(limits = c(0, 100))+
  rotate_x_text(60)

floro

#AMR resistance levels for sulfonamides by species

##Pooled prevalence for salmonella by antibiotic compounds- bargraphs
#aggregate all resistant isolates and all NIsolates by speciesag, compound and region
sulfSalRes = aggregate(no_isolates_resistant ~  antimicrobial_compound+species+country , data = AMR_sulf, FUN = sum)
sulfSalAll = aggregate(no_isolate ~  antimicrobial_compound+species+country , data = AMR_sulf, FUN = sum)

#and divide ResIso/NIsolates to return true mean
sulfSalMean = round((sulfSalRes$no_isolates_resistant /sulfSalAll$no_isolate)*100, digits = 0)
sulfSalMeandf = as.data.frame(cbind(sulfSalRes$antimicrobial_compound,  sulfSalRes$species, sulfSalMean, sulfSalAll$no_isolate, sulfSalRes$country))
colnames(sulfSalMeandf) = c("antimicrobial_compound", "Species","Mean","NIsolates", "country")
sulfSalMeandf$Mean = as.numeric(as.character(sulfSalMeandf$Mean))
sulfSalMeandf$NIsolates = as.numeric(as.character(sulfSalMeandf$NIsolates))

sulfSalMeandf$Mean <- as.numeric(as.character(sulfSalMeandf$Mean))
  
#group by species
sulfdf1 <- sulfSalMeandf %>%
  group_by(Species) %>%
  select(
    "country",
    "antimicrobial_compound",
    "Mean",
    "NIsolates")


##Aggregate mean resistance levels for sulfonamides by species
sulfdf2 <-aggregate(Mean~Species+antimicrobial_compound, data = sulfdf1, FUN=mean)


#changing percentage resistance into numeric class and rounding off

sulfdf2$Mean <- as.numeric(as.character(sulfdf2$Mean))

sulfdf2$Mean <- round(sulfdf2$Mean, digits = 0)

sulfdf2<- sulfdf2[which(sulfdf2$Mean>= 1),]

##species barplot
sulfspecies_plot<-ggplot(sulfdf2, aes(x = Species, y = Mean, fill = Species))+ 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA
  ) +
  
  scale_fill_brewer(palette = "Paired") +
  theme(legend.position = "top")+
  coord_cartesian(xlim = c(1.2, NA), clip = "off")

#theme_pubr
sulfpp1<- sulfspecies_plot+
  theme_pubr(
    base_size = 17,
    base_family = "",
    border = FALSE,
    margin = TRUE,
    legend = c("top"),
    x.text.angle = 0
  )

#labelling the plot

sulf<- sulfpp1 + 
  labs(
    x = "Host", 
    y = "Percent resistance", 
    fill = "Host", 
    color = "Host", 
    subtitle = "Sulfonamides"
  ) +
  scale_color_discrete(name = "Host") +
  scale_x_discrete(
    name = "Host", 
    limits = c("Cattle", "Chicken", "Goats", "Pigs", "Sheep", "Turkey", "Environment")
  ) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "none"
  ) +
  theme_Publication()

#call the plot

sulf<-sulf+ scale_y_continuous(limits = c(0, 100))+
  rotate_x_text(60)





cp.species_plot <- (tet | pen)/ (sulf | ceph)+
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")")&
  theme(plot.tag = element_text(face = "bold"))
print(cp.species_plot)


# plot.comb <- patchwork::patchworkGrob(cp.species_plot )
# 
# f.p<- grid.arrange(plot.comb, ncol = 1, nrow = 1,left = yleft, bottom = bottom)
# 
# f.p

# 
# pigmed<-aggregate(Mean~Species, data = sulfdf2, FUN = median)%>%
#   #reaodering the Mean column
#   arrange(desc(Mean))
# 
# 
# 
# 
# AMR_clean%>%
#   count(doi)

# 
# ggsave(filename = "tetracycline.tiff",
#               plot = tet,
#                device = "tiff",
#               width = 10,
#               height = 8,
#                units = "in",
#                dpi = 700,
#                compression = "lzw")
#        