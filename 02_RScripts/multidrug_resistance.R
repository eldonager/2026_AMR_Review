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


AMR_df<- AMR_clean %>%
  dplyr::select(country,region, antimicrobial_compound,antibiotic_class,
                species, pathogen, no_isolate,mdr_percentage)                                      



AMR_df1<-aggregate(mdr_percentage~antibiotic_class+antimicrobial_compound+species, 
          data = AMR_df,FUN = mean)


AMRdf2<- AMR_df1%>%
  na.omit


#Join 1st,2nd,3rd and 4th generation cephalosporins drugs to one class

AMRdf2$antibiotic_class[AMRdf2$antibiotic_class == "1st generation Cephalosporin"] <- "Cephalosporins"
AMRdf2$antibiotic_class[AMRdf2$antibiotic_class == "2nd generation Cephalosporin"] <- "Cephalosporins"
AMRdf2$antibiotic_class[AMRdf2$antibiotic_class == "3rd generation Cephalosporin"] <- "Cephalosporins"
AMRdf2$antibiotic_class[AMRdf2$antibiotic_class == "4th generation Cephalosporin"] <- "Cephalosporins"

#Join Amphenicols and Fluoroquinolones drugs to one class

AMRdf2$antibiotic_class[AMRdf2$antibiotic_class == "Amphenicols"] <- "Fluoroquinolones"


#Remove antibiotic class Lincosamides, Monobactam, Polymixins, Quinoxalines
# Filter the dataset for pigs
pigsAMRdf2 <- AMRdf2 %>%
  filter(species == "Pigs")
pigsAMRdf2 <- pigsAMRdf2  %>%
filter(antibiotic_class != "Macrolides") %>%
 filter(antibiotic_class != "Nitrofurans") %>%
  filter(antibiotic_class != "Quinoxalines") 
 
##pigs plot
pigsmdr_plot<-ggplot(pigsAMRdf2 , aes(x = antibiotic_class, y = mdr_percentage, fill = species))+ 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA
  ) +
  
  scale_fill_brewer(palette = "Paired") +
  theme(legend.position = "top")+
  coord_cartesian(xlim = c(1.2, NA), clip = "off")

pigspp1<- pigsmdr_plot+
  theme_pubr(
    base_size = 17,
    base_family = "",
    border = FALSE,
    margin = TRUE,
    legend = c("top"),
    x.text.angle = 60
  )

pigs.p<- pigspp1 + 
  labs(
    x = "Antibiotic class", 
    y = "MDR %", 
    fill = "antibiotic_class", 
    color = "antibiotic_classs", 
    subtitle = "Pigs"
  ) +
  theme_Publication()+
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "none"
  )


pigs.p<-pigs.p+ scale_y_continuous(limits = c(0, 100))

# Filter the dataset for cattle
cattleAMRdf2 <- AMRdf2 %>%
  filter(species == "Cattle")
cattleAMRdf2 <- cattleAMRdf2 %>%
  filter(antibiotic_class != "Macrolides") %>%
  filter(antibiotic_class != "Lincosamides") %>%
  filter(antibiotic_class != "Nitrofurans") %>%
  filter(antibiotic_class != "Polymyxins") 

##cattle plot
cattlemdr_plot<-ggplot(cattleAMRdf2 , aes(x = antibiotic_class, y = mdr_percentage, fill = species))+ 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA
  ) +
  
  scale_fill_brewer(palette = "Paired") +
  theme(legend.position = "top")+
  coord_cartesian(xlim = c(1.2, NA), clip = "off")

cattlepp2<- cattlemdr_plot+
  theme_pubr(
    base_size = 17,
    base_family = "",
    border = FALSE,
    margin = TRUE,
    legend = c("top"),
    x.text.angle = 60
  )

cattle.p<- cattlepp2 +
  labs(
    x = "Antibiotic class", 
    y = "MDR %", 
    fill = "antibiotic_class", 
    color = "antibiotic_classs", 
    subtitle = "Cattle"
  ) +
  theme_Publication()+
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "none"
  )

cattle.p<-cattle.p+ scale_y_continuous(limits = c(0, 100))



# Filter the dataset for chicken
chickAMRdf2 <- AMRdf2 %>%
  filter(species == "Chicken")

chickAMRdf2 <- chickAMRdf2 %>%
  filter(antibiotic_class != "Lincosamides") %>%
  filter(antibiotic_class != "Monobactam") %>%
  filter(antibiotic_class != "Polymyxins") %>%
  filter(antibiotic_class != "Quinoxalines")

##chicken plot
chickmdr_plot<-ggplot(chickAMRdf2 , aes(x = antibiotic_class, y = mdr_percentage, fill = species))+ 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA
  ) +
  
  scale_fill_brewer(palette = "Paired") +
  theme(legend.position = "top")+
  coord_cartesian(xlim = c(1.2, NA), clip = "off")

chickpp3<- chickmdr_plot+
  theme_pubr(
    base_size = 17,
    base_family = "",
    border = FALSE,
    margin = TRUE,
    legend = c("top"),
    x.text.angle = 60
  )

chick.p<- chickpp3 +
  labs(
    x = "Antibiotic class", 
    y = "MDR %", 
    fill = "antibiotic_class", 
    color = "antibiotic_classs", 
    subtitle = "Chicken"
  ) +
  theme_Publication()+
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "none"
  )

chick.p<-chick.p+ scale_y_continuous(limits = c(0, 100))


# Filter the dataset for environment
envAMRdf2 <- AMRdf2 %>%
  filter(species == "Environment")

envAMRdf2 <- envAMRdf2 %>%
  filter(antibiotic_class != "Lincosamides") %>%
  filter(antibiotic_class != "Monobactam") %>%
  filter(antibiotic_class != "Polymyxins") %>%
  filter(antibiotic_class != "Quinoxalines")

##environment plot
envmdr_plot<-ggplot(envAMRdf2 , aes(x = antibiotic_class, y = mdr_percentage, fill = species))+ 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA
  ) +
  
  scale_fill_brewer(palette = "Paired") +
  theme(legend.position = "top")+
  coord_cartesian(xlim = c(1.2, NA), clip = "off")

envpp4<- envmdr_plot+
  theme_pubr(
    base_size = 17,
    base_family = "",
    border = FALSE,
    margin = TRUE,
    legend = c("top"),
    x.text.angle = 60
  )

env.p<- envpp4 +
  labs(
    x = "Antibiotic class", 
    y = "MDR %", 
    fill = "antibiotic_class", 
    color = "antibiotic_classs", 
    subtitle = "Environment"
  ) +
  theme_Publication()+
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "none"
  )

env.p<-env.p+ scale_y_continuous(limits = c(0, 100))




#tiff("mdr_comb_plot.tiff", width=2800, height=2300, res=300)
mdr_comb_plot <- (chick.p | pigs.p)/ (cattle.p | env.p)+
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")")&
  theme(plot.tag = element_text(face = "bold"))
mdr_comb_plot 
# tiff("cp.species_plot.tiff", width=2600, height=2300, res=300)
#dev.off()
# plot.comb <- patchwork::patchworkGrob(mdr_comb_plot )
# 
# f.p<- grid.arrange(plot.comb, ncol = 1, nrow = 1,left = yleft, bottom = bottom)
# 
# f.p




chickmdr_data <- aggregate(mdr_percentage~antibiotic_class+species, 
          data = envAMRdf2,FUN = median)%>%
  #order mdr_percentage in descending order
  arrange(desc(mdr_percentage))





  