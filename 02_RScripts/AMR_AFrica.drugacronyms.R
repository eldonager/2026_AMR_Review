
library(tidyverse)
library(gt)

#Load data

library(kableExtra)
library(tinytex)
library(AER)
library(wesanderson)
library(modelsummary)
library(flextable)

#Load data

AMRdf2 <- read.csv("AMR_clean.csv")




#Join Amphenicols and Fluoroquinolones drugs to one class

AMRdf2$antibiotic_class[AMRdf2$antibiotic_class == "Amphenicols"] <- "Fluoroquinolones"


data3 <- AMRdf2 %>%
  select(antimicrobial,antimicrobial_compound, antibiotic_class,
         who_classification) %>%
  unique() %>%
  arrange(antimicrobial)

data4 <- with(data3,  data3[order(who_classification) , ])


my_table <- flextable(data4) %>% 
  autofit()

#aligning all colunms to center except antimicrobial
my_table <- my_table %>% 
  flextable::align(align = "center", j = c(2:4), part = "all") 

my_table <-  my_table %>%  
  fontsize(i = 1, size = 12, part = "header") %>%   # adjust font size of header
  bold(i = 1, bold = TRUE, part = "header") %>%
  set_header_labels(         # Rename the columns in original header row
    antimicrobial = "Antimicrobial", 
    antimicrobial_compound = "Antimicrobial compound",                  
    antimicrobial_class = "Antimicrobial class",
    who_classification = "WHO classification")

#write table as a csv file

#write.csv(data4, "AMR_Africa.drugacronyms.csv", row.names = FALSE)


