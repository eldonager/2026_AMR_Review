

library (tidyverse)
library(janitor)

AMR_sal <- read.csv("AMR_sal.csv", fileEncoding = "Latin1", check.names = F)

AMR_sal <- clean_names(AMR_sal)


AMR_sal %>% tabyl(country)


colnames(AMR_sal)

AMR_sal %>%
  count(doi)%>%
  View()

AMR_sal %>%
  count(title)%>%
  View()

AMR_sal %>%
  count(author,title,doi,antimicrobial)%>%
  View()

AMR_sal %>%
  count(antibiotic_class) %>%
  View()

AMR_sal %>%
  count(antimicrobial) %>%
  View()


AMR_sal<- AMR_sal %>%
  mutate(antimicrobial = 
           ifelse(antimicrobial == "sulfisoxazole", 
                  "Sulfisoxazole", antimicrobial))



AMR_sal<- AMR_sal %>%
  mutate(antimicrobial = 
           ifelse(antimicrobial == "sulfisoxazole", 
                  "Sulfisoxazole", antimicrobial))

AMR_sal %>%
  count(who_classification) %>%
  View()

###writing the new csv file
write.table(AMR_sal, file = "AMR_clean.csv", row.names = FALSE,
            sep = ",")
