
###loading packages
library(tidyverse)
library(tmap)
library(leaflet) 
library(sf) 
library(RColorBrewer)
library(htmltools) 
library(leafsync) 
library(kableExtra) 
library(rnaturalearth)
library(viridis)
library(mapview)


sf_use_s2(FALSE)
#> Spherical geometry (s2) switched off
##loading data
AMR_clean <- read.csv("AMR_clean.csv")

world <- ne_countries(type = 'countries', 
                      scale = 'small', 
                      returnclass = "sf") 


###using map view package to plot world AMR rates


df <- AMR_clean%>%
  group_by(country) %>%
  select(country, percent_resistant, x_coordinate,
         y_coordinate, species)%>%
  set_names(c("country", "percentage_resistance",
              "long", "lat", "species")) %>%
  na.omit()

df %>%
  count(country)%>%   
  view()
###changing percentage_resistance to numeric class
df$percentage_resistance <- as.numeric(as.character(df$percentage_resistance))

##Summarizing the data
df1 <- df %>%
  group_by(country) 
  # summarise(percentage_resistance = mean(percentage_resistance, na.rm = TRUE))

world_amr1 <- world %>%
  filter(continent == "Africa") %>%
  select(country = sovereignt,iso_a3, geometry) #%>%
 # left_join(df1, by = "country")

##Adding number of studies to the dataset
data1 <- AMR_clean %>% 
  group_by(country, author) %>% 
  filter(row_number()==1) %>% 
  dplyr:: select(country) %>%
  distinct() %>% as.data.frame() # drop duplicates.
data1

data2 <- data1 %>% 
  group_by( country) %>% 
  mutate(no_studies = n()) %>%
  dplyr:: select(country, no_studies) %>%
  distinct() %>% as.data.frame() # drop duplicates.
data2

#Joining the data to add number of studies
world_amr2 <- world_amr1 %>%
  left_join(data2, by = "country")%>%
  unique()

world_amr2 %>%
  mapview(zcol = "no_studies",
          col.regions = viridisLite::plasma) 


tm_shape(world_amr2) + 
  tm_polygons(col = "no_studies",  
              title = "number of studies",
              palette = "RdYlBu") +
  tm_text("no_studies", size = 0.7) 
  





# library(tmap)
# 
# tmap_mode("view")
# 
# map <- tm_shape(world_amr2) + 
#   tm_polygons(col = "no_studies",  
#               n = 1,
#               palette = "Blues") +
#   tm_text("iso_a3")  # Add country names
# 
# map <- map + tm_shape(world_amr2) +
#   tm_text("no_studies", size = "no_studies", root = 2)+
#   tm_text("iso_a3", size =5)# Add number of studies
# 
# map









