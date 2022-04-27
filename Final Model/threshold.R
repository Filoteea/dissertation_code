# 19/04/2022
# Filoteea Moldovan


# creating bar graphs for the occurrence of measurable data points 


# Libraries ----
library(ggplot2)  # data manipulation and visualization
library(viridis)    # colour-blind friendly colour scheme
library(gganimate)  # create animation
library(gifski)     # export animation as a GIF

# Import data ----
occurence <- read.csv("inputs/occurence.csv")

# Create palette
magic.palette <- c("#698B69", "#5D478B", "#5C5C5C", "#CD6090", "#EEC900", "#5F9EA0", "#6CA6CD")    # defining 7 colours
names(magic.palette) <- levels(occurence$Site)  

# Plot bar graph
(hist <- ggplot(occurence, aes(x = Scenario, y = Occurrence, fill = Site)) +
    geom_histogram(stat = "identity", position = "dodge") + 
    scale_y_continuous(limits = c(0, 100)) +
    scale_fill_manual(values = magic.palette,                       
                      name = "Measurement site") +                
    labs(title = "", 
         x = "", y = "Occurrence (%) \n") +  
    theme_bw() +
    theme(panel.grid = element_blank(), 
          axis.text = element_text(size = 12), 
          axis.text.x = element_text(angle = 45, hjust = 1), 
          axis.title = element_text(size = 12), 
          plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , "cm"), 
          legend.title = element_text(face = "bold"),
          legend.position = "top", 
          legend.box.background = element_rect(color = "grey", size = 0.3)))
