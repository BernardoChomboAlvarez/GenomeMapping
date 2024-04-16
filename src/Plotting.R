# Title: Plotting
# Author(s): Bernardo Chombo-Álvarez, Miguel Ángel Flores Varela
# Date: 11-04-2024
# Version: 1.0.0
# License: Apache-2.0
# 
# Description: Program for plotting all the outputs from the DataTreatment.R and the Mapping.R programs

## Libraries
library(dplyr)
library(ggplot2)
library(waffle)

## Set working directory
setwd("/home/chombo/Downloads/Proyecto_Genomica_Humana/CHOMBO/data")

## Import the data
chr_names <- read.csv("chr_ids.csv",header = FALSE,sep = ",")
names(chr_names) <- c("chr","chr_id","length")
chr_ids <- chr_names$chr_id

mutations <- read.csv("mutations_complete.csv",header = TRUE,sep = ",")
mutations <- mutations %>%
    distinct() %>%
    select(-c("position","gene","chr"))

## Extract the mutations per catgeory for each of the chromosomes
mut_arranged <- mutations %>%
    mutate(count = 1) %>%
    group_by(chr_id, feature) %>%
    summarize(count = sum(count)) %>%
    ungroup()

## Make a contigence table with the four features as columns and the counts for each one for each of the chromosomes
mut_contingence <- tidyr::pivot_wider(data = mut_arranged, names_from = feature, values_from = count)
mut_contingence <- replace(mut_contingence, is.na(mut_contingence), 0) %>%
    mutate(total = CDS+exon+intergenic+intron)

## Modify the contigence table so that now it has the percentages
mut_contingence.per <- mut_contingence %>%
    mutate(CDS = trunc((CDS/total)*100*10^2)/10^2,
           exon = trunc((exon/total)*100*10^2)/10^2,
           intergenic = trunc((intergenic/total)*100*10^2)/10^2,
           intron = trunc((intron/total)*100*10^2)/10^2)

## Clear the contingence table with percentages
mut_contingence.per <- merge(mut_contingence.per,chr_names, by = "chr_id") %>%
    select(-length) %>%
    mutate(id = paste0("Chr ",chr," (",total,")")) %>%
    select(-c("chr_id","chr","total"))

## Expand the table into two columns: id and feature
mut_expanded <- tidyr::pivot_longer(data = mut_contingence.per, 
                                    cols = -id, 
                                    names_to = "feature", 
                                    values_to = "count")

## Create a dataframe with the total counts for each feature
feature_counts <- as.data.frame(table(mutations$feature))
names(feature_counts) <- c("feature","count")
feature_counts <- feature_counts %>%
    mutate(percentage = (count/sum(count))*100) %>%
    mutate(percentage = trunc(percentage*10^2)/10^2)

## PLOTTING
## Bar plot of all the features
plot <- ggplot(feature_counts, aes(x = factor(feature), y = count, fill = factor(feature))) +
    geom_bar(stat = "identity") + 
    geom_text(aes(label = count), vjust = -0.5, size = 4) +
    scale_fill_manual(values = c("#8400EA", "#EA00DB", "#EA0066", "#EA8400", "#DBEA00")) +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.background = element_rect(colour = "white"),
        plot.background = element_rect(fill = "white", colour = "white"),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 1),
        strip.background = element_rect(fill = "white", colour = "white"),
        plot.title = element_text(size = 20)) +
    xlab("Features") +
    ylab("Frequency") +
    labs(title = "Mendelian disease-associated mutations") +
    annotate(geom="text", x=3.5, y=6e+05, label=paste0("Total = ",sum(feature_counts$count)),
             size = 5,color="black")

print(plot)

## Pie chart of all the general features
plot <- ggplot(feature_counts, aes(x = 1, y = count, fill = feature)) +
    geom_bar(stat = "identity") + 
    coord_polar(theta = "y") +
    geom_text(aes(x=1.6, label=paste0(percentage, "%")),
              position = position_stack(vjust=0.5)) +
    scale_fill_manual(values = c("#8400EA", "#EA00DB", "#EA0066", "#EA8400", "#DBEA00")) +
    theme_void() +
    guides(fill = guide_legend(title = "Features")) + 
    theme(legend.position = "right", 
          legend.background = element_rect(colour = "white"),
          plot.background = element_rect(fill = "white", colour = "white"),
          panel.background = element_rect(fill = "white", colour = "white", linewidth = 1),
          strip.background = element_rect(fill = "white", colour = "white"),
          plot.title = element_text(size = 20),
          axis.text.x = element_blank(),  # remove x axis labels
          axis.ticks.x = element_blank(),  # remove x axis ticks
          axis.text.y = element_blank(),  # remove y axis labels
          axis.ticks.y = element_blank()) +
    xlab("Features") +
    ylab("Frequency")

print(plot)

## Pie chart for each of the chromosomes
plot <- ggplot(mut_expanded, aes(x = "", y = count, fill = feature)) +
    geom_bar(stat = "identity") + 
    coord_polar("y", start = 0) +
    geom_text(aes(x=1.5, label=paste0(count, "%")),
              position = position_stack(vjust=0.5),size=2.5) +
    facet_wrap(vars(factor(id, levels = unique(id)))) +
    scale_fill_manual(values = c("#8400EA", "#EA00DB", "#EA0066", "#EA8400", "#DBEA00")) +
    theme_void() +
    guides(fill=guide_legend(title="Features")) + 
    theme(legend.position = "right", legend.background = element_rect(colour = "white"),
          plot.background = element_rect(fill = "white", colour = "white"),
          panel.background = element_rect(fill = "white", colour="white",linewidth = 1),
          strip.background = element_rect(fill="white",colour = "white"),
          plot.title = element_text(size = 20),
          axis.text.x = element_blank(), #remove x axis labels
          axis.ticks.x = element_blank(), #remove x axis ticks
          axis.text.y = element_blank(),  #remove y axis labels
          axis.ticks.y = element_blank()) +
    xlab("Features") +
    ylab("Frequency")

print(plot)


## Waffle plot for overview
x <- feature_counts$percentage
names(x) <- paste0(feature_counts$feature," (",feature_counts$percentage,"%)")
# Waffle chart
waffle(x, rows = 7,
       colors = c("#8400EA", "#EA00DB", "#EA0066", "#EA8400"),
       size = 1,legend_pos = "bottom")