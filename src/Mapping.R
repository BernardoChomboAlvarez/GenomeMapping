# Title: Mapping
# Author(s): Bernardo Chombo-Álvarez, Miguel Ángel Flores Varela
# Date: 11-04-2024
# Version: 1.0.0
# License: Apache-2.0
# 
# Description: Program for mapping all the clinvar mutations in the reference genome.
# Here we use the divide et vinci principle due to the large datasets we have. 

## Libraries
library(dplyr)
library(ggplot2)

## Set working directory
setwd("/home/chombo/Downloads/Proyecto_Genomica_Humana/CHOMBO/data")

## Import the data
## Chromosomes namnes
chr_names <- read.csv("chr_ids.csv",header = FALSE,sep = ",")
names(chr_names) <- c("chr","chr_id","length")
chr_ids <- chr_names$chr_id

## Mutations data
clinvar_df <- read.table("clinvar.cleaned.tsv",header = FALSE,sep = "\t")
names(clinvar_df) <- c("gene","chr","position","id")
clinvar_df <- merge(clinvar_df,chr_names, by = "chr", all.x = TRUE)
clinvar_df$chr_id <- as.character(clinvar_df$chr_id)
clinvar_df <- clinvar_df %>%
    select(-c("id","length"))

## Human Genome annotation data
grhc_complete <- read.csv("grhc38_complete.csv",header = TRUE,sep = ",")
grhc_complete$start <- as.numeric(grhc_complete$start)
grhc_complete$end <- as.numeric(grhc_complete$end)

## Characterize the mutations
## Divide the data for each of the chromosomes
mutations <- list()
for (chr in unique(grhc_complete$chr_id)) {
    ## Create a reference for each one of the chromosomes
    chr_ref <- grhc_complete %>%
        filter(chr_id == chr) %>%
        arrange(start) %>%
        distinct() %>%
        `rownames<-`(NULL) %>%
        tibble::rowid_to_column("index")
    
    ## Create a data subset for the mutations of each chromosomes
    chr_mut <- clinvar_df[clinvar_df$chr_id == chr,] %>%
        arrange(position) %>%
        `rownames<-`(NULL) %>%
        tibble::rowid_to_column("index")

    ## Map all the mutations in the reference genome data subset
    mapped <- data_frame()
    known_indexes <- c()
    unknown_indexes <- c()
    categories <- c("intergenic","CDS","intron","exon","gene") ## Each of this categories are hierarchically ordered
    for (category in categories) {
        ## Divide again the subset by each of the four categories
        reference <- chr_ref %>% filter(feature == category)
        
        ## Now itterate all the subset of mutations for each chromosome
        for (i in 1:length(chr_mut$position)) {
            ## Skip known indexes
            if (i %in% known_indexes) {
                next
            }
            
            ## Search index in the reference (only one result should be TRUE)
            logic_res <- with(reference, start <= chr_mut[i,"position"] & end >= chr_mut[i,"position"])
            
            ## Store the indexes of the TRUE result
            indexes <- which(logic_res)
            
            ## Continue if no index is found
            if (length(indexes) == 0) {
                next
            } else {
                ## Create a vector with the index and the category
                mapped <- rbind(mapped,c(i,category))
                known_indexes <- c(known_indexes,i)
            }
        }   
    }
    
    ## If there is any unknown index, compare it against the mutations indexes
    unknown_indexes <- chr_mut[!(chr_mut$index %in% known_indexes),"index"] %>% 
        unname()
    
    ## Create a dataframe with the unknown indexes if any
    if (length(unknown_indexes) != 0) {
        unmapped <- cbind(unknown_indexes,"UNMAPPED")
        unmapped <- as.data.frame(unmapped)
        names(unmapped) <- c("index","feature")
        unmapped$index <- as.numeric(unmapped$index)
        complete <- rbind(mapped,unmapped)
        complete <- complete[order(complete$index),]
    } else {
        ## Merge the known indexes
        complete <- mapped
        names(complete) <- c("index","feature")
        complete$index <- as.numeric(complete$index)
        complete <- complete[order(complete$index),]
    }
    
    ## Add the chromosome's id data
    complete <- merge(chr_mut,complete,by = "index") %>%
        select(-index) %>%
        `rownames<-`(NULL)
        
    ## Merge the data
    mutations[[chr]] <- complete
    
}

## Join all the data
mutations <- do.call(rbind, mutations) %>%
    `rownames<-`(NULL)
mutations <- merge(mutations,chr_names, by="chr_id")
mutations$chr <- paste0("Chr ",mutations$chr)
mutations <- mutations %>%
    select(-c("chr.y","length")) %>%
    mutate(chr = paste(chr,chr.x,sep = " ")) %>%
    select(-chr.x) %>%
    distinct()