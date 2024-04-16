# Title: DataTreatment
# Author(s): Bernardo Chombo-Álvarez, Miguel Ángel Flores Varela
# Date: 11-04-2024
# Version: 1.0.0
# License: Apache-2.0
# 
# Description: Script for defining intronic and intergenic regions all along the Human Genome.
# This can be adjusted for other genomes. The output features are: gene, exon, intron, intergenic and CDS,
# other features are deprecated in the program.

## Libraries
library(dplyr)

## Set working directory
setwd("/home/chombo/Downloads/Proyecto_Genomica_Humana/CHOMBO/data")

## Import the data
chr_names <- read.csv("chr_ids.csv",header = FALSE,sep = ",")
names(chr_names) <- c("chr","chr_id","stop")
chr_names <- chr_names %>%
    mutate(start = 1) %>%
    arrange("chr","chr_id","start","stop")
chr_ids <- chr_names$chr_id

## Human Genome GRHCh88.p14 assembly
grhc38_df <- read.table("GCF_000001405.40_GRCh38.p14_genomic.gff",skip = 9,header = FALSE,sep = "\t")
grhc38_df <- grhc38_df[,1:7]
names(grhc38_df) <- c("chr_id","source","feature","start","end","score","strand")
grhc38_df <- grhc38_df[!(grhc38_df$feature %in% c("region")),]
grhc38_df <- grhc38_df[(grhc38_df$chr_id %in% chr_ids),]

## Region search
full_data <- list()
for (i in 1:length(unique(chr_names$chr_id))) {
    ## Extract the data for only one chromosome
    chr <- chr_names[i,"chr_id"]
    chr_raw <- grhc38_df[grhc38_df$chr_id == chr,]
    
    ## Change the data type for "start" and "end" columns
    chr_raw$start <- as.numeric(chr_raw$start)
    chr_raw$end <- as.numeric(chr_raw$end)
    
    ## Reorder the dataframe and reset the index
    chr_raw <- chr_raw[order(chr_raw$start,decreasing = FALSE),]
    chr_raw <- chr_raw %>%
        select(chr_id,feature,start,end) %>%
        `rownames<-`(NULL)
    
    ## Create the object for intergenic regions
    intergenic_regions <- data_frame()
    
    ## Create a data subset only with the genes
    genes <- chr_raw[chr_raw$feature == "gene",]
    
    ## Define the limits of the chromosome coding region
    pos_start <- min(genes$start)
    pos_end <- max(genes$end)
    
    ## Create an intergenic region from the beginning of the chromosome to the start of the first gene
    intergenic_regions <- rbind(intergenic_regions,c(chr,"intergenic",1,pos_start - 1))
    
    ## Remove the data before the start and end coding regions
    ## Clean also other regions inside the genome
    ## IMPORTANT!! introns and exons: only the ones inside genes
    ## We are not counting the exons and introns from pseudogenes
    chr_raw <- chr_raw[!(chr_raw$start <= (pos_start - 1)) & !(chr_raw$end >= (pos_end + 1)),]  %>%
        `rownames<-`(NULL) %>%
        tibble::rowid_to_column("index")
    
    
    ## INTERGENIC REGION SEARCH
    ## There might be genes inside genes completely or partially so we search for the maximum gene length
    ## If the start position of the next gene is higher than the end of the previous one, and its end is higher than the previous' end and the max_end, then it's an intergenic region
    max_end <- genes[1,"end"]
    for (i in 2:length(genes$feature)) {
        ## If the previous gene's end is higher than the max_end, then we change it
        if (genes[i-1,"end"] >= max_end) {
            max_end <- genes[i-1,"end"]
        }
        
        ## Here we defint what it is an intergenic region
        if ((genes[i,"start"] > genes[i-1,"end"]) & (genes[i,"end"] > genes[i-1,"end"]) & (genes[i,"end"] >= max_end) & (genes[i,"start"] > max_end)) {
            intergenic <- c(chr, "intergenic", max_end + 1, genes[i,"start"] - 1)
            intergenic_regions <- rbind(intergenic_regions, intergenic)
            
            ## Additionally, we clean what's inside the chr_raw so that we get rid of everything that corresponds to an intergenic region
            pos_start <- max_end + 1
            pos_end <- genes[i, "start"] - 1
            chr_raw <- chr_raw[!((chr_raw$start >= pos_start) & (chr_raw$start <= pos_end) & (chr_raw$feature != "gene")),]
        }
    }
    
    ## We set the region corresponding to the end of the last gene and the end of the chromosome as an intergenic region
    end <- c(chr,"intergenic",max(genes$end) + 1,chr_names[chr_names$chr_id == chr,]$stop)
    intergenic_regions <- rbind(intergenic_regions,end)
    names(intergenic_regions) <- c("chr_id","feature","start","end")
    
    ## Merge the intergenic regions with the chromosome raw data
    chr_raw <- chr_raw %>% select(-index)
    chr_wintergen <- rbind(chr_raw,intergenic_regions)
    chr_wintergen$start <- as.numeric(chr_wintergen$start)
    chr_wintergen$end <- as.numeric(chr_wintergen$end)
    chr_wintergen <- chr_wintergen[order(chr_wintergen$start,decreasing = FALSE),]
    chr_wintergen <- chr_wintergen %>%
        select(chr_id,feature,start,end) %>%
        `rownames<-`(NULL)
    
    ## We then split the dataframe into different dataframes with CDS, exons and genes information
    cds <- chr_wintergen[chr_wintergen$feature == "CDS",]
    exons <- chr_wintergen[chr_wintergen$feature == "exon",]
    genes <- chr_wintergen[chr_wintergen$feature == "gene",] ## we repeat this so that we get no erros
    
    ## Now clean the merged dataframe so that we only stay with the exons and intergenic regions
    ## As every exon inside the data is now inside genes, we do this
    chr_wintergen <- chr_wintergen[(chr_wintergen$feature %in% c("intergenic","exon")),] %>%
        distinct()
    
    
    ## INTRONIC REGIONS
    ## With a similar approach for the intergenic regions, we extract the exons between intergenic regions
    ## We define an intron as the spaces between the first and last exons with the intergenic regions
    ## We look for exons that are not overlapped
    intronic_regions <- data_frame()
    for (i in 2:length(intergenic_regions$feature)) {
        ## Define the start position as the end of the previous intergenic region
        ## Define the end position as the start of the next intergenic region
        start <- as.numeric(intergenic_regions[i-1,"end"])
        end <- as.numeric(intergenic_regions[i,"start"])
        
        ## Extract the data, i.e. exons
        data <- exons[(exons$end < end & exons$start > start),] %>%
            `rownames<-`(NULL)
        
        ## In case there were no exons, we continue
        if (length(data$feature) == 0) {
            next
        }
        
        ## If there was a space between the first exon and the end of the previous intergenic region
        if (start - data[1,"start"] <= -2) {
            intron <- c(chr, "intron", start + 1, data[1,"start"] - 1)
            intronic_regions <- rbind(intronic_regions, intron)
        }
        
        ## If there was a space between the last exon and the beginning of the next intergenic region
        if (data[length(data$feature),"end"] - end <= -2) {
            intron <- c(chr, "intron", data[length(data$feature),"end"] + 1, end - 1)
            intronic_regions <- rbind(intronic_regions, intron)
        }
        
        ## If there was only one exon, we continue, else, the same exact logic for intergenic regions
        if (length(data$feature) == 1) {
            next
        } else {
            max_end <- data[1,"end"]
            for (i in 2:length(data$feature)) {
                ## Define the max_end
                if (data[i-1,"end"] >= max_end) {
                    max_end <- data[i-1,"end"]
                }
                
                ## Look for intronic regions
                if ((data[i,"start"] > data[i-1,"end"]) & (data[i,"end"] > data[i-1,"end"]) & (data[i,"end"] >= max_end) & (data[i,"start"] > max_end)) {
                    intron <- c(chr, "intron", max_end + 1, data[i,"start"] - 1)
                    intronic_regions <- rbind(intronic_regions, intron)
                }
            }
        }
    }
    names(intronic_regions) <- c("chr_id","feature","start","end")
    
    ## Merge the complete data for each chromosome
    complete <- rbind(genes,chr_wintergen,cds,exons,intronic_regions)
    complete$start <- as.numeric(complete$start)
    complete$end <- as.numeric(complete$end)
    complete <- complete[order(complete$start),] %>%
        `rownames<-`(NULL)
    
    ## Append them in a list
    full_data[[chr]] <- complete
}

## Transform the list
grhc_complete <- do.call(rbind, full_data) %>%
    `rownames<-`(NULL)

## Save the data
write.csv(grhc_complete,file = "grhc38_complete.csv",row.names = FALSE,col.names = FALSE,quote = FALSE)

## Save the data for each of the chromosomes
grhc_complete <- grhc_complete[order(grhc_complete$chr_id),]
for (chr in unique(grhc_complete$chr_id)) {
    data <- grhc_complete[grhc_complete$chr_id == chr,]
    data <- data[order(data$start),]
    write.csv(data,file = paste0(chr,"_data.csv"),row.names = FALSE,col.names = FALSE,quote = FALSE)

}
