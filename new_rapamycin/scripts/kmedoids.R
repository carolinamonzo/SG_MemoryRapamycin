suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(phangorn))
suppressPackageStartupMessages(library(DECIPHER))
suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(ANCOMBC))
suppressPackageStartupMessages(library(mvabund))
suppressPackageStartupMessages(library(corrr))
suppressPackageStartupMessages(library(cluster))
#library(microbial)
theme_set(theme_bw())
path <- "~/workspace/MPI/SG_MemoryRapamycin/new_rapamycin/"

new_day <- gsub("-", "", as.character(Sys.Date()))

setwd("~/workspace/MPI/SG_MemoryRapamycin/new_rapamycin/scripts/")
project_path <- "~/workspace/MPI/SG_MemoryRapamycin/new_rapamycin/"
new_day <- gsub("-", "", as.character(Sys.Date()))
metadata <- read.table("~/workspace/MPI/SG_MemoryRapamycin/new_rapamycin/metadata_rapaInter.csv", header = T, row.names = 1, stringsAsFactors = F, sep = ",")

metadata <- metadata %>% filter(Sex %in% c("m"))

df <- read.table("../analysis/norm-CLEAN_ASVs_counts_merged_20240530.tsv", sep = "\t", header = T, row.names = 1, stringsAsFactors = F)

# Remove OTUs that arent found in at least 290 samples (out of 300 sequenced)
df <- df[rowSums(is.na(df)) < 295,]

# Subset only those that are there in all 3 timepoints
subset_metadata <- metadata %>%
  #filter(Treatment_SG %in% c("Control", "Rapamycin conti (RC)"))
  filter(Treatment_SG %in% c("Control"))

subset_metadata$Age <- factor(subset_metadata$Age, levels = c("Young", "Mid", "Old"))
subset_metadata <- subset_metadata %>% arrange(Age)

# Sort the dataframe and keep only columns we are using.
df <- df[rownames(subset_metadata)]

df[is.na(df)] <- 0

# Get the list of ranges
ages_list <- list()
for (char in unique(subset_metadata$Age)){
  char_positions <- which(subset_metadata$Age == char)
  ages_list[[char]] <- char_positions
}

calculate_mean_per_list_df <- function(df, named_lists) {
  means <- list()
  for (name in names(named_lists)) {
    columns <- named_lists[[name]]
    means[[name]] <- rowMeans(df[, columns, drop = FALSE])
  }
  # Combine the list of means into a dataframe
  means_df <- do.call(cbind, means)
  # Add column names
  colnames(means_df) <- names(named_lists)
  return(means_df)
}

means_df <- calculate_mean_per_list_df(df, ages_list)


## Run DAISY
df_gower_dist <- cluster::daisy(means_df, metric = "gower")

sil_width <- c(NA)

for(i in 2:15){
  
  pam_fit <- pam(df_gower_dist,
                 diss = TRUE,
                 k = i)
  
  sil_width[i] <- pam_fit$silinfo$avg.width
  
}

sil_width <- data.frame(
  Position = 1:length(sil_width), 
  Value = sil_width)

ggplot(sil_width, aes(x = Position, y = Value)) + geom_point() + geom_line() +
  labs(x = "Number of medoids", y = "Pam average fit") + theme_classic()

ggsave(width = 3, height = 2, filename = "../analysis/plots/Kmedoids_scree_Control.pdf")

## Since the scree drops below 0.5 at 9 clusters, looking for that.

df_pam <- pam(df_gower_dist, diss = TRUE, k = 5)

pam_cluster <- as.data.frame(df_pam$clustering)
colnames(pam_cluster) <- c("pam_cluster")

pam_cluster <- as.data.frame(pam_cluster) %>% 
  mutate(clust = paste("clust_", pam_cluster,sep = ""))
pam_cluster$pam_cluster <- NULL

# Bring the cluster info for plotting
df_with_clust_info <- merge(means_df, pam_cluster, by = 0)

rownames(df_with_clust_info) <- df_with_clust_info$Row.names
df_with_clust_info$Row.names <- NULL

# visualise  each cluster 

mean_expression <- df_with_clust_info %>%
  group_by(clust) %>%
  summarise_all(mean)

dfa <- mean_expression %>% gather(key = "variable", value = "value", -c(1)) %>%  ###the last column is the cluster
  group_by(variable) %>%  
  mutate(row_num =  1:n()) 


df_with_clust_info %>% 
  gather(key = "variable" , value = "value", -c(ncol(df_with_clust_info))) %>%  ###the last column is the cluster
  group_by(variable) %>%  
  mutate(row_num =  1:n()) %>% 
  ggplot(aes(x =  variable , y = value , group = row_num)) +   
  geom_line(alpha = 1, color = "grey") + 
  geom_line(data = dfa, aes(x = variable, y= value, group = row_num), linewidth = 1, color = "black") +
  theme_bw() +  
  theme(legend.position = "none" , axis.text.x = element_text(angle = 45 , vjust = 0.4)) +
  facet_wrap(~clust, ncol = 5)
ggsave(width = 8, height = 3, "../analysis/plots/5medoids_Control.pdf")



###### Repeat everything with Rapa continuous

df <- read.table("../analysis/norm-CLEAN_ASVs_counts_merged_20240530.tsv", sep = "\t", header = T, row.names = 1, stringsAsFactors = F)

# Remove OTUs that arent found in at least 290 samples (out of 300 sequenced)
df <- df[rowSums(is.na(df)) < 295,]

# Subset only those that are there in all 3 timepoints
subset_metadata <- metadata %>%
  #filter(Treatment_SG %in% c("Control", "Rapamycin conti (RC)"))
  filter(Treatment_SG %in% c("Rapamycin conti (RC)"))

subset_metadata$Age <- factor(subset_metadata$Age, levels = c("Young", "Mid", "Old"))
subset_metadata <- subset_metadata %>% arrange(Age)

# Sort the dataframe and keep only columns we are using.
df <- df[rownames(subset_metadata)]

df[is.na(df)] <- 0

# Get the list of ranges
ages_list <- list()
for (char in unique(subset_metadata$Age)){
  char_positions <- which(subset_metadata$Age == char)
  ages_list[[char]] <- char_positions
}

calculate_mean_per_list_df <- function(df, named_lists) {
  means <- list()
  for (name in names(named_lists)) {
    columns <- named_lists[[name]]
    means[[name]] <- rowMeans(df[, columns, drop = FALSE])
  }
  # Combine the list of means into a dataframe
  means_df <- do.call(cbind, means)
  # Add column names
  colnames(means_df) <- names(named_lists)
  return(means_df)
}

means_df <- calculate_mean_per_list_df(df, ages_list)


## Run DAISY
df_gower_dist <- cluster::daisy(means_df, metric = "gower")

sil_width <- c(NA)

for(i in 2:15){
  
  pam_fit <- pam(df_gower_dist,
                 diss = TRUE,
                 k = i)
  
  sil_width[i] <- pam_fit$silinfo$avg.width
  
}

sil_width <- data.frame(
  Position = 1:length(sil_width), 
  Value = sil_width)

ggplot(sil_width, aes(x = Position, y = Value)) + geom_point() + geom_line() +
  labs(x = "Number of medoids", y = "Pam average fit") + theme_classic()

ggsave(width = 3, height = 2, filename = "../analysis/plots/Kmedoids_scree_RC.pdf")

## Since the scree drops below 0.5 at 5 clusters, looking for that.

df_pam <- pam(df_gower_dist, diss = TRUE, k = 5)

pam_cluster <- as.data.frame(df_pam$clustering)
colnames(pam_cluster) <- c("pam_cluster")

pam_cluster <- as.data.frame(pam_cluster) %>% 
  mutate(clust = paste("clust_", pam_cluster,sep = ""))
pam_cluster$pam_cluster <- NULL

# Bring the cluster info for plotting
df_with_clust_info <- merge(means_df, pam_cluster, by = 0)

# Sort by cluster
df_with_clust_info <- df_with_clust_info[order(df_with_clust_info$clust, decreasing = FALSE),]

rownames(df_with_clust_info) <- df_with_clust_info$Row.names
df_with_clust_info$Row.names <- NULL

# visualise  each cluster 

mean_expression <- df_with_clust_info %>%
  group_by(clust) %>%
  summarise_all(mean)

dfa <- mean_expression %>% gather(key = "variable", value = "value", -c(1)) %>%  ###the last column is the cluster
  group_by(variable) %>%  
  mutate(row_num =  1:n()) 


df_with_clust_info %>% 
  gather(key = "variable" , value = "value", -c(ncol(df_with_clust_info))) %>%  ###the last column is the cluster
  group_by(variable) %>%  
  mutate(row_num =  1:n()) %>% 
  ggplot(aes(x =  variable , y = value , group = row_num)) +   
  geom_line(alpha = 1, color = "grey") + 
  geom_line(data = dfa, aes(x = variable, y= value, group = row_num), linewidth = 1, color = "#c34e0d") +
  theme_bw() +  
  theme(legend.position = "none" , axis.text.x = element_text(angle = 45 , vjust = 0.4)) +
  facet_wrap(~clust, ncol = 5)
ggsave(width = 8, height = 3, "../analysis/plots/5medoids_RC.pdf")



###### Repeat everything with Rapa Early

df <- read.table("../analysis/norm-CLEAN_ASVs_counts_merged_20240530.tsv", sep = "\t", header = T, row.names = 1, stringsAsFactors = F)

# Remove OTUs that arent found in at least 290 samples (out of 300 sequenced)
df <- df[rowSums(is.na(df)) < 295,]

# Subset only those that are there in all 3 timepoints
subset_metadata1 <- metadata %>%
  #filter(Treatment_SG %in% c("Control", "Rapamycin conti (RC)"))
  filter(Treatment_SG %in% c("Rapa memory Early (RM E)"))

subset_metadata2 <- metadata %>% 
  filter(Treatment_SG %in% c("Rapamycin conti (RC)")) %>%
  filter(Age %in% c("Young"))

subset_metadata <- rbind(subset_metadata1, subset_metadata2)

subset_metadata$Age <- factor(subset_metadata$Age, levels = c("Young", "Mid", "Old"))
subset_metadata <- subset_metadata %>% arrange(Age)

# Sort the dataframe and keep only columns we are using.
df <- df[rownames(subset_metadata)]

df[is.na(df)] <- 0

# Get the list of ranges
ages_list <- list()
for (char in unique(subset_metadata$Age)){
  char_positions <- which(subset_metadata$Age == char)
  ages_list[[char]] <- char_positions
}

calculate_mean_per_list_df <- function(df, named_lists) {
  means <- list()
  for (name in names(named_lists)) {
    columns <- named_lists[[name]]
    means[[name]] <- rowMeans(df[, columns, drop = FALSE])
  }
  # Combine the list of means into a dataframe
  means_df <- do.call(cbind, means)
  # Add column names
  colnames(means_df) <- names(named_lists)
  return(means_df)
}

means_df <- calculate_mean_per_list_df(df, ages_list)


## Run DAISY
df_gower_dist <- cluster::daisy(means_df, metric = "gower")

sil_width <- c(NA)

for(i in 2:15){
  
  pam_fit <- pam(df_gower_dist,
                 diss = TRUE,
                 k = i)
  
  sil_width[i] <- pam_fit$silinfo$avg.width
  
}

sil_width <- data.frame(
  Position = 1:length(sil_width), 
  Value = sil_width)

ggplot(sil_width, aes(x = Position, y = Value)) + geom_point() + geom_line() +
  labs(x = "Number of medoids", y = "Pam average fit") + theme_classic()

ggsave(width = 3, height = 2, filename = "../analysis/plots/Kmedoids_scree_RE.pdf")

## Since the scree drops below 0.5 at 5 clusters, looking for that.

df_pam <- pam(df_gower_dist, diss = TRUE, k = 5)

pam_cluster <- as.data.frame(df_pam$clustering)
colnames(pam_cluster) <- c("pam_cluster")

pam_cluster <- as.data.frame(pam_cluster) %>% 
  mutate(clust = paste("clust_", pam_cluster,sep = ""))
pam_cluster$pam_cluster <- NULL

pam_cluster <- pam_cluster %>% mutate(across("clust", str_replace, "clust_3", "A"))
pam_cluster <- pam_cluster %>% mutate(across("clust", str_replace, "clust_4", "clust_3"))
pam_cluster <- pam_cluster %>% mutate(across("clust", str_replace, "A", "clust_4"))

# Bring the cluster info for plotting
df_with_clust_info <- merge(means_df, pam_cluster, by = 0)

# Sort by cluster
df_with_clust_info <- df_with_clust_info[order(df_with_clust_info$clust, decreasing = FALSE),]

rownames(df_with_clust_info) <- df_with_clust_info$Row.names
df_with_clust_info$Row.names <- NULL

# visualise  each cluster 

mean_expression <- df_with_clust_info %>%
  group_by(clust) %>%
  summarise_all(mean)

dfa <- mean_expression %>% gather(key = "variable", value = "value", -c(1)) %>%  ###the last column is the cluster
  group_by(variable) %>%  
  mutate(row_num =  1:n()) 


df_with_clust_info %>% 
  gather(key = "variable" , value = "value", -c(ncol(df_with_clust_info))) %>%  ###the last column is the cluster
  group_by(variable) %>%  
  mutate(row_num =  1:n()) %>% 
  ggplot(aes(x =  variable , y = value , group = row_num)) +   
  geom_line(alpha = 1, color = "grey") + 
  geom_line(data = dfa, aes(x = variable, y= value, group = row_num), linewidth = 1, color = "#5a9f68") +
  theme_bw() +  
  theme(legend.position = "none" , axis.text.x = element_text(angle = 45 , vjust = 0.4)) +
  facet_wrap(~clust, ncol = 5)
ggsave(width = 8, height = 3, "../analysis/plots/5medoids_RE.pdf")


###### Repeat everything with Rapa Late

df <- read.table("../analysis/norm-CLEAN_ASVs_counts_merged_20240530.tsv", sep = "\t", header = T, row.names = 1, stringsAsFactors = F)

# Remove OTUs that arent found in at least 290 samples (out of 300 sequenced)
df <- df[rowSums(is.na(df)) < 295,]

# Subset only those that are there in all 3 timepoints
subset_metadata1 <- metadata %>%
  #filter(Treatment_SG %in% c("Control", "Rapamycin conti (RC)"))
  filter(Treatment_SG %in% c("Rapa memory Late (RM L)"))

subset_metadata2 <- metadata %>% 
  filter(Treatment_SG %in% c("Control")) %>%
  filter(Age %in% c("Young", "Mid"))

subset_metadata <- rbind(subset_metadata1, subset_metadata2)

subset_metadata$Age <- factor(subset_metadata$Age, levels = c("Young", "Mid", "Old"))
subset_metadata <- subset_metadata %>% arrange(Age)

# Sort the dataframe and keep only columns we are using.
df <- df[rownames(subset_metadata)]

df[is.na(df)] <- 0

# Get the list of ranges
ages_list <- list()
for (char in unique(subset_metadata$Age)){
  char_positions <- which(subset_metadata$Age == char)
  ages_list[[char]] <- char_positions
}

calculate_mean_per_list_df <- function(df, named_lists) {
  means <- list()
  for (name in names(named_lists)) {
    columns <- named_lists[[name]]
    means[[name]] <- rowMeans(df[, columns, drop = FALSE])
  }
  # Combine the list of means into a dataframe
  means_df <- do.call(cbind, means)
  # Add column names
  colnames(means_df) <- names(named_lists)
  return(means_df)
}

means_df <- calculate_mean_per_list_df(df, ages_list)


## Run DAISY
df_gower_dist <- cluster::daisy(means_df, metric = "gower")

sil_width <- c(NA)

for(i in 2:15){
  
  pam_fit <- pam(df_gower_dist,
                 diss = TRUE,
                 k = i)
  
  sil_width[i] <- pam_fit$silinfo$avg.width
  
}

sil_width <- data.frame(
  Position = 1:length(sil_width), 
  Value = sil_width)

ggplot(sil_width, aes(x = Position, y = Value)) + geom_point() + geom_line() +
  labs(x = "Number of medoids", y = "Pam average fit") + theme_classic()

ggsave(width = 3, height = 2, filename = "../analysis/plots/Kmedoids_scree_RL.pdf")

## Since the scree drops below 0.5 at 5 clusters, looking for that.

df_pam <- pam(df_gower_dist, diss = TRUE, k = 5)

pam_cluster <- as.data.frame(df_pam$clustering)
colnames(pam_cluster) <- c("pam_cluster")

pam_cluster <- as.data.frame(pam_cluster) %>% 
  mutate(clust = paste("clust_", pam_cluster,sep = ""))
pam_cluster$pam_cluster <- NULL

# Bring the cluster info for plotting
df_with_clust_info <- merge(means_df, pam_cluster, by = 0)

# Sort by cluster
df_with_clust_info <- df_with_clust_info[order(df_with_clust_info$clust, decreasing = FALSE),]

rownames(df_with_clust_info) <- df_with_clust_info$Row.names
df_with_clust_info$Row.names <- NULL

# visualise  each cluster 

mean_expression <- df_with_clust_info %>%
  group_by(clust) %>%
  summarise_all(mean)

dfa <- mean_expression %>% gather(key = "variable", value = "value", -c(1)) %>%  ###the last column is the cluster
  group_by(variable) %>%  
  mutate(row_num =  1:n()) 


df_with_clust_info %>% 
  gather(key = "variable" , value = "value", -c(ncol(df_with_clust_info))) %>%  ###the last column is the cluster
  group_by(variable) %>%  
  mutate(row_num =  1:n()) %>% 
  ggplot(aes(x =  variable , y = value , group = row_num)) +   
  geom_line(alpha = 1, color = "grey") + 
  geom_line(data = dfa, aes(x = variable, y= value, group = row_num), linewidth = 1, color = "#4170b0") +
  theme_bw() +  
  theme(legend.position = "none" , axis.text.x = element_text(angle = 45 , vjust = 0.4)) +
  facet_wrap(~clust, ncol = 5)
ggsave(width = 8, height = 3, "../analysis/plots/5medoids_RL.pdf")



