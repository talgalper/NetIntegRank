# Run from the folder that holds SWISS-MODEL_Repository/ (unpacked from 9606_coords.tar.gz and 9606_meta.tar.gz);
# copies each protein's chosen model to structures_raw/<UniProt ID>.pdb. download_swissmodel_structures.sh runs it.
# (Only base R and progress are used; tidyverse was loaded here before but nothing in it was used.)
library(bio3d)

# read in metadata from SWISS-MODEL download
meta <- read.table("SWISS-MODEL_Repository/INDEX", sep = "\t", header = T)

# subset structure with highest qmeandisco_global score for each uniprot ID
meta_subset <- meta
meta_subset$uniprot_id <- sub("-\\d+$", "", meta_subset$UniProtKB_ac)

#meta_subset <- meta_subset[order(-meta_subset$qmeandisco_global), ]
#meta_subset <- meta_subset[!duplicated(meta_subset$uniprot_id), ]


# function to split uniprot ID into the three folder names SWISS-MODEL stores it under: characters 1-2, 3-4
# and the rest, so P04637 -> P0/46/37 and A0A075B734 -> A0/A0/75B734
split_into_folders <- function(input_string) {
  c(substr(input_string, 1, 2), substr(input_string, 3, 4), substr(input_string, 5, nchar(input_string)))
}

# create location ID to navigate SWISSMODEL directory
library(progress)
pb <- progress_bar$new(
  format = "  Splitting uniprot ID [:bar] :percent eta: :eta",
  total = length(meta_subset$uniprot_id), clear = FALSE)

loc_id <- list()
for (string in meta_subset$uniprot_id) {
  result <- split_into_folders(string)
  loc_id[[string]] <- result
  
  pb$tick()
}
rm(string, result)


## This code moves pdbs by coordinate ID from meta_subset
library(progress)
pb <- progress_bar$new(
  format = "  Copying files [:bar] :percent eta: :eta",
  total = length(unique(meta_subset$uniprot_id)), clear = FALSE)

# Loop over all Uniprot IDs
for (uniprot_id in unique(meta_subset$uniprot_id)) {
  # Filter meta_subset for the current Uniprot ID and seqid > 60
  subset_meta <- meta_subset[meta_subset$uniprot_id == uniprot_id & meta_subset$seqid >= 60, ]
  # If there are no structures with seqid >= 60, skip to the next Uniprot ID
  if (nrow(subset_meta) == 0 || all(is.na(subset_meta))) {  # added OR check recently
    pb$tick()
    next
  }
  # Find the entry with the highest qmeandisco_global score
  best_entry <- subset_meta[which.max(subset_meta$qmeandisco_global), ]
  # Extract the corresponding coordinate_id and directory
  coordinate_id <- best_entry$coordinate_id
  dir <- file.path("SWISS-MODEL_Repository", 
                   loc_id[[uniprot_id]][1], 
                   loc_id[[uniprot_id]][2], 
                   loc_id[[uniprot_id]][3], 
                   "swissmodel")
  # List all files in the directory
  files <- list.files(dir, full.names = TRUE)
  # Loop over the files and copy the one that matches the best coordinate_id
  for (file in files) {
    if (grepl(coordinate_id, basename(file))) {
      file.copy(file, file.path("structures_raw", paste0(uniprot_id, ".pdb")))
      break  # Once the file is copied, break out of the loop
    }
  }
  pb$tick()
}
rm(coordinate_id, dir, file, files, uniprot_id, subset_meta, pb)




