# Preparation----
source("Function.R")

pan_survival <- readRDS("Output/Steps/R01/pan_survival.Rds")
pan_anno_df <- readRDS(file = "Output/Steps/R01/pan_anno.Rds")
pathway_list <- readRDS("Output/Steps/R01/msigdb_selected_pathways.Rds")

with_seed(seed = 66, 
          pathway_list_sel <- sample(pathway_list, 60))

tumor_project <- unique(pan_anno_df$Project)
tumor_project_file <- paste0(tumor_project, ".Rds")

args <- commandArgs(T)
index1 <- as.numeric(args[1])
pathway <- pathway_list[[index1]]
genes <- pathway@geneIds

# Run----
all_cost_time <- system.time(
  tumor_list <- lapply(tumor_project, function(tumor){
  # tumor = tumor_project[3]
  tumor_file <- paste0(tumor, ".Rds")
  my_path <- "Output/Steps/R01/pan_RNA_exp/"
  test_exp <- readRDS(paste0(my_path, tumor_file))
  test_survival <- pan_anno_df[pan_anno_df$Project %in% tumor,]
  re <- tryCatch(directMode(test_exp,
                            test_survival,
                            genes, seeds = 60:61), error = function(x){ NA })
  return(re)
  })
)
tumor_list2 <- tumor_list[sapply(tumor_list, is.list)]
tumor_df <- do.call(rbind, tumor_list2)
tumor_df$Training <- as.numeric(tumor_df$Training)
tumor_df$Testing <- as.numeric(tumor_df$Testing)
tumor_df$AUC <- rowMeans(tumor_df[,1:2])
pathway_names <- pathway@setName
print("Finish Pathway!")

# Save
fwrite(tumor_df, file = paste0("Output/output_v1/", pathway_names, ".csv"))

time1 <- as.numeric(all_cost_time[3])
print(paste0("Time Consuming: ", time1))
