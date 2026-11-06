setwd("your data root")
getwd()
load("mRNA.rda")

ls()

class(total_clin_list)
length(total_clin_list)
names(total_clin_list)

out_dir <- "../clinical_csv"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

for (name in names(total_clin_list)) {
  
  clin <- total_clin_list[[name]]
  
  clin <- as.data.frame(clin)
  
  write.csv(
    clin,
    file = file.path(out_dir, paste0(name, "_clinical.csv")),
    row.names = TRUE
  )
}

class(total_expr_list)
length(total_expr_list)
names(total_expr_list)
out_dir <- "../expr_csv"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

for (name in names(total_expr_list)) {
  
  expr <- total_expr_list[[name]]
  
  expr <- as.data.frame(expr)
  
  write.csv(
    expr,
    file = file.path(out_dir, paste0(name, "_expr.csv")),
    row.names = TRUE
  )
}

# install.packages("data.table")
library(data.table)

clin_dir <- "../clinical_csv"

clin_files <- list.files(
  clin_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)

length(clin_files) 

clin_list <- lapply(clin_files, function(f) {
  
  df <- fread(f, data.table = FALSE)
  

  if (colnames(df)[1] == "") {
    rownames(df) <- df[[1]]
    df <- df[, -1, drop = FALSE]
  }
  
  # 记录来源（强烈推荐）
  df$source_file <- basename(f)
  
  return(df)
})

clinical_all <- rbindlist(
  clin_list,
  fill = TRUE
)
write.csv(
  clinical_all,
  file = "../clinical_all_merged.csv",
  row.names = FALSE
)
dim(clinical_all)
colnames(clinical_all)
head(clinical_all)


library(dplyr)
library(tidyr)


gene_list <- read.csv(
  "../common_external_gene_name.csv",
  stringsAsFactors = FALSE
)

gene_template <- data.frame(
  SYMBOL = gene_list$external_gene_name
)


expr_dir <- "../expr_csv"
out_dir  <- "../expr_csv_aligned"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


expr_files <- list.files(
  expr_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)


for (f in expr_files) {

  message("Processing: ", basename(f))


  expr_gene_final <- read.csv(
    f,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


  colnames(expr_gene_final)[1] <- "SYMBOL"


  expr_filtered <- gene_template %>%
    left_join(expr_gene_final, by = "SYMBOL") %>%
    mutate(across(where(is.numeric), ~replace_na(.x, 0))) %>%
    group_by(SYMBOL) %>%
    summarise(across(where(is.numeric), mean), .groups = "drop")


  stopifnot(!any(duplicated(expr_filtered$SYMBOL)))
  stopifnot(!any(is.na(expr_filtered)))


  expr_filtered_round <- expr_filtered
  expr_filtered_round[, -1] <- round(expr_filtered_round[, -1], 3)


  write.csv(
    expr_filtered_round,
    file = file.path(out_dir, basename(f)),
    row.names = FALSE
  )
}
