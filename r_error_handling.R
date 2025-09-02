# Error handling utilities for R scripts
# Source this file at the beginning of other R scripts for better error handling

# Function to safely check if required packages are installed
check_required_packages <- function(packages) {
  missing_packages <- packages[!packages %in% installed.packages()[,"Package"]]
  
  if (length(missing_packages) > 0) {
    cat("Missing required packages:", paste(missing_packages, collapse=", "), "\n")
    cat("Installing missing packages...\n")
    
    # Try to install from CRAN first, then Bioconductor
    for (pkg in missing_packages) {
      tryCatch({
        if (pkg %in% c("DESeq2", "tximport", "ComplexHeatmap", "EnhancedVolcano", 
                       "clusterProfiler", "pathview", "goseq", "fgsea", "enrichplot")) {
          BiocManager::install(pkg, ask = FALSE)
        } else {
          install.packages(pkg, dependencies = TRUE)
        }
        cat("Successfully installed:", pkg, "\n")
      }, error = function(e) {
        cat("Failed to install package:", pkg, "- Error:", e$message, "\n")
      })
    }
  }
}

# Function to safely load libraries with error handling
safe_library <- function(package_name) {
  tryCatch({
    library(package_name, character.only = TRUE)
    return(TRUE)
  }, error = function(e) {
    cat("Error loading package", package_name, ":", e$message, "\n")
    return(FALSE)
  })
}

# Function to check if a file exists and is readable
check_file_exists <- function(file_path, description = "File") {
  if (!file.exists(file_path)) {
    stop(paste(description, "not found:", file_path, 
               "\nPlease ensure the file exists and the path is correct."))
  }
  
  if (!file.access(file_path, mode = 4) == 0) {
    stop(paste(description, "exists but is not readable:", file_path))
  }
  
  cat("✅", description, "found and accessible:", file_path, "\n")
  return(TRUE)
}

# Function to safely create directories
safe_dir_create <- function(dir_path, description = "Directory") {
  tryCatch({
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
      cat("Created", description, ":", dir_path, "\n")
    } else {
      cat(description, "already exists:", dir_path, "\n")
    }
    return(TRUE)
  }, error = function(e) {
    cat("Error creating", description, dir_path, ":", e$message, "\n")
    return(FALSE)
  })
}

# Function to validate data frame columns
validate_columns <- function(df, required_cols, df_name = "data frame") {
  missing_cols <- required_cols[!required_cols %in% colnames(df)]
  
  if (length(missing_cols) > 0) {
    stop(paste("Required columns missing from", df_name, ":", 
               paste(missing_cols, collapse=", "), 
               "\nAvailable columns:", paste(colnames(df), collapse=", ")))
  }
  
  cat("✅ All required columns found in", df_name, "\n")
  return(TRUE)
}

cat("R error handling utilities loaded successfully.\n")