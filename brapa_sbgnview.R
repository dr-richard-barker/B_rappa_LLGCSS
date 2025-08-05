#-------------------------------------------------------------------------------
# SBGNview Pipeline for Brassica rapa
#
# This script adapts the SBGNview demonstration to create a workflow for
# B. rapa. It finds relevant pathways, downloads the corresponding maps,
# and plots your gene expression data onto all of them automatically.
#-------------------------------------------------------------------------------

# --- 1. SETUP ---
# Load the necessary library
if (!requireNamespace("SBGNview", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("SBGNview")
}
library(SBGNview)

# --- 2. FIND B. RAPA PATHWAYS ---
# We need to find pathways relevant to our organism. We will search for a
# classic Brassica pathway: "glucosinolate". You can change this keyword
# to anything relevant to your study (e.g., "flowering", "photosynthesis").
cat("Finding pathways related to 'glucosinolate'...\n")
brapa_pathways <- findPathways("glucosinolate")

# IMPORTANT CHECK: Stop if no pathways were found.
if (nrow(brapa_pathways) == 0) {
  stop("No pathways found for the keyword. Please try a different keyword.")
}

cat("Found the following pathways:\n")
print(head(brapa_pathways))

# Get the unique IDs of the pathways we found.
brapa_pathway_ids <- unique(brapa_pathways$pathway.id)

# --- 3. DOWNLOAD PATHWAY MAPS ---
# Download the SBGN map files (.sbgn) for all the pathways we found.
# These will be saved to a temporary directory by default.
cat("\nDownloading SBGN files for", length(brapa_pathway_ids), "pathway(s)...\n")
pathways_local_files <- downloadSbgnFile(brapa_pathway_ids)
print(pathways_local_files)


# --- 4. PREPARE YOUR B. RAPA GENE DATA ---
#
# In a real analysis, you would load your own data here.
# Your data should be a data frame or matrix with gene IDs in the first
# column and numeric data (like log2FoldChange) in the subsequent columns.
#
cat("\nLoading B. rapa gene expression data...\n")
DGE_output_table <- read.csv("Brapa_analysis/05-DESeq2_DGE/differential_expression.csv")
brapa_gene_data <- DGE_output_table[,c("KEGG_ID", "Log2fc_(scent)v(no_scent)")]
colnames(brapa_gene_data) <- c("ID", "logFC")
brapa_gene_data <- na.omit(brapa_gene_data)

cat("Using this example data:\n")
print(head(brapa_gene_data))

# --- 5. GENERATE ALL PATHWAY MAPS WITH YOUR DATA ---
#
# This is the core of the pipeline. We will loop through each pathway file
# we downloaded and plot our gene data onto it.
# The output files will be named automatically based on your output file
# prefix and the pathway ID (e.g., "brapa_glucosinolate_P00001.svg").
#
cat("\nGenerating pathway maps...\n")
# The main function call. It will process all pathways provided in 'input.sbgn'.
brapa_sbgnview_obj <- SBGNview(
    gene.data = brapa_gene_data,
    gene.id.type = "kegg", # Use 'kegg' because our example IDs are KEGG IDs
    input.sbgn = pathways_local_files, # Provide all the downloaded SBGN files
    output.file = "brapa_glucosinolate", # A prefix for your output files
    output.formats = c("png", "svg") # You can choose pdf, ps, png, or svg
)

cat("\n--- Workflow Complete ---\n")
cat("Pathway maps have been saved to your working directory with the prefix 'brapa_glucosinolate'.\n")
cat("For example: 'brapa_glucosinolate_R-BTA-75154.png'\n\n")


# --- 6. (Optional) HIGHLIGHTING SPECIFIC NODES ---
#
# If you want to customize a map, you can use the object created above.
# First, you need to find the IDs of the nodes you want to highlight within a
# specific pathway map. These are often not simple gene names.
#
# This example shows how to highlight nodes, but the IDs are placeholders.
# You would need to inspect a pathway map first to find the correct IDs.
#
cat("Demonstrating how to highlight nodes on the first pathway map...\n")

# Let's work with the first pathway object
first_pathway_obj <- brapa_sbgnview_obj +
                        # We use a placeholder ID here. You need to find real IDs from your map.
                        highlightNodes(node.set = c("a placeholder id"),
                                       stroke.width = 4,
                                       stroke.color = "green")

# Change the output file name for this specific highlighted map
outputFile(first_pathway_obj) <- "highlighted_map_example"

# Print the object to generate the file
# Note: this will re-render ONLY the first pathway map with highlighting
print(first_pathway_obj)

cat("A highlighted example map was saved as 'highlighted_map_example...png/svg'\n")


# --- 7. SESSION INFO ---
# Good practice to record the versions of packages used.
sessionInfo()
