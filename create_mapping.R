# Step 2: Load the library
library(biomaRt)

# Step 3: Connect to the Ensembl Plants database
# The 'host' argument is crucial for accessing plant genomes
plant_mart <- useMart(biomart = "plants_mart", host = "https://plants.ensembl.org")

# Step 4: Select the B. rapa dataset to use
# The dataset name for B. rapa is "braparapa_eg_gene"
brapa_dataset <- useDataset(dataset = "braparapa_eg_gene", mart = plant_mart)

# Step 5: Perform the query to get the mapping
# We retrieve the gene symbol ('external_gene_name') and the Ensembl gene ID ('ensembl_gene_id')
# The Ensembl ID is what KEGG uses as its identifier for B. rapa genes.
gene_map <- getBM(
  attributes = c('external_gene_name', 'ensembl_gene_id'),
  mart = brapa_dataset
)

# Step 6: Clean up and save the mapping file
# The result might contain rows with no gene symbol. Let's filter those out.
gene_map_clean <- gene_map[gene_map$external_gene_name != "" & !is.na(gene_map$external_gene_name), ]

# Rename columns for clarity. We call it 'kegg_id' because it's the identifier KEGG uses.
colnames(gene_map_clean) <- c("gene_symbol", "kegg_id")

# Display the first few lines of your new mapping table
print(head(gene_map_clean))

# Save the mapping file for future use
write.csv(gene_map_clean, "brapa_symbol_to_kegg_id_map.csv", row.names = FALSE)

print("Successfully created the mapping file: brapa_symbol_to_kegg_id_map.csv")
