# Load the biomaRt library
library(biomaRt)

# Set the Ensembl dataset for Brassica rapa
ensembl <- useEnsembl(biomart = "plants_mart",
                      dataset = "braparapa_eg_gene",
                      host = "https://plants.ensembl.org",
                      verbose = TRUE)

# Get all gene symbols and their corresponding Ensembl and KEGG IDs
gene_map <- getBM(attributes = c('external_gene_name', 'ensembl_gene_id', 'kegg_enzyme'),
                  mart = ensembl)

# Write the mapping to a CSV file
write.csv(gene_map, "brapa_symbol_to_kegg_id_map.csv", row.names = FALSE)
