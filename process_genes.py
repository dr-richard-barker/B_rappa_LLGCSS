import csv

with open('brapa_all_genes.tsv', 'r') as infile, open('brapa_symbol_to_kegg_id_map.csv', 'w', newline='') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(['kegg_id', 'symbol'])
    for line in infile:
        parts = line.strip().split('\t')
        if len(parts) > 1:
            kegg_id = parts[0]
            description = parts[3] # The description is in the 4th column
            symbol = ''
            if ';' in description:
                symbol = description.split(';')[0]
            else:
                symbol = description.split(' ')[0]
            writer.writerow([kegg_id, symbol])
