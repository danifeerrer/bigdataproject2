import pandas as pd

# Step 1

# Read nodes.tsv and edges.tsv into Pandas dataframes
nodes_df = pd.read_csv("nodes.tsv", sep='\t')
edges_df = pd.read_csv("edges.tsv", sep='\t')

# Step 2

# Map function
def mapper(row):
    if row['metaedge'] == 'CbG':
        compound_id = row['source'].split('::')[1]
        gene_id = row['target'].split('::')[1]
        return compound_id, gene_id
    else:
        return None

# Apply the Map function to each row of the dataframe
mapped = edges_df.apply(mapper, axis=1).dropna()

# Convert the mapped results to a dataframe
mapped_df = pd.DataFrame(mapped.tolist(), columns=['Compound', 'Gene'])

# Function to reduce the mapped results
def reduce_mapped(mapped_df):
    reduced = mapped_df.groupby('Compound').size()
    reduced_df = reduced.reset_index(name='Gene_Count')
    reduced_sorted = reduced_df.sort_values(by='Gene_Count', ascending=False)
    return reduced_sorted

# Reduce the mapped results
top_5_compounds = reduce_mapped(mapped_df).head(5)

# Display the top 5 compounds with the highest number of genes
print("Top 5 compounds with the highest number of genes:")
print(top_5_compounds)
