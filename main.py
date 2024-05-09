import pandas as pd
import sys

# Step 1

# Read nodes.tsv and edges.tsv into Pandas dataframes
nodes_df = pd.read_csv("nodes.tsv", sep='\t')
edges_df = pd.read_csv("edges.tsv", sep='\t')

# Step 2

# Map function
def mapper(row):
    if row['metaedge'] == 'CbG':
        compound_id = row['source'].split('::')[1][2:]
        gene_id = row['target'].split('::')[1]
        return compound_id, gene_id
    else:
        return None

# Apply the Map function to each row of the dataframe
mapped = edges_df.apply(mapper, axis=1).dropna()

# Convert the mapped results to a dataframe
mapped_df_task2 = pd.DataFrame(mapped.tolist(), columns=['Compound', 'Gene'])

# Function to reduce the mapped results
def reduce_mapped(mapped_df_task2):
    reduced = mapped_df_task2.groupby('Compound').size()
    reduced_df = reduced.reset_index(name='Gene_Count')
    reduced_sorted = reduced_df.sort_values(by='Gene_Count', ascending=False)
    return reduced_sorted

# Reduce the mapped results
reduced_mapped_df_task2 = reduce_mapped(mapped_df_task2)
top_5_compounds = reduced_mapped_df_task2.head(5)

# Display the top 5 compounds with the highest number of genes
print("Top 5 compounds with the highest number of genes:")
print(top_5_compounds)

# Step 3

# Map function for Task 3
def mapper_task3(row):
    if row['metaedge'] == 'DuG':
        disease_id = row['source'].split('::')[1].replace('DOID:', '')
        gene_id = row['target'].split('::')[1]
        return disease_id, gene_id
    else:
        return None

# Apply the Map function for Task 3 to each row of the dataframe
mapped_task3 = edges_df.apply(mapper_task3, axis=1).dropna()


# Convert the mapped results for Task 3 to a dataframe
mapped_df_task3 = pd.DataFrame(mapped_task3.tolist(), columns=['Disease', 'Gene'])


# Function to reduce the mapped results for Task 3
def reduce_mapped_task3(mapped_df_task3):
  reduced = mapped_df_task3.groupby('Disease').size()
  reduced_df = reduced.reset_index(name='Gene_Count')
  reduced_sorted = reduced_df.sort_values(by='Gene_Count', ascending=False)
  return reduced_sorted


# Reduce the mapped results for Task 3
reduced_mapped_task3 = reduce_mapped_task3(mapped_df_task3)
top_5_diseases = reduced_mapped_task3.head(5)

# Display the top 5 diseases with the highest number of upregulated genes
print("Top 5 diseases with the highest number of upregulated genes:")
print(top_5_diseases)

# Step 4
# Mid Square Hash Function for r=3
def mid_square_hash3(key):
    squared = key ** 2
    squared_str = str(squared)
    mid_index = len(squared_str) // 2
    start_index = mid_index - 1
    end_index = mid_index + 2
    hash_value = int(squared_str[start_index:end_index])
    return hash_value

# Mid Square Hash Function for r=4
def mid_square_hash4(key):
    squared = key ** 2
    squared_str = str(squared)
    mid_index = len(squared_str) // 2
    start_index = mid_index - 2
    end_index = mid_index + 2
    hash_value = int(squared_str[start_index:end_index])
    return hash_value

# Create Hash table function that calls hash function r=3
def create_hash_table_r3(mapped_df, is_compound):
    hash_table = []
    key_column = 'Compound' if is_compound else 'Disease'
    for key in mapped_df[key_column]:
        hash_value = mid_square_hash3(int(key))
        if len(hash_table) <= hash_value:
            hash_table.extend([[]] * (hash_value - len(hash_table) + 1))
        hash_table[hash_value].append(key)
    return hash_table

# Create Hash table function that calls hash function r=4
def create_hash_table_r4(mapped_df, is_compound):
    hash_table = []
    key_column = 'Compound' if is_compound else 'Disease'
    for key in mapped_df[key_column]:
        hash_value = mid_square_hash4(int(key))
        if len(hash_table) <= hash_value:
            hash_table.extend([[]] * (hash_value - len(hash_table) + 1))
        hash_table[hash_value].append(key)
    return hash_table

# Calculates the size of the linked lists in a given hash table
def calculate_linked_list_size(hash_table):
    linked_list_size = 0
    for bucket in hash_table:
        for item in bucket:
            linked_list_size += sys.getsizeof(item)
    return linked_list_size

# Generate hash tables for reduced_mapped_df_task2 and reduced_mapped_task3
hash_tables_task2_r3 = create_hash_table_r3(reduced_mapped_df_task2,True)
hash_tables_task3_r3 = create_hash_table_r3(reduced_mapped_task3, False)
hash_tables_task2_r4 = create_hash_table_r4(reduced_mapped_df_task2, True)
hash_tables_task3_r4 = create_hash_table_r4(reduced_mapped_task3, False)


# Calculate the total sizes of hash tables and linked lists to compare
total_storage_task2_r3 = sys.getsizeof(hash_tables_task2_r3) + calculate_linked_list_size(hash_tables_task2_r3)
total_storage_task3_r3 = sys.getsizeof(hash_tables_task3_r3) + calculate_linked_list_size(hash_tables_task3_r3)
total_storage_task2_r4 = sys.getsizeof(hash_tables_task2_r4) + calculate_linked_list_size(hash_tables_task2_r4)
total_storage_task3_r4 = sys.getsizeof(hash_tables_task3_r4) + calculate_linked_list_size(hash_tables_task3_r4)


print("Task 2 r=3:")
print("Hash Table Size:",sys.getsizeof(hash_tables_task2_r3))
print("list size:",calculate_linked_list_size(hash_tables_task2_r3))
print("total storage:", total_storage_task2_r3)
    
print("Task 2 r=4:")
print("Hash Table Size:",sys.getsizeof(hash_tables_task2_r4))
print("list size:",calculate_linked_list_size(hash_tables_task2_r4))
print("total storage:", total_storage_task2_r4)
    
print("Task 3 r=3:")
print("Hash Table Size:", sys.getsizeof(hash_tables_task3_r3))
print("list size:",calculate_linked_list_size(hash_tables_task3_r3))
print("total storage:", total_storage_task3_r3)
    
print("Task 3 r=4:")
print("Hash Table Size:",sys.getsizeof(hash_tables_task3_r4))
print("list size:", calculate_linked_list_size(hash_tables_task3_r4))
print("total storage:", total_storage_task3_r4)


if(total_storage_task2_r3>total_storage_task2_r4):
    print("For task 2, r=4 required less storage.")
else:
    print("For task 2, r=3 required less storage.")
    

if(total_storage_task3_r3>total_storage_task3_r4):
    print("For task 3, r=4 required less storage.")
else:
    print("For task 3, r=3 required less storage.")

# Folding method function
def folding_method(key, digit_size):
    key_str = str(key)
    chunks = [key_str[i:i + digit_size] for i in range(0, len(key_str), digit_size)]
    return sum(int(chunk) for chunk in chunks)

# Create hash table using folding method
def create_folding_hash_table(mapped_df, digit_size, is_compound):
    hash_table = []
    key_column = 'Compound' if is_compound else 'Disease'
    for key in mapped_df[key_column]:
        hash_value = folding_method(int(key), digit_size)
        if len(hash_table) <= hash_value:
            hash_table.extend([[]] * (hash_value - len(hash_table) + 1))
        hash_table[hash_value % 10].append(key)
    return hash_table

# Generate folding hash tables for reduced_mapped_df_task2 and reduced_mapped_task3
folding_hash_table_task2_digit2 = create_folding_hash_table(reduced_mapped_df_task2, 2, True)
folding_hash_table_task2_digit3 = create_folding_hash_table(reduced_mapped_df_task2, 3, True)
folding_hash_table_task3_digit2 = create_folding_hash_table(reduced_mapped_task3, 2, False)
folding_hash_table_task3_digit3 = create_folding_hash_table(reduced_mapped_task3, 3, False)

print("\n\nFoldiung Method:")
print("Task 2 digits=2:")
print("Hash Table Size:",sys.getsizeof(folding_hash_table_task2_digit2))
print("list size:",calculate_linked_list_size(folding_hash_table_task2_digit2))
    
print("Task 2 digits=3:")
print("Hash Table Size:",sys.getsizeof(folding_hash_table_task2_digit3))
print("list size:",calculate_linked_list_size(folding_hash_table_task2_digit3))

print("Task 3 digits=2:")
print("Hash Table Size:", sys.getsizeof(folding_hash_table_task3_digit2))
print("list size:",calculate_linked_list_size(folding_hash_table_task3_digit2))
    
print("Task 3 digits=3")
print("Hash Table Size:",sys.getsizeof(folding_hash_table_task3_digit3))
print("list size:", calculate_linked_list_size(folding_hash_table_task3_digit3))

'''
print("Hash Tables:")
print(hash_tables_task2_r3)
print("\n\n\n\n")
print(hash_tables_task2_r4)
print("\n\n\n\n")
print(hash_tables_task3_r3)
print("\n\n\n\n")
print(hash_tables_task3_r4)
print("\n\n\n\n")
print("Folding Hash Tables:")
print(folding_hash_table_task2_digit2)
print("\n\n\n\n")
print(folding_hash_table_task2_digit3)
print("\n\n\n\n")
print(folding_hash_table_task3_digit2)
print("\n\n\n\n")
print(folding_hash_table_task3_digit3)
print("\n\n\n\n")
'''