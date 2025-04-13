import os
import pandas as pd

# Get path to desktop folder from home directory
folder_path = os.path.join(os.path.expanduser("~"), "Desktop", "gene_sequences")
print("Looking for files in:", folder_path)

file_list = [os.path.join(folder_path, file) for file in os.listdir(folder_path) if file.endswith('.txt')]
print("Files found:", len(file_list))

data = []  # Create empty list
for file in file_list:
    try:
        temp_df = pd.read_csv(file, sep='\t')  
        data.append(temp_df)  # Append each dataframe to the list
    except Exception as e:
        print(f"Error reading file {file}: {e}")

# Only attempt to concatenate if data list is not empty
if data:
    combined_df = pd.concat(data, ignore_index=True)
    print(combined_df.head())
else:
    print("No data was successfully loaded.")