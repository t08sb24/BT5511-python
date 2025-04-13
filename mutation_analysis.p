combined_df.to_csv("combined_gene_data.csv", index=False)

numeric_columns = combined_df.select_dtypes(include=['float64', 'int64']).columns
print(combined_df[numeric_columns].describe())

print("Total number of genes:", combined_df['Gene'].nunique())
print("Dataset shape:", combined_df.shape)

print("\nMissing values per column:")
print(combined_df.isnull().sum())

def analyze_mutation(wt_seq, mut_seq):
    # If the sequences are exactly identical, then it's "no mutation".
    if wt_seq == mut_seq:
        return "no mutation"
    
    # If sequence lengths differ, it's an insertion or deletion.
    if len(wt_seq) > len(mut_seq):
        return "deletion"
    elif len(wt_seq) < len(mut_seq):
        return "insertion"
    else:
        # Both sequences have the same length but are not identical,
        # so we perform a base-by-base check.
        differences = sum(1 for wt_base, mut_base in zip(wt_seq, mut_seq) if wt_base != mut_base)
        # If there are any differences, it's a substitution.
        if differences > 0:
            return "substitution"
        else:
            # Although unlikely due to the earlier equality check, fallback to "no mutation".
            return "no mutation"

# Function to determine if mutation is in the promoter region.
def is_in_promoter(wt_seq, mut_seq):
    # Promoter is defined as the first 1000 nucleotides.
    wt_promoter = wt_seq[:1000]
    mut_promoter = mut_seq[:1000]
    return wt_promoter != mut_promoter

# Apply the functions to create new columns in your DataFrame.
combined_df['MutationType'] = combined_df.apply(lambda row: analyze_mutation(
    row['WildType.Sequence'], row['Mutant.Sequence']), axis=1)# Apply the functions to create new columns in your DataFrame.
combined_df['MutationType'] = combined_df.apply(lambda row: analyze_mutation(
    row['WildType.Sequence'], row['Mutant.Sequence']), axis=1)

combined_df['MutationInPromoter'] = combined_df.apply(lambda row: is_in_promoter(
    row['WildType.Sequence'], row['Mutant.Sequence']), axis=1)

# Count the mutation types.
print("\nMutation types:")
print(combined_df['MutationType'].value_counts())

# Count mutations in the promoter region vs. coding.
print("\nMutations in promoter region:")
print(combined_df['MutationInPromoter'].value_counts())



