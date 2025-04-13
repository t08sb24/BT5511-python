
import numpy as np
import pandas as pd

# -----------------------------
# 1. Compute Averages for Each Measurement
# -----------------------------
combined_df['Avg_mRNA_WT'] = combined_df[['mRNA.Expression.WT.Rep1', 
                                            'mRNA.Expression.WT.Rep2', 
                                            'mRNA.Expression.WT.Rep3']].mean(axis=1)

combined_df['Avg_mRNA_Mut'] = combined_df[['mRNA.Expression.Mut.Rep1', 
                                             'mRNA.Expression.Mut.Rep2', 
                                             'mRNA.Expression.Mut.Rep3']].mean(axis=1)

combined_df['Avg_Protein_WT'] = combined_df[['Protein.Expression.WT.Rep1', 
                                               'Protein.Expression.WT.Rep2', 
                                               'Protein.Expression.WT.Rep3']].mean(axis=1)

combined_df['Avg_Protein_Mut'] = combined_df[['Protein.Expression.Mut.Rep1', 
                                                'Protein.Expression.Mut.Rep2', 
                                                'Protein.Expression.Mut.Rep3']].mean(axis=1)

combined_df['Avg_Viability_WT'] = combined_df[['CellViability.WT.Rep1', 
                                                 'CellViability.WT.Rep2', 
                                                 'CellViability.WT.Rep3']].mean(axis=1)

combined_df['Avg_Viability_Mut'] = combined_df[['CellViability.Mut.Rep1', 
                                                  'CellViability.Mut.Rep2', 
                                                  'CellViability.Mut.Rep3']].mean(axis=1)


# -----------------------------
# 2. Insert Safeguard Code for Division by Zero and Negative Value Validation
# -----------------------------
def safe_divide(numerator, denominator, default=np.nan):
    """
    Safely divide two numbers.
    Returns numerator/denominator if denominator is non-zero;
    otherwise returns the default value (np.nan).
    """
    return numerator / denominator if denominator != 0 else default

# Recalculate fold changes using safe_divide.
combined_df['mRNA_Fold_Change'] = combined_df.apply(
    lambda row: safe_divide(row['Avg_mRNA_Mut'], row['Avg_mRNA_WT']), axis=1)

combined_df['Protein_Fold_Change'] = combined_df.apply(
    lambda row: safe_divide(row['Avg_Protein_Mut'], row['Avg_Protein_WT']), axis=1)

# Viability change is a simple subtraction; division safety is not required here.
combined_df['Viability_Change'] = combined_df['Avg_Viability_Mut'] - combined_df['Avg_Viability_WT']

# Validate for unexpected negative values in wild-type measurements.
negative_mRNA_wt = combined_df[combined_df['Avg_mRNA_WT'] < 0]
if not negative_mRNA_wt.empty:
    print("Warning: Negative mRNA wild-type values detected for genes:")
    print(negative_mRNA_wt[['Gene', 'Avg_mRNA_WT']])

negative_protein_wt = combined_df[combined_df['Avg_Protein_WT'] < 0]
if not negative_protein_wt.empty:
    print("Warning: Negative Protein wild-type values detected for genes:")
    print(negative_protein_wt[['Gene', 'Avg_Protein_WT']])

# -----------------------------
# 3. Compute Scores and Print Summary Statistics
# -----------------------------
combined_df['mRNA_Score'] = abs(combined_df['mRNA_Fold_Change'] - 1)
combined_df['Protein_Score'] = abs(combined_df['Protein_Fold_Change'] - 1)
combined_df['Viability_Score'] = abs(combined_df['Viability_Change'])
combined_df['Total_Score'] = (
    combined_df['mRNA_Score'] * 0.3 + 
    combined_df['Protein_Score'] * 0.3 + 
    combined_df['Viability_Score'] * 0.4
)

# Print summary stats of the fold changes and viability changes.
print("\nSummary of expression and viability changes:")
for col in ['mRNA_Fold_Change', 'Protein_Fold_Change', 'Viability_Change']:
    print(f"\n{col}:")
    print(combined_df[col].describe())


    combined_df['mRNA_Score'] = abs(combined_df['mRNA_Fold_Change'] - 1)
combined_df['Protein_Score'] = abs(combined_df['Protein_Fold_Change'] - 1)
combined_df['Viability_Score'] = abs(combined_df['Viability_Change'])

combined_df['Total_Score'] = (
    combined_df['mRNA_Score'] * 0.3 + 
    combined_df['Protein_Score'] * 0.3 + 
    combined_df['Viability_Score'] * 0.4
)

import difflib

def find_mutation_details(wt_seq, mut_seq):
    """Find the exact positions and changes in mutations"""
    if len(wt_seq) == len(mut_seq):  # Substitution
        positions = []
        changes = []
        for i, (wt_base, mut_base) in enumerate(zip(wt_seq, mut_seq)):
            if wt_base != mut_base:
                positions.append(i)
                changes.append(f"{wt_base}→{mut_base}")
        return positions, changes
    else:
        # For insertions/deletions, use difflib to find differences
        # This is simplified and may need refinement for complex indels
        d = difflib.SequenceMatcher(None, wt_seq, mut_seq)
        positions = []
        changes = []
        for tag, i1, i2, j1, j2 in d.get_opcodes():
            if tag == 'replace':
                positions.append(i1)
                changes.append(f"{wt_seq[i1:i2]}→{mut_seq[j1:j2]}")
            elif tag == 'delete':
                positions.append(i1)
                changes.append(f"{wt_seq[i1:i2]}→deletion")
            elif tag == 'insert':
                positions.append(i1)
                changes.append(f"insertion→{mut_seq[j1:j2]}")
        return positions, changescombined_df['MutationPositions'] = None
combined_df['MutationChanges'] = None 

for idx, row in combined_df.iterrows():
    positions, changes = find_mutation_details(row['WildType.Sequence'], row['Mutant.Sequence'])
    combined_df.at[idx, 'MutationPositions'] = str(positions)
    combined_df.at[idx, 'MutationChanges'] = str(changes)

top_genes = combined_df.sort_values('Total_Score', ascending=False).head(5)
print("\nTop 5 genes for prioritization:")
print(top_genes[['Gene', 'MutationType', 'MutationInPromoter', 
                'mRNA_Fold_Change', 'Protein_Fold_Change', 
                'Viability_Change', 'Total_Score']])
