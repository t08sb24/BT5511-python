import matplotlib.pyplot as plt
import seaborn as sns

# If needed, create a column for Mutation_Class based on MutationType
combined_df['Mutation_Class'] = combined_df['MutationType']

# Create a multi-panel figure (2x2 grid)
fig, axs = plt.subplots(2, 2, figsize=(14, 10))
axs = axs.flatten()  # Flatten to simplify indexing

# ----------------------------------------------------------------------------
# Panel 1: Boxplot of Cell Viability Change by Mutation Class and Promoter Status
# ----------------------------------------------------------------------------
sns.boxplot(
    data=combined_df,
    x='Mutation_Class',
    y='Viability_Change',
    hue='MutationInPromoter',   # Indicates if the mutation is in the promoter region
    palette='Set2',
    ax=axs[0]
)
axs[0].set_title("Viability Change by Mutation Type & Promoter Status")
axs[0].set_xlabel("Mutation Type")
axs[0].set_ylabel("Cell Viability Change")
axs[0].legend(title="In Promoter", loc="upper right")import matplotlib.pyplot as plt
import seaborn as sns

# If needed, create a column for Mutation_Class based on MutationType
combined_df['Mutation_Class'] = combined_df['MutationType']

# Create a multi-panel figure (2x2 grid)
fig, axs = plt.subplots(2, 2, figsize=(14, 10))
axs = axs.flatten()  # Flatten to simplify indexing

# ----------------------------------------------------------------------------
# Panel 1: Boxplot of Cell Viability Change by Mutation Class and Promoter Status
# ----------------------------------------------------------------------------
sns.boxplot(
    data=combined_df,
    x='Mutation_Class',
    y='Viability_Change',
    hue='MutationInPromoter',   # Indicates if the mutation is in the promoter region
    palette='Set2',
    ax=axs[0]
)
axs[0].set_title("Viability Change by Mutation Type & Promoter Status")
axs[0].set_xlabel("Mutation Type")
axs[0].set_ylabel("Cell Viability Change")
axs[0].legend(title="In Promoter", loc="upper right")# Hide the 4th subplot (unused)
# ----------------------------------------------------------------------------
axs[3].axis('off')

# Adjust layout
plt.tight_layout()
plt.show()