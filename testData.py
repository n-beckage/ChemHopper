from ChemTools import *
import csv
import pandas as pd
import time

# Read in csv
df = pd.read_csv('C_d1_ccFalse.csv')
df = df.sort_values(by='qed', ascending=False)

# get top 20% of molecules by qed, sorted top down
qed_cutoff = df['qed'].quantile(0.8)
df_top_qed = df[df['qed'] > qed_cutoff]


print(df_top_qed.head(10))

# top 8 molecules by qed, for table
top8 = df_top_qed.head(8)

# Loop through the dataframe and create a list of table rows for each molecule.
table_rows = []
for index, row in top8.iterrows():
    smiles = row['SMILE']
    props = [row[col] for col in top8.columns[1:]]  # Extract all molecular properties except the first column (smiles)
    table_row = [smiles]  # Combine all information into one row
    # creates an empty list but does the insertions
    [table_row.append(prop) for prop in props]
    table_rows.append(table_row)


### make figure and show plot
column_labels = ['SMILE string', 'NHD', 'NHA', 'MWT', 'MLP', 'MMR', 'NAT', 'PSA', 'QED', 'Molecule']

# Define the figure size and font size
fig, ax = plt.subplots(figsize=(10, 5))
font_size = 12

# Create the table
ax.table(cellText=table_rows, colLabels=column_labels, loc='center')
ax.axis('off')

# Adjust the font size of the table
plt.gca().get_xaxis().set_visible(False)
plt.gca().get_yaxis().set_visible(False)
plt.tight_layout()
plt.rc('font', size=font_size)

# Save the plot as a PDF file
plt.savefig('my_table.pdf', bbox_inches='tight')
