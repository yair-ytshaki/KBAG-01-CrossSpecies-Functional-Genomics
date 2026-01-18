import pandas as pd
import re

# 1. Load the original dataset
# Ensure the file path matches your environment
df = pd.read_csv('chr6_part1.csv.csv')

# 2. Define the 8 CC Founder strains + essential metadata columns
# Note: C57BL/6J is usually the 'Ref.' in Ensembl exports.
cc_founders_cols = [
    'Variant ID', 'Location', 'Class', 'Conseq. Type', 'Ref.',
    'A/J', '129S1/SvImJ', 'NOD/ShiLtJ', 'NZO/HILtJ', 
    'CAST/EiJ', 'PWK/PhJ', 'WSB/EiJ'
]

# Filter columns: keep only those that exist in the current file
df_filtered = df[[col for col in cc_founders_cols if col in df.columns]].copy()

# 3. Function to strip HTML tags and clean special characters
def clean_html_tags(text):
    if pd.isna(text) or not isinstance(text, str):
        return text
    # Remove HTML tags using Regex
    clean = re.sub('<.*?>', '', text)
    # Replace common HTML entities
    clean = clean.replace('&nbsp;', ' ')
    return clean.strip()

# Apply cleaning to every cell in the filtered dataframe
for col in df_filtered.columns:
    df_filtered[col] = df_filtered[col].apply(clean_html_tags)

# 4. Standardization (Optional but recommended for Week 1)
# In Ensembl, '|' often means 'Same as Reference'. 
# '.' often indicates missing data or was a centered placeholder in the UI.
# Depending on your downstream models (Evo 2/AlphaGenome), 
# you might want to replace these with the actual Ref allele or NaN.

# Save the processed file for your project
df_filtered.to_csv('chr6_CC_Founders_Cleaned.csv', index=False)

print(f"Preprocessing Complete. Remaining columns: {list(df_filtered.columns)}")