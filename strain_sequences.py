import pandas as pd
import os
import re

def generate_strain_sequences_v2(csv_path, fasta_path, output_dir="strains_output"):
    # 1. Load the FULL_MERGED CSV
    df = pd.read_csv(csv_path)
    
    # Define CC founder strains
    strains = ['A/J', '129S1/SvImJ', 'NOD/ShiLtJ', 'NZO/HILtJ', 'CAST/EiJ', 'PWK/PhJ', 'WSB/EiJ']
    
    # Helper to clean HTML remnants still in the "improved" CSV
    def clean_val(val):
        if pd.isna(val) or not isinstance(val, str): return val
        return re.sub('<.*?>', '', val).replace('&nbsp;', ' ').strip()

    # 2. Parse FASTA header for start coordinate (20239701)
    with open(fasta_path, 'r') as f:
        header = f.readline().strip()
        match = re.search(r'GRCm39:6:(\d+):', header)
        window_start_coord = int(match.group(1)) if match else 20239701
        ref_sequence = "".join(line.strip() for line in f).upper()
    
    original_ref_len = len(ref_sequence)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for strain in strains:
        # Work with list for mutable sequence
        strain_seq = list(ref_sequence)
        
        # Sort variants BACKWARDS (high to low) to keep indices stable
        df['pos_int'] = df['Location'].apply(lambda x: int(str(x).split(':')[1]))
        sorted_variants = df.sort_values(by='pos_int', ascending=False)
        
        for _, row in sorted_variants.iterrows():
            pos = row['pos_int']
            local_idx = pos - window_start_coord
            
            # Boundary check
            if local_idx < 0 or local_idx >= len(strain_seq):
                continue
            
            # Clean values for this specific variant
            allele = clean_val(row[strain])
            ref_val = clean_val(row['Ref.'])
            var_class = clean_val(row['Class']).lower()

            # Skip if same as reference
            if allele in ['.', '|'] or allele == ref_val:
                continue

            # --- Logic for different Classes ---
            
            if var_class == 'insertion' or (pd.isna(ref_val) and allele != '-'):
                # INSERTION: Insert between local_idx and local_idx + 1
                # We do NOT remove any reference bases (ref_len = 0)
                # Note: clean_val handles strings like 'T/TT/TTT' - 
                # ensure you choose the specific allele for the strain
                target_idx = local_idx + 1
                strain_seq[target_idx:target_idx] = list(allele)

            elif allele == '-':
                # DELETION: Remove the length of the reference sequence
                ref_len = len(ref_val) if pd.notna(ref_val) else 1
                del strain_seq[local_idx : local_idx + ref_len]

            else:
                # SNP or SEQUENCE ALTERATION:
                # Replace the exact length of ref_val with the allele string
                ref_len = len(ref_val) if pd.notna(ref_val) else 1
                strain_seq[local_idx : local_idx + ref_len] = list(allele)

        # 3. Final Trimming/Balancing to exactly 1MB
        final_str = "".join(strain_seq)
        diff = len(final_str) - original_ref_len
        
        # if diff > 0: # Trim symmetrically
        #     trim_start = diff // 2
        #     trim_end = diff - trim_start
        #     final_str = final_str[trim_start : len(final_str) - trim_end]
        # elif diff < 0: # Pad with 'N' if significantly shorter (rare for 1MB)
        #     final_str = final_str.ljust(original_ref_len, 'N')

        # 4. Save
        file_name = f"{strain.replace('/', '_')}_CHR6.txt"
        with open(os.path.join(output_dir, file_name), 'w') as out_f:
            out_f.write(final_str)
        
        print(f"Created {file_name} - Final Length: {len(final_str)}")

generate_strain_sequences_v2('chr6_FULL_MERGED.csv', 'chr6_ref.fasta')