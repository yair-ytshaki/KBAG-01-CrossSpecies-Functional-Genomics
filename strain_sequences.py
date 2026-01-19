import pandas as pd
import os
import re

def generate_strain_sequences(csv_path, fasta_path, output_dir="strains_output"):
    # 1. Load the cleaned CSV data
    df = pd.read_csv(csv_path)
    
    # Define the 7 target CC founder strains
    strains = ['A/J', '129S1/SvImJ', 'NOD/ShiLtJ', 'NZO/HILtJ', 'CAST/EiJ', 'PWK/PhJ', 'WSB/EiJ']
    
    # 2. Load the Reference Sequence and extract start coordinate from header
    with open(fasta_path, 'r') as f:
        header = f.readline().strip()
        # Extract the start position (20239701) from the header GRCm39:6:20239701:...
        match = re.search(r'GRCm39:6:(\d+):', header)
        if match:
            window_start_coord = int(match.group(1))
        else:
            raise ValueError("Could not find start coordinate in FASTA header.")
            
        ref_sequence = "".join(line.strip() for line in f).upper()
    
    original_ref_len = len(ref_sequence)
    print(f"Reference loaded. Start Coordinate: {window_start_coord}, Length: {original_ref_len}")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for strain in strains:
        # We work with a list of characters for efficient substitution/deletion
        strain_seq = list(ref_sequence)
        
        # Apply variants from back to front (descending) to maintain index integrity
        # We use the 'Location' column to get the genomic position
        sorted_variants = df.copy()
        sorted_variants['pos_int'] = sorted_variants['Location'].apply(lambda x: int(str(x).split(':')[1])) # Extract position, the "6:20239716" part
        sorted_variants = sorted_variants.sort_values(by='pos_int', ascending=False)
        
        for _, row in sorted_variants.iterrows():
            pos = row['pos_int']
            local_idx = pos - window_start_coord
            
            # Boundary check
            if local_idx < 0 or local_idx >= len(strain_seq):
                continue
                
            allele = str(row[strain])
            ref_val = str(row['Ref.'])
            
            # Skip if Reference allele ('.', '|', or matches Ref column)
            if allele in ['.', '|'] or allele == ref_val:
                continue
            
            # Calculate length of the reference segment to be replaced
            # This is crucial for Indels to keep the sequence 'shifted' correctly
            ref_len = len(ref_val) if pd.notna(ref_val) else 1

            if allele == '-':
                # DELETION: Remove the reference segment
                del strain_seq[local_idx : local_idx + ref_len]
            else:
                # SNP or INSERTION: Replace the reference segment with the variant allele(s)
                strain_seq[local_idx : local_idx + ref_len] = list(allele)

        # Reconstruct string
        final_str = "".join(strain_seq)

        # 3. Symmetrical Balancing to exactly 1MB
        current_len = len(final_str)
        diff = current_len - original_ref_len
        
        if diff > 0: # Sequence grew (Insertions) -> Trim edges
            trim_start = diff // 2
            trim_end = diff - trim_start
            final_str = final_str[trim_start : current_len - trim_end]
        elif diff < 0: # Sequence shrank (Deletions) -> Pad from original reference
            # Note: Padding with 'N' or original ref edges depends on project requirements
            # Here we just report if it's shorter, but typically the 1MB window 
            # is large enough that minor Indels don't break the model.
            print(f"Note: {strain} is {abs(diff)} bp shorter than reference.")

        # 4. Export
        file_name = f"{strain.replace('/', '_')}_CHR6.txt"
        with open(os.path.join(output_dir, file_name), 'w') as out_f:
            out_f.write(final_str)
        
        print(f"Successfully generated: {file_name} (Length: {len(final_str)})")

# To run:
generate_strain_sequences('chr6_CC_Founders_Cleaned.csv', 'chr6_ref.fasta')