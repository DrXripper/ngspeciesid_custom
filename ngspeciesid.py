import os
import subprocess
import tempfile
import shutil
import pandas as pd

# The top-level directory containing all sample folders
root_dir = os.getcwd()

# Walk through all subfolders
for dirpath, dirnames, filenames in os.walk(root_dir):
    if "sorted.fastq" in filenames:
        folder_name = os.path.basename(dirpath)
        fastq_path = os.path.join(dirpath, "sorted.fastq")

        print(f"🔹 Found sorted.fastq in folder: {folder_name}")

        # Create temp folder
        with tempfile.TemporaryDirectory() as tmpdir:
            print(f"    → Using temp folder: {tmpdir}")

            # Copy sorted.fastq into temp dir as 'sorted.fastq'
            temp_fastq_path = os.path.join(tmpdir, "sorted.fastq")
            shutil.copy(fastq_path, temp_fastq_path)

            # Also copy taxonomy.tsv and species_taxid.fasta into temp dir
            shutil.copy("taxonomy.tsv", os.path.join(tmpdir, "taxonomy.tsv"))
            shutil.copy("species_taxid.fasta", os.path.join(tmpdir, "species_taxid.fasta"))

            # Compose the bash pipeline as a multi-line string
            bash_script = """
            # (1) Convert FASTQ to FASTA
            seqkit fq2fa sorted.fastq -o sorted.fasta

            # (2) Run vsearch
            vsearch --usearch_global sorted.fasta \\
                --db species_taxid.fasta \\
                --id 0.9 \\
                --userout hits.txt \\
                --userfields query+target+id \\
                --threads 8

            # (3) Process vsearch hits and join with taxonomy.tsv
            awk '
                BEGIN {
                    FS=OFS="\\t"
                    # Read taxonomy.tsv into lookup table
                    while ((getline < "taxonomy.tsv") > 0) {
                        if (NR==1) {
                            header = $0
                            continue
                        }
                        tax[$1] = $0
                    }
                }
                {
                    split($2, parts, ":")
                    tid = parts[1]
                    identity = $3
                    if (tid in tax) {
                        print $1, identity, tax[tid]
                    } else {
                        print $1, identity, tid, "UNKNOWN"
                    }
                }
            ' hits.txt > classified_full_tmp.tsv

            # (4) Dynamically generate header
            {
                printf "centroid\\tidentity\\t%s\\n" "$(head -n1 taxonomy.tsv)"
                cat classified_full_tmp.tsv
            } > classified_full.tsv

            # (5) Clean up temp intermediate files
            rm sorted.fasta hits.txt classified_full_tmp.tsv
            """

            # Run the pipeline inside the temp dir
            subprocess.run(bash_script, shell=True, check=True, cwd=tmpdir)

            # -------------------------------
            # STEP 2: Compute abundances
            # -------------------------------
            classified_full_path = os.path.join(tmpdir, "classified_full.tsv")

            df = pd.read_csv(classified_full_path, sep="\t")

            # Replace spaces in species names for consistency
            species_norm = df["species"].str.replace(" ", "_", regex=False)

            # Compute abundance table
            abundance_table = (
                species_norm
                .value_counts(normalize=True)
                .rename_axis("species_norm")
                .reset_index(name="abundance")
            )

            # Map abundance back into main DataFrame
            df["species_norm"] = species_norm
            df = df.merge(abundance_table, on="species_norm", how="left")
            df.drop(columns=["species_norm"], inplace=True)

            # Save intermediate file
            intermediate_path = os.path.join(tmpdir, "classified_full_with_abundance.tsv")
            df.to_csv(intermediate_path, sep="\t", index=False)

            abundance_sum = df["abundance"].sum()
            print(f"    ✅ Total abundance sum: {abundance_sum:.10f}")

            # -------------------------------
            # STEP 3: Remove duplicate species
            # -------------------------------
            df = pd.read_csv(intermediate_path, sep="\t")

            # Drop duplicate species rows
            df_clean = df.drop_duplicates(subset=["species"])

            # Save cleaned TSV with the folder name as file name
            output_filename = f"{folder_name}.tsv"
            cleaned_path = os.path.join(tmpdir, output_filename)
            df_clean.to_csv(cleaned_path, sep="\t", index=False)

            total_abundance = df_clean["abundance"].sum()
            print(f"    ✅ Cleaned TSV ready: {output_filename}. Total abundance sum: {total_abundance:.10f}")

            # Determine target group folder (first 3 letters of folder name)
            group_prefix = folder_name[:3]
            group_folder = os.path.join(root_dir, group_prefix)

            # Create the group folder if it doesn't exist
            os.makedirs(group_folder, exist_ok=True)

            # Copy the TSV into the group folder
            target_path = os.path.join(group_folder, output_filename)
            shutil.copy(cleaned_path, target_path)

            print(f"    📂 TSV saved to: {target_path}")

print("🎉 All sorted.fastq files processed successfully!")
