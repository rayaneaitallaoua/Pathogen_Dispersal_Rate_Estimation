import os
import subprocess

from Bio import SeqIO

def write_genepop_settings(subdir):


    settings = f"""EstimationPloidy=Haploid
Mode=Batch
GenepopInputFile={os.path.basename(subdir)}_genepop.txt
MenuOptions=6.5
GeographicScale=Log
MinimalDistance=1
MantelPermutations=0"""
    with open(os.path.join(subdir, 'genepop.txt'), 'w') as f:
        f.write(settings)
    #print(f"genepop settings written to: {subdir.split('/')[-2]}/{subdir.split('/')[-1]}")
    return "success"

def convert_fasta_to_genepop(subdir_path):
    fasta_file = None
    for f in os.listdir(subdir_path) and subdir.startswith("dist"):
        if f.endswith('.fa'):
            fasta_file = f
            break
    if not fasta_file:
        print(f"No fasta file found in {subdir_path}")
        return

    fasta_path = os.path.join(subdir_path, fasta_file)
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if len(records) < 2:
        print(f"Not enough sequences in {fasta_path}")
        return

    output_lines = []
    seq_length = len(records[1].seq)
    output_lines.append(''.join([f"Pos_{i+1}\n" for i in range(seq_length)]))

    nt_code = {'A': '01', 'T': '02', 'C': '03', 'G': '04'}

    for rec in records[1:]:  # skip ancestral
        header_parts = rec.id.split("_")
        if len(header_parts) < 3:
            print(f"Invalid header format in {rec.id}")
            continue
        x, y = header_parts[9], header_parts[10]
        output_lines.append("pop")
        nt_seq = ' '.join([nt_code.get(nuc.upper(), "00") for nuc in rec.seq])
        output_lines.append(f"{x} {y} , {nt_seq}")

    output_file = os.path.join(subdir_path, f"{os.path.basename(subdir_path)}_genepop.txt")
    with open(output_file, 'w') as f:
        f.write('\n'.join(output_lines))
    #print(f"Genepop file written: {output_file.split('/')[-2]}/{output_file.split('/')[-1]}")
    return "success"

# --- SLURM integration ---
def submit_genepop_slurm_job(job_name="genepop_job"):
    """Submits a SLURM job for Genepop analysis in the current directory.
    The job will run 'genepop.txt'."""
    slurm_script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output=genepop_%j.out
#SBATCH --error=genepop_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=03:00:00

srun Genepop
"""

    with open("run_genepop.slurm", "w") as f:
        f.write(slurm_script)

    subprocess.run(["sbatch", "run_genepop.slurm"])

settings_success = 0
settings_fail = 0
conversion_success = 0
conversion_fail = 0
slurm_success = 0
slurm_fail = 0

for subdir in os.listdir(os.getcwd()):
    subdir_path = os.path.join(os.getcwd(), subdir)
    if os.path.isdir(subdir_path):
        try:
            write_genepop_settings(subdir_path)
            settings_success += 1
        except Exception as e:
            print(f"Settings failed for {subdir_path}: {e}")
            settings_fail += 1

        try:
            convert_fasta_to_genepop(subdir_path)
            conversion_success += 1
        except Exception as e:
            print(f"Conversion failed for {subdir_path}: {e}")
            conversion_fail += 1

        try:
            job_name = f"genepop_{os.path.basename(subdir_path)}"
            os.chdir(subdir_path)
            submit_genepop_slurm_job(job_name)
            os.chdir("..")
            slurm_success += 1
        except Exception as e:
            print(f"SLURM submission failed for {subdir_path}: {e}")
            slurm_fail += 1

print(f"{settings_success} settings written successfully, {settings_fail} failed.")
print(f"{conversion_success} genepop conversions succeeded, {conversion_fail} failed.")
print(f"{slurm_success} SLURM jobs submitted, {slurm_fail} failed.")