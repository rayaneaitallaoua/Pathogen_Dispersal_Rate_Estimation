from pathlib import Path
import os
import subprocess

#test

def run_beast(xml_file):
    """Submit BEAST job for a specific XML file"""
    base_name = Path(xml_file).stem
    job_name = f"BEAST_{base_name}"

    slurm_script = f"""#!/bin/bash
module load java-jdk/11.0.9.1

#SBATCH -c 1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=fast
#SBATCH --mem=4G
#SBATCH -t 24:00:00
#SBATCH -J {job_name}
#SBATCH -o {job_name}.%j.out
#SBATCH -e {job_name}.%j.err

srun beast {xml_file}
"""

    script_name = f"{job_name}_submit.sh"
    with open(script_name, "w") as f:
        f.write(slurm_script)

    subprocess.run(["chmod", "+x", script_name])
    result = subprocess.run(["sbatch", script_name], capture_output=True, text=True)
    print(f"Submitted {job_name}: {result.stdout.strip()}")


# Main execution
original_dir = Path.cwd()

# Process all directories containing XML files
for directory in Path('.').glob('*/'):  # Get all subdirectories
    if directory.is_dir():
        os.chdir(directory)

        # Find XML files in current directory
        xml_files = list(Path('.').glob('*.xml'))

        if xml_files:
            print(f"Processing {directory} with {len(xml_files)} XML files")
            for xml in xml_files:
                run_beast(xml_file=xml.name)
        else:
            print(f"Skipping {directory} - no XML files found")

        os.chdir(original_dir)
