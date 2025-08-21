# Experiment 2: Sampling Radii

## Goal
Simulate sequence data at multiple **sampling radii**, run **Genepop** (IBD), run **BEAST** (phylogenetic BD), and then compute and compare dispersal-rate estimates and IBD regressions across radii.

> **Repository rule:** all analysis scripts live in a **common** folder.  
> For this experiment, **copy the needed analysis scripts into the experiment directory** (next to the data) **before** running them.

---

## Directory layout (expected)
Set directory for the experiment with sub-directories for various sample size:
```
Experiment_1_samplig_radii
├── 4seqs   # sample size n=4
├── 20seqs  # sample size n=20
└── 50seqs  # sample size n=50
```
Then within each sub-directory, each simulation replicate lives in a directory named:
```
r{RADIUS}_sim{N}
```
(e.g. `r10_sim1`, `r20_sim15`, …).  
Many analysis scripts **detect** such directories via regex, so keep this naming or edit regex.

---

## Dependencies

### External software
- **GSpace** executable (path needed for `srun` in SLURM jobs).  
- **Genepop** executable (called as `Genepop` in SLURM).  
- **BEAST** (v1.10.x) accessible via `srun beast`.  
- **SLURM** job scheduler.

### Python packages
- **Required:** `biopython`, `pandas`, `matplotlib`, `scipy`  
- **Optional (for some plotting/stats scripts):** `seaborn`, `statsmodels`

---

## One-time preparation
1. **Copy analysis scripts** into each sample size subdirectory with the experiment:
   ```
   GSpace_simulations_per_radius.py
   genepop_analysis.py
   beast_xml_generation_per_radius.py
   beast_analysis_ifb_per_radius.py
   non_phylo_bd.py
   phylo_bd_beast.py
   reg_bd.py
   reg_ibd_genepop.py
   analyze_results_per_radius.py
   ```

2. In the following scripts, set:
   ```python
   MODE = "radius"
   ```
   - `non_phylo_bd.py`
   - `phylo_bd_beast.py`
   - `reg_bd.py`
   - `reg_ibd_genepop.py`


3. Set the GSpace path and sample size in `GSpace_simulations_per_radius.py`: 

    ```
    path_gspace = {path to gspace binary}
    nodes_sampled = {sample size (GSpace only samples one individual per node)}
    ```


4. **Check SLURM modules**:
   - `beast_analysis_ifb_per_radius.py` loads `java-jdk/11.0.9.1`.  
     Change this if your cluster uses a different Java module.  
   - `genepop_analysis.py` calls `srun Genepop`.  
     Ensure the `Genepop` binary is in your PATH.

---

## Step-by-step pipeline

### 1. Launch GSpace simulations (per radius)
```
GSpace_simulations_per_radius.py
```
- Generates circular sampling at radii `[10,20,30,40,50]` (defaults).  
- Creates subdirectories `r{radius}_sim{n}`.  
- Submits GSpace jobs via SLURM.

⚠️ Ensure `path_gspace` points to the correct binary.

---

### 2. Prepare and submit Genepop (IBD) jobs
```
 genepop_analysis.py
```
- Converts FASTA files to Genepop format (`*_genepop.txt`).  
- Writes Genepop settings file (`genepop.txt`).  
- Submits Genepop jobs with SLURM (`srun Genepop`).

⚠️ Requires `Genepop` in your PATH.

---

### 3. Generate BEAST XMLs
```
 beast_xml_generation_per_radius.py
```
- Scans for `*.fa` in each directory.  
- Creates BEAST v1.10.5 XML files for each dataset.

---

### 4. Submit BEAST jobs
```
 beast_analysis_ifb_per_radius.py
```
- Submits BEAST jobs for every XML file via SLURM.  
- Uses `srun beast {xml_file}`.

⚠️ Ensure the correct Java module is loaded on your cluster and the correct BEAGLE version is installed.

---

### 5. Aggregate & estimate diffusion rates

#### 5.a Non-phylogenetic BD
```
 non_phylo_bd.py
```
Computes:
- `var(X) / TMRCA`  
- Outputs: `non_phylo_bd_all.tsv`, `non_phylo_bd_grouped.tsv`

#### 5.b Phylogenetic BD (BEAST)
```
 phylo_bd_beast.py
```
Summarizes BEAST posterior diffusion rates.  
Outputs: `beast_diffusion_all.tsv`, `beast_diffusion_grouped.tsv`, `beast_diff_rate_per_{radius/maxDispDist}.png`

#### 5.c Regression (BD)
```
 reg_bd.py
```
Computes regression slopes: `geo_d² ~ gen_d`.  
Outputs: `regression_bd_slopes.tsv`, `regression_bd_summary.tsv`

#### 5.d Genepop IBD
```
 reg_ibd_genepop.py
```
Extracts slopes & confidence intervals from Genepop `.ISO` results.  
Outputs:
- `genepop_results.tsv`
- `genepop_inv_only.tsv`
- `genepop_slope_summary.tsv`
- `genepop_invSlope_summary.tsv`

---

### 6. Generate comparison plots using `plotting.py`

This script expects a **root** directory that contains **three subfolders**: `4seqs/`, `20seqs/`, `50seqs/`, each with grouped TSVs named:

- `non_phylo_bd_grouped.tsv`
- `beast_diffusion_grouped.tsv`
- `regression_bd_summary.tsv`
- `genepop_invSlope_summary.tsv`

By default it defines:

```python
import os

root = "."
```
Edit `root` to point to **your** archive path (or set it to `"."` if you run inside a folder containing `4seqs/`, `20seqs/`, `50seqs/`).  

**Run:**
```bash
python plotting.py
```
**Outputs:**  
- `4vs20vs50seqs_beast_quantiles_radii.png`
- `4vs20vs50seqs_BD_quantiles_radii.png`
- `4vs20vs50seqs_BD_reg_Quantiles_radii.png`
- `4vs20vs50seqs_BD_reg_noQuant_radii.png`
- `4vs20vs50seqs_IBD_reg_quantiles_radii.png`
- `4vs20vs50seqs_IBD_reg_noQuant_radii.png`

**Dependencies:** `pandas`, `matplotlib`.

---

## 7. Generate comparison statistics using `stats_for_various_radii.py`

This script also expects the same `root` structure with `4seqs/`, `20seqs/`, `50seqs/` and the same grouped TSVs as **plotting.py**. Edit `root` to your archive path.

**What it does:**
- Harmonizes column names and **fits linear models** of `estimation ~ radius` for each method/seq-count, writing results to:
  - `stats_for_radii.txt`
- Builds a **heatmap** of linear slopes (cells annotated with slope; `**` marking *p* < 0.05), saved to:
  - `slope_heatmap.png`

**Dependencies:** `pandas`, `numpy`, `matplotlib`, `statsmodels`, `seaborn`, `scipy`.

**Run:**
```bash
python stats_for_various_radii.py
```

**Notes:**
- The heatmap code fixes ordering of the sequence sizes to **4, 20, 50**:
  ```python
  slope_df["seqs"] = pd.Categorical(
      slope_df["seqs"],
      categories=["4 seqs", "20 seqs", "50 seqs"],
      ordered=True
  )
  ```
- Color scale is centered at 0 with `vmin=-0.05` and `vmax=0.05`; adjust if your slopes are outside that range.

---

## Notes & Warnings
- **Set `MODE = "radius"`** in: `non_phylo_bd.py`, `phylo_bd_beast.py`, `reg_bd.py`, `reg_ibd_genepop.py`.  
- **Edit `root`** in `plotting.py` and `stats_for_various_radii.py` to point to the directory containing `4seqs/`, `20seqs/`, `50seqs/`.  
- If you only have one experiment folder (no `4seqs/20seqs/50seqs` breakdown), **adapt the scripts or skip them**. They are designed for cross-sample-size comparisons.  
- Ensure external binaries (**GSpace**, **Genepop**, **BEAST**) are installed and accessible on your cluster (via modules or PATH).