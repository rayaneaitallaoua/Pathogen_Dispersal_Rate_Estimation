# Experiment 1: Varying Maximal Dispersal Distances

## Goal
Simulate sequence data under different **maximum dispersal distances**, run **Genepop** (IBD), run **BEAST** (phylogenetic BD), and then compute and compare dispersal-rate estimates and regressions across dispersal distance regimes.

> **Repository rule:** all analysis scripts live in a **common** folder.  
> For this experiment, **copy the needed analysis scripts into the experiment directory** (next to the data) **before** running them.

---

## Directory layout (expected)
Set up the directory for the experiment with sub-directories for various sample sizes:
```
Experiment_1_maxDispDist
├── 4seqs   # sample size n=4
├── 20seqs  # sample size n=20
└── 50seqs  # sample size n=50
```
Then within each sub-directory, each simulation replicate lives in a directory named:
```
dist{MAXDIST}_sim{N}
```
(e.g. `dist10_sim1`, `dist20_sim15`, …).  
Many analysis scripts detect such directories via regex, so keep this naming or edit regex.

---

## Dependencies

### External software (add to path)
- **GSpace**  
- **Genepop**  
- **BEAST** (v1.10.5).  
- **SLURM** job scheduler.

### Python packages
- **Required:** `biopython`, `pandas`, `matplotlib`, `scipy`, `seaborn`, `statsmodels`

---

## Preparation
1. **Copy analysis scripts** into each sample size subdirectory of the experiment:
   ```
   GSpace_simulations_per_MaxDist.py
   genepop_analysis.py
   beast_xml_generation_per_maxDist.py
   beast_analysis_ifb_per_maxDist.py
   analyze_results_per_maxDist.py
   ```

2. In the following scripts, set:

in the experiment root directory, copy the following scripts and set
   ```python
   MODE = "maxdist"
   ```
   - `non_phylo_bd.py`
   - `phylo_bd_beast.py`
   - `reg_bd.py`
   - `reg_ibd_genepop.py`

3. Set the sample size in `GSpace_simulations_per_MaxDist.py`:
   ```
   nodes_sampled = {sample size (GSpace only samples one individual per node)}
   ```
---

## Step-by-step pipeline

### 1. Launch GSpace simulations (per max dispersal distance)
in each sample size subdirectory, run the following script
```
GSpace_simulations_per_MaxDist.py
```
- Runs GSpace simulations at different max dispersal distances.  
- Creates subdirectories `dist{maxDist}_sim{n}`.  
- Submits GSpace jobs via SLURM.

>[!IMPORTANT]
> Requires `GSpace` in your PATH.

---

### 2. Prepare and submit Genepop (IBD) jobs
```
genepop_analysis.py
```
- Converts FASTA files to Genepop format (`*_genepop.txt`).  
- Writes Genepop settings file (`genepop.txt`).  
- Submits Genepop jobs with SLURM (`srun Genepop`).

>[!IMPORTANT]
> Requires `Genepop` in your PATH.

---

### 3. Generate BEAST XMLs
```
beast_xml_generation_per_maxDist.py
```
- Scans for `*.fa` in each directory.  
- Creates BEAST v1.10.5 XML files for each dataset.

---

### 4. Submit BEAST jobs
```
beast_analysis_ifb_per_maxDist.py
```
- Submits BEAST jobs for every XML file via SLURM.  
- Uses `srun beast {xml_file}`.

>[!IMPORTANT]
> - Ensure the correct Java module is loaded on your cluster  
> - Ensure the correct BEAGLE version is installed  
> - Requires `BEAST` in your PATH.

---

### 5. Aggregate & estimate diffusion rates

#### 5.a Non-phylogenetic BD
```
non_phylo_bd.py
```
Computes `var(X) / TMRCA`  
Outputs:
- `non_phylo_bd_all.tsv`
- `non_phylo_bd_grouped.tsv`

#### 5.b Phylogenetic BD (BEAST)
```
phylo_bd_beast.py
```
Summarizes BEAST posterior diffusion rates.  
Outputs:  
- `beast_diffusion_all.tsv`  
- `beast_diffusion_grouped.tsv`  
- `beast_diff_rate_per_{radius/maxDispDist}.png`

#### 5.c Regression (BD)
```
reg_bd.py
```
Computes regression slopes: `geo_d² ~ gen_d`.  
Outputs:  
- `regression_bd_slopes.tsv`  
- `regression_bd_summary.tsv`

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

### 6. Collate & summarize BEAST results
```
analyze_results_per_maxDist.py
```
- Aggregates results across all replicates.  
- Outputs grouped TSVs and figures.  
- Default SLURM helper script is provided:  
  ```
  run_analyze_results_per_maxDist.slurm
  ```

---

### 7. Generate comparison plots using `plotting.py`

This script expects a **root** directory that contains **three subfolders**: `4seqs/`, `20seqs/`, `50seqs/`, each with grouped TSVs named:

- `non_phylo_bd_grouped.tsv`  
- `beast_diffusion_grouped.tsv`  
- `regression_bd_summary.tsv`  
- `genepop_invSlope_summary.tsv`

By default it defines:
```
root = "."
```
Edit `root` to point to **your** archive path (or set it to `"."` if you run inside a folder containing `4seqs/`, `20seqs/`, `50seqs/`).  

**Outputs:**  
- `4vs20vs50seqs_beast_quantiles_maxDispDist.png`  
- `4vs20vs50seqs_BD_quantiles_maxDispDist.png`  
- `4vs20vs50seqs_BD_reg_Quantiles_maxDispDist.png`  
- `4vs20vs50seqs_BD_reg_noQuant_maxDispDist.png`  
- `4vs20vs50seqs_IBD_reg_quantiles_maxDispDist.png`  
- `4vs20vs50seqs_IBD_reg_noQuant_maxDispDist.png`

**Dependencies:** `pandas`, `matplotlib`.

---

### 8. Generate comparison statistics using `stats_for_disp_dist.py`

This script also expects the same `root` structure with `4seqs/`, `20seqs/`, `50seqs/` and the same grouped TSVs as **plotting.py**. Edit `root` to your archive path.

**What it does:**
- Harmonizes column names and fits linear models of `estimation ~ maxDispDist` for each method/seq-count, writing results to:
  - `stats_for_disp_dist.txt`
- Builds a heatmap of slopes (cells annotated with slope; `**` marking *p* < 0.05), saved to:
  - `slope_heatmap_max_disp_dist.png`

**Dependencies:** `pandas`, `numpy`, `matplotlib`, `statsmodels`, `seaborn`, `scipy`.

**Run:**
```bash
python stats_for_disp_dist.py
```

**Notes:**
- The heatmap code fixes ordering of the sequence sizes to **4, 20, 50**.  
- Color scale is centered at 0 with `vmin=-0.05` and `vmax=0.05`; adjust if your slopes are outside that range.

---

## Notes & Warnings
- **Set `MODE = "maxdist"`** in: `non_phylo_bd.py`, `phylo_bd_beast.py`, `reg_bd.py`, `reg_ibd_genepop.py`.  
- **Edit `root`** in `plotting.py` and `stats_for_disp_dist.py` to point to the directory containing `4seqs/`, `20seqs/`, `50seqs/`.  
- If you only have one experiment folder (no `4seqs/20seqs/50seqs` breakdown), adapt the scripts or skip the comparison scripts.  
- Ensure external binaries (**GSpace**, **Genepop**, **BEAST**) are installed and accessible on your cluster.  
- Non-standard Python packages required: **Biopython**, **statsmodels**, **seaborn**.
