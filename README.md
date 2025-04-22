# Estimating Pathogen Dispersal Rates Using Geo-Referenced Genetic Sequences

## Overview
This project is part of a 4-month M1 internship within the Bioinformatics Master's program at the University of Montpellier. The objective is to compare different probabilistic models describing the evolution of homologous genetic sequences over time and space. Specifically, it focuses on evaluating models that jointly infer evolutionary and migratory processes, as well as models that adopt a more generalized approach to migration dynamics.

Simulated datasets will be generated using **GSpace**, a tool designed to model genetic evolution in spatial and temporal contexts. These datasets will then be analyzed using two alternative modeling approaches: the **Relaxed Random Walk (RRW)** model and an **Individual-Based Model (IBM)**. The goal is to assess the ability of each model to accurately estimate pathogen dispersal rates from geo-referenced genetic data.

## Objectives
- Simulate geo-referenced genetic sequences using **GSpace**.
- Apply and compare the **RRW** and **IBM** models on the simulated data.
- Evaluate the strengths and limitations of each modeling approach in estimating dispersal rates.
- Provide a reproducible framework for future analyses of spatial genetic data.

## Methodology
1. **Data Simulation**  
   Use **GSpace** to generate synthetic datasets reflecting joint evolutionary and migratory processes.

2. **Model Application**  
   - Analyze the simulated data using the **RRW model**, which assumes a flexible random walk process for migration.
   - Apply the **IBM model**, which incorporates individual-based migration dynamics.

3. **Comparison & Evaluation**  
   - Compare model outputs in terms of dispersal rate estimation accuracy.
   - Assess computational efficiency and robustness under varying simulation parameters.

## Repository Structure
```
├── docs/                # Documentation and references
├── scripts/             # Analysis and simulation scripts
├── data/                # Simulated datasets (or links/instructions to generate them)
├── results/             # Outputs, figures, and comparative analyses
├── environment.yml      # Conda environment configuration
└── README.md            # Project description and guidelines
```

## Installation
It is recommended to use a Conda environment to manage dependencies.

```bash
git clone https://github.com/your-username/pathogen-dispersal-estimation.git
cd pathogen-dispersal-estimation
conda env create -f environment.yml
conda activate pathogen-dispersal
```

## Usage
1. **Simulate Data**  
   Run the GSpace simulation script:

   ```bash
   python scripts/run_gspace_simulation.py --config configs/simulation_config.yml
   ```

2. **Run Models**  
   Apply RRW and IBM analyses on the generated data:

   ```bash
   python scripts/run_rrw_analysis.py --input data/simulated_sequences.fasta
   python scripts/run_ibm_analysis.py --input data/simulated_sequences.fasta
   ```

3. **Compare Results**  
   Generate comparison reports and visualizations:

   ```bash
   python scripts/compare_models.py --results results/
   ```

## Requirements
- Python 3.x
- GSpace
- BEAST (for RRW model)
- Additional Python libraries (see `environment.yml`)

## References
1. Chalvet-Monfray K, et al. (2021). *GSpace: A framework for simulating genetic evolution in space and time.* [PubMed](https://pubmed.ncbi.nlm.nih.gov/33964130/)
2. Lemey P, et al. (2010). *Phylogeography takes a relaxed random walk in continuous space and time.* [PubMed](https://pubmed.ncbi.nlm.nih.gov/20203288/)
3. Preprint (2024). *An Individual-Based Modeling approach to pathogen dispersal.* [bioRxiv](https://www.biorxiv.org/content/10.1101/2024.06.06.597755v1.abstract)

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Contact
For questions or collaborations, please contact:  
**Rayane Ayoub Ait Allaoua**  
[GitHub](https://github.com/rayaneaitallaoua)
