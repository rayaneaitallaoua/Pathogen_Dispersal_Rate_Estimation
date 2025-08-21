# Estimating Pathogen Dispersal Rates Using Geo-Referenced Genetic Sequences

## Overview
This project is part of a 4-month M1 internship within the Bioinformatics Master's program at the University of Montpellier. The objective is to compare different probabilistic models describing the evolution of homologous genetic sequences over time and space. Specifically, it focuses on evaluating models that jointly infer evolutionary and migratory processes, as well as models that adopt a more generalized approach to migration dynamics.

Simulated datasets will be generated using **[GSpace](https://www1.montpellier.inra.fr/CBGP/software/gspace/index.html)**, a tool designed to model genetic evolution in spatial and temporal contexts. These datasets will then be analyzed using two alternative modeling approaches: the **Relaxed Random Walk (RRW)** model, **Isolation By Distance (IBD)** and **Brownian Diffusion (BD)**. The goal is to assess the ability of each model to accurately estimate pathogen dispersal rates from geo-referenced genetic data.

## Objectives
- Simulate geo-referenced genetic sequences using **GSpace**.
- Apply and compare the **RRW**, **IBD** and **BD** models on the simulated data.
- Evaluate the strengths and limitations of each modeling approach in estimating dispersal rates.
- Provide a reproducible framework for future analyses of spatial genetic data.

## Methodology
1. **Data Simulation**  
   Use **GSpace** to generate synthetic datasets reflecting joint evolutionary and migratory processes.

2. **Model Application**  
   - Analyze the simulated data using the **BD** model via regression and analytics approaches.
   - Analyze the simulated data using the **RRW model** via [BEAST](https://beast.community/index.html)
   - Analyze the simulated data using the **IBD** model using a regression approach via [Genepop](https://gitlab.mbb.univ-montp2.fr/francois/genepop)

3. **Comparison & Evaluation**  
   - Compare model outputs in terms of dispersal rate estimation accuracy.
   - Assess computational efficiency and robustness under varying simulation parameters.

## Repository Structure
```
├── dispersal_rate_per_maxDist/      # Experiment n°1: varying maximal dispersal distances
├── dispersal_rate_per_radius/       # Experiment n°2: varying the sampling radius
└── Analysis/                        # Analysis scripts
```

## Installation
It is recommended to use a Conda environment to manage dependencies.

```bash
git clone git@github.com:rayaneaitallaoua/Pathogen_Dispersal_Rate_Estimation.git
cd pathogen_dispersal_rate_estimation
conda env create -f environment.yml
conda phylo_disp
```

## Usage

The project was primarly ran on the [IFB core cluster](https://www.france-bioinformatique.fr/cluster-ifb-core/).

In each experiment directory are the following type of files:
   - GSpace simulation script.
   - BEAST configuration XML script.
   - BEAST Analysis script.
   - Genepop Analysis script
   - Results Analysis
   - A plots script
   - A statistical analysis script

> **_NOTE:_**  Plots and statistical analysis were done locally.

## Requirements
- Python 3.9
- Java-jdk 11.0.9.1
- GSpace 
- Genepop v4.8.3
- BEAST v10.5.0-beta5
- BEAGLE v4.0.1
- Additional Python libraries (see `environment.yml`)

## References
1. Virgoulay T, et al. (2021). *GSpace: an exact coalescence simulator of recombining genomes under isolation by distance.* [PubMed](https://pubmed.ncbi.nlm.nih.gov/33964130/)
2. Lemey P, et al. (2010). *Phylogeography takes a relaxed random walk in continuous space and time.* [PubMed](https://pubmed.ncbi.nlm.nih.gov/20203288/)
3. Rousset F, et al. (2008). *genepop’007: a complete re-implementation of the genepop software for Windows and Linux.* [PubMed](https://pubmed.ncbi.nlm.nih.gov/21585727/)
3. AIT ALLAOUA, R. (2024). “Projet bibliographique M1 : Estimer la vitesse de dispersion de pathogènes à partir de la comparaison de séquences génétiques géo-référencées”. [Rapport](https://github.com/rayaneaitallaoua/Pathogen_Dispersal_Rate_Estimation/blob/main/Projet_biblio__Rayane_2024.pdf)
4. AIT ALLAOUA, R. (2025). “Rapport de stage M1 : Estimer la vitesse de dispersion de pathogènes à partir de la comparaison de séquences génétiques géo-référencées”.[Rapport](https://github.com/rayaneaitallaoua/Pathogen_Dispersal_Rate_Estimation/blob/main/rapport_de_stage_Rayane.pdf)

## Contact
For questions or collaborations, please contact:  
**Rayane Ayoub Ait Allaoua**  
[GitHub](https://github.com/rayaneaitallaoua)