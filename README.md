# FIRST: Frontline Immunotherapy with Response-Guided Subsequent Treatment

This repository contains analysis code and supplemental materials for the **FIRST** project, an institutional cohort study examining real-world outcomes of frontline immunotherapy, including neoadjuvant and response-guided treatment strategies, in high-risk cutaneous squamous cell carcinoma (CSCC).

## Supplemental Materials and Interactive Workflows

Interactive supplemental materials and reproducible modeling workflows are available at:

**https://themillerlab.github.io/FIRST/**

### Core workflows

- **Response Modeling Workflow**  
  https://themillerlab.github.io/FIRST/generative_response.html

- **Event-Free Survival Workflow**  
  https://themillerlab.github.io/FIRST/generative_efs.html

These documents provide transparent, reproducible implementations of the Bayesian modeling framework used in the manuscript, including prior specification, generative assumptions, simulation-based validation, posterior predictive checks, and model diagnostics.

### Model-specific analyses

- **ORR Model: Linear Dose Specification**  
  https://themillerlab.github.io/FIRST/modeling_orr_linear.html

- **ORR Model: Log-Dose Specification**  
  https://themillerlab.github.io/FIRST/modeling_orr_log.html

- **ORR Model: Categorical Dose Specification**  
  https://themillerlab.github.io/FIRST/modeling_orr_categorical.html

- **Event-Free Survival Model**  
  https://themillerlab.github.io/FIRST/modeling_efs.html

These model-specific pages provide full code, posterior summaries, convergence diagnostics, posterior predictive checks, and figure-generation workflows corresponding to the primary and sensitivity analyses described in the manuscript.

## Project Overview

The overarching goals of FIRST are to:

- characterize response and survival outcomes following frontline immunotherapy,
- evaluate dose-response relationships and potential diminishing returns,
- apply Bayesian causal inference to address confounding and treatment selection,
- and explore response-guided treatment strategies relevant to real-world clinical decision-making.

All analyses emphasize reproducibility, transparency, and explicit modeling assumptions.

## Repository Structure

```text
FIRST/
├── docs/                  # GitHub Pages site and supplemental materials
├── scripts/
│   ├── DAGs/
│   ├── Generative_Models/
│   └── NLP/
└── README.md

```

## `scripts/` Overview

The `scripts/` directory contains the analytic codebase, organized by conceptual purpose rather than chronology.

### `scripts/DAGs/`

Code for defining, visualizing, and interrogating directed acyclic graphs (DAGs) used to guide causal modeling decisions.

Typical contents include:
- specification of assumed causal structures  
- identification of confounders, mediators, and colliders  
- sensitivity analyses around adjustment sets  

These DAGs inform model design across the project and serve as an explicit record of causal assumptions.

---

### `scripts/Generative_Models/`

Bayesian generative models used to estimate treatment effects and dose–response relationships.

This folder includes:
- logistic and survival models  
- dose–response formulations such as linear, diminishing-returns, and categorical threshold models  
- prior sensitivity analyses  
- posterior predictive checks  
- simulation-based validation workflows  

Models are written with an emphasis on interpretability and alignment with clinical questions rather than purely predictive performance.

---

### `scripts/NLP/`

Natural language processing workflows supporting cohort construction and outcome classification from clinical text.

These scripts include:
- text preprocessing and normalization  
- rule-based and model-assisted classification approaches  
- validation and error analysis  
- reproducible pipelines linking raw clinical text to analytic variables  

---

## Reproducibility and Style

- Analyses are primarily written in **R** using tidyverse-compatible workflows  
- Bayesian models rely on explicit priors and generative assumptions  
- Scripts are modular and designed to be run independently where possible  
- Intermediate outputs (figures, tables) are programmatically generated  
- Supplemental HTML documents are rendered from Quarto-based workflows  

---

## Data Availability

Patient-level data are **not publicly available** because of privacy and institutional review constraints.  
Scripts assume access to locally stored, IRB-approved datasets with appropriate permissions.

---

## Project Status

The FIRST project is active and evolving. Folder structure, modeling approaches, and scripts may change as analyses mature and manuscripts develop.
