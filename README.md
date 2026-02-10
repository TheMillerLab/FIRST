# FIRST: Frontline Immunotherapy with Response-Guided Subsequent Treatment

This repository contains analysis code for the **FIRST** project, an institutional cohort study examining real-world outcomes of frontline (including neoadjuvant) immunotherapy in high-risk cutaneous squamous cell carcinoma (CSCC).

The overarching goals of FIRST are to:
- Characterize response and survival outcomes following frontline immunotherapy  
- Evaluate **dose–response relationships** and diminishing returns  
- Use **Bayesian causal inference** to address confounding and treatment selection  
- Explore **response-guided treatment strategies** relevant to clinical decision-making  

All analyses emphasize reproducibility, transparency, and explicit modeling assumptions.

---

## Repository Structure

```text
FIRST/
├── scripts/
│   ├── DAGs/
│   ├── Generative_Models/
│   └── NLP/
└── README.md
```

---

## `scripts/` Overview

The `scripts/` directory contains all analytic code, organized by conceptual purpose rather than chronology.

### `scripts/DAGs/`
Code for defining, visualizing, and interrogating **Directed Acyclic Graphs (DAGs)** used to guide causal modeling decisions.

Typical contents include:
- Specification of assumed causal structures  
- Identification of confounders, mediators, and colliders  
- Sensitivity analyses around adjustment sets  

These DAGs inform model design across the project and serve as an explicit record of causal assumptions.

---

### `scripts/Generative_Models/`
Bayesian generative models used to estimate treatment effects and dose–response relationships.

This folder includes:
- Logistic and survival models  
- Dose–response formulations (e.g., linear, diminishing returns, categorical dose effects)  
- Prior sensitivity analyses (skeptical, weakly informative, informative)  
- Posterior predictive checks and model comparison  

Models are written with an emphasis on interpretability and alignment with clinical questions rather than purely predictive performance.

---

### `scripts/NLP/`
Natural language processing workflows supporting cohort construction and outcome classification from clinical text.

Includes:
- Text preprocessing and normalization  
- Rule-based and model-based classification approaches  
- Validation and error analysis  
- Reproducible pipelines linking raw text to analytic variables  

These scripts enable scalable abstraction of radiology, pathology, and clinical narratives relevant to FIRST outcomes.

---

## Reproducibility and Style

- Analyses are primarily written in **R**, using tidyverse-compatible workflows  
- Bayesian models rely on explicit priors and generative assumptions  
- Scripts are modular and designed to be run independently where possible  
- Intermediate outputs (figures, tables) are programmatically generated  

---

## Data Availability

Patient-level data are **not publicly available** due to privacy considerations.  
Scripts assume access to locally stored, IRB-approved datasets with appropriate permissions.

---

## Project Status

The FIRST project is active and evolving.  
Folder structure, modeling approaches, and scripts may change as analyses mature and manuscripts develop.
