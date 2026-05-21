# Frankneuro-bioinformatics
# Specialized Python Utilities for Neuro-Immune Assay Design and Sequence Processing

## Frank Neuro: Specialized Bioinformatics Toolkit
### Precision CNS R&D | Translational Neuro-Immune Strategy

This repository contains a suite of Python-driven bioinformatics tools developed by Frank Neuro to standardize and de-risk early-stage central nervous system (CNS) research. These tools are designed to bridge the gap between deep-bench neuro-immune intuition and the rigorous data standards required for IND-enabling studies.

## 🔬 About the Developer
**Matthew G. Frank, PhD** – *Specialized Advisor | Founder, Frank Neuro*

With over 100 publications and 11,000+ citations, my work in microglial priming, neuroinflammation, NLRP3 pathways, and alarmins serves as the scientific foundation for these tools. I developed this toolkit to provide the "Red Team" technical auditing required for high-fidelity CNS asset de-risking.

---

## 🛠️ Featured Tools

### 1. [🛡️ qPCR Diagnostic Suite (Live App)](https://frankneuro-qpcr.streamlit.app/)
An interactive data-auditing framework designed to optimize and validate Real-Time quantitative PCR (qPCR) kinetics using raw fluorescence curve principles. 
* **Kinetic Efficiency Overhaul:** Instead of assuming a baseline 100% amplification efficiency ($E = 2.0$), this application evaluates individual well kinetics using Ordinary Least Squares (OLS) linear regression models to compute a validated **Global Master Efficiency** parameter per gene target.
* **Variance Control Integration:** Minimizes mathematical bias across technical replicates by monitoring Standard Deviation and Coefficient of Variation ($CV\%$) directly within a localized 3-cycle log-linear window centered on the fractional $Cq$ point.
* **Plate Health Dashboard:** Evaluates multi-well kinetic consistencies across complete sample-target intersections, generating technical replicate jitter plots and auto-flagging high-variance replicate groups ($SD_{Cq} > 0.3$) that hint at pipetting inconsistencies or tissue-extracted enzyme inhibition.

### 2. [📍 Primer Designer (Live App)](https://frankneuro-primer.streamlit.app/)
A thermodynamic-first approach to DNA primer design, specifically optimized for PCR assays where target specificity in complex biological matrices is paramount.
* **Thermodynamic Accuracy:** Utilizes Nearest-Neighbor (NN) models to calculate melting temperatures ($T_m$), ensuring stability in complex neuro-inflammatory environments.
* **PCR Specificity:** Primers are built using structural index-handling models to target amplicons that span exon-exon boundaries, ensuring absolute exclusion of genomic DNA amplification.
* **Refined Strand Handling:** Features explicit user interface indicators clarifying that reported reverse sequences are reverse complements of the map to eliminate manufacturing order errors.
* **Production-Ready Layout:** Stripped of system-level shutdown routines to guarantee stable web hosting on public or institutional servers.

### 3. [🧬 FASTA Converter (Live App)](https://frankneuro-fasta.streamlit.app/)
A high-integrity utility for sanitizing and converting Genbank and raw sequence data into standardized FASTA formats for downstream bioinformatics pipelines. 
* **Data Integrity:** Automates the removal of artifacts, carriage returns, and formatting inconsistencies that often plague wet-lab sequence exports.
* **Pipeline Ready:** Outputs clean, header-standardized files ready for BLAST or alignment software. The converted sequence delineates exon-intron boundaries with exons in UPPERCASE and introns in lowercase, serving as a direct input for processing in the Primer Designer.
---

## 🚀 Installation & Usage

To use the Frank Neuro toolkit locally:

### 1. Clone the Repository
```bash
git clone [https://github.com/your-username/frankneuro-bio-tools.git](https://github.com/your-username/frankneuro-bio-tools.git)
cd frankneuro-bio-tools
```
**2. Install Dependencies**
Ensure you have Python installed, then install the data processing, visualization, and math libraries:
pip install streamlit pandas numpy plotly openpyxl

**3. Launch the Streamlit Applications**
Each tool functions as an independent, fully reactive web application. Launch any component instantly using the Streamlit execution layer:
streamlit run pcr.py
# OR
streamlit run primer.py
# OR
streamlit run fasta.py

🛡️ "Red Team" Technical Auditing
These tools represent the baseline of the technical oversight I provide. For Venture Capital firms and Biotech leadership, I offer:

Mechanism-of-Action (MoA) Validation: Using custom Python pipelines to audit preclinical data.

Biological "Noise" Identification: Identifying structural flaws in early-stage CNS targets.

Translational Strategy: Bridging molecular neuro-immunology with industrial R&D requirements.

📄 License & Citation
Repository License
This project is licensed under the MIT License—free for academic and industrial use with attribution.

### Academic Software Citation
If you utilize these diagnostic utilities, primer engines, or sequence formatting frameworks to audit molecular data for a peer-reviewed publication, please cite the respective software releases as follows:

*   **The Unified qPCR Module (`pcr.py`):**
    > **Frank, M. G.** (2026). *qPCR Diagnostic Suite: Global Kinetics & Local Precision Analytics for Relative Expression Optimization.* Open-source software framework for kinetic parameter evaluation. Available at: https://frankneuro-qpcr.streamlit.app/. Digital Object Identifier: https://doi.org/10.5281/zenodo.20329609

*   **The Baseline Toolkit (`primer.py` & `fasta.py`):**
    > **Frank, M. G.** (2026). *Frank Neuro Bioinformatics Toolkit: Specialized Python Utilities for Neuro-Immune Assay Design and Sequence Processing*. Digital Object Identifier:  https://doi.org/10.5281/zenodo.19372241

### Underlying Kinetic Reference
> Pfaffl, M. W. (2001). A new mathematical model for relative quantification in real-time RT–PCR. *Nucleic Acids Research*, *29*(9), e45–e45. https://doi.org/10.1093/nar/29.9.e45

📫 Contact
For strategic CNS advisory or to discuss a "Red Team" audit of your neuro-immune portfolio, contact:

Matthew G. Frank, PhD: matthew@frankneuro.com | frankneuro.com
