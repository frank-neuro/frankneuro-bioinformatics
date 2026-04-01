# Frankneuro-bioinformatics
Specialized Python utilities for neuro-immune assay design and sequence processing

Frank Neuro: Specialized Bioinformatics Toolkit 

Precision CNS R&D | Translational Neuro-Immune Strategy

This repository contains a suite of Python-driven bioinformatics tools developed by Frank Neuro to standardize and de-risk early-stage CNS research. These tools are designed to bridge the gap between deep-bench neuro-immune intuition and the rigorous data standards required for IND-enabling studies.

🔬 About the Developer

Matthew G. Frank, PhD Specialized Advisor | Founder, Frank Neuro 

With over 100 publications and 10,000+ citations, my work in microglial priming, neuroinflammation, NLRP3 pathways and alarmins serves as the foundation for these tools. I developed this toolkit to provide the "Red Team" technical auditing required for high-fidelity CNS asset de-risking.

🛠 Featured Tools. 

1. Primer Designer (primer.py). A thermodynamic-first approach to DNA primer design, specifically optimized for PCR assays where target specificity is paramount.
Thermodynamic Accuracy: Utilizes Nearest-Neighbor (NN) models to calculate melting temperatures (Tm), ensuring stability in complex neuro-inflammatory environments.
PCR Specificity: Primers are designed to generate amplicons that span exon/exon boundaries, thereby ensuring exclusion of genomic DNA amplification.
Specific Design: Engineered to minimize primer-dimer formation, non-specific binding in high-sensitivity qPCR/PCR workflows.
Streamlit Interface: Features a web-native UI for rapid iteration in the lab. 

3. Fasta Converter (fasta.py). A high-integrity utility for sanitizing and converting Genbank and raw sequence data into standardized FASTA formats for downstream bioinformatics pipelines. Data Integrity: Automates the removal of artifacts and formatting inconsistencies that often plague "wet-lab" sequence exports. Pipeline Ready: Outputs clean, header-standardized files ready for BLAST or alignment software. The converted sequence delineated exon/intron boundaries with exons in UPPERCASE and introns in lowercase, which serves as an input for processing in Primer Designer.

🚀 Installation & Usage

To use the Frank Neuro toolkit locally:

Clone the repository:

git clone https://github.com/your-username/frankneuro-bio-tools.git

cd frankneuro-bio-tools

Install dependencies:

pip install -r requirements.txt

Launch the Streamlit Apps:

streamlit run primer.py
# OR
streamlit run fasta.py

🛡️ "Red Team" Technical Auditing

These tools represent the baseline of the technical oversight I provide. For Venture Capital firms and Biotech leadership, I offer:

Mechanism-of-Action (MoA) Validation: Using custom Python pipelines to audit preclinical data.

Biological "Noise" Identification: Identifying structural flaws in early-stage CNS targets.

Translational Strategy: Bridging molecular neuro-immunology with industrial R&D requirements.

📄 License & Citation
License: This project is licensed under the MIT License—free for academic and industrial use with attribution.

📫 Contact
For strategic CNS advisory or to discuss a "Red Team" audit of your neuro-immune portfolio, contact:

Matthew G. Frank, PhD: matthew@frankneuro.com
